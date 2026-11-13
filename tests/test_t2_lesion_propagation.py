import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pytest
import SimpleITK as sitk

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))

from anatomy_review import file_record, save_new
from structural_multimodal import PARAMETERS, compose_pull, resample
from t2_lesion_propagation import SCHEMA, plan, run
from t2_stable import ENGINE_HASHES, PROFILE


@pytest.fixture
def frozen_run(tmp_path):
    def write_image(name, image):
        path = tmp_path / name
        sitk.WriteImage(image, str(path))
        return path

    source = sitk.GetImageFromArray(np.arange(6 * 8 * 8, dtype=np.float32).reshape((6, 8, 8)))
    source[4, 1, 2] = 0
    mask = sitk.Image(source.GetSize(), sitk.sitkUInt8)
    mask[4, 3, 2] = 1
    lesion = write_image("lesion.nii.gz", mask)
    mni_reference = sitk.Image(source)
    mni_reference.SetOrigin((0., 0., 1.))
    paths = {"source": write_image("source.nii.gz", source),
             "t1": write_image("t1.nii.gz", source),
             "mni": write_image("mni.nii.gz", mni_reference)}
    paths["source_mask"] = write_image("hd_bet.nii.gz", mask)
    paths["t1_mask"] = write_image("t1_mask.nii.gz", mask)
    for key in ("anatomy_review", "normalization_provenance", "mni_review"):
        paths[key] = tmp_path / f"{key}.json"
        save_new(paths[key], {"synthetic_fixture": key})
    t1_pull = sitk.AffineTransform([0., -1., 0., 1., 0., 0., 0., 0., 1.], [7., .4, 0.])
    mni_pull = sitk.TranslationTransform(3, [.4, 2., 0.])
    paths["mni_to_t1_pull"] = tmp_path / "mni_to_t1.h5"
    sitk.WriteTransform(mni_pull, str(paths["mni_to_t1_pull"]))
    identity = {"schema": "neo-multimodal-structural-1", "modality": "T2w",
                "files": {key: file_record(path) for key, path in paths.items()},
                "parameters": dict(PARAMETERS, N4_zero_policy=PROFILE["N4_zero_policy"]),
                "code_sha256": ENGINE_HASHES["structural_multimodal.py"],
                "runtime": {"SimpleITK": sitk.Version_VersionString()}}
    branch = "raw_independent_HDBET"
    parent = tmp_path / "structural"
    parent.mkdir()
    rigid = parent / f"{branch}_T1_to_source_pull.h5"
    sitk.WriteTransform(t1_pull, str(rigid))
    sitk.WriteImage(resample(source, source, t1_pull), str(parent / f"{branch}_T1.nii.gz"))
    sitk.WriteImage(resample(source, mni_reference, compose_pull(t1_pull, mni_pull)),
                    str(parent / f"{branch}_MNI_once.nii.gz"))
    fingerprint = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
    save_new(parent / "request.json", {"identity": identity, "fingerprint": fingerprint})
    save_new(parent / "receipt.json", {
        "status": "COMPUTED_REVIEW_PENDING", "branches": [{"branch": branch}],
        "fingerprint": fingerprint,
        "outputs": [{**file_record(path), "path": path.name} for path in parent.iterdir()],
    })
    return {"schema": SCHEMA, "structural_run": str(parent), "branch": branch,
            "lesion": str(lesion), "output": str(tmp_path / "propagated")}


def test_t2_propagation_composes_before_sampling_and_keeps_source_fov(frozen_run):
    prepared = plan(frozen_run)
    out = Path(prepared["output"])
    assert not out.exists()
    receipt = run(frozen_run)
    outputs = receipt["outputs"]
    assert receipt["native_lesion_voxels"] == 1
    assert set(outputs) == {"T1_lesion", "T1_source_fov", "MNI_lesion", "MNI_source_fov", "overlay"}
    t1 = sitk.ReadImage(outputs["T1_lesion"]["path"])
    mni = sitk.ReadImage(outputs["MNI_lesion"]["path"])
    t1_fov = sitk.ReadImage(outputs["T1_source_fov"]["path"])
    mni_fov = sitk.ReadImage(outputs["MNI_source_fov"]["path"])
    assert t1[3, 3, 2] == 1
    assert t1_fov[3, 3, 2] == 1
    assert mni.GetOrigin() == mni_fov.GetOrigin() == (0., 0., 1.)
    assert mni[2, 1, 1] == 1
    # Interpolating the intermediate T1 label would instead choose x=3.
    assert mni[3, 1, 1] == 0
    assert mni[0, 1, 1] == 0 and mni_fov[0, 1, 1] == 1
    assert mni[0, 6, 1] == 0 and mni_fov[0, 6, 1] == 0
    assert receipt["metrics"]["MNI"]["source_fov_fraction"] == pytest.approx((42 * 5) / (64 * 6))
    html = (out / "overlay.html").read_text(encoding="utf-8")
    payload = html.split('id="payload" type="application/json">', 1)[1].split("</script>", 1)[0]
    assert [space["name"] for space in json.loads(payload)] == ["原生 T2", "T1", "MNI"]
    with pytest.raises(FileExistsError):
        run(frozen_run)


def test_t2_propagation_refuses_a_changed_frozen_transform(frozen_run):
    prepared = plan(frozen_run)
    path = Path(prepared["inputs"]["t1_to_t2_pull"]["path"])
    with path.open("ab") as stream:
        stream.write(b"changed")
    with pytest.raises(ValueError, match="hash mismatch"):
        run(frozen_run)
    assert not Path(frozen_run["output"]).exists()


def test_t2_propagation_refuses_lesion_with_shifted_grid(frozen_run):
    path = frozen_run["lesion"]
    lesion = sitk.ReadImage(path)
    lesion.SetOrigin((1., 0., 0.))
    sitk.WriteImage(lesion, path)
    with pytest.raises(ValueError, match="same physical grid"):
        run(frozen_run)
    assert not Path(frozen_run["output"]).exists()
