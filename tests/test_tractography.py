import json
from hashlib import sha256
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest
import SimpleITK as sitk
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))

from tractography import workflow
from tractography.demo import demo
from spatial_contract import MNI2009C, MNI6, REVIEW_SCHEMA
from tractography.geometry import mapped_label
from tractography.screening import core_mask, ordered, slab_indices


@pytest.fixture(scope="module")
def example(tmp_path_factory):
    return {key: Path(path) for key, path in demo(tmp_path_factory.mktemp("tractography") / "example").items()}


def test_selected_lines_and_density_use_same_indices(example, tmp_path):
    results = example["results"]
    strict = "CST_PLIC_CP_Tang_pons"
    expanded = "CST_PLIC_plus2_CP_Tang_pons"
    result, run, root = workflow.result_run(results, "ST_iFOD2_L_M1", expanded)
    strict_ids = workflow.indices_for(root, run, strict)
    selected_ids = workflow.indices_for(root, run, expanded)
    assert len(strict_ids) == 5 and len(selected_ids) == 10
    assert set(strict_ids) < set(selected_ids)

    output = tmp_path / "selected"
    workflow.export(results, run["id"], expanded, output)
    np.testing.assert_array_equal(np.load(output / "indices.npy"), selected_ids)
    parent = nib.streamlines.load(run["parent"]["path"]).streamlines
    selected = nib.streamlines.load(output / "selected.tck").streamlines
    for i, actual in zip(selected_ids, selected, strict=True):
        np.testing.assert_array_equal(actual, parent[i])

    visits = nib.load(output / "visits.nii.gz")
    core = nib.load(output / "core_15pct.nii.gz")
    visit_values = np.asanyarray(visits.dataobj)
    assert visit_values.max() <= len(selected_ids)
    assert np.count_nonzero(visit_values) == run["gates"][expanded]["visited_voxels"]
    np.testing.assert_array_equal(core.dataobj, core_mask(
        visit_values, slab_indices(visits.shape, visits.affine, result["slab_mm"])))
    html = example["report"].read_text(encoding="utf-8")
    assert "CST_PLIC_plus2_CP_Tang_pons" in html
    assert html.count("data:image/png;base64,") == len(result["runs"])
    with pytest.raises(ValueError, match="Unavailable"):
        workflow.result_run(results, "Tensor_Det_L_M1", expanded)


def test_order_and_core_definition():
    assert not ordered({"a": np.array([3]), "b": np.array([1, 3])}, ("a", "b"))
    assert ordered({"a": np.array([3]), "b": np.array([1, 4])}, ("a", "b"))
    visits = np.array([100, 15, 14, 0, 20, 3, 2])
    slabs = np.array([0, 0, 0, 0, 1, 1, 1])
    np.testing.assert_array_equal(core_mask(visits, slabs), [1, 1, 0, 0, 1, 1, 0])


def test_native_atlas_mapping_uses_one_composed_pull(tmp_path, monkeypatch):
    image = sitk.Image([10, 10, 10], sitk.sitkUInt8)
    image.SetOrigin((1., -2., 0.))
    image.SetDirection((0., -1., 0., 1., 0., 0., 0., 0., 1.))
    ref = tmp_path / "b0.nii.gz"
    sitk.WriteImage(image, str(ref))
    image[3, 4, 1] = 1
    label = tmp_path / "atlas.nii.gz"
    sitk.WriteImage(image, str(label))
    middle = sitk.AffineTransform([0., -1., 0., 1., 0., 0., 0., 0., 1.], [-4., 5., 0.])
    bridge = sitk.TranslationTransform(3, [1., 0., 0.])
    paths = [tmp_path / "t1.h5", tmp_path / "bridge.h5"]
    for transform, path in zip((middle, bridge), paths, strict=True):
        sitk.WriteTransform(transform, str(path))
    ras = np.eye(4)
    ras[0, 3] = 2
    monkeypatch.setattr("spatial_contract.BRIDGE_SHA256", workflow.record(paths[1])["sha256"])
    visit = {"dataset": "fixture", "subject": "synthetic", "session": "01"}
    grid_roles = ("reference", "t1_reference", "mni2009c_reference", "atlas_reference")
    config = {"label": workflow.record(label), **{role: workflow.record(ref) for role in grid_roles},
              "source_space": MNI6, "visit": visit,
              "anatomy": {"dwi_to_t1_ras": ras.tolist(),
                          "mni2009c_to_t1_image_h5": workflow.record(paths[0]),
                          "mni6_to_mni2009c_image_h5": workflow.record(paths[1])}}
    directions = {"dwi_to_t1_ras": ["DWI", "T1w"], "mni2009c_to_t1_image_h5": ["T1w", MNI2009C],
                  "mni6_to_mni2009c_image_h5": [MNI2009C, MNI6]}
    points = [(-3., 1., 1.), (-4., 1., 1.), (-3., 2., 2.)]
    point_maps = (lambda x, y, z: (x - 2, y, z), lambda x, y, z: (-y - 4, x + 5, z),
                  lambda x, y, z: (x + 1, y, z))
    landmarks = tmp_path / "landmarks.json"
    grid_sha = {role: config[role]["sha256"] for role in grid_roles}
    landmarks.write_text(json.dumps({"visit": visit, "grids_sha256": grid_sha, "segments": {
        role: [dict(name=str(i), from_lps_mm=point, to_lps_mm=transform(*point))
               for i, point in enumerate(points)]
        for role, transform in zip(directions, point_maps, strict=True)}}))
    config["landmarks"] = workflow.record(landmarks)
    review = tmp_path / "review.json"
    decision = {"schema": REVIEW_SCHEMA, "decision": "APPROVED", "visit": visit,
                                  "scope": "cst_atlas_mapping", "point_directions": directions,
                                  "reviewer": "Synthetic fixture", "basis": "Known test coordinates",
                                  "reviewed_utc": "2000-01-01T00:00:00Z",
                                  "landmark_max_mm": {role: 0.01 for role in directions},
                                  "evidence": {"landmarks": config["landmarks"]},
                                  "input_sha256": {**grid_sha, "label": config["label"]["sha256"],
                                      "landmarks": config["landmarks"]["sha256"],
                                      "dwi_to_t1_ras": sha256(ras.astype("<f8").tobytes()).hexdigest(),
                                      **{role: config["anatomy"][role]["sha256"] for role in list(directions)[1:]}}}
    review.write_text(json.dumps(decision))
    config["spatial_review"] = workflow.record(review)
    target = tmp_path / "mapped.nii.gz"
    receipt = mapped_label(config, tmp_path, target)
    assert receipt["voxels"] == 1 and receipt["interpolations"] == 1
    assert all(row["max_mm"] == 0 for row in receipt["landmark_errors"].values())
    assert sitk.ReadImage(str(target))[2, 3, 1] == 1
    invalid = tmp_path / "invalid.nii.gz"
    wrong_grid = tmp_path / "wrong-grid.nii.gz"
    image.SetOrigin((2., -2., 0.))
    sitk.WriteImage(image, str(wrong_grid))
    valid_label = config["label"]
    config["label"] = workflow.record(wrong_grid)
    with pytest.raises(ValueError, match="same physical grid"):
        mapped_label(config, tmp_path, invalid, technical_replay=True)
    config["label"] = valid_label
    landmark_data = json.loads(landmarks.read_text())
    landmark_data["segments"]["dwi_to_t1_ras"][0]["to_lps_mm"][0] += 1.
    landmarks.write_text(json.dumps(landmark_data))
    config["landmarks"] = workflow.record(landmarks)
    decision["input_sha256"]["landmarks"] = config["landmarks"]["sha256"]
    decision["evidence"]["landmarks"] = config["landmarks"]
    review.write_text(json.dumps(decision))
    config["spatial_review"] = workflow.record(review)
    with pytest.raises(ValueError, match="Landmark error exceeds"):
        mapped_label(config, tmp_path, invalid)
    config["anatomy"]["mni6_to_mni2009c_image_h5"] = workflow.record(paths[0])
    with pytest.raises(ValueError, match="locked.*point map"):
        mapped_label(config, tmp_path, invalid)
    assert not invalid.exists()


def test_tracking_plan_and_replay_input_integrity(example, tmp_path):
    source = example["config"]
    config = workflow.absolute_records(workflow.read(source), source.parent)
    config["runs"] = [config["runs"][0]]
    run = config["runs"][0]
    run.update(tracker_input=run["parent"], rng_seed=1234, shell_b=1000)
    path = tmp_path / "local.json"
    path.write_text(json.dumps(config), encoding="utf-8")
    plan = workflow.tracking_plan(path, tmp_path / "plan")["jobs"][0]
    argv = plan["argv"]
    for key, value in (("-seeds", "32"), ("-select", "0"), ("-seed_cutoff", "0.10"), ("-cutoff", "0.10")):
        assert argv[argv.index(key) + 1] == value
    assert plan["rng_seed"] == 1234
    assert plan["shell_b"] == 1000
    run["seed_log"]["sha256"] = "0" * 64
    path.write_text(json.dumps(config), encoding="utf-8")
    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        workflow.screen(path, tmp_path / "invalid")
    assert not (tmp_path / "invalid" / "results.json").exists()


def test_native_preflight_header_and_shell(monkeypatch):
    ref = nib.Nifti1Image(np.zeros((5, 6, 7), dtype=np.uint8), np.diag([1., 2., 3., 1.]))
    report = "6 5 7 7\n2 1 3 1\n0 1 0 0\n1 0 0 0\n0 0 1 0\n0 0 0 1\n"
    monkeypatch.setattr(workflow, "mrinfo", lambda *args: report)
    assert workflow.check_tracker_grid(Path("model.mif"), ref) == 7

    dirs = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 0], [1, 0, 1], [0, 1, 1]], dtype=float)
    dirs /= np.linalg.norm(dirs, axis=1, keepdims=True)
    gradients = "0 0 0 0\n" + "\n".join(" ".join(map(str, [*v, 1000])) for v in dirs)
    monkeypatch.setattr(workflow, "mrinfo", lambda *args: gradients)
    workflow.check_gradients(Path("model.mif"), 7, 1000, 6)
    with pytest.raises(ValueError, match="declared shell"):
        workflow.check_gradients(Path("model.mif"), 7, 2000, 6)
