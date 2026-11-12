import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest
import SimpleITK as sitk

from tractography import workflow
from tractography.demo import demo
from tractography.geometry import MNI6, mapped_label
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


def test_native_atlas_mapping_uses_one_composed_pull(tmp_path):
    image = sitk.Image([10, 10, 10], sitk.sitkUInt8)
    ref = tmp_path / "b0.nii.gz"
    sitk.WriteImage(image, str(ref))
    image[4, 4, 1] = 1
    label = tmp_path / "atlas.nii.gz"
    sitk.WriteImage(image, str(label))
    middle = sitk.AffineTransform([0., -1., 0., 1., 0., 0., 0., 0., 1.], [8., 0., 0.])
    bridge = sitk.TranslationTransform(3, [0., 3., 0.])
    paths = [tmp_path / "t1.h5", tmp_path / "bridge.h5"]
    for transform, path in zip((middle, bridge), paths, strict=True):
        sitk.WriteTransform(transform, str(path))
    ras = np.eye(4)
    ras[0, 3] = 2
    config = {"label": workflow.record(label), "reference": workflow.record(ref), "source_space": MNI6,
              "anatomy": {"dwi_to_t1_ras": ras.tolist(),
                          "mni2009c_to_t1_image_h5": workflow.record(paths[0]),
                          "mni6_to_mni2009c_image_h5": workflow.record(paths[1])}}
    target = tmp_path / "mapped.nii.gz"
    receipt = mapped_label(config, tmp_path, target)
    assert receipt["voxels"] == 1 and receipt["interpolations"] == 1
    assert sitk.ReadImage(str(target))[3, 4, 1] == 1


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
