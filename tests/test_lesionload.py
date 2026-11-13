import json
from pathlib import Path
import sys

import nibabel as nib
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))

from lesionload import boyne
from lesionload.demo import demo
from lesionload.metrics import measure
from lesionload.workflow import run
from spatial_contract import digest


def test_demo_count_contract_and_output(tmp_path):
    info = demo(tmp_path / "demo")
    receipt = json.loads(Path(info["receipt"]).read_text(encoding="utf-8"))
    assert receipt["synthetic"] and set(receipt["methods"]) == {"SMATT", "HCPA", "Boyne2025"}
    method = receipt["methods"]["HCPA"]
    left, right, empty, union = method["summary"]
    assert [left[key] for key in ("atlas_voxels", "covered_voxels", "hits_voxels", "unknown_voxels")] == [15, 13, 4, 2]
    assert left["fraction_whole"] is None
    assert left["fraction_covered"] == pytest.approx(4 / 13)
    assert left["fraction_bounds"] == pytest.approx([4 / 15, 6 / 15])
    assert left["peak_fraction_bounds"] == pytest.approx([.3, .5])
    assert left["max_wll_bounds_voxel"] == [3., 5.]
    assert right["fraction_whole"] == 0
    assert empty["fraction_covered"] is None and empty["peak_rows"] == []
    assert union["atlas_voxels"] == 24
    outputs = method["outputs"]
    lesion = np.asanyarray(nib.load(outputs["lesion"]["path"]).dataobj)
    fov = np.asanyarray(nib.load(outputs["source_fov"]["path"]).dataobj)
    assert np.count_nonzero(lesion) == 4 and not np.any((lesion != 0) & (fov == 0))
    assert np.all(fov[4]) and not np.any(fov[5])
    assert "Left_CST,inside,4,4,3,0.625,0.625" in Path(outputs["fa_whole"]["path"]).read_text()
    html = Path(outputs["report"]["path"]).read_text(encoding="utf-8")
    payload = html.split('id="smatt-data">', 1)[1].split("</script>", 1)[0]
    region = json.loads(payload)["cases"]["current"]["regions"]["Left_CST"]
    assert region["slices"][1] == [1, 5, 5, .625, 4]
    assert region["slices"][2] == [3, 10, 8, .625, 8]


def test_peak_ties_keep_each_slice_denominator_and_cover_uses_union():
    tract = np.ones((4, 1, 3), dtype=bool)
    tract[2:, :, 1] = False
    lesion = np.zeros_like(tract)
    lesion[:2, :, 0] = True
    lesion[:1, :, 1] = True
    fov = np.ones_like(tract)
    fov[1:, :, 2] = False
    result, rows = measure(tract, lesion, fov, [2, 3, 4], hcpa=True)
    assert [(row["z_ras_mm"], row["atlas_voxels"]) for row in result["peak_rows"]] == [(2, 4), (3, 2)]
    assert result["peak_fraction_bounds"] == [.5, .75]
    assert result["max_wll_bounds_voxel"] == [2, 3]
    assert sum(row["covered_voxels"] for row in rows) == result["covered_voxels"]


def test_config_rejects_wrong_lesion_hash_before_creating_output(tmp_path):
    info = demo(tmp_path / "demo")
    config_path = Path(info["config"])
    config = json.loads(config_path.read_text(encoding="utf-8"))
    config["t2"]["lesion"]["sha256"] = "0" * 64
    config_path.write_text(json.dumps(config), encoding="utf-8")
    output = tmp_path / "invalid"
    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        run(config_path, output)
    assert not output.exists()


@pytest.mark.parametrize("change,message", [
    ("visit", "different visit"), ("decision", "no supported approval"),
    ("direction", "point directions"), ("lesion", "actual lesion"),
    ("transform", "actual t1_to_t2_pull"), ("fa", "actual fa"),
])
def test_spatial_review_rejects_mismatched_approved_inputs(tmp_path, change, message):
    info = demo(tmp_path / "demo")
    config_path = Path(info["config"])
    config = json.loads(config_path.read_text(encoding="utf-8"))
    key = "fa_spatial_review" if change == "fa" else "spatial_review"
    review_path = Path(config[key]["path"])
    review = json.loads(review_path.read_text(encoding="utf-8"))
    if change == "visit":
        review["visit"]["session"] = "02"
    elif change == "decision":
        review["decision"] = "PENDING"
    elif change == "direction":
        review["point_directions"]["mni_to_t1_pull"][1:] = ["T1w", "MNI152NLin2009cAsym"]
    elif change == "transform":
        review["input_sha256"]["t1_to_t2_pull"] = "0" * 64
    else:
        review["input_sha256"][change] = "0" * 64
    review_path.write_text(json.dumps(review), encoding="utf-8")
    config[key]["sha256"] = digest(review_path)
    config_path.write_text(json.dumps(config), encoding="utf-8")
    output = tmp_path / "invalid"
    with pytest.raises(ValueError, match=message):
        run(config_path, output)
    assert not output.exists()


def test_spatial_review_follows_same_bytes_after_relocation(tmp_path):
    info = demo(tmp_path / "demo")
    config_path = Path(info["config"])
    config = json.loads(config_path.read_text(encoding="utf-8"))
    old = Path(config["t2"]["lesion"]["path"])
    moved = tmp_path / "relocated_lesion.nii.gz"
    old.rename(moved)
    config["t2"]["lesion"]["path"] = str(moved)
    config_path.write_text(json.dumps(config), encoding="utf-8")
    output = tmp_path / "new"
    receipt = run(config_path, output)
    assert receipt["t2"]["input_review"]["decision"] == "APPROVED"
    assert receipt["t2"]["input_review"]["input_sha256"]["lesion"] == digest(moved)


def test_boyne_crosses_cells_between_vertices_and_deduplicates_visits(monkeypatch):
    monkeypatch.setattr(boyne, "STRIDES", np.array([25, 5, 1]))
    line = np.array([[0.0, 0., 0.], [4., 0., 0.], [0., 0., 0.]])
    visited = np.zeros(125, dtype=np.uint16)
    boyne.set_visits(visited, line, boyne.index_cells(line), np.array([2]), np.array([2], dtype=np.uint16))
    np.testing.assert_array_equal(np.flatnonzero(visited), np.arange(5) * 25)
    assert np.all(visited[visited != 0] == 2)
    first, last = np.array([[0., 0., 0.]]), np.array([[3., 3., 3.]])
    points = [tuple(cells[0]) for _, cells in boyne.traversed_cells(
        first, last, boyne.index_cells(first), boyne.index_cells(last))]
    assert points == [(1, 1, 1), (2, 2, 2)]
