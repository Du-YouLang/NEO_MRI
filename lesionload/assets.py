"""Load user-supplied atlases from a public source and SHA-256 lock."""

import hashlib
import json
from pathlib import Path
from zipfile import ZipFile

import nibabel as nib
import numpy as np

from anatomy_review import file_record
from . import boyne
from .spatial import MNI6, MNI2009C, binary, checked, image_grid, target_grid


LOCK_PATH = Path(__file__).with_name("reference_assets.lock.json")
HCPA_AFFINE = np.eye(4)
HCPA_AFFINE[:3, 3] = [-96, -132, -78]
BRIDGE_SOURCE = "tpl-MNI152NLin2009cAsym_from-MNI152NLin6Asym_mode-image_xfm.h5"
GEN6_BACKGROUND = "tpl-MNI152NLin6Asym_res-01_desc-brain_T1w.nii.gz"
METHODS = {"SMATT": MNI6, "HCPA": MNI2009C, "Boyne2025": MNI6}


def public_lock():
    return json.loads(LOCK_PATH.read_text(encoding="utf-8"))


def locked_file(directory, name):
    record = next(row for row in public_lock()["atlas_files"] if row["file"] == name)
    if isinstance(directory, dict):
        path = checked({"path": directory[name], "sha256": record["sha256"]}, Path.cwd())
    else:
        path = checked({"path": name, "sha256": record["sha256"]}, directory)
    if path.stat().st_size != record["bytes"]:
        raise ValueError(f"Reference size mismatch: {path}")
    return path


def load_smatt(directory, target):
    bank = {}
    with ZipFile(locked_file(directory, "SMATT.zip")) as archive:
        for key, entry in public_lock()["smatt_binary_members"].items():
            blob = archive.read(entry["member"])
            if hashlib.sha256(blob).hexdigest() != entry["member_sha256"]:
                raise ValueError(f"SMATT archive member hash mismatch: {key}")
            image = nib.Nifti1Image.from_bytes(blob)
            mask = binary(image)
            if (mask.shape != target.shape or not np.allclose(image.affine, target.affine, atol=1e-6, rtol=0)
                    or int(mask.sum()) != entry["voxel_count"]):
                raise ValueError(f"SMATT archive member grid or count mismatch: {key}")
            bank[key.removesuffix(".nii.gz")] = mask
    bank["Union"] = np.logical_or.reduce(list(bank.values()))
    return bank


def load_boyne(directory, target):
    bank = {}
    for side, name in boyne.SIDES.items():
        path = locked_file(directory, f"CRST_{side}_subject_count_1mm.nii.gz")
        image = nib.load(path)
        values = np.asanyarray(image.dataobj)
        if (values.shape != (*target.shape, len(boyne.ORIGINS))
                or not np.allclose(image.affine, target.affine, atol=1e-6, rtol=0)
                or values.dtype.kind not in "iu" or values.min() < 0 or values.max() > 50):
            raise ValueError(f"Invalid Boyne count map: {path}")
        if not np.array_equal(values[..., 0] > 0, np.any(values[..., 1:] > 0, axis=-1)):
            raise ValueError("Boyne whole support differs from the union of its origin channels")
        bank.update({f"{name}_{origin}": values[..., k] > 0 for k, origin in enumerate(boyne.ORIGINS)})
    bank["Union"] = bank["Left_Whole"] | bank["Right_Whole"]
    return bank


def real_atlas(directory, method):
    if method == "HCPA":
        target_path = locked_file(directory, "HCPA_CST_Left_MNI.nii")
        grid = target_grid(target_path)
        right_path = locked_file(directory, "HCPA_CST_Right_MNI.nii")
        right = image_grid(right_path, grid)
        if not np.array_equal(grid.affine, HCPA_AFFINE):
            raise ValueError("Unexpected HCPA 2009c geometry")
        bank = {"Left_CST": binary(grid), "Right_CST": binary(right)}
        return target_path, grid, bank, [file_record(path) for path in (target_path, right_path)]

    target_path = locked_file(directory, GEN6_BACKGROUND)
    grid = target_grid(target_path)
    # TemplateFlow background is distributed in RAS storage; the archived
    # tract banks are LAS. Their physical lattice is identical after x flip.
    target = grid.as_reoriented(
        nib.orientations.ornt_transform(nib.orientations.io_orientation(grid.affine),
                                        nib.orientations.axcodes2ornt(("L", "A", "S"))))
    if target.shape != boyne.SHAPE or not np.array_equal(target.affine, boyne.AFFINE):
        raise ValueError("Unexpected Gen6 atlas geometry")
    bank = load_smatt(directory, target) if method == "SMATT" else load_boyne(directory, target)
    names = ["SMATT.zip"] if method == "SMATT" else [
        f"CRST_{side}_subject_count_1mm.nii.gz" for side in boyne.SIDES]
    source_records = [file_record(target_path)]
    source_records.extend(file_record(locked_file(directory, name)) for name in names)
    return target_path, target, bank, source_records


def build_boyne_assets(directory, output):
    """Explicit offline rebuild from the two original TRKs and 18 ROI files."""
    paths = {f"{origin}_{side}.inDil2.nii.gz": locked_file(directory, f"{origin}_{side}.inDil2.nii.gz")
             for origin in boyne.ORIGINS[1:-1] for side in boyne.SIDES}
    lookup, sources = boyne.roi_lookup(paths)
    author_script = locked_file(directory, "AnalysisCode_2025-07-12.sh")
    output.mkdir(parents=True)
    rows = {}
    for side in boyne.SIDES:
        path = locked_file(directory, f"MNI_dsnFSL_IDall_CRST_{side}_manual.prune2.trk.gz")
        counts, evidence = boyne.voxelize(side, path, lookup)
        image = nib.Nifti1Image(np.moveaxis(counts, 0, -1), boyne.AFFINE)
        image.header.set_xyzt_units("mm")
        dest = output / f"CRST_{side}_subject_count_1mm.nii.gz"
        nib.save(image, dest)
        with dest.open("rb") as stream:
            digest = hashlib.file_digest(stream, "sha256").hexdigest()
        rows[side] = {"sha256": digest, "roi_inputs": sources[side],
                      "tractogram": file_record(path), "evidence": evidence}
    receipt = {"priority": list(boyne.ORIGINS), "method": "Boyne 2025 C DSN",
               "dlpfc": "Deferred; current branches replay the available masks",
               "author_script": file_record(author_script), "sides": rows}
    (output / "atlas_provenance.json").write_text(json.dumps(receipt, indent=2), encoding="utf-8")
    return receipt
