"""Propagate an original T2 lesion and its source FOV through a frozen registration."""

import argparse
import hashlib
import json
from itertools import product
from pathlib import Path

import nibabel as nib
import numpy as np

from anatomy_review import file_record, read, save_new, sha
from multimodal_review import export
from structural_multimodal import PARAMETERS, compose_pull, resample, sitk_module
from t2_stable import ENGINE_HASHES, PROFILE, verify_engine


SCHEMA = "neo-t2-lesion-propagation-1"
BRANCHES = ("raw_independent_HDBET", "n4_independent_HDBET", "raw_T1_only", "n4_T1_only")
REFERENCE_ROLES = ("source", "t1", "mni", "mni_to_t1_pull", "t1_mask", "anatomy_review",
                   "normalization_provenance", "mni_review")


def bound_file(path, expected):
    record = file_record(path)
    if record["sha256"] != expected:
        raise ValueError(f"Frozen input or output hash mismatch: {path}")
    return record


def plan(config):
    """Bind the selected completed structural run and its input grids by content."""
    if config.get("schema") != SCHEMA:
        raise ValueError(f"Expected schema {SCHEMA}")
    verify_engine()
    sitk = sitk_module()
    if sitk.Version_VersionString() != "2.5.6":
        raise ValueError("T2 lesion propagation v1 requires SimpleITK 2.5.6")
    branch = config["branch"]
    if branch not in BRANCHES:
        raise ValueError("Select an explicit T2 v1 branch; HD-BET is the primary brain domain")
    parent = Path(config["structural_run"]).resolve(strict=True)
    output = Path(config["output"]).resolve()
    if output.is_relative_to(parent) or output.exists():
        raise FileExistsError("Choose a new output directory outside the frozen structural run")

    receipt_path = parent / "receipt.json"
    receipt_record = file_record(receipt_path)
    receipt = read(receipt_path)
    if receipt["status"] != "COMPUTED_REVIEW_PENDING" or not any(
        row["branch"] == branch for row in receipt["branches"]
    ):
        raise ValueError("Selected branch is absent from a completed T2 structural receipt")
    saved = {record["path"]: record["sha256"] for record in receipt["outputs"]}
    request_record = bound_file(parent / "request.json", saved["request.json"])
    request = read(request_record["path"])
    identity = request["identity"]
    fingerprint = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
    if (fingerprint != request["fingerprint"] or fingerprint != receipt["fingerprint"]
            or identity["schema"] != "neo-multimodal-structural-1" or identity["modality"] != "T2w"
            or identity["code_sha256"] != ENGINE_HASHES["structural_multimodal.py"]
            or identity["parameters"] != dict(PARAMETERS, N4_zero_policy=PROFILE["N4_zero_policy"])
            or identity["runtime"]["SimpleITK"] != "2.5.6"):
        raise ValueError("Structural request or receipt does not match the frozen T2 v1 profile")

    roles = (*REFERENCE_ROLES, "source_mask") if branch.endswith("_HDBET") else REFERENCE_ROLES
    overrides = config.get("input_paths", {})
    if set(overrides) - set(roles):
        raise ValueError("An input path override must name a bound structural input")
    files = identity["files"]
    inputs = {key: bound_file(overrides.get(key, files[key]["path"]), files[key]["sha256"])
              for key in roles}
    for key, suffix in (("t1_to_t2_pull", "T1_to_source_pull.h5"),
                        ("t2_in_t1", "T1.nii.gz"), ("t2_in_mni", "MNI_once.nii.gz")):
        name = f"{branch}_{suffix}"
        inputs[key] = bound_file(parent / name, saved[name])
    inputs["lesion"] = file_record(config["lesion"])
    return {
        "schema": SCHEMA, "status": "PLANNED", "structural_profile": PROFILE["id"],
        "branch": branch, "structural_request": request_record, "structural_receipt": receipt_record,
        "structural_fingerprint": fingerprint, "inputs": inputs, "output": str(output),
        "pull_lps": {"T1": ["T1", "original T2"], "MNI": ["MNI", "T1", "original T2"]},
        "interpolation": "nearest_neighbor", "resampling_steps_per_output": 1,
        "clinical_validation": False,
    }


def check_nifti_grid(path):
    image = nib.load(path)
    if (len(image.shape) != 3 or not np.isfinite(image.affine).all()
            or abs(np.linalg.det(image.affine[:3, :3])) < 1e-8):
        raise ValueError(f"Expected a finite, nonsingular 3D NIfTI grid: {path}")
    qform, qcode = image.get_qform(coded=True)
    sform, scode = image.get_sform(coded=True)
    if qcode and scode:
        corners = np.asarray([(*index, 1) for index in product(*[(0, n - 1) for n in image.shape])])
        if np.max(np.linalg.norm((corners @ (qform - sform).T)[:, :3], axis=1)) > 0.1:
            raise ValueError(f"Conflicting qform/sform, review geometry before propagation: {path}")


def same_geometry(left, right):
    return left.IsSameImageGeometryAs(right, 1e-4, 1e-6)


def write_verified(sitk, image, target):
    sitk.WriteImage(image, str(target), True)
    readback = sitk.ReadImage(str(target))
    if (not same_geometry(image, readback)
            or not np.array_equal(sitk.GetArrayViewFromImage(image), sitk.GetArrayViewFromImage(readback))):
        raise ValueError(f"Written image failed geometry or voxel readback: {target}")


def run(config):
    request = plan(config)
    paths = {role: record["path"] for role, record in request["inputs"].items()}
    sitk = sitk_module()
    for key in ("source", "t1", "mni", "lesion"):
        check_nifti_grid(paths[key])
    source = sitk.ReadImage(paths["source"])
    lesion = sitk.ReadImage(paths["lesion"])
    if source.GetNumberOfComponentsPerPixel() != 1 or lesion.GetNumberOfComponentsPerPixel() != 1:
        raise ValueError("Original T2 and lesion must be scalar images")
    if not same_geometry(source, lesion):
        raise ValueError("Original T2 and lesion must share the same physical grid")
    lesion_values = sitk.GetArrayViewFromImage(lesion)
    if not np.isin(lesion_values, (0, 1)).all():
        raise ValueError("Original T2 lesion must be binary (0 or 1)")
    native_voxels = int(np.count_nonzero(lesion_values))
    fov = sitk.Image(source.GetSize(), sitk.sitkUInt8) + 1
    fov.CopyInformation(source)

    t1_pull = sitk.ReadTransform(paths["t1_to_t2_pull"])
    mni_pull = compose_pull(t1_pull, sitk.ReadTransform(paths["mni_to_t1_pull"]))
    out = Path(request["output"])
    out.mkdir(parents=True)
    save_new(out / "request.json", request)
    try:
        outputs = {}
        metrics = {}
        spaces = [{"name": "原生 T2", "background": paths["source"],
                   "frames": [{"name": "原始病灶", "image": paths["source"], "lesion": paths["lesion"]}]}]
        for space, pull in (("T1", t1_pull), ("MNI", mni_pull)):
            reference_path = paths[space.lower()]
            reference = sitk.ReadImage(reference_path)
            mapped_lesion = resample(lesion, reference, pull, label=True)
            mapped_fov = resample(fov, reference, pull, label=True)
            mask = sitk.GetArrayViewFromImage(mapped_lesion)
            support = sitk.GetArrayViewFromImage(mapped_fov)
            if np.any((mask != 0) & (support == 0)):
                raise ValueError(f"Lesion outside mapped source FOV in {space}")
            count = int(np.count_nonzero(mask))
            metrics[space] = {
                "lesion_voxels": count, "source_fov_voxels": int(np.count_nonzero(support)),
                "target_voxels": int(support.size), "source_fov_fraction": float(np.mean(support)),
                "lesion_empty_after_mapping": bool(native_voxels and not count),
            }
            for key, image in (("lesion", mapped_lesion), ("source_fov", mapped_fov)):
                target = out / f"{space}_{key}.nii.gz"
                write_verified(sitk, image, target)
                outputs[f"{space}_{key}"] = file_record(target)
            spaces.append({"name": space, "background": reference_path,
                           "frames": [{"name": request["branch"], "image": paths[f"t2_in_{space.lower()}"],
                                       "lesion": outputs[f"{space}_lesion"]["path"]}]})
        export(spaces, out / "overlay.html", "T2 病灶传播审核")
        outputs["overlay"] = file_record(out / "overlay.html")
        for record in (*request["inputs"].values(), request["structural_request"], request["structural_receipt"]):
            bound_file(record["path"], record["sha256"])
        receipt = {"schema": SCHEMA, "status": "COMPUTED_REVIEW_PENDING", "branch": request["branch"],
                   "request": file_record(out / "request.json"), "SimpleITK": sitk.Version_VersionString(),
                   "native_lesion_voxels": native_voxels, "metrics": metrics, "outputs": outputs,
                   "clinical_validation": False}
        save_new(out / "receipt.json", receipt)
        return receipt
    except Exception as error:
        save_new(out / "failure.json", {"status": "FAILED", "error": str(error)})
        raise


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("plan", "run"))
    parser.add_argument("--config", required=True)
    args = parser.parse_args(argv)
    config = read(args.config)
    result = (plan if args.action == "plan" else run)(config)
    print(json.dumps(result, ensure_ascii=True, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
