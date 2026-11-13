"""Use the reviewed T2 chain and strict native FA interpolation support."""

from pathlib import Path

import nibabel as nib
import numpy as np
from scipy.ndimage import map_coordinates

from anatomy_review import file_record
from structural_multimodal import compose_pull, resample, sitk_module
from t2_lesion_propagation import check_nifti_grid, plan, same_geometry


MNI2009C = "MNI152NLin2009cAsym"
MNI6 = "MNI152NLin6Asym"
RAS_LPS = np.array([-1., -1., 1.])


def checked(item, root):
    path = (root / item["path"]).resolve(strict=True)
    if file_record(path)["sha256"] != item["sha256"]:
        raise ValueError(f"File hash mismatch: {path}")
    return path


def image_grid(path, reference=None):
    check_nifti_grid(path)
    image = nib.load(path)
    if reference is not None and (image.shape != reference.shape or not np.allclose(
        image.affine, reference.affine, atol=1e-4, rtol=0
    )):
        raise ValueError(f"Physical grid mismatch: {path}")
    return image


def target_grid(path):
    image = image_grid(path)
    linear = image.affine[:3, :3]
    if not np.allclose(np.abs(linear), np.eye(3), atol=1e-5, rtol=0):
        raise ValueError(f"An axis-aligned 1 mm atlas grid is required: {path}")
    return image


def binary(image):
    values = np.asanyarray(image.dataobj)
    if not np.isin(values, (0, 1)).all():
        raise ValueError("Expected a binary mask with only 0 and 1")
    return values.astype(bool)


def transform_points(transform, points):
    return np.array([transform.TransformPoint(tuple(point)) for point in points]).reshape(-1, 3)


def t2_chain(config, base, stage):
    """Load the original mask and explicit T1/2009c point pulls."""
    sitk = sitk_module()
    if sitk.Version_VersionString() != "2.5.6":
        raise ValueError("This workflow requires SimpleITK 2.5.6")
    if "structural_run" in config:
        selected = dict(config)
        for key in ("structural_run", "lesion"):
            selected[key] = str((base / selected[key]).resolve(strict=True))
        selected["input_paths"] = {key: str((base / path).resolve(strict=True))
                                   for key, path in config.get("input_paths", {}).items()}
        selected["output"] = str(stage / "t2_plan_only")
        provenance = plan(selected)
        records = provenance["inputs"]
        paths = {key: Path(records[key]["path"]) for key in (
            "source", "lesion", "t1", "mni", "t1_to_t2_pull", "mni_to_t1_pull"
        )}
        if file_record(paths["lesion"])["sha256"] != config["lesion_sha256"]:
            raise ValueError("Original lesion hash mismatch")
    else:
        paths = {key: checked(config[key], base) for key in (
            "source", "lesion", "t1", "mni", "t1_to_t2_pull", "mni_to_t1_pull", "review"
        )}
        provenance = {"profile": "archived_explicit_chain",
                      "inputs": {key: file_record(path) for key, path in paths.items()}}
    source = image_grid(paths["source"])
    lesion_grid = image_grid(paths["lesion"], source)
    binary(lesion_grid)
    label = sitk.ReadImage(str(paths["lesion"]))
    source_image = sitk.ReadImage(str(paths["source"]))
    if not same_geometry(source_image, label):
        raise ValueError("T2 and lesion physical grids differ")
    for nifti, itk in ((source, source_image), (lesion_grid, label)):
        indexes = np.array([[0, 0, 0], *np.diag(np.array(nifti.shape) - 1)])
        expected = nib.affines.apply_affine(nifti.affine, indexes) * RAS_LPS
        observed = [itk.TransformIndexToPhysicalPoint(tuple(map(int, index))) for index in indexes]
        if not np.allclose(expected, observed, atol=.015, rtol=0):
            raise ValueError("ITK and NIfTI disagree on source physical geometry")
    fov = sitk.Image(source_image.GetSize(), sitk.sitkUInt8) + 1
    fov.CopyInformation(source_image)
    return paths, provenance, label, fov, sitk.ReadTransform(str(paths["t1_to_t2_pull"])), (
        sitk.ReadTransform(str(paths["mni_to_t1_pull"])))


def map_pair(label, fov, target, pull):
    sitk = sitk_module()
    sampler = sitk.ResampleImageFilter()
    sampler.SetReferenceImage(target)
    sampler.SetTransform(pull)
    sampler.SetInterpolator(sitk.sitkNearestNeighbor)
    sampler.SetOutputPixelType(sitk.sitkUInt8)
    # Match the audited sampling work units at nearest-neighbour boundaries.
    sampler.SetNumberOfWorkUnits(1)
    return tuple(sitk.GetArrayFromImage(sampler.Execute(image)).T.astype(bool)
                 for image in (label, fov))


def itk_grid(reference):
    """Build a blank SimpleITK grid from an axis-aligned RAS atlas header."""
    sitk = sitk_module()
    target = sitk.Image(reference.shape, sitk.sitkUInt8)
    target.SetOrigin(reference.affine[:3, 3] * RAS_LPS)
    target.SetDirection((RAS_LPS[:, None] * reference.affine[:3, :3]).ravel())
    return target


def fsl_coordinates(image):
    zooms = image.header.get_zooms()[:3]
    matrix = np.diag([*zooms, 1.])
    if np.linalg.det(image.affine[:3, :3]) > 0:
        matrix[0, 0] = -zooms[0]
        matrix[0, 3] = (image.shape[0] - 1) * zooms[0]
    return matrix


def native_fa(config, base, t1_path):
    paths = {key: checked(config[key], base) for key in (
        "b0", "fa", "md", "fit_mask", "flirt_dwi_to_t1", "flirt_inverse", "review", "reviewed_mask_in_t1"
    )}
    if file_record(t1_path)["sha256"] != config["t1_sha256"]:
        raise ValueError("FA registration belongs to another T1")
    b0, t1 = image_grid(paths["b0"]), image_grid(t1_path)
    fa, md = (np.asanyarray(image_grid(paths[key], b0).dataobj) for key in ("fa", "md"))
    fit = binary(image_grid(paths["fit_mask"], b0))
    valid = fit & np.isfinite(fa) & (fa >= 0) & (fa <= 1) & np.isfinite(md) & (md > 0)
    if int(valid.sum()) != config["valid_native_voxels"]:
        raise ValueError("Native valid FA support differs from the model receipt")
    matrix, reverse = (np.loadtxt(paths[key]) for key in ("flirt_dwi_to_t1", "flirt_inverse"))
    world = (t1.affine @ np.linalg.inv(fsl_coordinates(t1)) @ matrix
             @ fsl_coordinates(b0) @ np.linalg.inv(b0.affine))
    if (not np.allclose(matrix @ reverse, np.eye(4), atol=1e-6, rtol=0)
            or not np.allclose(world, config["world_forward_ras"], atol=1e-6, rtol=0)):
        raise ValueError("FLIRT transform differs from the reviewed physical registration")
    flip = np.diag([*RAS_LPS, 1.])
    lps = flip @ np.linalg.inv(world) @ flip
    sitk = sitk_module()
    pull = sitk.AffineTransform(lps[:3, :3].ravel(), lps[:3, 3])
    reference = sitk.ReadImage(str(t1_path))
    fit_image = sitk.ReadImage(str(paths["fit_mask"]))
    replay = sitk.GetArrayFromImage(resample(fit_image, reference, pull, label=True)).T > 0
    archived = binary(image_grid(paths["reviewed_mask_in_t1"], t1))
    denominator = archived.sum() + replay.sum()
    dice = float(2 * (archived & replay).sum() / denominator) if denominator else 0.
    if dice < .999:
        raise ValueError("FLIRT replay disagrees with the reviewed T1 mask")
    provenance = {"files": {key: file_record(path) for key, path in paths.items()},
                  "replay_t1_mask_dice": dice, "native_valid_voxels": int(valid.sum())}
    if config.get("model_receipt"):
        provenance["files"]["model_receipt"] = file_record(checked(config["model_receipt"], base))
    return b0.affine, fa, valid, pull, provenance


def sample_fa(model, template_to_t1, reference, union):
    affine, fa, valid, t1_to_dwi, _ = model
    points = np.argwhere(union)
    lps = nib.affines.apply_affine(reference.affine, points) * RAS_LPS
    dwi = transform_points(compose_pull(t1_to_dwi, template_to_t1), lps)
    coords = nib.affines.apply_affine(np.linalg.inv(affine), dwi * RAS_LPS)
    inside = ((coords >= 0) & (coords <= np.array(fa.shape) - 1)).all(axis=1)
    weights = map_coordinates(valid.astype(float), coords.T, order=1, prefilter=False)
    values = map_coordinates(np.where(valid, fa, 0).astype(float), coords.T, order=1, prefilter=False)
    sampled = np.full(reference.shape, np.nan)
    acquired = np.zeros(reference.shape, dtype=bool)
    keep = inside & (weights >= 1 - 1e-7)
    sampled[tuple(points[keep].T)] = values[keep]
    acquired[tuple(points[inside].T)] = True
    return sampled, acquired
