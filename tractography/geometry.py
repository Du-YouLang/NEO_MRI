"""Physical-space atlas mapping and validated native DWI grids."""
from pathlib import Path
import hashlib

import nibabel as nib
import numpy as np

MNI2009C = "MNI152NLin2009cAsym"
MNI6 = "MNI152NLin6Asym"


def digest(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def checked(record, base):
    path = (base / record["path"]).resolve(strict=True)
    if not path.is_file() or digest(path) != record["sha256"].lower():
        raise ValueError(f"Input SHA-256 mismatch: {path}")
    return path


def load_reference(path):
    ref = nib.as_closest_canonical(nib.load(path))
    affine = ref.affine
    if len(ref.shape) != 3 or not np.isfinite(affine).all():
        raise ValueError("A finite three-dimensional corrected DWI reference is required")
    axes = affine[:3, :3] / nib.affines.voxel_sizes(affine)
    if not np.allclose(axes.T @ axes, np.eye(3), atol=1e-5, rtol=0):
        raise ValueError("The DWI grid must be orthogonal; oblique rotation is supported")
    return ref


def load_mask(path, ref):
    image = nib.as_closest_canonical(nib.load(path))
    if image.shape != ref.shape or not np.allclose(image.affine, ref.affine, atol=1e-4, rtol=0):
        raise ValueError(f"Mask does not share the DWI physical grid: {path}")
    values = np.asanyarray(image.dataobj)
    if not np.isin(values, (0, 1)).all():
        raise ValueError(f"Expected a binary hard label: {path}")
    return values.astype(bool)


def atlas_pull(anatomy, base, source_space):
    """Map an output DWI LPS point to an input atlas LPS point."""
    import SimpleITK as sitk
    from structural_multimodal import compose_pull

    ras = np.asarray(anatomy["dwi_to_t1_ras"], dtype=float)
    if (ras.shape != (4, 4) or not np.isfinite(ras).all()
            or not np.allclose(ras[3], [0, 0, 0, 1], atol=1e-7, rtol=0)
            or not np.allclose(ras[:3, :3].T @ ras[:3, :3], np.eye(3), atol=1e-4, rtol=0)
            or not np.isclose(np.linalg.det(ras[:3, :3]), 1, atol=1e-4, rtol=0)):
        raise ValueError("DWI-to-T1 must be a rigid physical RAS matrix")
    flip = np.diag([-1, -1, 1, 1])
    lps = flip @ ras @ flip
    pose = sitk.AffineTransform(lps[:3, :3].ravel().tolist(), lps[:3, 3].tolist())
    inverse_image = sitk.ReadTransform(str(checked(anatomy["mni2009c_to_t1_image_h5"], base)))
    dwi_to_mni = compose_pull(inverse_image, pose)
    if source_space == MNI2009C:
        return dwi_to_mni
    if source_space == MNI6:
        bridge = sitk.ReadTransform(str(checked(anatomy["mni6_to_mni2009c_image_h5"], base)))
        return compose_pull(bridge, dwi_to_mni)
    raise ValueError(f"Unsupported atlas space: {source_space}")


def mapped_label(config, base, output):
    """Apply a composed pull to a binary atlas label in one interpolation."""
    import SimpleITK as sitk
    from structural_multimodal import resample

    source = sitk.ReadImage(str(checked(config["label"], base)))
    reference = sitk.ReadImage(str(checked(config["reference"], base)))
    if source.GetDimension() != 3 or reference.GetDimension() != 3:
        raise ValueError("Atlas and reference must be three-dimensional")
    if not np.isin(sitk.GetArrayViewFromImage(source), (0, 1)).all():
        raise ValueError("Threshold the source atlas to an explicit binary label first")
    if output.exists():
        raise FileExistsError(output)
    pull = atlas_pull(config["anatomy"], base, config["source_space"])
    result = resample(source, reference, pull, label=True)
    sitk.WriteImage(result, str(output))
    if not np.array_equal(sitk.GetArrayViewFromImage(sitk.ReadImage(str(output))),
                          sitk.GetArrayViewFromImage(result)):
        raise ValueError("Mapped label did not survive readback")
    return {"output": str(output), "sha256": digest(output),
            "voxels": int(sitk.GetArrayViewFromImage(result).sum()),
            "source_space": config["source_space"], "interpolations": 1}
