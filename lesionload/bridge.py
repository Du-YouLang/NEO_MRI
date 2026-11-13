"""Numerical inverse of the pinned TemplateFlow bridge on its target lattice."""

import nibabel as nib
import numpy as np
import SimpleITK as sitk
from scipy.ndimage import map_coordinates

from .spatial import RAS_LPS, transform_points

def check_bridge(forward, inverse, union, affine):
    zyx = np.argwhere(union)[::7]
    lps = nib.affines.apply_affine(affine, zyx[:, ::-1]) * RAS_LPS
    mapped = transform_points(forward, lps)
    error = np.linalg.norm(transform_points(inverse, mapped) - lps, axis=1)
    if not np.isfinite(error).all():
        raise ValueError("Nonfinite TemplateFlow roundtrip")
    return dict(samples=len(error), median_mm=float(np.median(error)),
                p95_mm=float(np.percentile(error, 95)), max_mm=float(error.max()))


def invert_template_bridge(published, target, union, ras_affine):
    """Invert A(x + u(x)) on the original SMATT output grid."""
    if (published.GetNumberOfTransforms() != 2
            or published.GetNthTransform(0).GetName() != "AffineTransform"
            or published.GetNthTransform(1).GetName() != "DisplacementFieldTransform"):
        raise ValueError("Unexpected published TemplateFlow transform structure")
    affine_inv = sitk.AffineTransform(published.GetNthTransform(0).GetInverse())
    matrix = np.array(affine_inv.GetMatrix()).reshape(3, 3)
    offset = affine_inv.TransformPoint((0., 0., 0.))
    displacement = sitk.DisplacementFieldTransform(published.GetNthTransform(1))
    field = displacement.GetDisplacementField()
    field_geometry = np.array(field.GetDirection()).reshape(3, 3) * field.GetSpacing()
    lps_to_index = np.linalg.inv(field_geometry)
    components = sitk.GetArrayViewFromImage(field)
    if not np.isfinite(components).all():
        raise ValueError("Nonfinite published displacement")

    nz, ny, nx = target.GetSize()[::-1]
    vector = np.empty((nz, ny, nx, 3), dtype=np.float64)
    tolerance = .0001
    maximum_iterations = 0
    maximum_residual = 0.
    for first in range(0, nz, 6):
        last = min(first + 6, nz)
        zyx = np.indices((last - first, ny, nx), dtype=np.float64).reshape(3, -1).T
        zyx[:, 0] += first
        output = nib.affines.apply_affine(ras_affine, zyx[:, ::-1]) * RAS_LPS
        base = output @ matrix.T + offset
        solution = base.copy()
        active = np.arange(len(output))
        for iteration in range(1, 81):
            x = solution[active]
            index = (x - field.GetOrigin()) @ lps_to_index.T
            inside = ((index >= 0) & (index <= np.array(field.GetSize()) - 1)).all(axis=1)
            if not inside.all():
                raise ValueError("Inverse solve left the sampled TemplateFlow displacement domain")
            shift = np.column_stack([map_coordinates(components[..., axis], index[:, ::-1].T,
                                                      order=1, prefilter=False) for axis in range(3)])
            delta = base[active] - x - shift
            residual = np.linalg.norm(delta, axis=1)
            moving = residual > tolerance
            maximum_residual = max(maximum_residual, float(residual[~moving].max(initial=0.)))
            active = active[moving]
            if not len(active):
                maximum_iterations = max(maximum_iterations, iteration)
                break
            solution[active] += .8 * delta[moving]
        else:
            raise ValueError(f"TemplateFlow inverse unconverged on {len(active)} voxels in slices {first}:{last}")
        vector[first:last] = (solution - output).reshape(last - first, ny, nx, 3)

    k, j, i = np.nonzero(union)
    if np.any((k == 0) | (k == nz - 1) | (j == 0) | (j == ny - 1) | (i == 0) | (i == nx - 1)):
        raise ValueError("SMATT union reaches a grid edge needed for Jacobian differences")
    jac = np.stack(((vector[k, j, i + 1] - vector[k, j, i - 1]) / 2,
                    (vector[k, j - 1, i] - vector[k, j + 1, i]) / 2,
                    (vector[k + 1, j, i] - vector[k - 1, j, i]) / 2), axis=-1)
    determinant = np.linalg.det(jac + np.eye(3))
    if not np.isfinite(determinant).all() or np.any(determinant <= 0):
        raise ValueError("Nonpositive numerical-inverse Jacobian in the SMATT union")
    image = sitk.GetImageFromArray(vector, isVector=True)
    image.CopyInformation(target)
    inverted = sitk.DisplacementFieldTransform(image)
    check = check_bridge(inverted, published, union, ras_affine)
    if check["max_mm"] > .001:
        raise ValueError(f"Numerical inverse failed independent ITK point check: {check}")
    return inverted, dict(**check, solved_voxels=nz * ny * nx, max_iterations=maximum_iterations,
                          max_solver_residual_mm=maximum_residual, min_atlas_jacobian=float(determinant.min()),
                          nonpositive_atlas_jacobians=0)
