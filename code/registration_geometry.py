"""Physical-coordinate deformation and surface diagnostics; no registration."""
import numpy as np
from scipy.ndimage import binary_erosion,distance_transform_edt

def jacobian_from_displacement(displacement_zyx_world,voxel_to_world):
    """det(I + du/dp), vector u in world components, p in physical mm.

    Array axes are z,y,x; voxel_to_world columns describe x,y,z index steps.
    This explicitly handles oblique/rotated reference grids, unlike treating
    vector components as voxel-axis components. Includes affine volume scale.
    """
    u=np.asarray(displacement_zyx_world);a=np.asarray(voxel_to_world,float)
    if u.ndim!=4 or u.shape[-1]!=3 or min(u.shape[:3])<2 or a.shape!=(3,3) or not np.isfinite(a).all() or abs(np.linalg.det(a))<1e-12 or not np.isfinite(u).all():raise ValueError('Finite 3D world-vector field and nonsingular physical grid required')
    gradient=np.empty(u.shape[:3]+(3,3),dtype=np.float32)
    for component in range(3):
        for index_axis,array_axis in enumerate([2,1,0]):gradient[...,component,index_axis]=np.gradient(u[...,component],axis=array_axis)
    world=np.einsum('...ij,jk->...ik',gradient,np.linalg.inv(a));world+=np.eye(3)
    return np.linalg.det(world)

def symmetric_surface_distances(a,b,spacing_zyx,evaluation_domain):
    """Extract real surfaces before restriction, avoiding halo-cut surfaces."""
    a=np.asarray(a,bool);b=np.asarray(b,bool);m=np.asarray(evaluation_domain,bool)
    if a.shape!=b.shape or a.shape!=m.shape or a.ndim!=3:raise ValueError('Same-grid 3D masks required')
    sa=(a&~binary_erosion(a))&m;sb=(b&~binary_erosion(b))&m
    if not sa.any() or not sb.any():return dict(status='INSUFFICIENT_SURFACE',mean_symmetric_mm=None,HD95_mm=None)
    distances=np.r_[distance_transform_edt(~sb,sampling=spacing_zyx)[sa],distance_transform_edt(~sa,sampling=spacing_zyx)[sb]]
    return dict(status='COMPUTED_DESCRIPTIVE_PROXY',mean_symmetric_mm=float(distances.mean()),HD95_mm=float(np.percentile(distances,95)),maximum_mm=float(distances.max()))

def dice(a,b,domain):
    a=np.asarray(a,bool)&domain;b=np.asarray(b,bool)&domain;den=a.sum()+b.sum()
    return float(2*(a&b).sum()/den) if den else None
