"""Small-array metrics, no registration execution or automatic branch winner."""
import numpy as np
def roundtrip_coordinates(original_mm,returned_mm):
    """Inputs already mapped through forward then inverse physical transforms."""
    x=np.asarray(original_mm,float);y=np.asarray(returned_mm,float)
    if x.shape!=y.shape or x.ndim!=2 or x.shape[1]!=3 or not x.size or not np.isfinite(x).all() or not np.isfinite(y).all():raise ValueError('Matched finite physical 3D coordinates required')
    d=np.linalg.norm(x-y,axis=1)
    return dict(mean_mm=float(d.mean()),maximum_mm=float(d.max()),values_mm=d.tolist(),anatomical_accuracy_proven=False)

def evaluate(original_lesion,roundtrip_lesion,source_landmarks_mm=None,target_landmarks_mm=None,jacobian=None,nonlesion_mask=None):
    arrays=[np.asarray(v) for v in [original_lesion,roundtrip_lesion]]
    if any(v.ndim!=3 or not np.isfinite(v).all() or not np.isin(v,[0,1]).all() for v in arrays):raise ValueError('Finite binary 3D lesion masks required')
    a,b=[v.astype(bool) for v in arrays]
    if a.shape!=b.shape or not a.any():raise ValueError('Matching grids and nonempty lesion required')
    out=dict(roundtrip_Dice=float(2*(a&b).sum()/(a.sum()+b.sum())),automatic_winner=False,
        interpretation='Inverse consistency alone does not establish anatomical accuracy.',anatomical_review='NOT_PROVIDED')
    if (source_landmarks_mm is None)!=(target_landmarks_mm is None):raise ValueError('Both sets of anatomical landmarks required')
    if source_landmarks_mm is not None:
        x=np.asarray(source_landmarks_mm,float);y=np.asarray(target_landmarks_mm,float)
        if x.shape!=y.shape or x.ndim!=2 or x.shape[1]!=3 or not x.size or not np.isfinite(x).all() or not np.isfinite(y).all():raise ValueError('Matched finite physical 3D landmarks required')
        # source landmarks already transformed into target physical coordinates.
        d=np.linalg.norm(x-y,axis=1);out['landmark_TRE_mm']=dict(mean=float(d.mean()),maximum=float(d.max()),values=d.tolist())
    if jacobian is not None:
        j=np.asarray(jacobian,float)
        if nonlesion_mask is None:raise ValueError('Explicit nonlesion evaluation domain required')
        m=np.asarray(nonlesion_mask,bool)
        if m.shape!=j.shape or not m.any() or not np.isfinite(j[m]).all():raise ValueError('Invalid Jacobian domain')
        v=j[m];out['jacobian']=dict(nonpositive_voxels=int((v<=0).sum()),minimum=float(v.min()),p01=float(np.percentile(v,1)),p99=float(np.percentile(v,99)),maximum=float(v.max()))
    return out
