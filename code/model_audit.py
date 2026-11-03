"""Design-only audit. Never reads BOLD voxels, fits activation, or edits designs."""
import numpy as np

def audit_design(design,contrast,task_columns):
    x=np.asarray(design,dtype=float);c=np.asarray(contrast,dtype=float)
    if x.ndim!=2 or c.shape!=(x.shape[1],) or not np.isfinite(x).all() or not np.isfinite(c).all():raise ValueError('Finite design and matching contrast required')
    names=list(map(str,design.columns));n,p=x.shape
    if not n or not p or len(set(names))!=p:raise ValueError('Nonempty unique design columns required')
    norms=np.linalg.norm(x,axis=0);safe=np.where(norms>0,norms,1.)
    a=x/safe;s=np.linalg.svd(a,compute_uv=False);rank=int(np.linalg.matrix_rank(a));cs=c/safe
    projection=np.linalg.pinv(a)@a
    estimable=bool(np.allclose(projection@cs,cs,rtol=1e-7,atol=1e-9))
    variance=float(np.linalg.norm(cs@np.linalg.pinv(a))**2) if estimable else None
    details=[]
    for name in task_columns:
        if name not in names:raise ValueError('Missing task column '+name)
        j=names.index(name);v=x[:,j];center=v-v.mean();total=float(center@center)
        other=np.delete(a,j,axis=1);res=v-other@np.linalg.lstsq(other,v,rcond=None)[0];rss=float(res@res)
        correlations=[]
        for k in range(p):
            if k==j:continue
            z=x[:,k]-x[:,k].mean();den=np.linalg.norm(center)*np.linalg.norm(z)
            if den>1e-12:correlations.append((names[k],float(center@z/den)))
        correlations.sort(key=lambda z:abs(z[1]),reverse=True)
        details.append(dict(column=name,VIF=total/rss if total>1e-12 and rss>total*1e-12 else None,
            VIF_undefined_or_infinite=total<=1e-12 or rss<=total*1e-12,largest_correlations=correlations[:5]))
    return dict(observations=n,columns=p,rank=rank,nominal_residual_dof=n-rank,full_column_rank=rank==p,
        condition_number_unit_norm=float(s[0]/s[-1]) if rank==p else None,
        contrast_estimable=estimable,unit_noise_OLS_contrast_variance=variance,task_columns=details,
        column_groups={k:sum(test(v) for v in names) for k,test in {
            'spike':lambda v:v.startswith('spike_'),'cosine':lambda v:v.startswith('cosine'),
            'compcor':lambda v:v.startswith('a_comp_cor_'),'constant':lambda v:v=='constant'}.items()},
        interpretation='Design geometry only; OLS unit independent noise reference, not fitted AR1 uncertainty, power, overfitting diagnosis, or independent timepoint count. No automatic regressor removal.')

def audit_case(c):
    from pathlib import Path
    import nibabel as nib
    from staged_glm import input_path,prepare_design
    from staged_workflow import StageError,sha
    records=[]
    for r in c['analysis']['runs']:
        im=nib.load(input_path(c,r,'bold'));tr=float(r['tr_seconds'])
        if tr<=0 or not np.isclose(im.header.get_zooms()[3],tr,rtol=1e-5,atol=1e-6):raise StageError('Nominal TR / NIfTI header TR conflict')
        try:design,contrast,diag,_,_=prepare_design(c,r,im.shape[3],tr)
        except StageError as exc:
            records.append(dict(run=r['id'],TR_seconds=tr,status='BLOCKED_DESIGN_REVIEW_REQUIRED',reason=str(exc)))
            continue
        a=audit_design(design,contrast,list(c['analysis']['statistics']['contrast']))
        a.update(run=r['id'],TR_seconds=tr,design_diagnostics=diag,status='AUDITED_DESIGN_ONLY');records.append(a)
    return dict(runs=records,fit_allowed=bool(records) and all(r['status']=='AUDITED_DESIGN_ONLY' for r in records),
        audit_engine_sha256=sha(Path(__file__)),BOLD_voxels_read=False,GLM_fits=0,preprocessing_invocations=0)
