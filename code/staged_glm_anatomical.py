"""Versioned anatomical-intersection engine; original baseline engine is unchanged."""
from pathlib import Path
import gc,importlib.metadata as md,json,warnings
import nibabel as nib
import numpy as np
import pandas as pd
from scipy import stats
from nilearn.glm.first_level import FirstLevelModel,make_first_level_design_matrix
from nilearn.glm import compute_fixed_effects,threshold_stats_img
from staged_workflow import derivatives,resolve,write,StageError,task_gate,THREADS

MOTION=['trans_x','trans_y','trans_z','rot_x','rot_y','rot_z']
MOTION24=MOTION+[x+'_derivative1' for x in MOTION]+[x+'_power2' for x in MOTION]+[x+'_derivative1_power2' for x in MOTION]

def input_path(c,r,key):
    p=(derivatives(c)/r[key]).resolve()
    if not p.is_relative_to(derivatives(c).resolve()):raise StageError('Run input lies outside registered derivatives')
    if not p.is_file():raise StageError('Missing '+key+': '+str(p))
    return p

def prepare_design(c,r,n,tr):
    stem=Path(r['bold']).name.split('_space-')[0]
    if not Path(r['confounds']).name.startswith(stem+'_desc-confounds_timeseries') or Path(r['confounds_json']).name!=Path(r['confounds']).with_suffix('.json').name:
        raise StageError('BOLD/confounds belong to different run names')
    s=c['analysis']['statistics'];df=pd.read_csv(input_path(c,r,'confounds'),sep='\t')
    meta=json.loads(input_path(c,r,'confounds_json').read_text(encoding='utf-8'))
    if len(df)!=n:raise StageError('BOLD / confounds frame mismatch')
    events=pd.read_csv(resolve(r['events']),sep='\t')
    if not {'onset','duration','trial_type'}<=set(events):raise StageError('Events require onset, duration, trial_type')
    if not len(events) or not np.isfinite(events[['onset','duration']].to_numpy()).all():raise StageError('Empty/nonfinite events')
    if (events.onset<0).any() or (events.duration<0).any() or (events.onset+events.duration).max()>n*tr+1e-4:
        raise StageError('Events extend outside the saved BOLD duration; no trimming/rescaling is allowed')
    cc=[x for x in df if x.startswith('a_comp_cor_') and meta.get(x,{}).get('Mask')=='combined' and meta.get(x,{}).get('Retained',True)]
    count=s['compcor_components']
    if not isinstance(count,int) or count<0 or len(cc)<count:raise StageError('Insufficient actual combined-mask aCompCor components')
    columns=MOTION24+cc[:count]+[x for x in df if x.startswith('cosine')]
    nuisance=df[columns].copy()
    for name in columns:
        bad=nuisance[name].isna()
        if bad.any() and not (bad.sum()==1 and bad.iloc[0] and '_derivative1' in name):
            raise StageError('Unexpected missing confound: '+name)
    nuisance=nuisance.fillna(0.)
    for name in columns:
        sd=float(nuisance[name].std(ddof=0))
        if not np.isfinite(sd) or sd<=1e-12:raise StageError('Constant/invalid confound; no automatic removal: '+name)
        nuisance[name]=(nuisance[name]-nuisance[name].mean())/sd
    spike=(df.framewise_displacement.fillna(0)>s['fd_threshold_mm'])|(df.std_dvars.fillna(0)>s['dvars_threshold'])
    ns=[x for x in df if x.startswith('non_steady_state_outlier')]
    if ns:spike|=df[ns].fillna(0).sum(axis=1)>0
    indices=np.flatnonzero(spike.to_numpy())
    for k in indices:nuisance[f'spike_{k:03d}']=(np.arange(n)==k).astype(float)
    design=make_first_level_design_matrix(np.arange(n)*tr,events,hrf_model=s['hrf_model'],drift_model=None,
        add_regs=nuisance.to_numpy(),add_reg_names=list(nuisance.columns))
    x=design.to_numpy();rank=int(np.linalg.matrix_rank(x));singular=np.linalg.svd(x,compute_uv=False)
    if rank!=x.shape[1] or n<=rank:raise StageError(f'Non-estimable design: {n} observations, {x.shape[1]} columns, rank {rank}')
    condition=float(singular[0]/singular[-1])
    if condition>=1e8:raise StageError('Ill-conditioned design; substantive model review required')
    weights=s['contrast'];missing=set(weights)-set(design)
    if missing:raise StageError('Contrast references missing conditions: '+str(missing))
    contrast=np.array([weights.get(col,0.) for col in design.columns],float)
    if not np.any(contrast):raise StageError('All-zero contrast')
    diag=dict(run=r['id'],volumes=n,TR_seconds=tr,columns=x.shape[1],rank=rank,nominal_residual_dof=n-rank,
        condition_number=condition,combined_compcor=cc[:count],cosine_columns=[x for x in columns if x.startswith('cosine')],
        spikes=indices.tolist(),event_count=len(events),event_end_seconds=float((events.onset+events.duration).max()))
    return design,contrast,diag,events,nuisance

def save_maps(folder,maps):
    folder.mkdir(parents=True,exist_ok=True)
    names={'effect_size':'effect','effect_variance':'variance','stat':'t','z_score':'z','p_value':'one_sided_p'}
    for name,im in maps.items():im.to_filename(folder/(names[name]+'.nii.gz'))

def inference(folder,zimg,mask,s):
    binary=mask.get_fdata()>.5;z=zimg.get_fdata()[binary]
    if not np.isfinite(z).all():raise StageError('Nonfinite Z values inside analysis mask')
    q=s['fdr_q']
    if not 0<q<1:raise StageError('FDR q must be between zero and one')
    with warnings.catch_warnings(record=True) as ws:
        thimg,threshold=threshold_stats_img(zimg,mask_img=mask,alpha=q,height_control='fdr',two_sided=True,cluster_threshold=0)
    thimg.to_filename(folder/'z_fdr.nii.gz')
    p=np.minimum(1.,2*stats.norm.sf(np.abs(z)));ix=np.argsort(p);adj=np.minimum.accumulate((p[ix]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]
    corrected=np.empty(len(p));corrected[ix]=np.minimum(1.,adj)
    if not np.array_equal(corrected<q,np.abs(thimg.get_fdata()[binary])>0):raise StageError('Independent FDR q-map mismatch')
    data=np.ones(binary.shape,np.float32);data[binary]=corrected
    nib.Nifti1Image(data,mask.affine).to_filename(folder/'q_fdr.nii.gz')
    info=dict(fdr_q=q,two_sided=True,analysis_voxels=int(binary.sum()),
        z_threshold=float(threshold) if np.isfinite(threshold) else None,
        positive_voxels=int(np.sum(thimg.get_fdata()[binary]>0)),negative_voxels=int(np.sum(thimg.get_fdata()[binary]<0)),
        peak_positive_z=float(z.max()),peak_negative_z=float(z.min()),warnings=[str(x.message) for x in ws])
    write(folder/'statistics.json',info);return info

def fit(c,dest):
    task_gate(c);a=c['analysis'];s=a['statistics'];runs=a['runs']
    if not runs:raise StageError('No task runs configured')
    if s['noise_model'] not in ['ar1','ols']:raise StageError('Unsupported noise model')
    masks=[nib.load(input_path(c,r,'mask')) for r in runs];ref=masks[0]
    for m in masks:
        if m.shape!=ref.shape or not np.allclose(m.affine,ref.affine):raise StageError('Run brain masks have incompatible grids')
    binary=np.logical_and.reduce([m.get_fdata()>.5 for m in masks])
    anatomical=nib.load(input_path(c,{'mask':a['anatomical_mask']},'mask'))
    if anatomical.shape!=ref.shape or not np.allclose(anatomical.affine,ref.affine):
        raise StageError('Anatomical/functional mask grid mismatch; no implicit resampling')
    anatomy=anatomical.get_fdata()>.5
    functional=binary.copy();binary &= anatomy
    write(dest/'mask_summary.json',dict(kind='FUNCTIONAL_INTERSECT_ANATOMICAL',
        functional_voxels=int(functional.sum()),anatomical_voxels=int(anatomy.sum()),
        intersection_voxels=int(binary.sum()),removed_from_functional=int((functional & ~anatomy).sum()),
        note='Statistical domain restriction; does not repair geometric registration.'))
    minimum=32 if a['timing']['status']=='engineering_subset' else 10000
    if binary.sum()<minimum:raise StageError('Insufficient common-mask coverage')
    mask=nib.Nifti1Image(binary.astype(np.uint8),ref.affine);mask.to_filename(dest/'common_mask.nii.gz')
    effects=[];variances=[];dofs=[];records=[]
    for r in runs:
        bold=input_path(c,r,'bold');im=nib.load(bold);header_tr=float(im.header.get_zooms()[3]);tr=float(r['tr_seconds']);n=im.shape[3]
        if tr<=0 or not np.isclose(header_tr,tr,rtol=1e-5,atol=1e-6):raise StageError('Nominal TR / NIfTI header TR conflict; no inference/rescaling is allowed')
        if im.shape[:3]!=ref.shape or not np.allclose(im.affine,ref.affine):raise StageError('BOLD/mask grid mismatch')
        design,contrast,diag,events,nuisance=prepare_design(c,r,n,tr)
        folder=dest/r['id'];folder.mkdir()
        design.to_csv(folder/'design.tsv',sep='\t',index_label='time_seconds');events.to_csv(folder/'events.tsv',sep='\t',index=False)
        nuisance.to_csv(folder/'nuisance.tsv',sep='\t',index=False);write(folder/'diagnostics.json',diag)
        # One job avoids nested pools; BLAS/OMP/ITK are fixed to four threads by
        # the entry point before scientific libraries are imported.
        model=FirstLevelModel(t_r=tr,slice_time_ref=0.,hrf_model=s['hrf_model'],drift_model=None,
            mask_img=mask,smoothing_fwhm=s['smoothing_fwhm_mm'],noise_model=s['noise_model'],
            standardize=False,signal_scaling=0,n_jobs=1,minimize_memory=True,random_state=20261026,verbose=0)
        model.fit(str(bold),design_matrices=design)
        maps=model.compute_contrast(contrast,stat_type='t',output_type='all');save_maps(folder,maps)
        stats_info=inference(folder,maps['z_score'],mask,s)
        effects.append(maps['effect_size']);variances.append(maps['effect_variance']);dofs.append(diag['nominal_residual_dof'])
        records.append(dict(id=r['id'],diagnostics=diag,statistics=stats_info));del model,maps;gc.collect()
    combined=dest/'combined';combined.mkdir()
    if len(runs)>1:
        eff,var,t,z=compute_fixed_effects(effects,variances,mask=mask,precision_weighted=True,dofs=dofs,return_z_score=True)
        values=np.stack([v.get_fdata()[binary] for v in variances]);betas=np.stack([v.get_fdata()[binary] for v in effects])
        if not np.isfinite(values).all() or (values<=0).any():raise StageError('Nonpositive/invalid contrast variance')
        w=1/values
        if not np.allclose(eff.get_fdata()[binary],(w*betas).sum(axis=0)/w.sum(axis=0),rtol=1e-5,atol=1e-6):raise StageError('Fixed-effects independent check failed')
        if not np.allclose(var.get_fdata()[binary],1/w.sum(axis=0),rtol=1e-5,atol=1e-6):raise StageError('Fixed-effects variance check failed')
        for tag,img in [('effect',eff),('variance',var),('t',t),('z',z)]:img.to_filename(combined/(tag+'.nii.gz'))
        summary_stats=inference(combined,z,mask,s)
    else:
        import shutil
        for p in (dest/runs[0]['id']).glob('*.nii.gz'):shutil.copy2(p,combined/p.name)
        summary_stats=records[0]['statistics']
    result=dict(status='COMPLETED_NUMERICALLY_REVIEW_PENDING',case_id=c['case_id'],runs=records,combined=summary_stats,
        statistics=s,timing=a['timing'],tool_threads=THREADS,n_jobs=1,preprocessing_rerun=False,
        scientific_scope='ENGINEERING_SUBSET_NOT_CLINICAL' if a['timing']['status']=='engineering_subset' else 'Single-subject task model; no clinical eligibility inference',
        software={x:md.version(x) for x in ['nilearn','nibabel','numpy','scipy','pandas','scikit-learn']})
    write(dest/'summary.json',result);return result
