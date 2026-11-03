"""Reproducible tiny synthetic engineering example; no patient inputs."""
from pathlib import Path
import json, sys
import numpy as np
import pandas as pd
import nibabel as nib
import staged_workflow as w

def create_case(tmp_path):
    project=tmp_path/'source';base=project/'derivatives/fmriprep';func=base/'sub-test/func';func.mkdir(parents=True)
    n=100;shape=(6,6,6);aff=np.diag([2.,2.,2.,1.]);rng=np.random.default_rng(883)
    stem='sub-test_task-test_run-01';prefix=stem+'_space-MNI152NLin2009cAsym_res-2'
    image=nib.Nifti1Image((100+rng.normal(size=shape+(n,))).astype(np.float32),aff);image.header.set_zooms((2,2,2,2))
    bold=func/(prefix+'_desc-preproc_bold.nii.gz');image.to_filename(bold)
    mask=func/(prefix+'_desc-brain_mask.nii.gz');nib.Nifti1Image(np.ones(shape,np.uint8),aff).to_filename(mask)
    df=pd.DataFrame();motion=['trans_x','trans_y','trans_z','rot_x','rot_y','rot_z']
    for k in motion:
        x=rng.normal(size=n);d=np.r_[np.nan,np.diff(x)]
        df[k]=x;df[k+'_derivative1']=d;df[k+'_power2']=x*x;df[k+'_derivative1_power2']=d*d
    meta={}
    for j in range(5):name=f'a_comp_cor_{j:02d}';df[name]=rng.normal(size=n);meta[name]={'Mask':'combined','Retained':True}
    df['framewise_displacement']=np.r_[np.nan,np.zeros(n-1)];df['std_dvars']=np.r_[np.nan,np.ones(n-1)]
    conf=func/(stem+'_desc-confounds_timeseries.tsv');df.to_csv(conf,sep='\t',index=False);w.write(conf.with_suffix('.json'),meta)
    events=tmp_path/'events.tsv';pd.DataFrame({'onset':[10,50,90,130,170],'duration':[10]*5,'trial_type':['movement']*5}).to_csv(events,sep='\t',index=False)
    proof=project/'subset_provenance.json';w.write(proof,dict(kind='SYNTHETIC_ENGINEERING_INPUTS',clinical=False,RepetitionTime=2.))
    anat=base/'T1.nii.gz';nib.Nifti1Image(np.ones(shape,np.float32),aff).to_filename(anat)
    rel=lambda p:p.relative_to(base).as_posix()
    c=dict(schema='neo-bold-stages-1',case_id='synthetic',workspace='results/state',preprocessing=dict(mode='engineering_subset',project='source'),
        analysis=dict(timing=dict(status='engineering_subset',evidence='source/subset_provenance.json',interpretation='Synthetic test, not human activation'),anatomy=rel(anat),
            runs=[dict(id='test_1',bold=rel(bold),mask=rel(mask),confounds=rel(conf),confounds_json=rel(conf.with_suffix('.json')),events='events.tsv',tr_seconds=2.,tr_evidence='source/subset_provenance.json')],
            statistics=dict(hrf_model='spm + derivative',smoothing_fwhm_mm=0.,noise_model='ols',compcor_components=5,fd_threshold_mm=.5,dvars_threshold=1.5,contrast={'movement':1.},fdr_q=.05)),report=dict(title='Synthetic tests'))
    p=tmp_path/'config.json';w.write(p,c);c=w.load(p);return c

def main(argv=None):
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args(argv);root=args.output.resolve()
    if root.exists():raise SystemExit('Choose a new output directory; no previous result is overwritten')
    root.mkdir(parents=True);w.REPO=root;c=create_case(root)
    actions=[w.dispatch(c,stage) for stage in ['preprocess','audit','glm','report']]
    actions.append(w.dispatch(c,'glm'));actions.append(w.dispatch(c,'report'))
    record=dict(scope='SYNTHETIC_6x6x6_100_FRAMES_NOT_PATIENT_OR_FMRIPREP_VALIDATION',actions=actions)
    w.write(root/'demo_receipt.json',record);print(json.dumps(record,ensure_ascii=False));return 0
