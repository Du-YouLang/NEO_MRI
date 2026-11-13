"""Small generated task BOLD dataset for the tutorial and CPU tests."""

import nibabel as nib
import numpy as np
import pandas as pd
from nilearn.glm.first_level import compute_regressor

import bold_workflow as w
from bold_ab import LAGS, RATIOS
from bold_models import MOTION


def make_dataset(tmp_path):
    project = tmp_path/"source"
    derivatives = project/"derivatives/fmriprep"
    func = derivatives/"sub-synthetic/func"
    func.mkdir(parents=True)
    n = 110
    rng = np.random.default_rng(71)
    shape = (4,4,4)
    affine = np.diag([2.,2.,2.,1.])
    mask = np.zeros(shape,np.uint8)
    mask[1:3,1:3,1:3] = 1
    onsets = np.array([14.,65.,122.,171.])
    durations = np.full(len(onsets),10.)
    reg,_ = compute_regressor(np.vstack((onsets+3.,durations*1.15,np.ones(len(onsets)))),
                              "spm",np.arange(n)*2,oversampling=50)
    w.write(project/"subset_provenance.json",{"kind":"SYNTHETIC_ENGINEERING_INPUTS","RepetitionTime":2.})
    t1 = derivatives/"synthetic_T1.nii.gz"
    t1_mask = derivatives/"synthetic_T1_mask.nii.gz"
    nib.save(nib.Nifti1Image(np.arange(64,dtype=np.float32).reshape(shape)+2,affine),t1)
    nib.save(nib.Nifti1Image(np.ones(shape,np.uint8),affine),t1_mask)
    runs = []
    for i,task in enumerate(("MI","PM"),1):
        stem = f"sub-synthetic_task-test_run-{i:02d}"
        prefix = stem+"_space-MNI152NLin2009cAsym_res-2"
        bold = func/(prefix+"_desc-preproc_bold.nii.gz")
        img = nib.Nifti1Image((100+reg[:,0]*4+rng.normal(size=shape+(n,))*.35).astype(np.float32),affine)
        img.header.set_zooms((2,2,2,2))
        img.header.set_xyzt_units("mm","sec")
        nib.save(img,bold)
        mask_path = func/(prefix+"_desc-brain_mask.nii.gz")
        nib.save(nib.Nifti1Image(mask,affine),mask_path)
        sidecar = func/(prefix+"_desc-preproc_bold.json")
        w.write(sidecar,dict(RepetitionTime=2.,SliceTimingCorrected=False))
        data = {}
        for motion in MOTION:
            a = rng.normal(size=n)
            diff = np.r_[np.nan,np.diff(a)]
            data.update({motion:a,motion+"_derivative1":diff,motion+"_power2":a*a,motion+"_derivative1_power2":diff*diff})
        data.update(framewise_displacement=np.r_[np.nan,np.zeros(n-1)],std_dvars=np.r_[np.nan,np.ones(n-1)])
        confounds = func/(stem+"_desc-confounds_timeseries.tsv")
        pd.DataFrame(data).to_csv(confounds,sep="\t",index=False)
        w.write(confounds.with_suffix(".json"),{})
        events = tmp_path/f"events_{i}.tsv"
        pd.DataFrame(dict(onset=onsets,duration=durations,trial_type=task)).to_csv(events,sep="\t",index=False)
        rel = lambda path:path.relative_to(derivatives).as_posix()
        runs.append(dict(id=f"run_{i}",task_id=task,estimand_id=task,reference_space="MNI152NLin2009cAsym",
                         tr_seconds=2.,stc="off",sdc="off",bold=rel(bold),mask=rel(mask_path),bold_json=rel(sidecar),
                         confounds=rel(confounds),confounds_json=rel(confounds.with_suffix(".json")),
                         viewer_anatomy=rel(t1),viewer_anatomy_mask=rel(t1_mask),events=events.name,
                         tr_evidence="source/subset_provenance.json",dynamic_condition=task,contrast={task:1}))
    config = dict(schema="neo-bold-stages-2",case_id="synthetic",data_root="data",workspace="data/synthetic/s01/results/task_fMRI/ab",
        preprocessing=dict(project="source",mode="engineering_subset"),
        analysis=dict(runs=runs,kernels_mm=[0,2,3],fir_delays_tr=[0,1,2,3],
                      timing=dict(status="engineering_subset",evidence="source/subset_provenance.json",interpretation="Synthetic fixture"),
                      double_scanner=dict(lags_seconds=list(LAGS),duration_ratios=list(RATIOS),selection="raw_r2"),
                      double_scanner_ab=dict(seed=31,device="cpu"),
                      statistics=dict(hrf_model="spm + derivative",smoothing_fwhm_mm=2,noise_model="ar1",compcor_components=0,
                                      fd_threshold_mm=.5,dvars_threshold=1.5,contrast={"MI":1},fdr_q=.05),
                      comparisons=[dict(id="mi_pm",left="run_1",right="run_2",model="DoubleScanner_AB")]))
    path = tmp_path/"config.json"
    w.write(path,config)
    return path
