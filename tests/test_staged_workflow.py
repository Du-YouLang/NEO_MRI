"""Meaningful stage/cache/timing guards and independent OLS numeric reference."""
from pathlib import Path
import copy,json,sys
import numpy as np
import pandas as pd
import nibabel as nib
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'code'))
import staged_workflow as w

@pytest.fixture
def case(tmp_path,monkeypatch):
    monkeypatch.setattr(w,'REPO',tmp_path)
    from synthetic_demo import create_case
    return create_case(tmp_path)

def forbid_preprocessing(monkeypatch):
    import staged_preprocess
    def fail(*a,**k):raise AssertionError('GLM/report attempted to launch preprocessing')
    monkeypatch.setattr(staged_preprocess,'launch_new',fail)

def test_no_implicit_preprocessing(case,monkeypatch):
    forbid_preprocessing(monkeypatch)
    with pytest.raises(w.StageError,match='explicit preprocess'):w.dispatch(case,'glm')

def test_cache_parameters_and_reports_are_separate(case,monkeypatch):
    forbid_preprocessing(monkeypatch);w.dispatch(case,'preprocess')
    first=w.dispatch(case,'glm');second=w.dispatch(case,'glm')
    assert first['glm_fits']==1 and second['glm_fits']==0 and first['key']==second['key']
    cosmetic=copy.deepcopy(case);cosmetic['report']['title']='new title'
    assert w.dispatch(cosmetic,'glm')['action']=='REUSED_IDENTICAL_STATISTICS'
    changed=copy.deepcopy(case);changed['analysis']['statistics']['compcor_components']=4
    newer=w.dispatch(changed,'glm');assert newer['key']!=first['key'] and newer['fmriprep_invocations']==0
    import staged_glm
    monkeypatch.setattr(staged_glm,'fit',lambda *a:(_ for _ in ()).throw(AssertionError('Report fitted GLM')))
    r=w.dispatch(case,'report');assert r['glm_fits']==0 and r['fmriprep_invocations']==0
    assert w.dispatch(case,'report')['action']=='REUSED_IDENTICAL_REPORT'

def test_ols_effect_and_variance_match_independent_linear_algebra(case):
    w.dispatch(case,'preprocess');out=Path(w.dispatch(case,'glm')['output'])/'test_1'
    x=pd.read_csv(out/'design.tsv',sep='\t',index_col=0).to_numpy();names=list(pd.read_csv(out/'design.tsv',sep='\t',index_col=0))
    image=nib.load(w.derivatives(case)/case['analysis']['runs'][0]['bold']);y=image.get_fdata().reshape(-1,100).T
    scaled=100*(y/y.mean(axis=0)-1);pinv=np.linalg.pinv(x);coef=pinv@scaled;resid=scaled-x@coef
    c=np.array([1 if n=='movement' else 0 for n in names]);effect=c@coef
    variance=((resid*resid).sum(axis=0)/(len(x)-np.linalg.matrix_rank(x)))*(c@pinv@pinv.T@c)
    np.testing.assert_allclose(nib.load(out/'effect.nii.gz').get_fdata().ravel(),effect,rtol=2e-4,atol=2e-5)
    np.testing.assert_allclose(nib.load(out/'variance.nii.gz').get_fdata().ravel(),variance,rtol=2e-4,atol=2e-5)

def test_paused_patient_never_fits_or_creates_glm(case,monkeypatch):
    w.dispatch(case,'preprocess');case['analysis']['timing'].update(status='paused',reason='600 vs 1200 seconds unresolved')
    forbid_preprocessing(monkeypatch)
    with pytest.raises(w.TimingBlocked):w.dispatch(case,'glm')
    assert not (w.workspace(case)/'glm').exists()
    r=w.dispatch(case,'report');assert r['glm_fits']==0

def test_patient_shared_events_require_case_and_input_bound_authorization(case):
    w.dispatch(case,'preprocess')
    case['preprocessing']['mode']='existing'
    case['report']['participant_group']='chronic_stroke'
    case['analysis']['statistics']['contrast']={'image':1.,'abstract':-1.}
    evidence=w.resolve('patient_authorization.json')
    case['analysis']['timing'].update(status='conditional_shared_events_user_authorized',evidence=str(evidence))
    r=case['analysis']['runs'][0]
    auth=dict(status='USER_AUTHORIZED_EXPLORATORY_SHARED_EVENTS',case_id=case['case_id'],
        preprocessing_fingerprint=w.read(w.workspace(case)/'preprocessing/receipt.json')['fingerprint'],
        independent_stimulus_log_verified=False,user_answer='Accept publisher events as an exploratory assumption.',
        event_sha256_by_run={r['id']:w.sha(w.resolve(r['events']))},
        TR_source_sha256_by_run={r['id']:w.sha(w.resolve(r['tr_evidence']))})
    w.write(evidence,auth);assert w.task_gate(case)['status']=='conditional_shared_events_user_authorized'
    for field,value in [('case_id','another_patient'),('preprocessing_fingerprint','wrong')]:
        bad=copy.deepcopy(auth);bad[field]=value;w.write(evidence,bad)
        with pytest.raises(w.TimingBlocked):w.task_gate(case)
    bad=copy.deepcopy(auth);bad['event_sha256_by_run'][r['id']]='wrong';w.write(evidence,bad)
    with pytest.raises(w.TimingBlocked):w.task_gate(case)
    w.write(evidence,auth);case['analysis']['statistics']['contrast']={'movement':1.}
    with pytest.raises(w.TimingBlocked,match='image > abstract'):w.task_gate(case)

def test_events_cannot_be_halved_or_trimmed_implicitly(case):
    w.dispatch(case,'preprocess');p=w.resolve('events.tsv');p.write_text('onset\tduration\ttrial_type\n1000\t10\tmovement\n')
    with pytest.raises(w.StageError,match='outside'):w.dispatch(case,'glm')

def test_tampered_input_and_changed_preprocessing_are_blocked(case):
    w.dispatch(case,'preprocess');p=w.derivatives(case)/case['analysis']['runs'][0]['confounds'];p.write_text(p.read_text()+'\n')
    with pytest.raises(w.StageError,match='changed'):w.dispatch(case,'glm')

def test_changed_preprocessing_mode_cannot_reuse_old_receipt(case):
    w.dispatch(case,'preprocess');case['preprocessing']['mode']='existing'
    with pytest.raises(w.StageError,match='request changed'):w.dispatch(case,'glm')



def test_loaded_blas_and_environment_are_limited_to_four():
    import os
    from threadpoolctl import threadpool_info
    w.configure_threads()
    for name in ['OMP_NUM_THREADS','ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:
        assert os.environ[name]=='4'
    assert all(p['num_threads']<=4 for p in threadpool_info())

def test_full_hash_detects_changes_even_with_preserved_size_and_timestamp(case):
    import os
    w.dispatch(case,'preprocess');p=w.derivatives(case)/case['analysis']['runs'][0]['confounds'];s=p.stat()
    data=bytearray(p.read_bytes());data[-2]=ord('9') if data[-2]!=ord('9') else ord('8');p.write_bytes(data)
    os.utime(p,ns=(s.st_atime_ns,s.st_mtime_ns))
    with pytest.raises(w.StageError,match='Hash changed'):w.dispatch(case,'glm',full=True)






def test_report_on_unfitted_parameters_does_not_show_stale_fit(case,monkeypatch):
    w.dispatch(case,'preprocess');w.dispatch(case,'glm')
    case['analysis']['statistics']['compcor_components']=4
    import staged_glm
    monkeypatch.setattr(staged_glm,'fit',lambda *a:(_ for _ in ()).throw(AssertionError('Report fitted GLM')))
    r=w.dispatch(case,'report');page=Path(r['output']).read_text(encoding='utf-8')
    assert 'has not been fitted' in page and '<figure>' not in page

def test_true_TR_conflict_is_not_treated_as_float_rounding(case):
    p=w.derivatives(case)/case['analysis']['runs'][0]['bold'];im=nib.load(p);data=im.get_fdata(dtype=np.float32)
    replacement=nib.Nifti1Image(data,im.affine);replacement.header.set_zooms((2,2,2,3));replacement.to_filename(p)
    w.dispatch(case,'preprocess')
    with pytest.raises(w.StageError,match='TR conflict'):w.dispatch(case,'glm')

def test_design_audit_never_fits_or_reads_bold_voxels(case,monkeypatch):
    w.dispatch(case,'preprocess')
    expected=w.digest(dict(mode='engineering_subset',project=str(w.source_root(case))))
    assert w.preprocessing_request(case)==expected
    import staged_glm
    monkeypatch.setattr(staged_glm,'fit',lambda *a:pytest.fail('Audit attempted GLM fit'))
    monkeypatch.setattr(nib.Nifti1Image,'get_fdata',lambda *a,**k:pytest.fail('Audit read image voxels'))
    forbid_preprocessing(monkeypatch)
    result=w.dispatch(case,'audit');audit=w.read(result['output'])
    assert audit['BOLD_voxels_read'] is False and result['glm_fits']==0
    assert audit['runs'][0]['contrast_estimable']
    assert not (w.workspace(case)/'glm').exists()

def test_invalid_design_audit_saves_reason_without_running_fit(case,monkeypatch):
    path=w.derivatives(case)/case['analysis']['runs'][0]['confounds']
    df=pd.read_csv(path,sep='\t');df['trans_x']=0.;df.to_csv(path,sep='\t',index=False)
    w.dispatch(case,'preprocess')
    import staged_glm
    monkeypatch.setattr(staged_glm,'fit',lambda *a:pytest.fail('Invalid audit reached voxel fit'))
    result=w.dispatch(case,'audit');audit=w.read(result['output'])
    assert audit['fit_allowed'] is False and 'Constant/invalid confound' in audit['runs'][0]['reason']
    with pytest.raises(w.StageError,match='audit blocked'):w.dispatch(case,'glm')
