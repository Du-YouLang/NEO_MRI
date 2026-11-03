from pathlib import Path
import copy,json,sys,subprocess
import numpy as np
import pandas as pd
import nibabel as nib
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'code'))
import scientific_protocol as s
from model_audit import audit_design
from registration_audit import evaluate,roundtrip_coordinates

def profile(**extra):
    p=dict(schema='neo-mri-protocol-1',execution_enabled=False,normalization='standard',sdc='off',stc='off',input_ids={'t1':'a','bold':'b','lesion':'c','fieldmaps':'d'})
    p.update(extra);return p

def test_legacy_arguments_and_sdc():
    assert s.arguments(None)==['--ignore','slicetiming','fieldmaps','flair']
    args=s.arguments(profile(sdc='syn',bold_metadata='bold.json'))
    assert args[args.index('--use-syn-sdc')+1]=='error'
    assert args[-2:]==['--bold2anat-init','t1w']
    args=s.arguments(profile(sdc='pepolar',bold_metadata='bold.json',fieldmaps=['fmap_epi.json']))
    assert 'fieldmaps' not in args
    with pytest.raises(ValueError):s.arguments(profile(sdc='pepolar'))
    with pytest.raises(ValueError):s.arguments(profile(stc='on'))
    with pytest.raises(ValueError):s.arguments(profile(normalization='cfm'))

def test_dependency_boundaries():
    a=s.plan(profile())['dependency_keys'];b=s.plan(profile(normalization='cfm',t1w='t1.nii.gz',lesion_roi='lesion_roi.nii.gz'))['dependency_keys']
    assert a['t1_reference']==b['t1_reference'] and a['native_bold_confounds']==b['native_bold_confounds']
    assert a['normalization']!=b['normalization'] and a['mni_resampling']!=b['mni_resampling']
    d=s.plan(profile(sdc='syn',bold_metadata='bold.json'))['dependency_keys']
    assert a['normalization']==d['normalization'] and a['native_bold_confounds']!=d['native_bold_confounds']

def test_cfm_real_small_headers_and_binary(tmp_path):
    b=tmp_path/'bids';(b/'sub-s/anat').mkdir(parents=True);(b/'.bidsignore').write_text('*lesion_roi.nii.gz')
    t='sub-s/anat/sub-s_T1w.nii.gz';r='sub-s/anat/sub-s_label-lesion_roi.nii.gz'
    nib.Nifti1Image(np.ones((4,4,4)),np.eye(4)).to_filename(b/t)
    mask=np.zeros((4,4,4));mask[1,1,1]=1;nib.Nifti1Image(mask,np.eye(4)).to_filename(b/r)
    p=profile(normalization='cfm',t1w=t,lesion_roi=r);assert len(s.validate_inputs(tmp_path,p)['inputs'])==2
    moved=np.eye(4);moved[0,3]=2;nib.Nifti1Image(mask,moved).to_filename(b/r)
    with pytest.raises(ValueError,match='grid'):s.validate_inputs(tmp_path,p)
    with pytest.raises(ValueError,match='ROI'):s.validate_inputs(tmp_path,profile(t1w=t))

def test_method_review_does_not_accept_requested_as_observed():
    p=profile(sdc='syn',bold_metadata='bold.json')
    with pytest.raises(ValueError):s.require_review(p,'hash',{})
    r=dict(status='REVIEWED_MATCHES_REQUEST',preprocessing_fingerprint='hash',observed_normalization='standard',observed_sdc='off',reviewer='synthetic',evidence_files=[dict(path='log',sha256='0'*64)])
    with pytest.raises(ValueError,match='fallback'):s.require_review(p,'hash',r)
    r['observed_sdc']='syn';s.require_review(p,'hash',r)

def test_roundtrip_alone_cannot_select_alignment():
    mask=np.zeros((5,5,5));mask[2,2,2]=1
    r=evaluate(mask,mask,[[10,0,0]],[[0,0,0]],np.ones(mask.shape),mask==0)
    assert r['roundtrip_Dice']==1 and r['landmark_TRE_mm']['mean']==10 and r['automatic_winner'] is False
    assert r['jacobian']['nonpositive_voxels']==0
    assert roundtrip_coordinates([[0,0,0]],[[0,0,2]])['maximum_mm']==2
    with pytest.raises(ValueError,match='binary'):evaluate(mask*0.5,mask)

def test_more_repeated_independent_observations_halves_ols_variance():
    x=pd.DataFrame({'task':np.r_[np.zeros(20),np.ones(20)],'constant':np.ones(40)})
    a=audit_design(x,[1,0],['task']);b=audit_design(pd.concat([x,x],ignore_index=True),[1,0],['task'])
    assert b['unit_noise_OLS_contrast_variance']==pytest.approx(a['unit_noise_OLS_contrast_variance']/2)
    assert b['nominal_residual_dof']==78

def test_confound_collinearity_and_nonestimable_contrast():
    rng=np.random.default_rng(8);x=rng.normal(size=40)
    d=pd.DataFrame({'task':x,'motion':x+0.001*rng.normal(size=40),'constant':np.ones(40)})
    a=audit_design(d,[1,0,0],['task']);assert a['task_columns'][0]['VIF']>1000
    d['motion']=x;a=audit_design(d,[1,0,0],['task']);assert not a['contrast_estimable'] and not a['full_column_rank']
    assert a['unit_noise_OLS_contrast_variance'] is None


def test_planner_no_subprocess_or_patient_reads(monkeypatch):
    monkeypatch.setattr(subprocess,'run',lambda *a,**k:pytest.fail('Subprocess called'))
    monkeypatch.setattr(Path,'open',lambda *a,**k:pytest.fail('File opened'))
    assert s.plan(profile())['status']=='PLAN_ONLY_NO_RUNTIME_PROBED'

def test_plan_entry_never_dispatches_or_configures_runtime(monkeypatch,tmp_path,capsys):
    import neo_mri_cli as entry
    monkeypatch.setattr(entry,'dispatch',lambda *a,**k:pytest.fail('Dispatched scientific stage'))
    f=tmp_path/'profile.json';f.write_text(json.dumps(profile()))
    assert entry.main(['plan','--config',str(f)])==0
    assert json.loads(capsys.readouterr().out)['status']=='PLAN_ONLY_NO_RUNTIME_PROBED'


def test_measured_sdc_association_and_content_identity(tmp_path):
    b=tmp_path/'bids';(b/'sub-s/anat').mkdir(parents=True);(b/'sub-s/fmap').mkdir();(b/'sub-s/func').mkdir()
    t='sub-s/anat/sub-s_T1w.nii.gz';f='sub-s/fmap/sub-s_dir-PA_epi.json';m='sub-s/func/sub-s_task-test_bold.json'
    nib.Nifti1Image(np.ones((4,4,4)),np.eye(4)).to_filename(b/t)
    ni=(b/f).with_suffix('.nii.gz');nib.Nifti1Image(np.ones((4,4,4)),np.eye(4)).to_filename(ni)
    (b/m).write_text(json.dumps(dict(PhaseEncodingDirection='j',TotalReadoutTime=.03,B0FieldSource='pair')))
    field=dict(PhaseEncodingDirection='j-',TotalReadoutTime=.03,B0FieldIdentifier='pair')
    (b/f).write_text(json.dumps(field));p=profile(t1w=t,sdc='pepolar',bold_metadata=m,fieldmaps=[f])
    first=s.validate_inputs(tmp_path,p)
    nib.Nifti1Image(np.full((4,4,4),2.),np.eye(4)).to_filename(ni)
    assert s.identity(tmp_path,p)!=first
    field['B0FieldIdentifier']='another_pair';(b/f).write_text(json.dumps(field))
    with pytest.raises(ValueError,match='linked'):s.validate_inputs(tmp_path,p)
    field.update(B0FieldIdentifier='pair',PhaseEncodingDirection='j');(b/f).write_text(json.dumps(field))
    with pytest.raises(ValueError,match='opposite'):s.validate_inputs(tmp_path,p)
