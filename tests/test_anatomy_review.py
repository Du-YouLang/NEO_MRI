from pathlib import Path
import json,sys
import numpy as np,nibabel as nib,pytest
from PIL import Image
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'code'))
from anatomy_review import comparison,decision,require_approved,manual_export,manual_import,export_precomputed,verify_applied_mask,read,sha,save_new,file_record
from registration_blink import render

@pytest.fixture
def anatomy(tmp_path):
    shape=(12,14,16);g=np.indices(shape);data=(g[0]+2*g[1]+3*g[2]+50).astype(np.float32)
    aff=np.diag([1,1,1,1]);aff[:3,3]=[-6,-7,-8]
    t1=tmp_path/'sub-demo_T1w.nii.gz';mask=tmp_path/'mask.nii.gz';m=np.zeros(shape,np.uint8);m[2:10,2:12,2:14]=1
    nib.Nifti1Image(data,aff).to_filename(t1);nib.Nifti1Image(m,aff).to_filename(mask)
    comparison(dict(case_id='synthetic',t1=str(t1),output=str(tmp_path/'comparison'),candidates=[dict(id='test',label='SYNTHETIC',mask=str(mask))]))
    return dict(root=tmp_path,t1=t1,mask=mask,manifest=tmp_path/'comparison/comparison.json',review=tmp_path/'approval.json')

def approve(a):decision(a['manifest'],'test','APPROVED','SYNTHETIC_TEST_ONLY','Synthetic reference only, no patient approved',True,a['review'])

def test_no_automatic_approval_and_explicit_human_confirmation(anatomy):
    a=anatomy;assert read(a['manifest'])['human_decision']=='PENDING'
    with pytest.raises(ValueError,match='confirmation'):decision(a['manifest'],'test','APPROVED','test','notes',False,a['review'])
    decision(a['manifest'],'test','NEEDS_MANUAL_EDIT','test','posterior border',False,a['review'])
    with pytest.raises(ValueError,match='approval'):require_approved(a['review'])

def test_approval_binds_input_mask_and_comparison(anatomy):
    a=anatomy;approve(a);require_approved(a['review'])
    original=a['mask'].read_bytes();a['mask'].write_bytes(original+b'x')
    with pytest.raises(ValueError,match='changed'):require_approved(a['review'])


def test_explicit_user_approval_does_not_claim_full_slice_review(anatomy):
    a=anatomy
    decision(a['manifest'],'test','APPROVED','SYNTHETIC_USER','Explicit authorization fixture',False,a['review'],user_statement='SYNTHETIC explicit approval')
    r=require_approved(a['review'])
    assert r['reviewed_all_slices'] is False
    assert r['approval_basis']=='EXPLICIT_USER_APPROVAL'

def test_manual_roundtrip_preserves_original_and_requires_new_review(anatomy):
    a=anatomy;before=sha(a['mask']);bundle=manual_export(a['manifest'],'test',a['root']/'edit')
    mask=nib.load(a['mask']);m=mask.get_fdata().astype(np.uint8);m[1,3,3]=1
    corrected=a['root']/'edit/mask_corrected.nii.gz';nib.Nifti1Image(m,mask.affine).to_filename(corrected)
    manual_import(bundle,corrected,a['root']/'manual_review')
    assert sha(a['mask'])==before
    assert read(a['root']/'manual_review/comparison.json')['human_decision']=='PENDING'
    moved=mask.affine.copy();moved[0,3]+=1;nib.Nifti1Image(m,moved).to_filename(corrected)
    with pytest.raises(ValueError,match='grid'):manual_import(bundle,corrected,a['root']/'bad')

def test_export_reuses_only_approved_native_t1_and_mask(anatomy):
    a=anatomy;approve(a);d=Path(export_precomputed(a['review'],a['t1'],a['root']/'export'))
    e=read(d/'reviewed_anatomy.json');assert e['mask']['sha256']==sha(a['mask'])
    assert not list(d.rglob('*dseg*')) and not list(d.rglob('*.h5'))

def blink_config(a):
    template=a['t1'];image=nib.load(template);x=image.get_fdata();shift=np.roll(x,2,axis=0)
    moving=a['root']/'normalized.nii.gz';nib.Nifti1Image(shift,image.affine).to_filename(moving)
    transform=a['root']/'synthetic_transform.txt';transform.write_text('SYNTHETIC shift x=2 vox; not actual ANTs')
    prov=a['root']/'normalization.json';save_new(prov,dict(t1_sha256=sha(a['t1']),mask_sha256=sha(a['mask']),normalized_t1_sha256=sha(moving),template_sha256=sha(template),transforms=[file_record(transform)]))
    return dict(anatomy_review=str(a['review']),normalization_provenance=str(prov),normalized_t1=str(moving),template=str(template),space='SYNTHETIC_GRID_NOT_MNI',cut_coords_world_mm=[0,0,0],output=str(a['root']/'blink'),scope='SYNTHETIC_TEST_ONLY')

def test_blink_is_two_frames_and_never_estimates_transform(anatomy):
    a=anatomy;approve(a);c=blink_config(a);render(c)
    im=Image.open(a['root']/'blink/registration_blink.gif');assert im.n_frames==2
    assert read(a['root']/'blink/blink_receipt.json')['registration_fits']==0
    assert 'SYNTHETIC' in (a['root']/'blink/blink.html').read_text(encoding='utf-8')

def test_blink_display_enlarges_small_slices_with_physical_aspect():
    from registration_blink import display_tile
    tile=display_tile(np.zeros((10,20),np.uint8),[1,2])
    assert tile.size==(310,310)
    tile=display_tile(np.zeros((10,20),np.uint8),[1,1])
    assert tile.size==(310,155)


def test_blink_blocked_before_approval_and_on_wrong_source(anatomy):
    a=anatomy;decision(a['manifest'],'test','REJECTED','test','missing cortex',True,a['review']);c=blink_config(a)
    with pytest.raises(ValueError,match='approval'):render(c)

def test_blink_refuses_old_normalization_for_new_mask(anatomy):
    a=anatomy;approve(a);c=blink_config(a);p=read(c['normalization_provenance']);p['mask_sha256']='wrong';Path(c['normalization_provenance']).write_text(json.dumps(p))
    with pytest.raises(ValueError,match='bound'):render(c)

def test_blink_refuses_geometry_mismatch(anatomy):
    a=anatomy;approve(a);c=blink_config(a);p=Path(c['normalized_t1']);im=nib.load(p);aff=im.affine.copy();aff[0,3]+=2;nib.Nifti1Image(im.get_fdata(),aff).to_filename(p)
    prov=read(c['normalization_provenance']);prov['normalized_t1_sha256']=sha(p);Path(c['normalization_provenance']).write_text(json.dumps(prov))
    with pytest.raises(ValueError,match='grid'):render(c)

def test_candidate_command_is_cpu_bounded_and_plan_never_launches(anatomy,monkeypatch):
    import subprocess
    from brain_extract import plan
    a=anatomy;script=a['root']/'mri_synthstrip';script.write_text('SYNTHETIC stub, not executed')
    weights=a['root']/'model.pt';weights.write_text('SYNTHETIC stub, not a model')
    config=dict(method='synthstrip',t1=str(a['t1']),output=str(a['root']/'candidate'),runtime=dict(python=sys.executable,script=str(script),model=str(weights)),timeout_seconds=1800)
    monkeypatch.setattr(subprocess,'Popen',lambda *args,**kwargs:pytest.fail('Plan launched process'))
    r=plan(config);assert r['gpu'] is False and r['threads']==4 and '-g' not in r['command']
    config['timeout_seconds']=1801
    with pytest.raises(ValueError,match='1800'):plan(config)


def test_missing_method_selects_hdbet_without_launching_or_approving(anatomy,monkeypatch):
    import subprocess
    from brain_extract import plan,selection_policy,main
    model=anatomy['root']/'models';model.mkdir();(model/'weights.pth').write_bytes(b'SYNTHETIC NOT A MODEL')
    monkeypatch.setattr(subprocess,'Popen',lambda *a,**k:pytest.fail('Plan launched a process'))
    c=dict(t1=str(anatomy['t1']),output=str(anatomy['root']/'hd'),runtime=dict(python=sys.executable,model_dir=str(model)))
    r=plan(c)
    assert r['method']=='hd_bet_cpu_no_tta' and not r['selection']['explicit_method']
    assert r['gpu'] is False and r['threads']==4 and r['timeout_seconds']==1800
    assert selection_policy()['human_mask_review_required'] and not selection_policy()['automatic_fallback']
    assert main(['policy'])==0
    assert not Path(c['output']).exists()
    assert read(anatomy['manifest'])['human_decision']=='PENDING'


def test_missing_hdbet_model_does_not_fall_back_to_available_synthstrip(anatomy):
    from brain_extract import plan
    c=dict(t1=str(anatomy['t1']),output=str(anatomy['root']/'hd'),
           runtimes=dict(synthstrip=dict(python=sys.executable,script='available',model='available')))
    with pytest.raises(ValueError,match='no automatic fallback'):plan(c)
    c['runtimes']['hd_bet_cpu_no_tta']=dict(python=sys.executable,model_dir=str(anatomy['root']/'absent'))
    with pytest.raises(ValueError,match='model directory'):plan(c)


@pytest.mark.parametrize('method,template_id',[('ants_oasis','OASIS30ANTs'),('ants_mni2009c','MNI152NLin2009cAsym')])
def test_named_ants_choices_preserve_explicit_template_identity(anatomy,method,template_id):
    from brain_extract import plan
    root=anatomy['root'];script=root/'antsBrainExtraction.sh';script.write_text('SYNTHETIC stub')
    runtime=dict(script=str(script),template_id=template_id)
    for key,suffix in [('template','T1w'),('prior','probseg'),('registration_mask','mask')]:
        p=root/('tpl-'+template_id+'_'+suffix+'.nii.gz');p.write_bytes(b'SYNTHETIC resource');runtime[key]=str(p)
    c=dict(method=method,t1=str(anatomy['t1']),output=str(root/'ants'),runtime=runtime)
    r=plan(c)
    assert r['method']==method and r['selection']['explicit_method']
    assert runtime['template'] in r['command'] and r['selection']['automatic_fallback'] is False
    runtime['template_id']='wrong'
    with pytest.raises(ValueError,match='template identity'):plan(c)
    runtime['template_id']=template_id;runtime['prior']=runtime['prior'].replace('tpl-'+template_id,'tpl-Other')
    with pytest.raises(ValueError,match='template identity'):plan(c)


def test_viewer_prefers_hdbet_without_changing_reference_or_approval(anatomy):
    a=anatomy;c=dict(case_id='synthetic',t1=str(a['t1']),output=str(a['root']/'preferred'),
        candidates=[dict(id='baseline',label='Reference',mask=str(a['mask'])),
                    dict(id='hd_bet_cpu_no_tta',label='HD-BET',mask=str(a['mask']))])
    comparison(c);r=read(Path(c['output'])/'comparison.json')
    assert r['initial_candidate']=='hd_bet_cpu_no_tta' and r['candidates'][0]['id']=='baseline'
    assert r['human_decision']=='PENDING' and not r['clinical_validation']
    assert '<option value="1" selected>HD-BET</option>' in (Path(c['output'])/'review.html').read_text(encoding='utf-8')
    c.update(output=str(a['root']/'explicit'),initial_candidate='baseline');comparison(c)
    assert read(Path(c['output'])/'comparison.json')['initial_candidate']=='baseline'
    c.update(output=str(a['root']/'invalid'),initial_candidate='absent')
    with pytest.raises(ValueError,match='Unknown initial'):comparison(c)

def test_postrun_mask_check_rejects_silent_skullstrip_fallback(anatomy):
    import shutil
    a=anatomy;approve(a);der=a['root']/'derivatives';folder=der/'sub-demo/anat';folder.mkdir(parents=True)
    mask=folder/'sub-demo_desc-brain_mask.nii.gz';shutil.copyfile(a['mask'],mask)
    assert verify_applied_mask(a['review'],der)['identical_native_binary_voxels']
    im=nib.load(mask);m=im.get_fdata();m[0,0,0]=1;nib.Nifti1Image(m,im.affine).to_filename(mask)
    with pytest.raises(ValueError,match='differs'):verify_applied_mask(a['review'],der)


def test_reused_inputs_require_observed_node_and_voxel_evidence(anatomy):
    a=anatomy;approve(a);folder=a['root']/'runtime';folder.mkdir()
    image=nib.load(a['t1']);mask=nib.load(a['mask']).get_fdata()>.5
    brain=folder/'brain.nii.gz';nib.Nifti1Image(image.get_fdata()*mask,image.affine).to_filename(brain)
    nodes=[]
    for name in ['apply','trunc','registration']:
        p=folder/(name+'.pklz');p.write_bytes(b'SYNTHETIC node receipt, never executed');nodes.append(file_record(p))
    evidence=folder/'evidence.json'
    e=dict(schema='neo-mask-runtime-evidence-1',review_sha256=sha(a['review']),normalization_input_chain_verified=True,
        expected_mask_sha256=sha(a['mask']),node_result_files=nodes,observed_mask=file_record(a['mask']),observed_t1=file_record(a['t1']),masked_brain=file_record(brain))
    evidence.write_text(json.dumps(e))
    assert verify_applied_mask(a['review'],folder/'no_native_reexport',runtime_evidence=evidence)['identical_native_binary_voxels']
    bad=nib.load(brain).get_fdata();bad[0,0,0]=999;nib.Nifti1Image(bad,image.affine).to_filename(brain)
    e['masked_brain']=file_record(brain);evidence.write_text(json.dumps(e))
    with pytest.raises(ValueError,match='application'):verify_applied_mask(a['review'],folder/'no_native_reexport',runtime_evidence=evidence)
