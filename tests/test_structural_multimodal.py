from pathlib import Path
import sys,json
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'code'))
s=pytest.importorskip('SimpleITK')
from structural_multimodal import compose_pull,resample,n4,check_cache,rigid
from anatomy_review import sha

def test_pull_composition_is_ordered_in_physical_coordinates():
    a=s.Euler3DTransform();a.SetRotation(.1,-.2,.15);a.SetTranslation([2,3,1])
    b=s.AffineTransform(3);b.Scale([1.1,.9,1.05]);b.Translate([-3,1,5])
    c=compose_pull(a,b)
    for p in [(0,0,0),(10,-5,17),(-20,31,2)]:
        assert np.allclose(c.TransformPoint(p),a.TransformPoint(b.TransformPoint(p)))
    assert not np.allclose(c.TransformPoint((10,2,5)),b.TransformPoint(a.TransformPoint((10,2,5))))

def test_one_interpolation_preserves_signal_that_double_interpolation_blurs():
    a=np.zeros((24,24,24),np.float32);a[5:19,5:19,::2]=1
    im=s.GetImageFromArray(a);im.SetSpacing([.8,1.1,1.5]);im.SetOrigin([-4,2,1])
    rot=np.array([[0,-1,0],[1,0,0],[0,0,1]],float);im.SetDirection(rot.ravel())
    # Two half-voxel shifts compose to an exact integer voxel translation.
    step=rot@np.array([.4,0,0]);t=s.TranslationTransform(3,step.tolist())
    once=resample(im,im,compose_pull(t,t));twice=resample(resample(im,im,t),im,t)
    x=s.GetArrayFromImage(once)[7:17,7:17,3:19];y=s.GetArrayFromImage(twice)[7:17,7:17,3:19]
    assert np.var(x)>.2 and np.var(y)<.01
    label=s.Cast(im>.5,s.sitkUInt8)
    assert set(np.unique(s.GetArrayFromImage(resample(label,label,compose_pull(t,t),True))))<={0,1}

def test_cache_never_accepts_filename_alone_or_changed_outputs(tmp_path):
    request=dict(output=str(tmp_path/'v1'),fingerprint='first');out=Path(request['output'])
    assert check_cache(request) is None
    out.mkdir()
    with pytest.raises(FileExistsError,match='Incomplete'):check_cache(request)
    f=out/'map.txt';f.write_text('original')
    r=dict(status='COMPUTED_REVIEW_PENDING',fingerprint='first',outputs=[dict(path=f.name,sha256=sha(f))])
    (out/'receipt.json').write_text(json.dumps(r));assert check_cache(request)==r
    with pytest.raises(ValueError,match='different input'):check_cache(dict(request,fingerprint='other-T1-or-transform'))
    f.write_text('tampered')
    with pytest.raises(ValueError,match='changed'):check_cache(request)

def test_n4_reduces_smooth_synthetic_bias_without_moving_voxels():
    g=np.indices((32,32,32));mask=((g-15.5)**2).sum(axis=0)<13**2
    bias=np.exp(.6*(g[2]-15.5)/32);a=(100*bias).astype('float32')
    im=s.GetImageFromArray(a);im.SetSpacing([2,2,2]);m=s.GetImageFromArray(mask.astype('uint8'));m.CopyInformation(im)
    corrected,field=n4(im,m);b=s.GetArrayFromImage(corrected);f=s.GetArrayFromImage(field)
    assert np.std(b[mask])/np.mean(b[mask]) < np.std(a[mask])/np.mean(a[mask])*.6
    assert corrected.GetSize()==im.GetSize() and corrected.GetOrigin()==im.GetOrigin()
    assert np.all(f>0) and np.allclose(b*f,a,rtol=1e-5)


def test_n4_zero_exclusion_is_explicit_bounded_and_does_not_offset_data():
    from structural_multimodal import n4_mask_audit
    g=np.indices((24,24,24));m=((g-11.5)**2).sum(axis=0)<10**2
    a=np.full(m.shape,100,dtype='float32');points=np.argwhere(m);p=tuple(points[10]);a[p]=0
    im=s.GetImageFromArray(a);mask=s.GetImageFromArray(m.astype('uint8'))
    assert n4_mask_audit(im,mask)['zero_voxels']==1
    with pytest.raises(ValueError,match='zero'):n4(im,mask)
    corrected,field=n4(im,mask,'exclude_zero_up_to_1pct')
    assert s.GetArrayFromImage(corrected)[p]==0
    for point in points[:100]:a[tuple(point)]=0
    im=s.GetImageFromArray(a)
    with pytest.raises(ValueError,match='at most 1%'):n4(im,mask,'exclude_zero_up_to_1pct')

def test_rigid_recovers_synthetic_translation():
    z,y,x=np.indices((40,40,40));a=(np.exp(-((x-15)**2+(y-21)**2+(z-17)**2)/30)*100+
        np.exp(-((x-26)**2+(y-13)**2+(z-25)**2)/15)*70).astype('float32')
    fixed=s.GetImageFromArray(a);known=s.TranslationTransform(3,[2,-1,1])
    moving=resample(fixed,fixed,known.GetInverse());mask=s.Cast(fixed>1,s.sitkUInt8)
    transform,info=rigid(fixed,moving,mask)
    points=[(15.,21.,17.),(26.,13.,25.)]
    assert max(np.linalg.norm(np.array(transform.TransformPoint(p))-known.TransformPoint(p)) for p in points)<1.
    assert info['valid_points']>0


def test_plan_and_small_pipeline_bind_approvals_inputs_and_inverse(tmp_path):
    import nibabel as nib
    from anatomy_review import comparison,decision,file_record
    from structural_multimodal import plan,run
    g=np.indices((24,24,24));m=((g-11.5)**2).sum(axis=0)<10**2
    a=(30+g[0]*2+g[1]+np.exp(-((g-10)**2).sum(axis=0)/15)*50).astype('float32')
    t1=tmp_path/'t1.nii.gz';mask=tmp_path/'mask.nii.gz';source=tmp_path/'t2.nii.gz'
    for p,x in [(t1,a),(source,a),(mask,m.astype('uint8'))]:nib.Nifti1Image(x,np.eye(4)).to_filename(p)
    comparison(dict(case_id='SYNTHETIC',t1=str(t1),output=str(tmp_path/'native'),candidates=[dict(id='test',label='Synthetic',mask=str(mask))]))
    review=tmp_path/'native_approval.json';decision(tmp_path/'native/comparison.json','test','APPROVED','TEST_ONLY','Synthetic fixture',False,review,user_statement='Synthetic authorization only')
    xf=tmp_path/'forward.h5';xi=tmp_path/'inverse.h5'
    for f in [xf,xi]:s.WriteTransform(s.Transform(3,s.sitkIdentity),str(f))
    prov=tmp_path/'provenance.json';prov.write_text(json.dumps(dict(t1_sha256=sha(t1),mask_sha256=sha(mask),normalized_t1_sha256=sha(t1),transforms=[file_record(xf),file_record(xi)])))
    approval=tmp_path/'mni_approval.json';approval.write_text(json.dumps(dict(decision='APPROVED',approval_basis='EXPLICIT_USER_APPROVAL',user_statement='Synthetic only',branches=[dict(decision='APPROVED',normalization_provenance=file_record(prov))])))
    c=dict(modality='T2w',source=str(source),t1=str(t1),mni=str(t1),anatomy_review=str(review),normalization_provenance=str(prov),mni_review=str(approval),mni_to_t1_pull=str(xf),t1_to_mni_pull=str(xi),output=str(tmp_path/'result'),execution_enabled=True)
    first=plan(c)['fingerprint'];r=run(c)
    assert r['status']=='COMPUTED_REVIEW_PENDING' and len(r['branches'])==2
    assert run(c)==r  # No new fitting when the receipt and all outputs match.
    changed=a.copy();changed[2,2,2]+=1;nib.Nifti1Image(changed,np.eye(4)).to_filename(source)
    assert plan(c)['fingerprint']!=first
    with pytest.raises(ValueError,match='different input'):run(c)
    affine=np.eye(4);affine[0,3]=5;im=nib.load(source);im.set_qform(affine,1);im.set_sform(np.eye(4),1);im.to_filename(source)
    with pytest.raises(ValueError,match='qform/sform'):plan(c)


def test_structural_viewer_refuses_implicit_resampling(tmp_path):
    import nibabel as nib
    from multimodal_review import export
    a=np.arange(12**3,dtype='float32').reshape(12,12,12)+1
    p=tmp_path/'a.nii.gz';q=tmp_path/'b.nii.gz';nib.Nifti1Image(a,np.eye(4)).to_filename(p)
    aff=np.eye(4);aff[1,3]=1;nib.Nifti1Image(a,aff).to_filename(q)
    with pytest.raises(ValueError,match='same grid'):export([dict(name='Synthetic',background=str(p),frames=[dict(name='bad',image=str(q))])],tmp_path/'bad.html','Test')
