"""Versioned T2/FLAIR experiment: native N4, rigid MI and explicit pull transforms.

SimpleITK is an optional, isolated dependency. Never starts WSL or services.
The approved T1/MNI pair is reused; no new T1 normalization or lesion segmentation.
"""
from pathlib import Path
import argparse,datetime,hashlib,json,os,sys,time
import numpy as np
from anatomy_review import read,sha,file_record,save_new,require_approved,checked_mask

PARAMETERS=dict(threads=4,registration='SimpleITK Mattes MI rigid, fixed T1 -> moving source LPS pull',
    bins=50,sampling=.2,seed=20261104,iterations=200,shrink_factors=[4,2,1],
    registration_smoothing_mm=[2,1,0],N4_iterations=[40,30,20],N4_shrink=3,
    interpolation='linear intensities / nearest-neighbor labels',
    no_final_image_smoothing=True,automatic_approval=False)

def sitk_module():
    import SimpleITK as sitk
    sitk.ProcessObject.SetGlobalDefaultNumberOfThreads(4)
    return sitk

def plan(c):
    if c.get('modality') not in ['T2w','FLAIR']:raise ValueError('Explicit T2w or FLAIR modality required')
    if set(c.get('parameters',{}))-{'threads','N4_zero_policy'} or c.get('parameters',{}).get('threads',4)!=4:
        raise ValueError('Unreviewed parameter override; new experiment protocol required')
    zero_policy=c.get('parameters',{}).get('N4_zero_policy','stop')
    if zero_policy not in ['stop','exclude_zero_up_to_1pct']:raise ValueError('Explicit supported N4 zero policy required')
    review=require_approved(c['anatomy_review'],c['t1'])
    prov=read(c['normalization_provenance']);approval=read(c['mni_review'])
    if approval.get('decision')!='APPROVED' or approval.get('approval_basis')!='EXPLICIT_USER_APPROVAL' or not approval.get('user_statement') or not any(
        b.get('decision')=='APPROVED' and b['normalization_provenance']['sha256']==sha(c['normalization_provenance'])
        for b in approval.get('branches',[])):
        raise ValueError('Explicit approval of this normalization required')
    if prov['t1_sha256']!=sha(c['t1']) or prov['mask_sha256']!=review['mask']['sha256'] or prov['normalized_t1_sha256']!=sha(c['mni']):
        raise ValueError('Normalization does not match selected approved T1/mask/MNI')
    if prov['transforms'][0]['sha256']!=sha(c['mni_to_t1_pull']):raise ValueError('Image-forward pull transform identity mismatch')
    source=Path(c['source']).resolve();out=Path(c['output']).resolve()
    if out==source.parent or out in source.parents:raise ValueError('Separate result directory required')
    import nibabel as nib
    image=nib.load(source)
    if len(image.shape)!=3:raise ValueError('Only explicit 3D scalar inputs; split identified Dixon/echo components first')
    if not np.isfinite(image.affine).all() or abs(np.linalg.det(image.affine[:3,:3]))<1e-8:raise ValueError('Invalid source physical geometry')
    q,qcode=image.get_qform(coded=True);s,scode=image.get_sform(coded=True)
    if qcode and scode:
        corners=np.array([[i,j,k,1] for i in [0,image.shape[0]-1] for j in [0,image.shape[1]-1] for k in [0,image.shape[2]-1]])
        if np.max(np.linalg.norm((corners@(q-s).T)[:,:3],axis=1))>.1:raise ValueError('Conflicting qform/sform; manual geometry review required')
    for key in ['source_mask','lesion']:
        if c.get(key):checked_mask(source,c[key])
    keys=['source','t1','mni','mni_to_t1_pull','anatomy_review','normalization_provenance','mni_review']
    if c.get('t1_to_mni_pull'):
        if prov['transforms'][1]['sha256']!=sha(c['t1_to_mni_pull']):raise ValueError('Saved inverse normalization identity mismatch')
        keys.append('t1_to_mni_pull')
    elif c.get('lesion'):raise ValueError('Saved inverse normalization required for lesion roundtrip')
    keys += [k for k in ['source_mask','lesion'] if c.get(k)]
    files={k:file_record(c[k]) for k in keys};files['t1_mask']=review['mask']
    identity=dict(schema='neo-multimodal-structural-1',modality=c['modality'],files=files,
        parameters=dict(PARAMETERS,N4_zero_policy=zero_policy),code_sha256=sha(__file__),runtime=dict(SimpleITK=sitk_module().Version_VersionString(),numpy=np.__version__,python=sys.version))
    fingerprint=hashlib.sha256(json.dumps(identity,sort_keys=True).encode()).hexdigest()
    return dict(identity=identity,fingerprint=fingerprint,output=str(out),clinical_validation=False)

def resample(source,reference,pull,label=False):
    s=sitk_module()
    return s.Resample(source,reference,pull,s.sitkNearestNeighbor if label else s.sitkLinear,0,s.sitkUInt8 if label else s.sitkFloat32)

def compose_pull(fixed_to_source,mni_to_fixed):
    """Source sample = fixed_to_source(mni_to_fixed(output MNI point))."""
    s=sitk_module();t=s.CompositeTransform(3)
    t.AddTransform(fixed_to_source);t.AddTransform(mni_to_fixed)
    t.FlattenTransform()
    return t

def rigid(fixed,moving,fixed_mask,moving_mask=None,initial=None):
    s=sitk_module();r=s.ImageRegistrationMethod()
    r.SetMetricAsMattesMutualInformation(PARAMETERS['bins'])
    r.SetMetricSamplingStrategy(r.RANDOM);r.SetMetricSamplingPercentage(PARAMETERS['sampling'],PARAMETERS['seed'])
    r.SetMetricFixedMask(fixed_mask)
    if moving_mask is not None:r.SetMetricMovingMask(moving_mask)
    r.SetInterpolator(s.sitkLinear)
    r.SetOptimizerAsRegularStepGradientDescent(learningRate=1.,minStep=.001,
        numberOfIterations=PARAMETERS['iterations'],relaxationFactor=.5,gradientMagnitudeTolerance=1e-6)
    r.SetOptimizerScalesFromPhysicalShift()
    r.SetShrinkFactorsPerLevel(PARAMETERS['shrink_factors']);r.SetSmoothingSigmasPerLevel(PARAMETERS['registration_smoothing_mm'])
    r.SmoothingSigmasAreSpecifiedInPhysicalUnitsOn()
    start=initial if initial is not None else s.CenteredTransformInitializer(fixed,moving,s.Euler3DTransform(),s.CenteredTransformInitializerFilter.GEOMETRY)
    r.SetInitialTransform(start,inPlace=False)
    transform=r.Execute(s.Cast(fixed,s.sitkFloat32),s.Cast(moving,s.sitkFloat32))
    if not np.isfinite(transform.GetParameters()).all():raise ValueError('Nonfinite rigid transform')
    if hasattr(transform,'FlattenTransform'):transform.FlattenTransform()
    return transform,dict(metric=float(r.GetMetricValue()),stop=r.GetOptimizerStopConditionDescription(),
        valid_points=int(r.GetMetricNumberOfValidPoints()),parameters=list(transform.GetParameters()))

def n4_mask_audit(source,mask):
    s=sitk_module();a=s.GetArrayFromImage(source);m=s.GetArrayFromImage(mask)>0;n=int(m.sum())
    return dict(ROI_voxels=n,positive_voxels=int((a[m]>0).sum()),zero_voxels=int((a[m]==0).sum()),negative_voxels=int((a[m]<0).sum()),
        zero_fraction=float((a[m]==0).sum()/n) if n else None,
        note='Counts in actual preliminary T1-propagated N4 ROI, not independent HD-BET mask')

def n4(source,mask,zero_policy='stop'):
    s=sitk_module();a=s.GetArrayFromImage(source);m=s.GetArrayFromImage(mask)>0
    audit=n4_mask_audit(source,mask)
    if audit['positive_voxels']<1000 or audit['negative_voxels']:raise ValueError('N4 ROI has insufficient positive voxels or negative values; inspect saved audit')
    if audit['zero_voxels']:
        if zero_policy!='exclude_zero_up_to_1pct' or audit['zero_fraction']>.01:raise ValueError('N4 ROI contains zero values; inspect saved audit, explicit zero policy required (at most 1%)')
        effective=s.GetImageFromArray((m&(a>0)).astype('uint8'));effective.CopyInformation(mask);mask=effective
    filt=s.N4BiasFieldCorrectionImageFilter();filt.SetMaximumNumberOfIterations(PARAMETERS['N4_iterations'])
    filt.Execute(s.Shrink(source,[PARAMETERS['N4_shrink']]*3),s.Shrink(mask,[PARAMETERS['N4_shrink']]*3))
    field=s.Exp(filt.GetLogBiasFieldAsImage(source));corrected=source/field
    return corrected,field

def check_cache(request):
    out=Path(request['output']);receipt=out/'receipt.json'
    if not out.exists():return None
    if not receipt.exists():raise FileExistsError('Incomplete output protected; no automatic rerun')
    r=read(receipt)
    if r.get('status')!='COMPUTED_REVIEW_PENDING' or r['fingerprint']!=request['fingerprint']:
        raise ValueError('Existing output belongs to different input/parameters/code or failed run; use a new version')
    for f in r['outputs']:
        if sha(out/f['path'])!=f['sha256']:raise ValueError('Cached output changed: '+f['path'])
    return r

def run(c):
    request=plan(c);cached=check_cache(request)
    if cached:return cached
    if not c.get('execution_enabled'):raise ValueError('Explicit execution_enabled required')
    s=sitk_module();out=Path(request['output']);out.mkdir(parents=True)
    save_new(out/'request.json',request);start=time.monotonic()
    try:
        source=s.ReadImage(c['source'],s.sitkFloat32);t1=s.ReadImage(c['t1'],s.sitkFloat32);mni=s.ReadImage(c['mni'],s.sitkFloat32)
        t1_mask=s.ReadImage(request['identity']['files']['t1_mask']['path'],s.sitkUInt8)
        for im in [source,t1,mni]:
            if im.GetDimension()!=3 or not np.isfinite(s.GetArrayViewFromImage(im)).all():raise ValueError('Finite 3D images required')
        initial,initial_info=rigid(t1,source,t1_mask)
        source_roi=resample(t1_mask,source,initial.GetInverse(),True)
        s.WriteTransform(initial,str(out/'initial_T1_to_source_pull.h5'))
        save_new(out/'initial_registration.json',initial_info)
        s.WriteImage(source_roi,str(out/'T1_brain_in_source.nii.gz'),True)
        save_new(out/'N4_mask_audit.json',n4_mask_audit(source,source_roi))
        corrected,field=n4(source,source_roi,request['identity']['parameters']['N4_zero_policy'])
        s.WriteImage(corrected,str(out/'N4.nii.gz'),True);s.WriteImage(field,str(out/'N4_bias_field.nii.gz'),True)
        source_mask=s.ReadImage(c['source_mask'],s.sitkUInt8) if c.get('source_mask') else None
        mni_pull=s.ReadTransform(c['mni_to_t1_pull']);rows=[]
        lesion=s.ReadImage(c['lesion'],s.sitkUInt8) if c.get('lesion') else None
        mask_domain=s.GetArrayFromImage(t1_mask)>0
        for correction,image in [('raw',source),('n4',corrected)]:
            for mask_mode in (['T1_only','independent_HDBET'] if source_mask is not None else ['T1_only']):
                name=correction+'_'+mask_mode
                transform,info=(initial,initial_info) if name=='raw_T1_only' else rigid(t1,image,t1_mask,source_mask if mask_mode=='independent_HDBET' else None,initial)
                s.WriteTransform(transform,str(out/(name+'_T1_to_source_pull.h5')))
                composed=compose_pull(transform,mni_pull);s.WriteTransform(composed,str(out/(name+'_MNI_to_source_pull.h5')))
                in_t1=resample(image,t1,transform);once=resample(image,mni,composed);twice=resample(in_t1,mni,mni_pull)
                for suffix,im in [('T1',in_t1),('MNI_once',once),('MNI_twice',twice)]:s.WriteImage(im,str(out/(name+'_'+suffix+'.nii.gz')),True)
                diff=s.GetArrayFromImage(once)-s.GetArrayFromImage(twice)
                row=dict(branch=name,registration=info,single_vs_double_RMSE_all_FOV=float(np.sqrt(np.mean(diff.astype(float)**2))),
                    interpretation='Difference from interpolation only, not accuracy against ground truth')
                if source_mask is not None:
                    moved=resample(source_mask,t1,transform,True);s.WriteImage(moved,str(out/(name+'_HDBET_in_T1.nii.gz')),True)
                    b=s.GetArrayFromImage(moved)>0;row['brain_Dice_proxy']=float(2*(b&mask_domain).sum()/(b.sum()+mask_domain.sum()))
                if lesion is not None:
                    l1=resample(lesion,t1,transform,True);l_once=resample(lesion,mni,composed,True);l_twice=resample(l1,mni,mni_pull,True)
                    back=resample(l_once,lesion,composed.GetInverse() if composed.IsLinear() else compose_inverse(transform,c),True)
                    for suffix,im in [('lesion_T1',l1),('lesion_MNI_once',l_once),('lesion_MNI_twice',l_twice)]:s.WriteImage(im,str(out/(name+'_'+suffix+'.nii.gz')),True)
                    a=s.GetArrayFromImage(lesion)>0;b=s.GetArrayFromImage(back)>0;x=s.GetArrayFromImage(l_once)>0;y=s.GetArrayFromImage(l_twice)>0
                    row.update(lesion_roundtrip_Dice=float(2*(a&b).sum()/(a.sum()+b.sum())),
                        lesion_native_ml=float(a.sum()*np.prod(lesion.GetSpacing())/1000),
                        lesion_MNI_once_ml=float(x.sum()*np.prod(mni.GetSpacing())/1000),lesion_MNI_twice_ml=float(y.sum()*np.prod(mni.GetSpacing())/1000),
                        lesion_single_double_disagreement_voxels=int(np.count_nonzero(x!=y)),roundtrip_is_not_anatomical_accuracy=True)
                    if source_mask is not None:row['native_lesion_retained_by_HDBET_fraction']=float((a&(s.GetArrayFromImage(source_mask)>0)).sum()/a.sum())
                rows.append(row)
        record=dict(status='COMPUTED_REVIEW_PENDING',fingerprint=request['fingerprint'],seconds=time.monotonic()-start,
            branches=rows,clinical_validation=False,SimpleITK=s.Version_VersionString(),numpy=np.__version__,python=sys.version,
            note='Registration uses internal multiresolution smoothing; no smoothing of final images. N4 mask uses common preliminary rigid T1 ROI; no lesion exclusion. Independent HD-BET fixed across raw/N4.')
        for f in request['identity']['files'].values():
            if sha(f['path'])!=f['sha256']:raise ValueError('Input changed during calculation; stop acceptance')
        record['outputs']=[dict(path=p.relative_to(out).as_posix(),sha256=sha(p)) for p in sorted(out.rglob('*')) if p.is_file()]
        save_new(out/'receipt.json',record);return record
    except Exception as exc:
        save_new(out/'failure.json',dict(status='FAILED_NO_AUTOMATIC_RETRY',error=str(exc),seconds=time.monotonic()-start));raise

def compose_inverse(t1_to_source,c):
    s=sitk_module()
    if not c.get('t1_to_mni_pull'):raise ValueError('Saved inverse normalization required for lesion roundtrip')
    t=s.CompositeTransform(3);t.AddTransform(s.ReadTransform(c['t1_to_mni_pull']));t.AddTransform(t1_to_source.GetInverse());t.FlattenTransform();return t

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=['plan','run']);p.add_argument('--config',required=True);a=p.parse_args(argv)
    result=plan(read(a.config)) if a.action=='plan' else run(read(a.config));print(json.dumps(result,ensure_ascii=False,indent=2));return 0

if __name__=='__main__':sys.exit(main())
