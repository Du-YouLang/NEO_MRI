"""Explicit CFM/SDC policy, pure planning and input checks. Never starts a runtime."""
from pathlib import Path
import hashlib,json

def key(value):return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def validate(p):
    allowed={'schema','execution_enabled','normalization','sdc','t1w','lesion_roi','bold_metadata','fieldmaps','input_ids','stc'}
    if set(p)-allowed:raise ValueError('Unknown scientific protocol keys: '+str(set(p)-allowed))
    if p.get('schema')!='neo-mri-protocol-1':raise ValueError('Unknown scientific protocol schema')
    if not isinstance(p.get('execution_enabled'),bool):raise ValueError('Explicit execution_enabled boolean required')
    if p.get('normalization') not in ['standard','cfm']:raise ValueError('Select standard or cfm')
    if p.get('sdc') not in ['off','pepolar','fieldmap','syn']:raise ValueError('Unknown SDC method')
    if p.get('stc','off')!='off':raise ValueError('STC changes are deferred in this revision')
    if p['normalization']=='cfm' and not p.get('lesion_roi'):raise ValueError('CFM requires a T1-grid lesion ROI')
    if p['normalization']=='cfm' and not p.get('t1w'):raise ValueError('CFM requires selected T1 input')
    if p['normalization']=='standard' and p.get('lesion_roi'):raise ValueError('Standard branch must not ingest a lesion ROI')
    if p['sdc'] in ['pepolar','fieldmap'] and not p.get('fieldmaps'):raise ValueError('Requested measured SDC requires linked fieldmap inputs')
    if p['sdc'] in ['off','syn'] and p.get('fieldmaps'):raise ValueError('Off/SyN branches must not declare measured fieldmaps')
    if p['sdc']!='off' and not p.get('bold_metadata'):raise ValueError('SDC requires BOLD metadata')
    return p

def arguments(p):
    """Only arguments; real inputs must separately pass validate_inputs."""
    if p is None:return ['--ignore','slicetiming','fieldmaps','flair']
    validate(p);args=['--ignore','slicetiming','flair']
    if p['sdc'] in ['off','syn']:args+=['fieldmaps']
    if p['sdc']=='syn':args+=['--use-syn-sdc','error']
    args+=['--bold2anat-init','t1w']
    # CFM is a BIDS lesion _roi input, NOT an invented --lesion-mask option.
    return args

def inputs(p):
    return [p[k] for k in ['t1w','lesion_roi','bold_metadata'] if p.get(k)]+list(p.get('fieldmaps',[]))

def identity(root,p):
    validate(p);root=Path(root).resolve();rows=[]
    paths=inputs(p)+[str(Path(x).with_suffix('.nii.gz')).replace('\\','/') for x in p.get('fieldmaps',[])]
    for name in paths:
        f=(root/'bids'/name).resolve()
        if not f.is_relative_to(root/'bids'):raise ValueError('Protocol input outside BIDS')
        with f.open('rb') as stream:h=hashlib.file_digest(stream,'sha256').hexdigest()
        rows.append(dict(path=name,sha256=h))
    return dict(protocol=p,inputs=rows)

def validate_inputs(root,p):
    """Explicit preprocessing only; planning does not call this or open images."""
    validate(p);root=Path(root);bids=root/'bids';identity(root,p)
    if not p.get('t1w') or not p['t1w'].endswith('_T1w.nii.gz'):raise ValueError('Selected input T1 required')
    t1_inputs=list(bids.glob('sub-*/**/*_T1w.nii.gz'))
    if len(t1_inputs)!=1 or t1_inputs[0].resolve()!=(bids/p['t1w']).resolve():raise ValueError('New branch BIDS must contain exactly the selected T1')
    rois=list(bids.glob('sub-*/**/*_roi.nii.gz'))
    if p['normalization']=='standard' and rois:raise ValueError('Standard BIDS contains ROI; split branch inputs')
    if p['normalization']=='cfm':
        import nibabel as nib
        import numpy as np
        roi=bids/p['lesion_roi']
        if not roi.name.endswith('_roi.nii.gz') or len(rois)!=1 or rois[0].resolve()!=roi.resolve():raise ValueError('Exactly one declared _roi required for this case branch')
        t1=nib.load(bids/p['t1w']);mask=nib.load(roi);a=mask.get_fdata()
        if mask.shape!=t1.shape or not np.allclose(mask.affine,t1.affine,rtol=0,atol=1e-5):raise ValueError('CFM mask must match selected input T1 grid')
        if not np.isfinite(a).all() or not set(np.unique(a))<={0,1} or not a.any() or a.all():raise ValueError('CFM lesion must be a nonempty binary proper subset')
        ignored=(bids/'.bidsignore').read_text() if (bids/'.bidsignore').exists() else ''
        if '*lesion_roi.nii.gz' not in ignored:raise ValueError('Document lesion ROI in .bidsignore')
    if p['sdc']!='off':
        meta=json.loads((bids/p['bold_metadata']).read_text())
        if meta.get('PhaseEncodingDirection') not in ['i','i-','j','j-','k','k-']:raise ValueError('Missing/invalid phase encoding')
        if p['sdc'] in ['pepolar','fieldmap']:
            if float(meta.get('TotalReadoutTime',0))<=0:raise ValueError('Measured SDC requires verified readout time')
            if not meta.get('B0FieldSource'):raise ValueError('Explicit B0FieldSource association required')
            sources=meta['B0FieldSource'];sources=[sources] if isinstance(sources,str) else sources
            field=[]
            for name in p['fieldmaps']:
                if not name.endswith('.json'):raise ValueError('fieldmaps lists JSON metadata with accompanying NIfTI')
                path=bids/name
                if not path.with_suffix('.nii.gz').is_file():raise ValueError('Fieldmap NIfTI missing')
                m=json.loads(path.read_text());ids=m.get('B0FieldIdentifier',[]);ids=[ids] if isinstance(ids,str) else ids
                if not set(ids)&set(sources):raise ValueError('Fieldmap not linked to requested BOLD')
                field.append((path,m))
            if p['sdc']=='pepolar':
                pe=meta['PhaseEncodingDirection'];opposite=pe[:-1] if pe.endswith('-') else pe+'-'
                if not any(x.name.endswith('_epi.json') and m.get('PhaseEncodingDirection')==opposite and float(m.get('TotalReadoutTime',0))>0 for x,m in field):raise ValueError('PEPOLAR requires associated opposite-encoding EPI and readout metadata')
            elif not any(x.name.endswith(('_phasediff.json','_fieldmap.json','_phase1.json')) for x,m in field):raise ValueError('No recognized measured fieldmap metadata')
    return identity(root,p)

def plan(p):
    validate(p);ids=p.get('input_ids',{});nodes={}
    def add(name,value):nodes[name]=key(value);return nodes[name]
    t1=add('t1_reference',ids.get('t1','UNBOUND'))
    norm=add('normalization',dict(t1=t1,method=p['normalization'],lesion=ids.get('lesion','UNBOUND') if p['normalization']=='cfm' else None))
    sdc=add('functional_geometry',dict(bold=ids.get('bold','UNBOUND'),t1=t1,sdc=p['sdc'],fieldmaps=ids.get('fieldmaps','UNBOUND') if p['sdc'] in ['pepolar','fieldmap'] else None))
    add('native_bold_confounds',dict(geometry=sdc,t1=t1))
    add('mni_resampling',dict(geometry=sdc,normalization=norm))
    required_ids={'t1','bold'}|({'lesion'} if p['normalization']=='cfm' else set())|({'fieldmaps'} if p['sdc'] in ['pepolar','fieldmap'] else set())
    return dict(status='PLAN_ONLY_NO_RUNTIME_PROBED',execution_enabled=p['execution_enabled'],arguments=arguments(p),
        dependency_keys=nodes,identities_bound=required_ids<=set(ids) and all(ids[k]!='UNBOUND' for k in required_ids),
        required_inputs=inputs(p),qc=['lesion_roundtrip_Dice','roundtrip_coordinate_error_mm','nonlesion_anatomical_correspondence','deformation_Jacobian','human_review'],
        automatic_winner=False,interpolation='compose physical transforms; resample once when feasible',
        limitation='Declared dependency keys are planning identities, not verified cache receipts or proof that a branch ran.')

def require_review(protocol,fingerprint,review):
    if review.get('status')!='REVIEWED_MATCHES_REQUEST' or review.get('preprocessing_fingerprint')!=fingerprint:
        raise ValueError('New protocol requires input-bound method review before GLM')
    if review.get('observed_normalization')!=protocol['normalization'] or review.get('observed_sdc')!=protocol['sdc']:
        raise ValueError('Observed method differs from requested method; no silent fallback')
    if not review.get('evidence_files') or not review.get('reviewer'):raise ValueError('Method review requires reviewer and evidence files')
    if any(not isinstance(x,dict) or set(x)!={'path','sha256'} for x in review['evidence_files']):raise ValueError('Evidence files require path and SHA256')
