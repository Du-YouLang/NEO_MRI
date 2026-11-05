"""Explicit human gate for native T1 brain masks. Never estimates registration."""
from pathlib import Path
import argparse,base64,datetime,gzip,hashlib,html,json,shutil
import nibabel as nib
import numpy as np
from scipy.ndimage import label

def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def sha(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def save_new(p,obj):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('x',encoding='utf-8') as f:json.dump(obj,f,ensure_ascii=False,indent=2,allow_nan=False)
def file_record(p):
    p=Path(p).resolve();return dict(path=str(p),sha256=sha(p))
def check_file(record):
    p=Path(record['path'])
    if not p.is_file() or sha(p)!=record['sha256']:raise ValueError('Reviewed file missing or changed: '+str(p))
    return p
def checked_mask(t1,mask):
    image=nib.load(t1);seg=nib.load(mask)
    if len(image.shape)!=3 or seg.shape!=image.shape or not np.allclose(image.affine,seg.affine,rtol=0,atol=1e-4):
        raise ValueError('Native T1/mask grid mismatch; never repair this with implicit scaling/resampling')
    data=seg.get_fdata(dtype=np.float32)
    if not np.isfinite(data).all() or not set(np.unique(data))<={0,1} or not data.any() or data.all():
        raise ValueError('Mask must be a finite nonempty binary proper subset')
    return image,seg,data.astype(bool)
def pack(a):return base64.b64encode(gzip.compress(np.ascontiguousarray(a).tobytes(),mtime=0)).decode()

def comparison(config):
    """One native full-head volume, multiple masks, common display window."""
    t1=Path(config['t1']).resolve();image=nib.as_closest_canonical(nib.load(t1))
    data=image.get_fdata(dtype=np.float32)
    if not np.isfinite(data).all():raise ValueError('Nonfinite T1')
    masks=[];records=[];ids=set()
    for item in config['candidates']:
        if item['id'] in ids:raise ValueError('Unique candidate IDs required')
        ids.add(item['id']);_,seg,mask=checked_mask(t1,item['mask']);canonical=nib.as_closest_canonical(seg).get_fdata()>.5
        _,n=label(mask,structure=np.ones((3,3,3)))
        masks.append(dict(id=item['id'],label=item['label'],mask=pack(canonical.astype(np.uint8))))
        records.append(dict(id=item['id'],label=item['label'],mask=file_record(item['mask']),
            volume_ml=float(mask.sum()*abs(np.linalg.det(seg.affine[:3,:3]))/1000),components_26=int(n),human_decision='PENDING',
            computation=item.get('computation',{})))
    if not masks:raise ValueError('At least one candidate required')
    initial=config.get('initial_candidate')
    if initial is None:initial='hd_bet_cpu_no_tta' if 'hd_bet_cpu_no_tta' in ids else masks[0]['id']
    if initial not in ids:raise ValueError('Unknown initial candidate; do not silently select another method')
    reference=nib.as_closest_canonical(nib.load(records[0]['mask']['path'])).get_fdata()>.5
    low,high=np.percentile(data[reference],[1,99])
    if high<=low:raise ValueError('Invalid intensity window')
    scaled=np.clip((data-low)/(high-low)*255,0,255).astype(np.uint8)
    payload=dict(shape=list(image.shape),affine=image.affine.tolist(),image=pack(scaled),masks=masks,
        crosshair=np.round(np.argwhere(reference).mean(axis=0)).astype(int).tolist(),spacing=list(map(float,image.header.get_zooms())))
    dest=Path(config['output']).resolve()
    if dest.exists():raise FileExistsError('New comparison version required')
    dest.mkdir(parents=True)
    manifest=dict(schema='neo-t1-mask-comparison-1',case_id=config['case_id'],t1=file_record(t1),candidates=records,
        human_decision='PENDING',clinical_validation=False,initial_candidate=initial,display_window=[float(low),float(high)],
        note='Full native spatial grid; only display intensity quantized. Candidate disagreement is not accuracy. No automatic winner.')
    save_new(dest/'comparison.json',manifest)
    from anatomy_viewer_assets import SCRIPT
    title=html.escape(config.get('title',config['case_id']))
    options=''.join('<option value="'+str(i)+'"'+(' selected' if m['id']==initial else '')+'>'+html.escape(m['label'])+'</option>' for i,m in enumerate(masks))
    page='''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><style>
body{font:17px/1.8 system-ui,"Microsoft YaHei";margin:20px;background:#edf3f7;color:#173e55}main{max-width:1400px;margin:auto}select,button{font:inherit}#panels{display:flex;flex-wrap:wrap;gap:16px}.panel{flex:1;min-width:310px}canvas{background:black;width:100%;image-rendering:pixelated;cursor:crosshair}input[type=range]{width:98%}.note{padding:18px;background:white;border-left:5px solid #e1a32c}table{border-collapse:collapse}td,th{padding:8px;border:1px solid #aabccc}
</style><main><h1>'''+title+'''：脑提取人工审核</h1><p class="note">逐层看完整头部红线，再看“排除部分”有无连续脑实质；红色表示当前候选掩膜，不是病灶。尚未人工通过，不能继续正式MNI验收。显示以同一原生T1为背景，不自动配准。</p>
<label>候选 <select id="candidate">'''+options+'''</select></label> <label>模式 <select id="mode"><option value="full">完整头部</option><option value="keep">仅保留部分</option><option value="drop">仅排除部分</option><option value="diff">与第一项差异（新增绿/移除紫）</option></select></label>
<label><input id="outline" type="checkbox" checked>红线</label> <label><input id="shade" type="checkbox">红色填充</label> <label>透明度 <input id="alpha" type="range" min="0" max="1" step="0.05" value="0.35"></label><div id="message">本地解压中…</div><div id="panels"></div>
<p>左右为神经学显示：横轴左L→右R；矢状横轴后P→前A。原生斜切仅重排轴，不假装MNI坐标。使用浏览器放大查看薄皮层。审核重点：额/颞极、顶部硬膜/静脉窦、后颅窝/小脑、脑干、病灶旁和M1邻近。</p>
<h2>如何记录决定</h2><p>本页不会因滑过所有切片就自动标记通过。请记录候选、问题所在坐标/索引和意见；确认后用neo-mri-mask decision写入带输入/掩膜SHA的明确决定。需要手修则先导出ITK-SNAP编辑副本，再回读核查；手修结果也必须重新审核。</p><table><tr><th>候选</th><th>掩膜体积mL</th><th>连通块</th></tr>'''
    for r in records:page+=f'<tr><td>{html.escape(r["label"])}</td><td>{r["volume_ml"]:.1f}</td><td>{r["components_26"]}</td></tr>'
    page+='</table><p>体积和连通性不等于准确性；无人工真值，不提供“漏掉0%皮层”。完整输入及模型记录见 comparison.json。</p></main><script id="payload" type="application/json">'+json.dumps(payload,separators=(',',':')).replace('<','\\u003c')+'</script><script>'+SCRIPT+'</script></html>'
    (dest/'review.html').write_text(page,encoding='utf-8');return str(dest/'review.html')

def decision(manifest,candidate,verdict,reviewer,notes,reviewed_all_slices,output,*,user_statement=None):
    m=read(manifest);t1=check_file(m['t1']);r=next((v for v in m['candidates'] if v['id']==candidate),None)
    if r is None:raise ValueError('Unknown candidate')
    mask=check_file(r['mask']);checked_mask(t1,mask)
    if verdict not in ['APPROVED','REJECTED','NEEDS_MANUAL_EDIT']:raise ValueError('Explicit review decision required')
    if not reviewer.strip() or not notes.strip():raise ValueError('Reviewer and actual review notes required')
    if verdict=='APPROVED' and not reviewed_all_slices and not (isinstance(user_statement,str) and user_statement.strip()):raise ValueError('Approval requires explicit human confirmation of full slice review or explicit user approval statement')
    record=dict(schema='neo-t1-mask-review-1',case_id=m['case_id'],t1=m['t1'],mask=r['mask'],candidate=candidate,
        decision=verdict,reviewer=reviewer,notes=notes,reviewed_all_slices=bool(reviewed_all_slices),
        comparison=file_record(manifest),reviewed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
    if user_statement:
        record.update(approval_basis='EXPLICIT_USER_APPROVAL',user_statement=user_statement)
    save_new(output,record);return str(output)

def require_approved(path,expected_t1=None):
    r=read(path)
    explicit=r.get('approval_basis')=='EXPLICIT_USER_APPROVAL' and isinstance(r.get('user_statement'),str) and bool(r['user_statement'].strip())
    if r.get('schema')!='neo-t1-mask-review-1' or r.get('decision')!='APPROVED' or (r.get('reviewed_all_slices') is not True and not explicit) or not r.get('reviewer') or not r.get('notes'):
        raise ValueError('T1 brain-mask human approval missing; stop before normalization')
    t1=check_file(r['t1']);mask=check_file(r['mask']);check_file(r['comparison']);checked_mask(t1,mask)
    if expected_t1 is not None and sha(expected_t1)!=r['t1']['sha256']:raise ValueError('Approval belongs to another T1')
    return r

def manual_export(manifest,candidate,output):
    m=read(manifest);check_file(m['t1']);r=next(v for v in m['candidates'] if v['id']==candidate);source=check_file(r['mask']);checked_mask(m['t1']['path'],source)
    output=Path(output)
    if output.exists():raise FileExistsError('Choose a fresh manual edit directory')
    output.mkdir(parents=True);shutil.copy2(source,output/'mask_to_edit.nii.gz')
    save_new(output/'manual_edit.json',dict(schema='neo-t1-manual-edit-1',t1=m['t1'],parent_mask=r['mask'],comparison=file_record(manifest),
        output_expected='mask_corrected.nii.gz',instructions='ITK-SNAP: load T1 as main image, mask_to_edit as segmentation; paint label1/erase0; save as mask_corrected.nii.gz, never overwrite original candidate.'))
    return str(output/'manual_edit.json')

def manual_import(bundle,edited,output):
    b=read(bundle);check_file(b['t1']);check_file(b['parent_mask']);checked_mask(b['t1']['path'],edited)
    return comparison(dict(case_id=read(check_file(b['comparison']))['case_id'],t1=b['t1']['path'],output=str(output),
        title='手工修正后重新审核',candidates=[dict(id='manual_corrected',label='ITK-SNAP手工修正（待审核）',mask=str(Path(edited).resolve()),computation=dict(parent_mask=b['parent_mask'],manual=True))]))

def export_precomputed(review,raw_t1,output):
    """Export only reviewed preproc+mask; never carry old tissue maps/transforms."""
    r=require_approved(review);checked_mask(raw_t1,r['mask']['path']);raw=Path(raw_t1)
    if not raw.name.endswith('_T1w.nii.gz') or not raw.name.startswith('sub-'):raise ValueError('BIDS T1w filename required')
    # Existing single-input T1 geometry must match the reviewed native image.
    stem=raw.name[:-len('_T1w.nii.gz')];entities=dict(v.split('-',1) for v in stem.split('_') if '-' in v)
    dest=Path(output)
    if dest.exists():raise FileExistsError('Fresh reviewed derivative export required')
    anat=dest/('sub-'+entities['sub'])
    if 'ses' in entities:anat/='ses-'+entities['ses']
    anat/='anat';anat.mkdir(parents=True)
    t1=anat/(stem+'_desc-preproc_T1w.nii.gz');mask=anat/(stem+'_desc-brain_mask.nii.gz')
    shutil.copy2(r['t1']['path'],t1);shutil.copy2(r['mask']['path'],mask)
    save_new(t1.with_name(t1.name.replace('.nii.gz','.json')),dict(SkullStripped=False))
    save_new(mask.with_name(mask.name.replace('.nii.gz','.json')),dict(Type='Brain'))
    save_new(dest/'dataset_description.json',dict(Name='Human-reviewed NEO anatomical inputs',BIDSVersion='1.9.0',DatasetType='derivative',GeneratedBy=[dict(Name='NEO_MRI',Version='0.1.0.dev2')]))
    save_new(dest/'reviewed_anatomy.json',dict(schema='neo-reviewed-anatomy-export-1',review=file_record(review),raw_t1=file_record(raw),t1=file_record(t1),mask=file_record(mask)))
    return str(dest)

def verify_applied_mask(review,derivatives,runtime_evidence=None):
    """Post-run proof: require exact native mask voxels, not only a CLI flag."""
    r=require_approved(review);files=[p for p in Path(derivatives).glob('sub-*/**/*desc-brain_mask.nii.gz') if p.parent.name=='anat' and '_space-' not in p.name]
    if not files and runtime_evidence is not None:
        e=read(runtime_evidence)
        if e.get('schema')!='neo-mask-runtime-evidence-1' or e.get('review_sha256')!=sha(review) or e.get('normalization_input_chain_verified') is not True or e.get('expected_mask_sha256')!=r['mask']['sha256']:
            raise ValueError('Runtime mask evidence not bound to this approval')
        for item in e.get('node_result_files',[]):check_file(item)
        if len(e.get('node_result_files',[]))!=3:raise ValueError('Actual mask/truncation/registration node records required')
        mask=check_file(e['observed_mask']);t1=check_file(e['observed_t1']);brain=check_file(e['masked_brain'])
        _,_,actual=checked_mask(r['t1']['path'],mask);expected=nib.load(r['mask']['path']).get_fdata()>.5
        if not np.array_equal(actual,expected):raise ValueError('Observed runtime mask differs from approved mask')
        source=nib.load(r['t1']['path']);observed=nib.load(t1);b=nib.load(brain)
        if observed.shape!=source.shape or b.shape!=source.shape or not np.allclose(source.affine,observed.affine,rtol=0,atol=1e-4) or not np.allclose(source.affine,b.affine,rtol=0,atol=1e-4):raise ValueError('Runtime T1/brain grid differs')
        x=source.get_fdata(dtype=np.float32)
        if not np.array_equal(x,observed.get_fdata(dtype=np.float32)) or not np.allclose(x*expected,b.get_fdata(dtype=np.float32),rtol=1e-6,atol=1e-5):raise ValueError('Runtime mask application not reproduced')
        return dict(status='APPROVED_MASK_REUSED_VERIFIED',review=file_record(review),approved_mask=r['mask'],observed_mask=e['observed_mask'],
            runtime_evidence=file_record(runtime_evidence),evidence_mode='ACTUAL_NIPYPE_INPUT_CHAIN_NOT_REEXPORTED_NATIVE_DERIVATIVE',identical_native_binary_voxels=True,clinical_validation=False)
    if len(files)!=1:raise ValueError('Exactly one native anatomical mask output required')
    _,actual,data=checked_mask(r['t1']['path'],files[0]);expected=nib.load(r['mask']['path']).get_fdata()>.5
    if not np.array_equal(data,expected):raise ValueError('fMRIPrep native mask differs from approved mask; stop, no implicit replacement')
    return dict(status='APPROVED_MASK_REUSED_VERIFIED',review=file_record(review),approved_mask=r['mask'],output_mask=file_record(files[0]),
        identical_native_binary_voxels=True,clinical_validation=False)

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='action',required=True)
    q=sub.add_parser('compare');q.add_argument('--config',required=True)
    q=sub.add_parser('decision');q.add_argument('--manifest',required=True);q.add_argument('--candidate',required=True);q.add_argument('--verdict',required=True,choices=['APPROVED','REJECTED','NEEDS_MANUAL_EDIT']);q.add_argument('--reviewer',required=True);q.add_argument('--notes',required=True);q.add_argument('--reviewed-all-slices',action='store_true');q.add_argument('--output',required=True)
    q=sub.add_parser('manual-export');q.add_argument('--manifest',required=True);q.add_argument('--candidate',required=True);q.add_argument('--output',required=True)
    q=sub.add_parser('manual-import');q.add_argument('--bundle',required=True);q.add_argument('--edited',required=True);q.add_argument('--output',required=True)
    q=sub.add_parser('export-reviewed');q.add_argument('--review',required=True);q.add_argument('--raw-t1',required=True);q.add_argument('--output',required=True)
    a=p.parse_args(argv)
    if a.action=='compare':result=comparison(read(a.config))
    elif a.action=='decision':result=decision(a.manifest,a.candidate,a.verdict,a.reviewer,a.notes,a.reviewed_all_slices,a.output)
    elif a.action=='manual-export':result=manual_export(a.manifest,a.candidate,a.output)
    elif a.action=='manual-import':result=manual_import(a.bundle,a.edited,a.output)
    else:result=export_precomputed(a.review,a.raw_t1,a.output)
    print(result);return 0
