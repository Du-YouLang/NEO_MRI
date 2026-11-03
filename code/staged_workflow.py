"""Independent preprocessing, statistical and reporting stages.

Only the explicit preprocess stage can reach a container launcher. Completed
source experiments are referenced read-only; all state/results live elsewhere.
"""
from pathlib import Path
from contextlib import contextmanager
import datetime,hashlib,json,os,re,time,sys

REPO=Path(os.environ.get('NEO_MRI_ROOT',Path.cwd())).resolve()
CODE=Path(__file__).resolve().parent
THREADS=4

class StageError(RuntimeError):pass
class TimingBlocked(StageError):pass

def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def read(path):return json.loads(Path(path).read_text(encoding='utf-8'))
def digest(obj):return hashlib.sha256(json.dumps(obj,sort_keys=True,separators=(',',':'),ensure_ascii=False,allow_nan=False).encode()).hexdigest()
def sha(path):
    with Path(path).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def write(path,obj):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_name(path.name+'.tmp')
    tmp.write_text(json.dumps(obj,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8');tmp.replace(path)
def logical_resolve(path):
    p=Path(path);return p.resolve() if p.is_absolute() else (REPO/p).resolve()

def resolve(path):
    """Relative data paths use the explicit --root, never a sibling repository."""
    return logical_resolve(path)
def configure_threads():
    global _THREAD_LIMITER
    for name in ['OMP_NUM_THREADS','ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:
        os.environ[name]=str(THREADS)
    # Also constrain libraries that an interactive caller imported beforehand.
    from threadpoolctl import threadpool_limits
    _THREAD_LIMITER=threadpool_limits(limits=THREADS)
def load(path):
    path=resolve(path)
    c=read(path)
    if c.get('template_only'):raise StageError('Copy this template, fill case/BIDS/runtime/timing fields, then remove template_only; it is not a runnable analysis')
    if c.get('schema')!='neo-bold-stages-1':raise StageError('Unsupported staged workflow configuration')
    if not re.fullmatch(r'[A-Za-z0-9_-]+',c['case_id']):raise StageError('Invalid case_id')
    state=resolve(c['workspace']);source=resolve(c['preprocessing']['project'])
    if state==source or source in state.parents or state in source.parents:
        raise StageError('Stage workspace must be separate from the preserved experiment')
    original_state=logical_resolve(c['workspace'])
    approved=any(original_state.is_relative_to(resolve(p)) for p in c.get('output_roots',[]))
    if not original_state.is_relative_to(REPO/'results') and not approved:
        raise StageError('Stage workspace must be under repository results/ or an explicitly approved individual result root')
    ids=[r['id'] for r in c['analysis']['runs']]
    if len(ids)!=len(set(ids)) or any(not re.fullmatch(r'[A-Za-z0-9_-]+',x) for x in ids):raise StageError('Run IDs must be unique safe names')
    allowed={'hrf_model','smoothing_fwhm_mm','noise_model','compcor_components','fd_threshold_mm','dvars_threshold','contrast','fdr_q'}
    if set(c['analysis']['statistics'])!=allowed:raise StageError('Unknown or missing statistical parameter; no silently ignored options')
    if c['analysis'].get('engine_variant','baseline') not in ['baseline','anatomical_intersection_v1']:
        raise StageError('Unknown GLM engine variant')
    c['_config_path']=str(Path(path).resolve());return c
def workspace(c):return resolve(c['workspace'])
def source_root(c):return resolve(c['preprocessing']['project'])
def derivatives(c):return source_root(c)/'derivatives/fmriprep'
def preprocessing_request(c):
    # Path relocation is engineering state. Keep the ORIGINAL configured path
    # in the request identity so existing science does not become a new model.
    p=c['preprocessing'];r=dict(mode=p['mode'],project=str(logical_resolve(p['project'])))
    if p['mode']=='new':r.update(runtime_config_sha256=sha(resolve(p['runtime_config'])),bids_filter=p['bids_filter'],
        expected_image_id=p['expected_image_id'])
    if p['mode']=='new':
        protocol=read(resolve(p['runtime_config'])).get('scientific_protocol')
        if protocol:
            from scientific_protocol import identity
            r['scientific_protocol_inputs']=identity(source_root(c),protocol)
            r['scientific_protocol_engine_sha256']=sha(CODE/'scientific_protocol.py')
    return digest(r)

@contextmanager
def lock(c):
    p=workspace(c)/'.stage.lock';p.parent.mkdir(parents=True,exist_ok=True)
    try:f=p.open('x',encoding='utf-8')
    except FileExistsError:raise StageError(f'Another stage owns {p}; inspect its PID before removing a stale lock')
    try:
        with f:f.write(json.dumps(dict(pid=os.getpid(),started_UTC=now())))
        yield
    finally:p.unlink(missing_ok=True)

def descriptor(p,base):
    s=p.stat();return dict(path=p.relative_to(base).as_posix(),bytes=s.st_size,mtime_ns=s.st_mtime_ns,sha256=sha(p))
def file_identity(rows):return [dict(path=r['path'],bytes=r['bytes'],sha256=r['sha256']) for r in rows]
def check_files(base,rows,full=False):
    """Fast path uses pinned size/mtime; changed stamps require SHA verification.

    --verify-full ignores stamps. No mismatched input is silently adopted.
    """
    checked=0
    for r in rows:
        p=base/r['path']
        if not p.is_file():raise StageError('Required file missing: '+str(p))
        s=p.stat()
        if s.st_size!=r['bytes']:raise StageError('Size changed: '+str(p))
        if full or s.st_mtime_ns!=r['mtime_ns']:
            if sha(p)!=r['sha256']:raise StageError('Hash changed: '+str(p)+'; preserve this run and register a new version')
            checked+=1
    return checked

def preprocessed_files(c):
    base=derivatives(c)
    files=sorted(p for p in base.rglob('*') if p.is_file() and
        (p.name.endswith(('.nii.gz','.h5','.tsv','.json'))) and 'logs' not in p.relative_to(base).parts)
    if not files:raise StageError('No preprocessed scientific outputs found: '+str(base))
    return files

def preprocessing(c,full=False):
    """Register existing outputs once, or explicitly compute a NEW version."""
    start=time.perf_counter();state=workspace(c);receipt=state/'preprocessing/receipt.json'
    if receipt.is_file():
        r=read(receipt)
        if resolve(r['source_project'])!=source_root(c):raise StageError('Preprocessing source changed; use a new stage workspace')
        if r['request_key']!=preprocessing_request(c):raise StageError('Preprocessing request changed; use a new output version/workspace')
        checked=check_files(derivatives(c),r['files'],full)
        check_files(source_root(c),r['proof_files'],full)
        if r.get('BIDS_files'):check_files(source_root(c)/'bids',r['BIDS_files'],full)
        return dict(stage='preprocess',action='REUSED_VERIFIED_OUTPUTS',fingerprint=r['fingerprint'],
                    files_rehashed=checked,fmriprep_invocations=0,elapsed_seconds=time.perf_counter()-start)
    mode=c['preprocessing']['mode']
    invocations=0
    if mode=='new':
        from staged_preprocess import launch_new
        launch_new(c);invocations=1
    elif mode not in ['existing','engineering_subset']:raise StageError('Unknown preprocessing mode')
    project=source_root(c)
    if mode=='engineering_subset':
        proof=read(project/'subset_provenance.json')
        if proof.get('kind') not in ['CROPPED_PREPROCESSED_ENGINEERING_INPUTS','SYNTHETIC_ENGINEERING_INPUTS']:raise StageError('Missing engineering input provenance')
        source_kind=proof['kind']+'_NOT_A_NEW_FMRIPREP_RUN';actual_threads=None
        proof_paths=[project/'subset_provenance.json']
    else:
        proof=read(project/'logs/preprocessing_status.json')
        if proof.get('exit_code')!=0:raise StageError('Existing preprocessing has no successful completion record')
        run=read(project/'logs/fmriprep_run_configuration.json')
        args=run['arguments'];actual_threads=int(args[args.index('--omp-nthreads')+1])
        if mode=='new' and actual_threads!=THREADS:raise StageError('New preprocessing must use four threads')
        source_kind='NEW_COMPLETED_FMRIPREP' if mode=='new' else 'EXISTING_COMPLETED_FMRIPREP';proof=dict(status=proof,run_configuration=run)
        proof_paths=[project/'logs/preprocessing_status.json',project/'logs/fmriprep_run_configuration.json']
    rows=[descriptor(p,derivatives(c)) for p in preprocessed_files(c)]
    fingerprint=digest(dict(files=file_identity(rows),source_kind=source_kind,proof=proof))
    record=dict(schema='preprocessing-receipt-1',case_id=c['case_id'],source_project=str(project),
        source_kind=source_kind,historical_actual_omp_threads=actual_threads,new_run_thread_policy=THREADS,
        fingerprint=fingerprint,request_key=preprocessing_request(c),files=rows,
        proof_files=[descriptor(p,project) for p in proof_paths],registered_UTC=now(),registration_hash_policy='Full SHA256 at first registration',
        BIDS_files=[descriptor(p,project/'bids') for p in sorted((project/'bids').rglob('*')) if p.is_file()] if mode=='new' else [],
        historical_results_modified=False)
    if mode=='new':
        protocol=read(resolve(c['preprocessing']['runtime_config'])).get('scientific_protocol')
        if protocol:record['scientific_protocol']=protocol
    write(receipt,record)
    return dict(stage='preprocess',action='REGISTERED_COMPLETED_OUTPUTS' if not invocations else 'NEW_PREPROCESSING_COMPLETED',
        fingerprint=fingerprint,files_rehashed=len(rows),fmriprep_invocations=invocations,elapsed_seconds=time.perf_counter()-start)

def require_preprocessing(c,full=False):
    p=workspace(c)/'preprocessing/receipt.json'
    if not p.is_file():raise StageError('Run the explicit preprocess stage once to register/compute preprocessing; GLM never launches it')
    r=read(p)
    if resolve(r['source_project'])!=source_root(c):raise StageError('Preprocessing source does not match registered source')
    if r['request_key']!=preprocessing_request(c):raise StageError('Preprocessing request changed; use a new workspace')
    check_files(derivatives(c),r['files'],full);check_files(source_root(c),r['proof_files'],full)
    if r.get('scientific_protocol'):
        from scientific_protocol import require_review
        review=source_root(c)/'logs/scientific_protocol_review.json'
        if not review.is_file():raise StageError('New branch method evidence review missing; GLM remains blocked')
        require_review(r['scientific_protocol'],r['fingerprint'],read(review))
        for evidence in read(review)['evidence_files']:
            path=(source_root(c)/evidence['path']).resolve()
            if not path.is_relative_to(source_root(c)) or sha(path)!=evidence['sha256']:raise StageError('Method evidence file missing, outside project, or changed')
    return r

def task_gate(c):
    a=c['analysis'];gate=a['timing']
    if gate['status']=='paused':raise TimingBlocked(gate['reason'])
    if gate['status'] not in ['confirmed','conditional_user_authorized','conditional_shared_events_user_authorized','engineering_subset']:
        raise TimingBlocked('Task timing has not been authorized or confirmed')
    evidence=resolve(gate['evidence'])
    if not evidence.is_file():raise TimingBlocked('Timing evidence file missing')
    if gate['status']=='confirmed':
        clearance=read(evidence)
        if clearance.get('status')!='CONFIRMED_SCAN_SPECIFIC_TIMING' or clearance.get('case_id')!=c['case_id']:
            raise TimingBlocked('Confirmed timing requires a scan-specific clearance record for this case')
        registered=read(workspace(c)/'preprocessing/receipt.json')
        if clearance.get('preprocessing_fingerprint')!=registered['fingerprint']:
            raise TimingBlocked('Timing clearance is not bound to the registered preprocessing inputs')
        for run in a['runs']:
            if clearance.get('event_sha256_by_run',{}).get(run['id'])!=sha(resolve(run['events'])):
                raise TimingBlocked('Timing clearance is not bound to this run/events file')
    if gate['status']=='conditional_user_authorized' and c.get('report',{}).get('participant_group')!='healthy_control':
        raise TimingBlocked('Conditional healthy timing cannot release the paused stroke-patient GLM')
    if gate['status']=='conditional_shared_events_user_authorized':
        clearance=read(evidence)
        if clearance.get('status')!='USER_AUTHORIZED_EXPLORATORY_SHARED_EVENTS' or clearance.get('case_id')!=c['case_id']:
            raise TimingBlocked('Shared-event exploration requires case-specific user authorization')
        registered=read(workspace(c)/'preprocessing/receipt.json')
        if clearance.get('preprocessing_fingerprint')!=registered['fingerprint']:
            raise TimingBlocked('Shared-event authorization is not bound to these preprocessed inputs')
        if clearance.get('independent_stimulus_log_verified') is not False or not clearance.get('user_answer'):
            raise TimingBlocked('Shared-event evidence must preserve the unverified log status and actual user answer')
        if c.get('report',{}).get('participant_group')!='chronic_stroke':
            raise TimingBlocked('This shared-event authorization schema is for explicitly approved stroke exploration')
        if a['statistics']['contrast']!={'image':1.0,'abstract':-1.0}:
            raise TimingBlocked('Authorized shared-event contrast is image > abstract; review a different scientific contrast')
        for run in a['runs']:
            if clearance.get('event_sha256_by_run',{}).get(run['id'])!=sha(resolve(run['events'])):
                raise TimingBlocked('Shared-event authorization does not match this run/events file')
            if clearance.get('TR_source_sha256_by_run',{}).get(run['id'])!=sha(resolve(run['tr_evidence'])):
                raise TimingBlocked('Shared-event authorization does not match this run/TR evidence')
    if c['preprocessing']['mode']=='engineering_subset' and gate['status']!='engineering_subset':
        raise TimingBlocked('Subset must be labelled engineering, not clinical/confirmed')
    return dict(status=gate['status'],evidence_sha256=sha(evidence),interpretation=gate['interpretation'])

def statistical_request(c,pre):
    timing=task_gate(c);a=c['analysis'];files=[]
    for run in a['runs']:
        p=resolve(run['events'])
        if not p.is_file():raise TimingBlocked('Run-specific events missing: '+str(p))
        tr_file=resolve(run['tr_evidence']);tr_metadata=read(tr_file)
        if float(tr_metadata['RepetitionTime'])!=float(run['tr_seconds']):raise StageError('Configured nominal TR differs from its source JSON')
        files.append(dict(run=run['id'],sha256=sha(p),TR_evidence_sha256=sha(tr_file)))
    spec=dict(preprocessing_fingerprint=pre['fingerprint'],runs=a['runs'],statistics=a['statistics'],
        timing=timing,event_contents=files,engine_sha256=sha(CODE/'staged_glm.py'),thread_policy=THREADS)
    spec['release_policy_sources']={name:sha(CODE/name) for name in ['staged_workflow.py','model_audit.py','scientific_protocol.py']}
    if a.get('engine_variant','baseline')=='anatomical_intersection_v1':
        relative=a['anatomical_mask']
        path=(derivatives(c)/relative).resolve()
        registered={r['path']:r for r in pre['files']}
        if not path.is_relative_to(derivatives(c).resolve()) or relative not in registered:
            raise StageError('Anatomical mask must belong to registered preprocessing')
        spec.update(engine_variant='anatomical_intersection_v1',
            base_engine_sha256=spec['engine_sha256'],engine_sha256=sha(CODE/'staged_glm_anatomical.py'),
            anatomical_mask=dict(path=relative,sha256=sha(path)))
    return digest(spec),spec

def completed(folder,full=False):
    p=folder/'receipt.json'
    if not p.is_file():return None
    r=read(p)
    if r.get('status')!='COMPLETED_NUMERICALLY_REVIEW_PENDING':raise StageError('Output receipt is not completed')
    check_files(folder,r['files'],full);return r

def glm(c,full=False):
    start=time.perf_counter();pre=require_preprocessing(c,full);key,spec=statistical_request(c,pre)
    dest=workspace(c)/'glm'/key
    if completed(dest,full):return dict(stage='glm',action='REUSED_IDENTICAL_STATISTICS',key=key,output=str(dest),
        fmriprep_invocations=0,glm_fits=0,elapsed_seconds=time.perf_counter()-start)
    if dest.exists():raise StageError('Incomplete statistical directory preserved: '+str(dest)+'; inspect it before an explicit new attempt')
    dest.mkdir(parents=True);write(dest/'request.json',spec)
    try:
        from model_audit import audit_case
        model_audit=audit_case(c);write(dest/'model_audit.json',model_audit)
        if not model_audit['fit_allowed']:
            raise StageError('Model audit blocked fitting: '+'; '.join(r.get('reason','') for r in model_audit['runs'] if r['status']=='BLOCKED_DESIGN_REVIEW_REQUIRED'))
        if c['analysis'].get('engine_variant','baseline')=='anatomical_intersection_v1':
            from staged_glm_anatomical import fit
        else:
            from staged_glm import fit
        summary=fit(c,dest)
        rows=[descriptor(p,dest) for p in sorted(dest.rglob('*')) if p.is_file()]
        write(dest/'receipt.json',dict(status='COMPLETED_NUMERICALLY_REVIEW_PENDING',key=key,files=rows,
            preprocessing_fingerprint=pre['fingerprint'],finished_UTC=now(),clinical_validation=False))
    except Exception as exc:
        write(dest/'failure.json',dict(error=str(exc),UTC=now()));raise
    return dict(stage='glm',action='FITTED_NEW_STATISTICAL_VERSION',key=key,output=str(dest),
        fmriprep_invocations=0,glm_fits=len(c['analysis']['runs']),elapsed_seconds=time.perf_counter()-start)

def report(c,full=False):
    start=time.perf_counter();pre=require_preprocessing(c,full);stat=None;reason=None;key=None
    try:
        key,spec=statistical_request(c,pre);folder=workspace(c)/'glm'/key
        if completed(folder,full):stat=folder
        else:reason='Current statistical request has not been fitted; report stage does not launch GLM.'
    except TimingBlocked as exc:reason='TASK_GLM_PAUSED: '+str(exc)
    request=dict(preprocessing_fingerprint=pre['fingerprint'],statistics_key=key,reason=reason,
        context=c.get('report',{}),statistics_configuration=c['analysis']['statistics'],renderer_sha256=sha(CODE/'staged_report.py'))
    report_key=digest(request);dest=workspace(c)/'reports'/report_key
    if (dest/'report_receipt.json').is_file():
        check_files(dest,read(dest/'report_receipt.json')['files'],full)
        return dict(stage='report',action='REUSED_IDENTICAL_REPORT',output=str(dest/'report.html'),
            fmriprep_invocations=0,glm_fits=0,elapsed_seconds=time.perf_counter()-start)
    if dest.exists():raise StageError('Incomplete report directory preserved: '+str(dest))
    from staged_report import render
    dest.mkdir(parents=True);render(c,pre,stat,reason,dest)
    rows=[descriptor(p,dest) for p in sorted(dest.rglob('*')) if p.is_file()]
    write(dest/'report_receipt.json',dict(status='STATIC_REPORT_GENERATED_INTERACTION_NOT_VERIFIED',files=rows,
        request=request,finished_UTC=now(),glm_fits=0,fmriprep_invocations=0))
    return dict(stage='report',action='RENDERED_FROM_SAVED_OUTPUTS',output=str(dest/'report.html'),
        fmriprep_invocations=0,glm_fits=0,elapsed_seconds=time.perf_counter()-start)

def status(c):
    pre=workspace(c)/'preprocessing/receipt.json'
    return dict(case_id=c['case_id'],workspace=str(workspace(c)),preprocessing_registered=pre.is_file(),
        timing_status=c['analysis']['timing']['status'],new_tool_threads=THREADS,
        glm_versions=len(list((workspace(c)/'glm').glob('*/receipt.json'))),
        reports=len(list((workspace(c)/'reports').glob('*/report_receipt.json'))))

def dispatch(c,stage,full=False):
    configure_threads()
    with lock(c):
        if stage=='status':return status(c)
        try:r={'preprocess':preprocessing,'glm':glm,'report':report,'audit':audit}[stage](c,full)
        except Exception as exc:
            write(workspace(c)/'last_action.json',dict(stage=stage,status='BLOCKED' if isinstance(exc,TimingBlocked) else 'FAILED',error=str(exc),UTC=now()));raise
        write(workspace(c)/'last_action.json',dict(**r,UTC=now()));return r

def audit(c,full=False):
    pre=require_preprocessing(c,full);key,_=statistical_request(c,pre)
    from model_audit import audit_case
    result=audit_case(c);folder=workspace(c)/'audits'/key
    write(folder/'model_audit.json',result)
    return dict(stage='audit',output=str(folder/'model_audit.json'),glm_fits=0,fmriprep_invocations=0)
