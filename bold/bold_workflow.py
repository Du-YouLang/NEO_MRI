"""Portable lifecycle for BOLD v2. Existing completed derivatives only."""
from pathlib import Path
from contextlib import contextmanager
import datetime, hashlib, json, os, re, time

CODE = Path(__file__).resolve().parent
THREADS = 4


class StageError(RuntimeError):
    pass


class TimingBlocked(StageError):
    pass


def resolve(path):
    # load normalizes every declared portable path before any stage executes.
    return Path(path).resolve()


def data_root(config):
    return resolve(config["data_root"])


def workspace(config):
    return resolve(config["workspace"])


def source_root(config):
    return resolve(config["preprocessing"]["project"])


def derivatives(config):
    return source_root(config)/"derivatives/fmriprep"


def preprocessing_request(config):
    return digest(dict(mode=config["preprocessing"]["mode"],project=str(source_root(config))))


def load(path):
    from bold_runs import SCHEMA, validate
    config_path = resolve(path)
    c = json.loads(config_path.read_text(encoding="utf-8-sig"))
    if c.get("template_only"):
        raise StageError("Copy and complete the example before use")
    if c.get("schema") != SCHEMA or not re.fullmatch(r"[A-Za-z0-9_-]+",c["case_id"]):
        raise StageError("Expected a version 2 configuration and a safe case identifier")

    def pin(record,key):
        if record.get(key):
            record[key]=str((config_path.parent/record[key]).resolve())

    pin(c,"data_root")
    pin(c,"workspace")
    pin(c["preprocessing"],"project")
    state,source=workspace(c),source_root(c)
    if state==source or state in source.parents or source in state.parents:
        raise StageError("Choose an output root separate from the preserved preprocessing project")
    if not state.is_relative_to(data_root(c)):
        raise StageError("Choose an output under the declared data_root")
    analysis=c["analysis"]
    pin(analysis["timing"],"evidence")
    for run in analysis["runs"]:
        for key in ("events","tr_evidence","method_review"):
            pin(run,key)
        pin(run.get("timing",{}),"evidence")
    known={"hrf_model","smoothing_fwhm_mm","noise_model","compcor_components","fd_threshold_mm","dvars_threshold","contrast","fdr_q"}
    if set(analysis["statistics"])!=known:
        raise StageError("Missing or unknown statistics parameter")
    if analysis.get("engine_variant","baseline") not in ("baseline","anatomical_intersection_v1"):
        raise StageError("Unknown GLM engine variant")
    validate(c)
    c["_config_path"]=str(config_path)
    return c


def dispatch(c,stage,full=False):
    from bold_runs import inventory
    if stage=="inventory":
        return inventory(c)
    import bold_stage
    functions={"preprocess":preprocessing,"audit":bold_stage.audit,"glm":bold_stage.glm,
               "compare":bold_stage.compare,"report":bold_stage.report}
    if stage=="status":
        base=workspace(c)
        return dict(preprocessing_registered=(base/"preprocessing/receipt.json").is_file(),
                    fits=len(list((base/"glm").glob("*/receipt.json"))),
                    reports=len(list((base/"reports").glob("*/report_receipt.json"))))
    configure_threads()
    with lock(c):
        try:
            result=functions[stage](c,full)
        except Exception as exc:
            write(workspace(c)/"last_action.json",dict(stage=stage,action="BLOCKED",reason=str(exc),UTC=now()))
            raise
        write(workspace(c)/"last_action.json",dict(**result,UTC=now()))
        return result


def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def read(path):return json.loads(Path(path).read_text(encoding='utf-8'))
def digest(obj):return hashlib.sha256(json.dumps(obj,sort_keys=True,separators=(',',':'),ensure_ascii=False,allow_nan=False).encode()).hexdigest()
def sha(path):
    with Path(path).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def write(path,obj):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_name(path.name+'.tmp')
    tmp.write_text(json.dumps(obj,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8');tmp.replace(path)
def configure_threads():
    global _THREAD_LIMITER
    for name in ['OMP_NUM_THREADS','ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:
        os.environ[name]=str(THREADS)
    # Also constrain libraries that an interactive caller imported beforehand.
    from threadpoolctl import threadpool_limits
    _THREAD_LIMITER=threadpool_limits(limits=THREADS)
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
        (p.name.endswith(('.nii.gz','.h5','.tsv','.json','.html'))) and 'logs' not in p.relative_to(base).parts)
    if not files:raise StageError('No preprocessed scientific outputs found: '+str(base))
    return files

def preprocessing(c,full=False):
    """Register already completed outputs. This stage never launches fMRIPrep."""
    start=time.perf_counter();state=workspace(c);receipt=state/'preprocessing/receipt.json'
    if receipt.is_file():
        r=read(receipt)
        if r['case_id']!=c['case_id'] or resolve(r['source_project'])!=source_root(c):
            raise StageError('Preprocessing case or source changed; use a new stage workspace')
        if r['request_key']!=preprocessing_request(c):raise StageError('Preprocessing request changed; use a new output version/workspace')
        checked=check_files(derivatives(c),r['files'],full)
        check_files(source_root(c),r['proof_files'],full)
        if r.get('BIDS_files'):check_files(source_root(c)/'bids',r['BIDS_files'],full)
        return dict(stage='preprocess',action='REUSED_VERIFIED_OUTPUTS',fingerprint=r['fingerprint'],
                    files_rehashed=checked,fmriprep_invocations=0,elapsed_seconds=time.perf_counter()-start)
    mode=c['preprocessing']['mode']
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
        source_kind='EXISTING_COMPLETED_FMRIPREP';proof=dict(status=proof,run_configuration=run)
        proof_paths=[project/'logs/preprocessing_status.json',project/'logs/fmriprep_run_configuration.json']
    rows=[descriptor(p,derivatives(c)) for p in preprocessed_files(c)]
    fingerprint=digest(dict(files=file_identity(rows),source_kind=source_kind,proof=proof))
    record=dict(schema='preprocessing-receipt-1',case_id=c['case_id'],source_project=str(project),
        source_kind=source_kind,historical_actual_omp_threads=actual_threads,new_run_thread_policy=THREADS,
        fingerprint=fingerprint,request_key=preprocessing_request(c),files=rows,
        proof_files=[descriptor(p,project) for p in proof_paths],registered_UTC=now(),registration_hash_policy='Full SHA256 at first registration',
        historical_results_modified=False)
    write(receipt,record)
    return dict(stage='preprocess',action='REGISTERED_COMPLETED_OUTPUTS',
        fingerprint=fingerprint,files_rehashed=len(rows),fmriprep_invocations=0,elapsed_seconds=time.perf_counter()-start)

def require_preprocessing(c,full=False):
    p=workspace(c)/'preprocessing/receipt.json'
    if not p.is_file():raise StageError('Run the explicit preprocess stage once to register completed preprocessing')
    r=read(p)
    if r['case_id']!=c['case_id'] or resolve(r['source_project'])!=source_root(c):
        raise StageError('Preprocessing case or source does not match the registered input')
    if r['request_key']!=preprocessing_request(c):raise StageError('Preprocessing request changed; use a new workspace')
    check_files(derivatives(c),r['files'],full);check_files(source_root(c),r['proof_files'],full)
    if r.get('scientific_protocol'):
        raise StageError("Use a new BOLD workspace to register derivatives under the current reviewed method evidence")
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
    if gate['status']=='conditional_user_authorized':
        clearance=read(evidence)
        registered=read(workspace(c)/'preprocessing/receipt.json')
        if (clearance.get('status')!='USER_AUTHORIZED_CONDITIONAL_TIMING'
            or clearance.get('case_id')!=c['case_id']
            or clearance.get('preprocessing_fingerprint')!=registered['fingerprint']
            or clearance.get('independent_stimulus_log_verified') is not False
            or not clearance.get('user_answer')):
            raise TimingBlocked('Conditional healthy timing requires case-specific written authorization')
        for run in a['runs']:
            if clearance.get('event_sha256_by_run',{}).get(run['id'])!=sha(resolve(run['events'])):
                raise TimingBlocked('Conditional timing does not identify this run/events')
            if clearance.get('TR_source_sha256_by_run',{}).get(run['id'])!=sha(resolve(run['tr_evidence'])):
                raise TimingBlocked('Conditional timing does not identify the TR source')
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

def completed(folder,full=False):
    p=folder/'receipt.json'
    if not p.is_file():return None
    r=read(p)
    if r.get('status')!='COMPLETED_NUMERICALLY_REVIEW_PENDING':raise StageError('Output receipt is not completed')
    check_files(folder,r['files'],full);return r
