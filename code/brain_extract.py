"""Bounded local Linux brain-mask candidates. Dependencies/models installed separately."""
from pathlib import Path
import argparse,os,signal,subprocess,time
from anatomy_review import read,save_new,file_record,checked_mask

DEFAULT_METHOD='hd_bet_cpu_no_tta'
METHODS=('hd_bet_cpu_no_tta','synthstrip','ants_oasis','ants_mni2009c','ants_template')

def selection_policy():
    return dict(schema='neo-brain-extraction-policy-1',revision='2026-11-04',default_method=DEFAULT_METHOD,
        optional_methods=['synthstrip','ants_oasis','ants_mni2009c'],legacy_explicit_method='ants_template',
        automatic_fallback=False,human_mask_review_required=True,automatic_MNI_approval=False,
        hd_bet_profile=dict(device='cpu',tta=False,threads=4),
        templates=dict(ants_oasis='OASIS30ANTs',ants_mni2009c='MNI152NLin2009cAsym'),
        note='User-selected default for future inputs. Preserve explicit existing method choices; never select by largest mask or highest Dice. ANTs extraction templates do not change the target MNI space.')

def plan(c):
    t1=Path(c['t1']).resolve();out=Path(c['output']).resolve();method=c.get('method',DEFAULT_METHOD)
    if method not in METHODS:raise ValueError('Unknown candidate method')
    runtime=c.get('runtimes',{}).get(method,c.get('runtime'))
    if not isinstance(runtime,dict):raise ValueError('Runtime for selected method '+method+' must be configured; no automatic fallback')
    if c.get('threads',4)!=4 or not 1<=c.get('timeout_seconds',1800)<=1800:raise ValueError('Four threads and at most 1800 seconds required')
    if out==t1.parent or out in t1.parents:raise ValueError('Separate candidate output directory required')
    mask=out/'mask.nii.gz';artifacts=[]
    if method=='synthstrip':
        script=Path(runtime['script']).resolve();model=Path(runtime['model']).resolve();artifacts=[script,model]
        argv=[runtime['python'],str(script),'-i',str(t1),'-m',str(mask),'--model',str(model),'-t','4']
    elif method=='hd_bet_cpu_no_tta':
        model=Path(runtime['model_dir']).resolve()
        artifacts=sorted(p for p in model.rglob('*') if p.is_file())
        if not artifacts:raise ValueError('Preinstalled HD-BET model directory required')
        argv=[runtime['python'],str(Path(__file__).with_name('hdbet_cpu.py')),'--input',str(t1),'--output',str(mask),'--model-dir',str(model)]
    else:
        script=Path(runtime['script']).resolve();template=Path(runtime['template']).resolve();prior=Path(runtime['prior']).resolve();regmask=Path(runtime['registration_mask']).resolve()
        artifacts=[script,template,prior,regmask];mask=out/'ants_BrainExtractionMask.nii.gz'
        if method in ['ants_oasis','ants_mni2009c']:
            template_id=selection_policy()['templates'][method]
            if runtime.get('template_id')!=template_id or not all(p.name.startswith('tpl-'+template_id+'_') for p in [template,prior,regmask]):
                raise ValueError('Selected ANTs template identity and all three resources must match '+template_id)
        argv=[str(script),'-d','3','-a',str(t1),'-e',str(template),'-m',str(prior),'-f',str(regmask),'-u','0','-o',str(out/'ants_')]
    return dict(schema='neo-t1-extraction-request-1',method=method,input=file_record(t1),artifacts=[file_record(p) for p in artifacts],
        command=argv,mask_path=str(mask),output=str(out),threads=4,gpu=False,timeout_seconds=c.get('timeout_seconds',1800),
        selection=dict(policy_revision=selection_policy()['revision'],default_method=DEFAULT_METHOD,explicit_method='method' in c,automatic_fallback=False),
        methodological_note='HD-BET noTTA is an explicit speed profile; ANTs script repeats N4 and is not a pure one-factor comparison against historical fMRIPrep.')

def run(c):
    if os.name!='posix':raise ValueError('Use a Linux runtime; this command never starts WSL or services')
    request=plan(c);out=Path(request['output'])
    if out.exists():raise FileExistsError('Candidate version already exists; no automatic rerun/overwrite')
    if not c.get('execution_enabled'):raise ValueError('Explicit execution_enabled required')
    out.mkdir(parents=True);save_new(out/'request.json',request)
    env=os.environ.copy()
    for name in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS','NUMEXPR_NUM_THREADS']:env[name]='4'
    env['CUDA_VISIBLE_DEVICES']='-1';start=time.monotonic()
    with (out/'console.log').open('w') as log:
        process=subprocess.Popen(request['command'],stdout=log,stderr=subprocess.STDOUT,env=env,start_new_session=True)
        status='COMPLETED_UNREVIEWED'
        try:code=process.wait(timeout=request['timeout_seconds'])
        except subprocess.TimeoutExpired:
            status='TIME_BUDGET_EXCEEDED_NOT_COMPLETED';os.killpg(process.pid,signal.SIGTERM)
            try:code=process.wait(timeout=15)
            except subprocess.TimeoutExpired:os.killpg(process.pid,signal.SIGKILL);code=process.wait()
        if code and status=='COMPLETED_UNREVIEWED':status='FAILED'
    record=dict(status=status,exit_code=code,elapsed_seconds=time.monotonic()-start,request=request,human_decision='PENDING')
    if status=='COMPLETED_UNREVIEWED':
        try:
            checked_mask(c['t1'],request['mask_path']);record['mask']=file_record(request['mask_path'])
        except (ValueError,OSError) as exc:record.update(status='INVALID_OR_MISSING_MASK',error=str(exc))
    save_new(out/'receipt.json',record);return record

def main(argv=None):
    import json
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=['policy','plan','run']);p.add_argument('--config');a=p.parse_args(argv)
    if a.action=='policy':print(json.dumps(selection_policy(),ensure_ascii=False,indent=2));return 0
    if not a.config:p.error('--config is required for plan/run')
    result=plan(read(a.config)) if a.action=='plan' else run(read(a.config));print(json.dumps(result,ensure_ascii=False));return 0 if result.get('status','COMPLETED_UNREVIEWED')=='COMPLETED_UNREVIEWED' else 2
