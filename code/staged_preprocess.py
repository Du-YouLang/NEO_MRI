"""Explicit local Linux Docker backend. Planning never calls a process.

No daemon startup, WSL discovery, image pull, service repair, or automatic retry.
This portable backend is NEW and has not had a real-image container acceptance run.
"""
from pathlib import Path
import json, os, re, subprocess, time
from staged_workflow import read, write, resolve, source_root, StageError, sha
from scientific_protocol import arguments, validate_inputs

IMAGE='nipreps/fmriprep:25.2.7'

def local_linux():
    return os.name=='posix'

def build_plan(c):
    p=c['preprocessing'];runtime=read(resolve(p['runtime_config']))
    required={'container_image','omp_nthreads','nprocs','memory_mb','container_memory_mb',
              'container_cpus','license_file','work_dir','templateflow_dir','output_spaces','participant_label','scientific_protocol'}
    if set(runtime)-required-{'anat_only','execution_enabled'} or required-set(runtime):
        raise StageError('Unknown or missing portable runtime options')
    if runtime['container_image']!=IMAGE:raise StageError('Only the pinned fMRIPrep version is supported')
    if any(type(runtime[k]) is not int for k in ['omp_nthreads','nprocs','memory_mb','container_memory_mb','container_cpus']):
        raise StageError('Thread, CPU and memory budgets must be integers')
    if runtime['omp_nthreads']!=4 or not 4<=runtime['nprocs']<=runtime['container_cpus']:
        raise StageError('Require four tool threads and 4 <= nprocs <= container CPU limit')
    if not 1000<=runtime['memory_mb']<runtime['container_memory_mb']:raise StageError('Invalid memory budget')
    digest=p.get('expected_image_id','')
    if not re.fullmatch(r'sha256:[0-9a-f]{64}',digest):raise StageError('Require locally inspected immutable Docker image ID')
    if not re.fullmatch(r'[A-Za-z0-9]+',runtime['participant_label']):raise StageError('Invalid participant label')
    root=source_root(c);work=resolve(runtime['work_dir']);templates=resolve(runtime['templateflow_dir'])
    if work==root or work.is_relative_to(root) or root.is_relative_to(work):raise StageError('Work directory must be separate from input/output project')
    protocol=runtime['scientific_protocol'];scientific=arguments(protocol)
    spaces=runtime['output_spaces']
    if not spaces or not all(isinstance(s,str) and re.fullmatch(r'[A-Za-z0-9_:.-]+',s) for s in spaces):raise StageError('Explicit valid output spaces required')
    output=root/'derivatives/fmriprep';license_file=resolve(runtime['license_file'])
    def mount(source,target,ro=False):
        if ',' in str(source):raise StageError('Comma in Docker mount path is unsupported')
        return ['--mount',f'type=bind,source={source},target={target}'+(',readonly' if ro else '')]
    cmd=['docker','run','--name','neo-mri-exclusive','--pull=never',
         '--cpus',str(runtime['container_cpus']),'--memory',str(runtime['container_memory_mb'])+'m',
         '--shm-size','4g','--user',f'{os.getuid()}:{os.getgid()}' if hasattr(os,'getuid') else '0:0']
    for name,value in [('OMP_NUM_THREADS','4'),('ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS','4'),('OPENBLAS_NUM_THREADS','4'),('MKL_NUM_THREADS','4'),('TEMPLATEFLOW_HOME','/templateflow')]:
        cmd+=['-e',name+'='+value]
    cmd+=mount(root/'bids','/data',True)+mount(output,'/out')+mount(work,'/work')+mount(templates,'/templateflow')+mount(license_file,'/license.txt',True)
    cmd += [digest,'/data','/out','participant','--participant-label',runtime['participant_label'],
            '--fs-license-file','/license.txt','--fs-no-reconall','--nprocs',str(runtime['nprocs']),
            '--omp-nthreads','4','--mem',str(runtime['memory_mb']),'-w','/work',
            '--random-seed','20261026','--output-spaces',*spaces,'--bids-filter-file','/data/neo_bids_filter.json',*scientific]
    if runtime.get('anat_only'):cmd+=['--anat-only']
    return dict(status='COMMAND_PLAN_ONLY_NO_DOCKER_PROBED',arguments=cmd,runtime=runtime,
                image_id=digest,project=str(root),work=str(work),templateflow=str(templates),license_path=str(license_file),
                new_backend_real_data_verified=False)

def launch_new(c):
    plan=build_plan(c);runtime=plan['runtime'];root=Path(plan['project'])
    if not runtime.get('execution_enabled') or not runtime['scientific_protocol'].get('execution_enabled'):
        raise StageError('Portable preprocessing is disabled; explicitly review and enable both execution flags')
    if not local_linux():raise StageError('Run portable preprocessing from a local Linux/WSL shell; no WSL is started automatically')
    if (root/'logs').exists() or (root/'derivatives').exists():raise StageError('Fresh project required; retain any previous outputs/logs before a separately reviewed retry')
    if not (root/'bids/dataset_description.json').is_file():raise StageError('Reviewed BIDS dataset required')
    validate_inputs(root,runtime['scientific_protocol'])
    if not Path(plan['license_path']).is_file():raise StageError('Local license file is missing')
    work=Path(plan['work']);templates=Path(plan['templateflow'])
    if work.exists() and any(work.iterdir()):raise StageError('Fresh case-specific work directory required in this preview')
    if not templates.is_dir() or not any(templates.iterdir()):raise StageError('Pre-populate and review TemplateFlow cache before execution')
    # Only query the existing daemon. Never pull or start/restart services.
    for image in [plan['image_id'],IMAGE]:
        probe=subprocess.run(['docker','image','inspect',image,'--format','{{.Id}}'],capture_output=True,text=True,timeout=20)
        if probe.returncode or probe.stdout.strip()!=plan['image_id']:raise StageError('Pinned local image/tag not available or mismatched; no image pull attempted')
    expected_filter=c['preprocessing']['bids_filter'];filter_path=root/'bids/neo_bids_filter.json'
    if not filter_path.is_file() or read(filter_path)!=expected_filter:
        raise StageError('Provide matching neo_bids_filter.json in the reviewed BIDS staging copy; original inputs are never edited')
    logs=root/'logs';logs.mkdir();out=root/'derivatives/fmriprep';out.mkdir(parents=True);work.mkdir(parents=True,exist_ok=True)
    write(logs/'fmriprep_run_configuration.json',dict(arguments=plan['arguments'],image_id=plan['image_id'],runtime=runtime,
        bids_files=[dict(path=p.relative_to(root/'bids').as_posix(),sha256=sha(p)) for p in sorted((root/'bids').rglob('*')) if p.is_file()]))
    started=time.monotonic()
    with (logs/'fmriprep_console.log').open('w',encoding='utf-8') as console:
        process=subprocess.Popen(plan['arguments'],stdout=console,stderr=subprocess.STDOUT)
        write(logs/'process.json',dict(pid=process.pid,container_name='neo-mri-exclusive',automatic_retry=False))
        code=process.wait()
    write(logs/'preprocessing_status.json',dict(exit_code=code,elapsed_seconds=time.monotonic()-started,
        actual_container_verification='EXIT_CODE_ONLY_REQUIRES_IMAGE_AND_METHOD_REVIEW'))
    if code:raise StageError('Preprocessing failed; container and outputs retained. No automatic retry or service changes.')
