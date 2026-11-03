from pathlib import Path
import copy, importlib.util, json, subprocess, sys
import nibabel as nib
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'code'))
import staged_workflow as w
from staged_preprocess import build_plan, launch_new
from test_staged_workflow import case

def config(tmp_path,monkeypatch):
    monkeypatch.setattr(w,'REPO',tmp_path)
    protocol=dict(schema='neo-mri-protocol-1',execution_enabled=False,normalization='standard',sdc='off',stc='off',t1w='sub-demo/anat/sub-demo_T1w.nii.gz')
    runtime=dict(container_image='nipreps/fmriprep:25.2.7',execution_enabled=False,omp_nthreads=4,nprocs=8,memory_mb=8000,
        container_memory_mb=9000,container_cpus=8,license_file='secrets/license.txt',work_dir='work/run1',templateflow_dir='templates',
        output_spaces=['T1w','MNI152NLin2009cAsym:res-2'],participant_label='demo',scientific_protocol=protocol)
    w.write(tmp_path/'runtime.json',runtime)
    return dict(case_id='demo',workspace='results/demo',preprocessing=dict(mode='new',project='project',runtime_config='runtime.json',
         expected_image_id='sha256:'+'a'*64,bids_filter={}))

def test_plan_never_queries_or_starts_runtime(tmp_path,monkeypatch):
    c=config(tmp_path,monkeypatch)
    monkeypatch.setattr(subprocess,'run',lambda *a,**k:pytest.fail('Docker probed'))
    monkeypatch.setattr(subprocess,'Popen',lambda *a,**k:pytest.fail('Process started'))
    p=build_plan(c);cmd=p['arguments']
    assert '--pull=never' in cmd and '--sloppy' not in cmd and '--low-mem' not in cmd
    assert cmd[cmd.index('--omp-nthreads')+1]=='4'
    assert '/work' in cmd and not (tmp_path/'work').exists()
    with pytest.raises(w.StageError,match='disabled'):launch_new(c)

@pytest.mark.parametrize('threads',[1,3,12])
def test_thread_change_rejected(tmp_path,monkeypatch,threads):
    c=config(tmp_path,monkeypatch);r=w.read(tmp_path/'runtime.json');r['omp_nthreads']=threads;w.write(tmp_path/'runtime.json',r)
    with pytest.raises(w.StageError,match='four'):build_plan(c)

def test_input_work_overlap_rejected(tmp_path,monkeypatch):
    c=config(tmp_path,monkeypatch);r=w.read(tmp_path/'runtime.json');r['work_dir']='project/work';w.write(tmp_path/'runtime.json',r)
    with pytest.raises(w.StageError,match='separate'):build_plan(c)

def test_source_and_installed_entry_have_same_main():
    import neo_mri_cli
    p=Path(__file__).resolve().parents[1]/'run.py'
    spec=importlib.util.spec_from_file_location('entry',p);entry=importlib.util.module_from_spec(spec);spec.loader.exec_module(entry)
    assert entry.main is neo_mri_cli.main

def test_anatomical_intersection_is_an_explicit_new_version(case):
    w.dispatch(case,'preprocess');baseline=w.dispatch(case,'glm')
    anatomy=w.derivatives(case)/'T1.nii.gz';im=nib.load(anatomy)
    mask=np.ones(im.shape,np.uint8);mask[0]=0
    anatomy_mask=w.derivatives(case)/'anatomy_mask.nii.gz';nib.Nifti1Image(mask,im.affine).to_filename(anatomy_mask)
    # A new preprocessing registration is required to introduce a new mask;
    # never edit or reuse the existing immutable receipt.
    changed=copy.deepcopy(case);changed['workspace']='results/intersection';changed['analysis'].update(engine_variant='anatomical_intersection_v1',anatomical_mask='anatomy_mask.nii.gz')
    w.dispatch(changed,'preprocess');new=w.dispatch(changed,'glm')
    summary=w.read(Path(new['output'])/'mask_summary.json')
    assert summary['removed_from_functional']==36 and summary['intersection_voxels']==180
    assert new['key']!=baseline['key'] and new['fmriprep_invocations']==0

@pytest.mark.parametrize('returncode',[0,1])
def test_launcher_records_actual_mocked_exit_without_starting_services(tmp_path,monkeypatch,returncode):
    import staged_preprocess as backend
    c=config(tmp_path,monkeypatch);r=w.read(tmp_path/'runtime.json')
    r['execution_enabled']=r['scientific_protocol']['execution_enabled']=True
    w.write(tmp_path/'runtime.json',r)
    root=tmp_path/'project';bids=root/'bids';t1=bids/r['scientific_protocol']['t1w'];t1.parent.mkdir(parents=True)
    nib.Nifti1Image(np.ones((4,4,4)),np.eye(4)).to_filename(t1)
    w.write(bids/'dataset_description.json',{'Name':'Synthetic launcher test','BIDSVersion':'1.9.0'})
    w.write(bids/'neo_bids_filter.json',{})
    templates=tmp_path/'templates';templates.mkdir();(templates/'synthetic-cache-marker').write_text('not real templates')
    license_path=tmp_path/'secrets/license.txt';license_path.parent.mkdir();license_path.write_text('SYNTHETIC PLACEHOLDER; NOT A VALID LICENSE')
    monkeypatch.setattr(backend,'local_linux',lambda:True)
    calls=[]
    def inspect(command,**kwargs):
        calls.append(command);assert command[:3]==['docker','image','inspect']
        return subprocess.CompletedProcess(command,0,c['preprocessing']['expected_image_id']+'\n','')
    class FakeProcess:
        pid=99999
        def __init__(self,command,**kwargs):
            calls.append(command);assert command[:2]==['docker','run']
        def wait(self):return returncode
    monkeypatch.setattr(subprocess,'run',inspect);monkeypatch.setattr(subprocess,'Popen',FakeProcess)
    if returncode:
        with pytest.raises(w.StageError,match='failed'):launch_new(c)
    else:launch_new(c)
    assert w.read(root/'logs/preprocessing_status.json')['exit_code']==returncode
    assert len(calls)==3 and (root/'logs/fmriprep_run_configuration.json').is_file()
    assert all('wsl.exe' not in command and 'systemctl' not in command and command[:2]!=['docker','pull'] for command in calls)
    assert '--pull=never' in calls[-1]
    with pytest.raises(w.StageError,match='Fresh project'):launch_new(c)

def test_explicit_subject_output_root_and_source_protection(case):
    changed=copy.deepcopy(case);changed['workspace']='data/demo/demo_s01/results/task_fMRI/v1'
    path=w.REPO/'subject.json';w.write(path,changed)
    with pytest.raises(w.StageError,match='approved'):w.load(path)
    changed['output_roots']=['data/demo/demo_s01/results'];w.write(path,changed)
    assert w.load(path)['workspace']==changed['workspace']
    changed['workspace']='source/unsafe';changed['output_roots']=['source'];w.write(path,changed)
    with pytest.raises(w.StageError,match='separate'):w.load(path)

def test_viewer_refuses_grid_mismatch_and_writes_anatomy(tmp_path):
    from image_viewer import export
    data=np.indices((9,9,9)).sum(axis=0).astype(np.float32)
    bg=tmp_path/'t1.nii.gz';nib.Nifti1Image(data,np.eye(4)).to_filename(bg)
    bad=tmp_path/'bad.nii.gz';aff=np.eye(4);aff[0,3]=1;nib.Nifti1Image(data,aff).to_filename(bad)
    with pytest.raises(ValueError,match='grid'):export(str(bg),tmp_path/'bad.html',str(bad),'native T1')
    p=export(str(bg),tmp_path/'anatomy.html',space='synthetic native T1')
    assert 'synthetic native T1' in p.read_text(encoding='utf-8')
    with pytest.raises(FileExistsError):export(str(bg),p)
