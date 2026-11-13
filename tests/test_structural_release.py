from pathlib import Path
import sys,json,shutil,pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'code'))
from test_anatomy_review import anatomy,approve
from anatomy_review import export_precomputed
from structural_context import load,StageError
from t1_preprocess import build_plan

def test_t1_plan_is_anatomy_only_and_never_launches(anatomy,monkeypatch):
    import subprocess
    monkeypatch.setattr(subprocess,'Popen',lambda *a,**k:pytest.fail('Plan launched a process'))
    a=anatomy;approve(a);root=a['root'];project=root/'project';raw=project/'bids/sub-demo/anat/sub-demo_T1w.nii.gz';raw.parent.mkdir(parents=True);shutil.copy2(a['t1'],raw)
    (project/'bids/.bidsignore').write_text('neo_bids_filter.json\n')
    exported=export_precomputed(a['review'],raw,root/'reviewed')
    example=Path(__file__).resolve().parents[1]/'configs/examples/t1_runtime.example.json'
    runtime=json.loads(example.read_text(encoding='utf-8'));runtime.update(work_dir=str(root/'work'),license_file=str(root/'license.txt'),templateflow_dir=str(root/'templates'))
    (root/'runtime.json').write_text(json.dumps(runtime))
    c=dict(schema='neo-t1-anatomy-1',anatomy_review=str(a['review']),preprocessing=dict(project=str(project),runtime_config=str(root/'runtime.json'),reviewed_anatomy_dir=str(exported),expected_image_id='sha256:'+'a'*64,bids_filter={}))
    cfg=root/'case.json';cfg.write_text(json.dumps(c));plan=build_plan(load(cfg))
    assert '--anat-only' in plan['arguments'] and '--fs-no-reconall' in plan['arguments']
    assert not plan['new_backend_real_data_verified']
    runtime['anat_only']=False;(root/'runtime.json').write_text(json.dumps(runtime))
    with pytest.raises(StageError,match='anat_only'):build_plan(c)
    runtime['anat_only']=True
    shutil.copy2(raw,project/'external_T1w.nii.gz')
    runtime['scientific_protocol']['t1w']='../external_T1w.nii.gz'
    (root/'runtime.json').write_text(json.dumps(runtime))
    with pytest.raises(StageError,match='outside BIDS'):build_plan(c)

def test_release_has_no_functional_entry_or_module():
    import tomllib
    root=Path(__file__).resolve().parents[1];p=tomllib.loads((root/'pyproject.toml').read_text(encoding='utf-8'))
    modules=p['tool']['setuptools']['py-modules'];scripts=p['project']['scripts']
    assert scripts['neo-mri']=='structural_cli:main'
    assert not any(x in modules for x in ['staged_glm','staged_glm_anatomical','staged_workflow','staged_report'])
    assert 'neo-mri-demo' not in scripts
