from pathlib import Path
import sys,pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'code'))
from t2_stable import prepare_config,verify_engine,ENGINE_HASHES,main,selection

def test_frozen_engine_identity_and_readonly_policy():
    assert verify_engine()=='T2_structural_v1'
    assert main(['policy'])==0

def test_different_source_cannot_silently_become_frozen_v1(tmp_path):
    import shutil
    root=Path(__file__).resolve().parents[1]/'code'
    for name in ENGINE_HASHES:shutil.copy2(root/name,tmp_path/name)
    p=tmp_path/'structural_multimodal.py';p.write_text(p.read_text(encoding='utf-8')+'\n# unvalidated change\n',encoding='utf-8')
    with pytest.raises(ValueError,match='source changed'):verify_engine(tmp_path)

def test_frozen_defaults_are_explicit_and_never_mutate_caller():
    c=dict(modality='T2w',source_mask='independent_mask.nii.gz',parameters={});r=prepare_config(c)
    assert c['parameters']=={} and r['parameters']==dict(threads=4,N4_zero_policy='exclude_zero_up_to_1pct')
    with pytest.raises(ValueError,match='parameter mismatch'):prepare_config(dict(modality='T2w',source_mask='mask.nii.gz',parameters=dict(threads=8)))
    with pytest.raises(ValueError,match='FLAIR'):prepare_config(dict(modality='FLAIR'))

def test_hdbet_is_default_and_missing_mask_never_falls_back():
    with pytest.raises(ValueError,match='no automatic T1 fallback'):prepare_config(dict(modality='T2w'))
    c=prepare_config(dict(modality='T2w',source_mask='independent_mask.nii.gz'))
    assert selection(c)['primary_branches']==['raw_independent_HDBET','n4_independent_HDBET']
    original=dict(modality='T2w',brain_domain='t1_domain',source_mask='unused.nii.gz')
    optional=prepare_config(original)
    assert 'source_mask' not in optional and original['source_mask']=='unused.nii.gz'
    assert selection(optional)['primary_branches']==['raw_T1_only','n4_T1_only']
    with pytest.raises(ValueError,match='brain_domain'):prepare_config(dict(modality='T2w',brain_domain='guess'))
