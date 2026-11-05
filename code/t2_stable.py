"""Frozen T2 structural v1 interface; original validated engine is unchanged."""
from pathlib import Path
import argparse,copy,hashlib,json

# ENGINE_HASHES is populated from the reviewed source, using UTF-8/LF text.
ENGINE_HASHES = {'structural_multimodal.py': 'bb54e616bd1028dd55f3a7a928ffba820ac178c6809f753d8fa39a80d4b7282b', 'anatomy_review.py': '7f2a1d4ab0634703d95d1b58ab6a9a9936d167fd16fabe601998c423443b24f7', 'hdbet_cpu.py': '1a01b11dea3682a2530f6022220a0420ebf0fcb40e4116b394412e33b6d60177', 'multimodal_review.py': '0966a6acc8881e0c0cf77d0fe79a028c01b14c2a6512671318edf2ac5a1702be', 'anatomy_viewer_assets.py': '910173a1dfa6ead91b004e27d404757f8c3acc970295a61e1ee7e785b73a241f', 'state_io.py': 'a11b580cc0b99f417c092184bdee5c4d584da7e7e7df881d060615539e32b5bf'}

PROFILE=dict(id='T2_structural_v1',approval_date='2026-11-04',modality='T2w',threads=4,
    N4_zero_policy='exclude_zero_up_to_1pct',registration='SimpleITK 2.5.6 Mattes MI rigid 6DOF',
    intensity_options=['raw','n4'],brain_domain_options=['approved_T1_domain','independent_HDBET_additional_domain'],
    default_brain_domain='hd_bet',optional_brain_domain='t1_domain',automatic_fallback=False,
    automatic_method_winner=False,new_input_review_required=True,clinical_validation=False,
    resampling='Composed pull, linear intensity / nearest-neighbor labels; two-step output retained for comparison',
    no_new_T1_MNI_fit=True,no_new_BOLD_GLM=True)

def source_digest(path):
    return hashlib.sha256(Path(path).read_text(encoding='utf-8').replace('\r\n','\n').encode('utf-8')).hexdigest()

def verify_engine(directory=None):
    root=Path(directory) if directory else Path(__file__).resolve().parent
    if not ENGINE_HASHES:raise ValueError('Frozen engine manifest missing')
    for name,expected in ENGINE_HASHES.items():
        if source_digest(root/name)!=expected:raise ValueError('Frozen T2 v1 source changed: '+name+'; validate and release a new version')
    return PROFILE['id']

def prepare_config(config):
    c=copy.deepcopy(config)
    if c.get('modality')!='T2w':raise ValueError('T2 v1 approval is for T2w; FLAIR remains a separate reviewed route')
    mode=c.get('brain_domain','hd_bet')
    if mode not in ['hd_bet','t1_domain']:raise ValueError('Explicit hd_bet or t1_domain brain_domain required')
    if mode=='hd_bet' and not c.get('source_mask'):raise ValueError('Default HD-BET route requires independent source_mask; no automatic T1 fallback')
    if mode=='t1_domain':c.pop('source_mask',None)
    parameters=c.setdefault('parameters',{})
    parameters.setdefault('threads',4);parameters.setdefault('N4_zero_policy',PROFILE['N4_zero_policy'])
    if parameters['threads']!=4 or parameters['N4_zero_policy']!=PROFILE['N4_zero_policy']:
        raise ValueError('Frozen T2 v1 parameter mismatch; use the exploratory route for different settings')
    return c

def selection(config):
    mode=config.get('brain_domain','hd_bet');suffix='independent_HDBET' if mode=='hd_bet' else 'T1_only'
    return dict(brain_domain=mode,primary_branches=['raw_'+suffix,'n4_'+suffix],automatic_fallback=False,
        note='All registrations retain the approved fixed T1 brain domain; HD-BET adds the moving-image domain. Comparison outputs remain available in the HD-BET route.')

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=['policy','plan','run']);p.add_argument('--config');a=p.parse_args(argv)
    verify_engine()
    if a.action=='policy':result=dict(PROFILE,source_hash_algorithm='sha256 UTF-8/LF',engine=ENGINE_HASHES)
    else:
        if not a.config:p.error('--config required for plan/run')
        from structural_multimodal import plan,run
        c=prepare_config(json.loads(Path(a.config).read_text(encoding='utf-8')))
        result=plan(c) if a.action=='plan' else run(c)
        result=dict(result,T2_v1_selection=selection(c))
    print(json.dumps(result,ensure_ascii=False,indent=2));return 0
