"""T1/T2 structural release entry; no DWI, BOLD, GLM or task-model commands."""
import argparse,json
def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('stage',choices=['t1-plan','t1-run','t2-policy','t2-plan','t2-run']);p.add_argument('--config');a=p.parse_args(argv)
    if a.stage.startswith('t2-'):
        from t2_stable import main as t2
        return t2([a.stage[3:]]+(['--config',a.config] if a.config else []))
    if not a.config:p.error('--config required')
    from structural_context import load
    from t1_preprocess import build_plan,launch_new
    c=load(a.config)
    if a.stage=='t1-plan':print(json.dumps(build_plan(c),ensure_ascii=False,indent=2))
    else:launch_new(c)
    return 0
