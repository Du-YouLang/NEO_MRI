"""Task BOLD entry. Registers completed derivatives; never launches preprocessing."""
from pathlib import Path
import argparse,json,sys
sys.path.insert(0,str(Path(__file__).resolve().parent/'bold'))
from bold_workflow import load,dispatch,StageError,TimingBlocked

def main(argv=None):
    p=argparse.ArgumentParser(description='Independent MRI preprocessing / GLM / report; fixed four tool threads')
    p.add_argument('stage',choices=['preprocess','glm','report','status','audit','inventory','compare'])
    p.add_argument('--config',required=True,help='Path to a completed BOLD v2 configuration')
    p.add_argument('--verify-full',action='store_true',help='Rehash every referenced input/output, ignoring cached file stamps')
    a=p.parse_args(argv)
    try:result=dispatch(load(a.config),a.stage,a.verify_full)
    except TimingBlocked as exc:print(json.dumps(dict(status='TASK_GLM_PAUSED',reason=str(exc)),ensure_ascii=False));return 3
    except (StageError,OSError,ValueError,KeyError) as exc:print(json.dumps(dict(status='STOPPED',reason=str(exc)),ensure_ascii=False));return 2
    print(json.dumps(result,ensure_ascii=False));return 0
if __name__=='__main__':raise SystemExit(main())
