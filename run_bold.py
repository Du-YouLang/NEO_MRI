"""Task BOLD entry. Registers completed derivatives; never launches preprocessing."""
from pathlib import Path
import argparse,json,sys
sys.path.insert(0,str(Path(__file__).resolve().parent/'bold'))
from bold_workflow import load,dispatch,StageError,TimingBlocked

def main(argv=None):
    p=argparse.ArgumentParser(description='Register completed task fMRI derivatives, fit models and generate reports')
    p.add_argument('stage',choices=['preprocess','glm','report','status','audit','inventory','compare','demo'])
    p.add_argument('--config',help='Path to a completed BOLD v2 configuration')
    p.add_argument('--output',type=Path,help='New directory for generated inputs (demo only)')
    p.add_argument('--verify-full',action='store_true',help='Rehash every referenced input/output, ignoring cached file stamps')
    a=p.parse_args(argv)
    if a.stage=='demo':
        if a.output is None or a.config is not None or a.verify_full:
            p.error('demo requires --output and accepts neither --config nor --verify-full')
        from bold_demo import make_dataset

        output = a.output.resolve()
        try:
            output.mkdir(parents=True)
            config = make_dataset(output)
        except (OSError,ValueError) as exc:
            print(json.dumps(dict(status='STOPPED',reason=str(exc))))
            return 2
        print(json.dumps(dict(status='SYNTHETIC_INPUTS_CREATED',config=str(config))))
        return 0
    if a.config is None or a.output is not None:
        p.error('this stage requires --config and does not accept --output')
    try:result=dispatch(load(a.config),a.stage,a.verify_full)
    except TimingBlocked as exc:print(json.dumps(dict(status='TASK_GLM_PAUSED',reason=str(exc)),ensure_ascii=False));return 3
    except (StageError,OSError,ValueError,KeyError) as exc:print(json.dumps(dict(status='STOPPED',reason=str(exc)),ensure_ascii=False));return 2
    print(json.dumps(result,ensure_ascii=False));return 0
if __name__=='__main__':raise SystemExit(main())
