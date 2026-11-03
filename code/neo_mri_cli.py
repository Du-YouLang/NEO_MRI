"""Installed and source-tree entry; no implicit execution of other stages."""
from pathlib import Path
import argparse, json, sys
import staged_workflow as workflow
from staged_workflow import load, dispatch, StageError, TimingBlocked

def main(argv=None):
    parser=argparse.ArgumentParser(description='NEO MRI research preview: independent stages')
    parser.add_argument('stage',choices=['preprocess','glm','report','status','plan','audit','container-plan'])
    parser.add_argument('--config',required=True)
    parser.add_argument('--root',type=Path,default=None,help='Base of all relative data paths; default current directory')
    parser.add_argument('--verify-full',action='store_true')
    args=parser.parse_args(argv)
    if args.root is not None:workflow.REPO=args.root.resolve()
    try:
        if args.stage=='plan':
            from scientific_protocol import plan
            result=plan(workflow.read(workflow.resolve(args.config)))
        elif args.stage=='container-plan':
            from staged_preprocess import build_plan
            result=build_plan(load(args.config))
        else:result=dispatch(load(args.config),args.stage,args.verify_full)
    except TimingBlocked as exc:
        print(json.dumps(dict(status='TASK_GLM_PAUSED',reason=str(exc))));return 3
    except (StageError,OSError,ValueError,KeyError) as exc:
        print(json.dumps(dict(status='STOPPED',reason=str(exc))));return 2
    print(json.dumps(result,ensure_ascii=False));return 0

if __name__=='__main__':raise SystemExit(main())
