"""Independent manual lesion and offline SynthStroke entry point."""
from pathlib import Path
import argparse
import json

from lesion_state import StageError


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage', choices=['prepare', 'edit', 'import', 'review', 'export', 'status', 'synthstroke'])
    parser.add_argument('--config', required=True)
    parser.add_argument('--root', type=Path, help='Base of relative paths; default configuration file directory')
    parser.add_argument('--mask', type=Path)
    parser.add_argument('--label', type=int)
    parser.add_argument('--reorient-to-reference', action='store_true')
    parser.add_argument('--version')
    parser.add_argument('--reviewer')
    parser.add_argument('--review-note')
    parser.add_argument('--confirm-empty', action='store_true')
    parser.add_argument('--purpose', choices=['analysis'])
    parser.add_argument('--export-dir', type=Path)
    parser.add_argument('--editor-exe', type=Path)
    parser.add_argument('--dry-run', action='store_true', help='Only edit and synthstroke support planning without execution')
    args = parser.parse_args(argv)
    try:
        if args.dry_run and args.stage not in ['edit', 'synthstroke']:
            raise StageError(f'{args.stage} does not support --dry-run; no operation was performed')
        root = args.root.resolve() if args.root else Path(args.config).resolve().parent
        config = (root / args.config).resolve() if args.root else Path(args.config).resolve()
        if args.stage == 'synthstroke':
            import lesion_synthstroke
            result = lesion_synthstroke.run(lesion_synthstroke.load_config(config, root), args.dry_run)
        else:
            import lesion_workflow
            for name in ['mask', 'export_dir', 'editor_exe']:
                value = getattr(args, name)
                if value is not None:
                    setattr(args, name, (root / value).resolve())
            result = lesion_workflow.dispatch(lesion_workflow.load_config(config, root), 'lesion-' + args.stage, args)
    except (StageError, OSError, ValueError, KeyError) as exc:
        print(json.dumps(dict(status='STOPPED', reason=str(exc)), ensure_ascii=False))
        return 2
    print(json.dumps(result, ensure_ascii=False))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
