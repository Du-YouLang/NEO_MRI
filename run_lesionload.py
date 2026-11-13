"""CLI for the SMATT, HCPA and Boyne 2025 research workflows."""

import argparse
import json
from pathlib import Path
import sys


sys.path.insert(0, str(Path(__file__).resolve().parent / "code"))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("demo", "run", "build-boyne"))
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--config", type=Path)
    parser.add_argument("--assets", type=Path)
    args = parser.parse_args()
    if args.stage == "run":
        if not args.config:
            parser.error("run requires --config")
        from lesionload.workflow import run
        result = run(args.config, args.output)
    elif args.stage == "build-boyne":
        if not args.assets:
            parser.error("build-boyne requires --assets")
        from lesionload.assets import build_boyne_assets
        result = build_boyne_assets(args.assets.resolve(), args.output.resolve())
    else:
        from lesionload.demo import demo
        result = demo(args.output)
    print(json.dumps(result, ensure_ascii=True, allow_nan=False))


if __name__ == "__main__":
    main()
