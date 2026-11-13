"""CST/CRP commands for a checkout of NEO_MRI."""
import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "code"))

from tractography import workflow
from tractography.geometry import mapped_label


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("demo", "fit-plan", "fit", "track-plan", "track",
                                          "map-label", "audit-label", "screen", "export", "report"))
    parser.add_argument("--config", type=Path)
    parser.add_argument("--results", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--run")
    parser.add_argument("--gate")
    args = parser.parse_args()
    if args.stage in {"fit-plan", "fit", "track-plan", "track", "screen", "map-label", "audit-label"} and not args.config:
        parser.error("--config is required")
    if args.stage in {"export", "report"} and not args.results:
        parser.error("--results is required")
    if args.stage == "demo":
        from tractography.demo import demo
        result = demo(args.output)
    elif args.stage in {"map-label", "audit-label"}:
        config = args.config.resolve(strict=True)
        result = mapped_label(workflow.read(config), config.parent, args.output,
                              technical_replay=args.stage == "audit-label")
    elif args.stage == "export":
        if not args.run or not args.gate:
            parser.error("export requires --run and --gate")
        result = workflow.export(args.results, args.run, args.gate, args.output)
    elif args.stage == "report":
        result = workflow.report(args.results, args.output)
    elif args.stage == "screen":
        result = workflow.screen(args.config, args.output)
        result = {"results": str(args.output.resolve()), "runs": len(result["runs"])}
    else:
        phases = {"fit-plan": workflow.fit_plan, "fit": workflow.fit,
                  "track-plan": workflow.tracking_plan, "track": workflow.track}
        result = phases[args.stage](args.config, args.output.resolve())
    print(json.dumps(result, ensure_ascii=True))


if __name__ == "__main__":
    main()
