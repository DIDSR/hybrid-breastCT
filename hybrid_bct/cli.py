import argparse
from .config import load_config
from .pipeline import run_hybrid_simulation


def main() -> None:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)

    run_parser = subparsers.add_parser("run")
    run_parser.add_argument("--config", required=True)
    run_parser.add_argument("--scan-id", type=int, required=True)
    run_parser.add_argument("--cluster-diameter-mm", type=float, required=True)
    run_parser.add_argument("--num-calcs", type=int, required=True)
    run_parser.add_argument("--calc-diameter-mm", type=float, required=True)

    args = parser.parse_args()

    if args.command == "run":
        cfg = load_config(args.config)
        run_hybrid_simulation(
            cfg=cfg,
            scan_id=args.scan_id,
            cluster_diameter_mm=args.cluster_diameter_mm,
            num_calcs=args.num_calcs,
            calc_diameter_mm=args.calc_diameter_mm,
        )
