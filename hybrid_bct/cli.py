import argparse

from .config import load_config
from .pipeline import run_hybrid_simulation
from .systems.doheny import DohenySystem


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Hybrid breast CT virtual imaging trial framework"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    run_parser = subparsers.add_parser("run", help="Run hybrid simulation")
    run_parser.add_argument("--config", required=True, help="Path to YAML config file")
    run_parser.add_argument("--scan-id", type=int, required=True, help="Patient scan ID")
    run_parser.add_argument(
        "--cluster-diameter-mm",
        type=float,
        required=True,
        help="Cluster diameter in mm",
    )
    run_parser.add_argument(
        "--num-calcs",
        type=int,
        required=True,
        help="Number of calcifications in the cluster",
    )
    run_parser.add_argument(
        "--calc-diameter-mm",
        type=float,
        required=True,
        help="Calcification diameter in mm",
    )

    validate_parser = subparsers.add_parser("validate", help="Validate inputs only")
    validate_parser.add_argument(
        "--config", required=True, help="Path to YAML config file"
    )

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

    elif args.command == "validate":
        cfg = load_config(args.config)

        if cfg.get("system", "").lower() == "doheny":
            system = DohenySystem(cfg)
        else:
            raise ValueError(f"Unsupported system: {cfg.get('system', '')}")

        system.validate()
        print("Input validation passed.")
        print(system.summary())
        
if __name__ == "__main__":
    main()
