import argparse
from hybrid_bct.config import load_config
from hybrid_bct.systems.doheny import DohenySystem


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    args = parser.parse_args()

    cfg = load_config(args.config)

    if cfg.get("system", "").lower() == "doheny":
        system = DohenySystem(cfg)
    else:
        raise ValueError("Unsupported system")

    system.validate()
    print("Input validation passed.")
    print(system.summary())


if __name__ == "__main__":
    main()
