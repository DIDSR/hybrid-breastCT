import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from hybrid_bct.config import load_config
from hybrid_bct.systems.doheny import DohenySystem


def main():
    parser = argparse.ArgumentParser(description="Validate hybrid-bct inputs")
    parser.add_argument("--config", required=True, help="Path to YAML config file")
    args = parser.parse_args()

    cfg = load_config(args.config)

    if cfg.get("system", "").lower() == "doheny":
        system = DohenySystem(cfg)
    else:
        raise ValueError("Unsupported system")

    system.validate()

    repo_root = Path(cfg["_repo_root"])
    files_cfg = cfg["files"]
    paths_cfg = cfg["paths"]

    print("Basic config validation passed.")
    print(system.summary())

    # Check scanlog can be opened
    scanlog_path = Path(files_cfg["scanlog"])
    if not scanlog_path.is_absolute():
        scanlog_path = repo_root / scanlog_path
    df = pd.read_excel(scanlog_path, sheet_name="log")
    print(f"Scanlog loaded successfully: {scanlog_path}")
    print(f"Scanlog rows: {len(df)}")

    # Check detector MTF can be opened
    detector_mtf_path = Path(files_cfg["detector_mtf"])
    if not detector_mtf_path.is_absolute():
        detector_mtf_path = repo_root / detector_mtf_path

    try:
        mtf_data = np.loadtxt(detector_mtf_path, delimiter=",")
    except Exception:
        mtf_data = np.loadtxt(detector_mtf_path)

    print(f"Detector MTF loaded successfully: {detector_mtf_path}")
    print(f"Detector MTF shape: {mtf_data.shape}")

    # Check VOI center directory exists
    voi_center_dir = Path(paths_cfg["voi_center_dir"]).expanduser()
    if not voi_center_dir.exists():
        raise FileNotFoundError(f"VOI center directory not found: {voi_center_dir}")
    print(f"VOI center directory exists: {voi_center_dir}")

    print("Input validation passed.")


if __name__ == "__main__":
    main()
