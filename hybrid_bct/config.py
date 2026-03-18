from pathlib import Path
import yaml


def load_config(config_path: str | Path) -> dict:
    config_path = Path(config_path).resolve()
    with open(config_path, "r") as f:
        cfg = yaml.safe_load(f)

    cfg["_config_path"] = config_path
    cfg["_repo_root"] = config_path.parent.parent
    return cfg
