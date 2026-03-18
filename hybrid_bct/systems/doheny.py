from pathlib import Path
from .base import BaseSystem


class DohenySystem(BaseSystem):
    def __init__(self, cfg: dict):
        super().__init__(cfg)

    def validate(self) -> None:
        files = self.cfg.get("files", {})
        required = ["spectrum", "detector_mtf", "scanlog"]
        for key in required:
            if key not in files:
                raise ValueError(f"Missing files.{key} in config")

        repo_root = Path(self.cfg["_repo_root"])
        for key in required:
            p = Path(files[key])
            if not p.is_absolute():
                p = repo_root / p
            if not p.exists():
                raise FileNotFoundError(f"Missing {key}: {p}")

    def summary(self) -> dict:
        return {
            "system": "doheny",
            "spectrum": self.cfg["files"]["spectrum"],
            "detector_mtf": self.cfg["files"]["detector_mtf"],
            "scanlog": self.cfg["files"]["scanlog"],
        }
