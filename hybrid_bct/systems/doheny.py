from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import tigre

from .base import BaseSystem


def _resolve_path(path_str: str | Path, base_dir: str | Path | None = None) -> Path:
    path = Path(path_str).expanduser()
    if not path.is_absolute() and base_dir is not None:
        path = Path(base_dir) / path
    return path.resolve()


def _is_float(value: str) -> bool:
    try:
        float(value)
        return True
    except Exception:
        return False


class DohenySystem(BaseSystem):
    def __init__(self, cfg: dict):
        super().__init__(cfg)

    def validate(self) -> None:
        files = self.cfg.get("files", {})
        paths = self.cfg.get("paths", {})
    
        required_files = ["spectrum", "detector_mtf", "scanlog"]
        for key in required_files:
            if key not in files:
                raise ValueError(f"Missing files.{key} in config")
    
        required_paths = ["patient_data_dir", "output_dir", "voi_center_dir"]
        for key in required_paths:
            if key not in paths:
                raise ValueError(f"Missing paths.{key} in config")
    
        repo_root = Path(self.cfg["_repo_root"])
    
        for key in required_files:
            p = Path(files[key])
            if not p.is_absolute():
                p = repo_root / p
            if not p.exists():
                raise FileNotFoundError(f"Missing {key}: {p}")
    
        material_files = files.get("material_files", {})
        required_materials = ["calc", "adipose", "glandular", "csI"]
    
        for key in required_materials:
            if key not in material_files:
                raise ValueError(f"Missing files.material_files.{key} in config")
    
            p = Path(material_files[key])
            if not p.is_absolute():
                p = repo_root / p
            if not p.exists():
                raise FileNotFoundError(f"Missing material file {key}: {p}")

    def summary(self) -> dict:
        repo_root = self.cfg.get("_repo_root", Path.cwd())
        return {
            "system": "doheny",
            "patient_data_dir": str(Path(self.cfg["paths"]["patient_data_dir"]).expanduser()),
            "output_dir": str(Path(self.cfg["paths"]["output_dir"]).expanduser()),
            "voi_center_dir": str(Path(self.cfg["paths"]["voi_center_dir"]).expanduser()),
            "spectrum": str(_resolve_path(self.cfg["files"]["spectrum"], repo_root)),
            "detector_mtf": str(_resolve_path(self.cfg["files"]["detector_mtf"], repo_root)),
            "scanlog": str(_resolve_path(self.cfg["files"]["scanlog"], repo_root)),
        }


def load_scanlog_doheny(cfg: dict) -> tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    """
    Load the Doheny example scan log.

    Expected config keys:
      cfg["files"]["scanlog"]
      cfg["_repo_root"]

    Returns
    -------
    df : pd.DataFrame
        Raw scanlog DataFrame.
    VscanID : np.ndarray
        Scan ID column.
    VcropLocs : np.ndarray
        Crop coordinates from columns C:H in the old spreadsheet layout.
    """
    repo_root = cfg.get("_repo_root", Path.cwd())
    scanlog_path = _resolve_path(cfg["files"]["scanlog"], repo_root)

    df = pd.read_excel(scanlog_path, sheet_name="log")
    VscanID = df.iloc[:, 1].values
    VcropLocs = df.iloc[:, 2:8].values

    return df, VscanID, VcropLocs


def get_scanlog_index_doheny(scan_id: int, VscanID: np.ndarray) -> int:
    matches = np.where(VscanID == scan_id)[0]
    if len(matches) == 0:
        raise ValueError(f"scanID {scan_id} not found in Doheny scan log.")
    return int(matches[0])


def load_seg_doheny(
    scan_id: int,
    patient_data_dir: str | Path,
) -> tuple[np.ndarray, list[int], float, dict[str, int]]:
    """
    Load Doheny segmentation volume.

    This intentionally preserves the old Doheny folder layout:
      patient_data_dir/Doheny_CT/CTD####/CTi_Corrected/...
      patient_data_dir/Doheny_CT/CTD####/SEG/...

    Returns
    -------
    seg_volume : np.ndarray
    Nxyz : list[int]
    original_vx_um : float
    labels : dict[str, int]
    """
    root_dir = Path(patient_data_dir).expanduser().resolve()
    print(f"Loading seg {scan_id:04d}... ", end="", flush=True)

    header_file = root_dir / f"Doheny_CT/CTD{scan_id:04d}/CTi_Corrected/CT{scan_id:04d}_01.0001"
    if not header_file.exists():
        raise FileNotFoundError(f"Could not find reconstructed header file: {header_file}")

    with open(header_file, "rb") as file:
        header = np.fromfile(file, dtype="float32", count=20)

    Nz = int(header[14])
    matrixsize = int(header[0])
    original_vx_um = float(header[9]) * 1000.0

    Nxyz = [Nz, matrixsize, matrixsize]
    seg_volume = np.zeros(Nxyz, dtype=np.uint8)

    for islice in range(1, Nz + 1):
        slice_file = root_dir / f"Doheny_CT/CTD{scan_id:04d}/SEG/SEG{scan_id:04d}_01.{islice:04d}"
        if not slice_file.exists():
            raise FileNotFoundError(f"Missing segmentation slice: {slice_file}")

        with open(slice_file, "rb") as file:
            imslice = np.fromfile(file, dtype="uint8", count=matrixsize * matrixsize)
            imslice = imslice.reshape(matrixsize, matrixsize)
            seg_volume[islice - 1, :, :] = imslice

    # Preserve old remapping behavior
    seg_volume[seg_volume == 4] = 1
    seg_volume[seg_volume == 3] = 2

    labels = {
        "air": 0,
        "adipose": 1,
        "glandular": 2,
        "skin": 5,
    }

    print("done.", flush=True)
    return seg_volume, Nxyz, original_vx_um, labels


def load_projections_and_geometry_doheny(
    scan_id: int,
    patient_data_dir: str | Path,
    xsiz: int,
    ysiz: int,
    dexel_mm: float,
    vertical_offset_mm: float,
) -> tuple[np.ndarray, Any, np.ndarray]:
    """
    Load Doheny projection data and TIGRE geometry.

    This is the Doheny example implementation and intentionally keeps the old
    directory assumptions. Adjust this function for other systems.

    Expected folders under patient_data_dir:
      Doheny_CT/CTD####/SIN/
      Doheny_CT/CTD####/XXX/
      Doheny_CT/CTD####/CAL/
    """
    root_dir = Path(patient_data_dir).expanduser().resolve()
    print(f"\nLoading metadata and projection data for scan {scan_id}", flush=True)

    xxx_file = root_dir / f"Doheny_CT/CTD{scan_id:04d}/XXX/CT{scan_id:04d}.xxx"
    if not xxx_file.exists():
        raise FileNotFoundError(f"Could not find XXX geometry file: {xxx_file}")

    with open(xxx_file, "r") as file:
        lines = file.readlines()

    dsd_mm = None
    dso_mm = None
    num_proj = None
    angles_deg = []

    for line in lines:
        lower = line.lower().strip()

        if "source to detector" in lower or "sdd" in lower or "dsd" in lower:
            tokens = lower.replace("=", " ").replace(",", " ").split()
            vals = [t for t in tokens if _is_float(t)]
            if vals:
                dsd_mm = float(vals[-1])

        elif "source to object" in lower or "sod" in lower or "dso" in lower:
            tokens = lower.replace("=", " ").replace(",", " ").split()
            vals = [t for t in tokens if _is_float(t)]
            if vals:
                dso_mm = float(vals[-1])

        elif "number of projections" in lower or "num projections" in lower:
            tokens = lower.replace("=", " ").replace(",", " ").split()
            vals = [t for t in tokens if t.isdigit()]
            if vals:
                num_proj = int(vals[-1])

        elif "angle" in lower:
            tokens = lower.replace("=", " ").replace(",", " ").split()
            vals = [t for t in tokens if _is_float(t)]
            if vals:
                angles_deg.append(float(vals[-1]))

    sin_dir = root_dir / f"Doheny_CT/CTD{scan_id:04d}/SIN"
    if not sin_dir.exists():
        raise FileNotFoundError(f"Could not find SIN projection directory: {sin_dir}")

    sin_files = sorted(sin_dir.glob(f"CT{scan_id:04d}_01.*"))
    if len(sin_files) == 0:
        raise FileNotFoundError(f"No SIN projection files found in: {sin_dir}")

    if num_proj is None:
        num_proj = len(sin_files)

    prjstack = np.zeros((num_proj, ysiz, xsiz), dtype=np.float32)

    for iproj, proj_file in enumerate(sin_files[:num_proj]):
        with open(proj_file, "rb") as file:
            proj = np.fromfile(file, dtype="float32", count=xsiz * ysiz)
            proj = proj.reshape(ysiz, xsiz)
            prjstack[iproj, :, :] = proj

    if len(angles_deg) >= num_proj:
        ang = np.deg2rad(np.asarray(angles_deg[:num_proj], dtype=np.float32))
    else:
        ang = np.linspace(0, 2 * np.pi, num_proj, endpoint=False).astype(np.float32)

    if dsd_mm is None:
        dsd_mm = 878.0
    if dso_mm is None:
        dso_mm = 511.0

    geo = tigre.geometry()
    geo.DSD = float(dsd_mm)
    geo.DSO = float(dso_mm)

    geo.nDetector = np.array((ysiz, xsiz), dtype=np.int32)
    geo.dDetector = np.array((dexel_mm, dexel_mm), dtype=np.float32)
    geo.sDetector = geo.nDetector * geo.dDetector

    # These defaults are typically overwritten later
    geo.nVoxel = np.array((512, 512, 512), dtype=np.int32)
    geo.dVoxel = np.array((0.1, 0.1, 0.1), dtype=np.float32)
    geo.sVoxel = geo.nVoxel * geo.dVoxel

    geo.offOrigin = np.array((0.0, 0.0, 0.0), dtype=np.float32)
    geo.offDetector = np.array((vertical_offset_mm, 0.0), dtype=np.float32)
    geo.rotDetector = np.array((0.0, 0.0, 0.0), dtype=np.float32)
    geo.mode = "cone"

    print("done.", flush=True)
    return prjstack, geo, ang


def summarize_doheny_inputs(cfg: dict) -> dict[str, str]:
    repo_root = cfg.get("_repo_root", Path.cwd())
    patient_data_dir = cfg["paths"]["patient_data_dir"]

    return {
        "system": "doheny",
        "patient_data_dir": str(Path(patient_data_dir).expanduser()),
        "scanlog": str(_resolve_path(cfg["files"]["scanlog"], repo_root)),
        "spectrum": str(_resolve_path(cfg["files"]["spectrum"], repo_root)),
        "detector_mtf": str(_resolve_path(cfg["files"]["detector_mtf"], repo_root)),
    }
