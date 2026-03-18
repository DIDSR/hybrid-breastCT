from __future__ import annotations

from pathlib import Path
from typing import Any
import json


def _format_float_list(values) -> str:
    return ", ".join([f"{float(v):.6g}" for v in values])


def _ensure_dir(path: str | Path) -> Path:
    p = Path(path).expanduser().resolve()
    p.mkdir(parents=True, exist_ok=True)
    return p


def build_patch_folder_stem(
    calc_diameter_mm: float,
    cluster_diameter_mm: float,
    num_calcs: int,
    calc_density: float,
    recon_alg: str,
    folder_suffix: str,
) -> str:
    """
    Match the old naming convention closely enough for compatibility.
    """
    return (
        f"{calc_diameter_mm:.2f}mmCalc_"
        f"{cluster_diameter_mm:.1f}mmCluster_"
        f"{num_calcs:02d}Calcs_"
        f"CaOx_{calc_density:1.2f}_"
        f"{recon_alg}_{folder_suffix}"
    )


def write_metadata(
    output_patch_root: str | Path,
    cfg: dict,
    run_params: dict,
    density: dict,
    geo: Any,
    ang,
    recon_alg: str,
    folder_suffix: str,
    kernel: str,
) -> tuple[Path, Path]:
    """
    Write metadata text files for both MIP and VOI patch directories.

    Expected run_params keys:
      scan_id
      effective_energy_keV
      energy_bin_size_keV
      csI_thickness_mm
      vertical_offset_mm
      new_vx_um
      calc_shape
      calc_diameter_mm
      num_calcs
      cluster_diameter_mm
      clusterCenteredFLAG
      num_SPvois_perbreast
      num_SAvois_perbreast
      voi_size_mm
      voi_size_vx
      flagHU
      mu_water

    Returns
    -------
    tuple[Path, Path]
        Paths to written metadata files for MIP and VOI directories.
    """
    output_patch_root = _ensure_dir(output_patch_root)
    repo_root = cfg.get("_repo_root", Path.cwd())

    files_cfg = cfg.get("files", {})
    paths_cfg = cfg.get("paths", {})

    spectrum_path = _resolve_optional(files_cfg.get("spectrum"), repo_root)
    detector_mtf_path = _resolve_optional(files_cfg.get("detector_mtf"), repo_root)
    scanlog_path = _resolve_optional(files_cfg.get("scanlog"), repo_root)

    material_files = files_cfg.get("material_files", {})

    stem = build_patch_folder_stem(
        calc_diameter_mm=run_params["calc_diameter_mm"],
        cluster_diameter_mm=run_params["cluster_diameter_mm"],
        num_calcs=run_params["num_calcs"],
        calc_density=density["calc"],
        recon_alg=recon_alg,
        folder_suffix=folder_suffix,
    )

    mip_dir = output_patch_root / f"{stem}_MIP"
    voi_dir = output_patch_root / f"{stem}_VOI"
    mip_dir.mkdir(parents=True, exist_ok=True)
    voi_dir.mkdir(parents=True, exist_ok=True)

    metadata_txt = _render_metadata_text(
        cfg=cfg,
        paths_cfg=paths_cfg,
        spectrum_path=spectrum_path,
        detector_mtf_path=detector_mtf_path,
        scanlog_path=scanlog_path,
        material_files=material_files,
        run_params=run_params,
        density=density,
        geo=geo,
        ang=ang,
        recon_alg=recon_alg,
        folder_suffix=folder_suffix,
        kernel=kernel,
    )

    metadata_json = _render_metadata_json(
        cfg=cfg,
        paths_cfg=paths_cfg,
        spectrum_path=spectrum_path,
        detector_mtf_path=detector_mtf_path,
        scanlog_path=scanlog_path,
        material_files=material_files,
        run_params=run_params,
        density=density,
        geo=geo,
        ang=ang,
        recon_alg=recon_alg,
        folder_suffix=folder_suffix,
        kernel=kernel,
    )

    mip_txt = mip_dir / "METADATA_ALLPATCHES.txt"
    voi_txt = voi_dir / "METADATA_ALLPATCHES.txt"
    mip_json = mip_dir / "metadata.json"
    voi_json = voi_dir / "metadata.json"

    mip_txt.write_text(metadata_txt)
    voi_txt.write_text(metadata_txt)
    mip_json.write_text(json.dumps(metadata_json, indent=2))
    voi_json.write_text(json.dumps(metadata_json, indent=2))

    return mip_txt, voi_txt


def _render_metadata_text(
    cfg: dict,
    paths_cfg: dict,
    spectrum_path: str | None,
    detector_mtf_path: str | None,
    scanlog_path: str | None,
    material_files: dict,
    run_params: dict,
    density: dict,
    geo: Any,
    ang,
    recon_alg: str,
    folder_suffix: str,
    kernel: str,
) -> str:
    lines = []
    lines.append("METADATA FOR ALL PATCHES IN THIS DIRECTORY:\n")

    lines.append("[INPUT FILES]")
    if spectrum_path is not None:
        lines.append(f"Energy spectrum: {spectrum_path}")
    if detector_mtf_path is not None:
        lines.append(f"Detector MTF: {detector_mtf_path}")
    if scanlog_path is not None:
        lines.append(f"Scan log: {scanlog_path}")

    if material_files:
        lines.append(f"Calc material file: {material_files.get('calc', '')}")
        lines.append(f"Adipose material file: {material_files.get('adipose', '')}")
        lines.append(f"Glandular material file: {material_files.get('glandular', '')}")
        lines.append(f"CsI material file: {material_files.get('csI', '')}")

    lines.append("")
    lines.append("[PATHS]")
    if "patient_data_dir" in paths_cfg:
        lines.append(f"Patient data directory: {paths_cfg['patient_data_dir']}")
    if "output_dir" in paths_cfg:
        lines.append(f"Output directory: {paths_cfg['output_dir']}")
    if "voi_center_dir" in paths_cfg:
        lines.append(f"VOI center directory: {paths_cfg['voi_center_dir']}")

    lines.append("")
    lines.append("[MATERIAL DENSITIES]")
    lines.append(f"Calc density: {density['calc']:1.4f}")
    lines.append(f"Adipose density: {density['adipose']:1.4f}")
    lines.append(f"Glandular density: {density['glandular']:1.4f}")
    lines.append(f"CsI density: {density['csI']:1.4f}")

    lines.append("")
    lines.append("[SIMULATION PARAMETERS]")
    lines.append(f"scanID: {run_params['scan_id']}")
    lines.append(f"Effective energy (keV): {run_params['effective_energy_keV']}")
    lines.append(f"Energy bin size (keV): {run_params['energy_bin_size_keV']}")
    lines.append(f"CsI thickness (mm): {run_params['csI_thickness_mm']}")
    lines.append(f"Additional vertical offset of object (mm): {run_params['vertical_offset_mm']}")
    lines.append(f"Object voxel size (um): {run_params['new_vx_um']}")
    lines.append(f"Calc shape: {run_params['calc_shape']}")
    lines.append(f"Calc diameter (mm): {run_params['calc_diameter_mm']}")
    lines.append(f"Number of calcs: {run_params['num_calcs']}")
    lines.append(f"Cluster diameter (mm): {run_params['cluster_diameter_mm']}")
    lines.append(f"Cluster centered in VOI flag: {run_params['clusterCenteredFLAG']}")
    lines.append(f"Signal-present VOIs requested per breast: {run_params['num_SPvois_perbreast']}")
    lines.append(f"Signal-absent VOIs requested per breast: {run_params['num_SAvois_perbreast']}")
    lines.append(f"VOI size (mm): {run_params['voi_size_mm']}")
    lines.append(f"VOI size (voxels): {run_params['voi_size_vx']}")

    lines.append("")
    lines.append("[RECONSTRUCTION PARAMETERS]")
    lines.append(f"Reconstruction algorithm: {recon_alg}")
    lines.append(f"Folder suffix: {folder_suffix}")
    lines.append(f"Kernel: {kernel}")
    lines.append(f"Output in HU flag: {run_params['flagHU']}")
    lines.append(f"mu_water: {run_params['mu_water']}")

    lines.append("")
    lines.append("[GEOMETRY]")
    lines.append(f"DSD (mm): {getattr(geo, 'DSD', 'NA')}")
    lines.append(f"DSO (mm): {getattr(geo, 'DSO', 'NA')}")
    lines.append(f"Detector pixels: {getattr(geo, 'nDetector', 'NA')}")
    lines.append(f"Detector pixel size (mm): {getattr(geo, 'dDetector', 'NA')}")
    lines.append(f"Detector size (mm): {getattr(geo, 'sDetector', 'NA')}")
    lines.append(f"Voxel counts: {getattr(geo, 'nVoxel', 'NA')}")
    lines.append(f"Voxel size (mm): {getattr(geo, 'dVoxel', 'NA')}")
    lines.append(f"Volume size (mm): {getattr(geo, 'sVoxel', 'NA')}")
    lines.append(f"Detector offsets: {getattr(geo, 'offDetector', 'NA')}")
    lines.append(f"Origin offsets: {getattr(geo, 'offOrigin', 'NA')}")
    lines.append(f"Number of projection angles: {len(ang)}")

    try:
        lines.append(f"Projection angles (rad): {_format_float_list(ang)}")
    except Exception:
        lines.append("Projection angles (rad): unavailable")

    lines.append("")
    lines.append("[CONFIG]")
    lines.append(f"System: {cfg.get('system', 'unknown')}")
    lines.append(f"Config path: {cfg.get('_config_path', 'unknown')}")

    return "\n".join(lines) + "\n"


def _render_metadata_json(
    cfg: dict,
    paths_cfg: dict,
    spectrum_path: str | None,
    detector_mtf_path: str | None,
    scanlog_path: str | None,
    material_files: dict,
    run_params: dict,
    density: dict,
    geo: Any,
    ang,
    recon_alg: str,
    folder_suffix: str,
    kernel: str,
) -> dict:
    return {
        "system": cfg.get("system", "unknown"),
        "config_path": str(cfg.get("_config_path", "unknown")),
        "input_files": {
            "spectrum": spectrum_path,
            "detector_mtf": detector_mtf_path,
            "scanlog": scanlog_path,
            "material_files": material_files,
        },
        "paths": {
            "patient_data_dir": paths_cfg.get("patient_data_dir"),
            "output_dir": paths_cfg.get("output_dir"),
            "voi_center_dir": paths_cfg.get("voi_center_dir"),
        },
        "density": density,
        "run_params": run_params,
        "reconstruction": {
            "algorithm": recon_alg,
            "folder_suffix": folder_suffix,
            "kernel": kernel,
        },
        "geometry": {
            "DSD": _safe_list_or_scalar(getattr(geo, "DSD", None)),
            "DSO": _safe_list_or_scalar(getattr(geo, "DSO", None)),
            "nDetector": _safe_list_or_scalar(getattr(geo, "nDetector", None)),
            "dDetector": _safe_list_or_scalar(getattr(geo, "dDetector", None)),
            "sDetector": _safe_list_or_scalar(getattr(geo, "sDetector", None)),
            "nVoxel": _safe_list_or_scalar(getattr(geo, "nVoxel", None)),
            "dVoxel": _safe_list_or_scalar(getattr(geo, "dVoxel", None)),
            "sVoxel": _safe_list_or_scalar(getattr(geo, "sVoxel", None)),
            "offDetector": _safe_list_or_scalar(getattr(geo, "offDetector", None)),
            "offOrigin": _safe_list_or_scalar(getattr(geo, "offOrigin", None)),
            "angles_rad": _safe_list_or_scalar(ang),
            "num_angles": len(ang),
        },
    }


def _resolve_optional(path_str: str | None, base_dir: str | Path | None = None) -> str | None:
    if path_str is None:
        return None
    path = Path(path_str).expanduser()
    if not path.is_absolute() and base_dir is not None:
        path = Path(base_dir) / path
    return str(path.resolve())


def _safe_list_or_scalar(value):
    if value is None:
        return None
    try:
        if hasattr(value, "tolist"):
            return value.tolist()
        return value
    except Exception:
        return str(value)
