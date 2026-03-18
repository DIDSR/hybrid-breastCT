from __future__ import annotations

from pathlib import Path

from hybrid_bct.systems.doheny import (
    DohenySystem,
    load_scanlog_doheny,
    get_scanlog_index_doheny,
    load_seg_doheny,
    load_projections_and_geometry_doheny,
)
from hybrid_bct.io.spectra import read_energy_spectrum
from hybrid_bct.io.materials import read_material_file
from hybrid_bct.metadata import write_metadata
from hybrid_bct.simulation.calc_models import fxn_generate_calc
from hybrid_bct.simulation.insertion import fxn_insert_calc_cluster_new

from hybrid_bct.simulation.volume import (
    fxn_crop_volume,
    fxn_upsample_volume_in_sections,
)

def _resolve_cfg_path(cfg: dict, path_str: str) -> Path:
    p = Path(path_str).expanduser()
    if not p.is_absolute():
        p = Path(cfg["_repo_root"]) / p
    return p.resolve()

def run_hybrid_simulation(
    cfg: dict,
    scan_id: int,
    cluster_diameter_mm: float,
    num_calcs: int,
    calc_diameter_mm: float,
) -> None:
    system_name = cfg.get("system", "").lower()

    if system_name != "doheny":
        raise ValueError(f"Unsupported system: {system_name}")

    system = DohenySystem(cfg)
    system.validate()

    paths_cfg = cfg["paths"]
    files_cfg = cfg["files"]

    patient_data_dir = paths_cfg["patient_data_dir"]
    output_dir = paths_cfg["output_dir"]
    voi_center_dir = paths_cfg["voi_center_dir"]

    # simulation / detector / recon params
    sim_cfg = cfg.get("simulation", {})
    recon_cfg = cfg.get("reconstruction", {})
    density = cfg.get("density", {})

    xsiz = sim_cfg["detector"]["xsiz"]
    ysiz = sim_cfg["detector"]["ysiz"]
    dexel_mm = sim_cfg["detector"]["dexel_mm"]
    vertical_offset_mm = sim_cfg["detector"]["vertical_offset_mm"]

    energy_bin_size_keV = sim_cfg["energy_bin_size_keV"]
    csI_thickness_mm = sim_cfg["csI_thickness_mm"]
    new_vx_um = sim_cfg["new_vx_um"]

    # load scanlog
    df_scanlog, VscanID, VcropLocs = load_scanlog_doheny(cfg)
    iscan = get_scanlog_index_doheny(scan_id, VscanID)

    # load segmentation
    seg_volume, Nxyz, original_vx_um, labels = load_seg_doheny(
        scan_id=scan_id,
        patient_data_dir=patient_data_dir,
    )

    VcropLocs_scan = VcropLocs[iscan:iscan+1, :]
    FLAGchestwall = 0
    
    seg_volume_cropped, shift_mm_3D = fxn_crop_volume(
        iscan=0,
        volume=seg_volume,
        original_vx_um=original_vx_um,
        Nxyz=Nxyz,
        VcropLocs=VcropLocs_scan,
        FLAGchestwall=FLAGchestwall,
    )
    
    seg_volume_HR = fxn_upsample_volume_in_sections(
        seg_volume=seg_volume_cropped,
        original_vx_um=original_vx_um,
        new_vx_um=new_vx_um,
    )

    print(f"Cropped segmentation shape: {seg_volume_cropped.shape}")
    print(f"Upsampled segmentation shape: {seg_volume_HR.shape}")

    calc_shape = sim_cfg.get("calc_shape", "sphere")
    voxel_size_mm = 0.001 * new_vx_um
    
    calc = fxn_generate_calc(
        voxel_size_mm=voxel_size_mm,
        default_size_calc_mm=calc_diameter_mm,
        shape=calc_shape,
        saveFLAG=0,
    )

    print(f"Generated calc shape: {calc.shape}")

    # load projections and geometry
    prjstack, geo, ang = load_projections_and_geometry_doheny(
        scan_id=scan_id,
        patient_data_dir=patient_data_dir,
        xsiz=xsiz,
        ysiz=ysiz,
        dexel_mm=dexel_mm,
        vertical_offset_mm=vertical_offset_mm,
    )

    # load spectrum
    spectrum_path = _resolve_cfg_path(cfg, files_cfg["spectrum"])
    energy_keV, fluence = read_energy_spectrum(spectrum_path)

    # load material files
    material_files = files_cfg.get("material_files", {})
    calc_energy, calc_mu = read_material_file(_resolve_cfg_path(cfg, material_files["calc"]))
    adipose_energy, adipose_mu = read_material_file(_resolve_cfg_path(cfg, material_files["adipose"]))
    gland_energy, gland_mu = read_material_file(_resolve_cfg_path(cfg, material_files["glandular"]))
    csi_energy, csi_mu = read_material_file(_resolve_cfg_path(cfg, material_files["csI"]))

    print("Configuration validated.")
    print(f"Loaded scanlog row index: {iscan}")
    print(f"Segmentation shape: {seg_volume.shape}")
    print(f"Projection stack shape: {prjstack.shape}")
    print(f"Spectrum bins: {len(energy_keV)}")
    print(f"Output directory: {output_dir}")
    print(f"VOI center directory: {voi_center_dir}")

    # ---------------------------------------------------------
    # TODO: add actual hybrid simulation workflow here
    # - crop volume
    # - upsample segmentation
    # - generate calc cluster
    # - insert calc cluster
    # - ray trace
    # - blur projections
    # - combine with patient projections
    # - reconstruct
    # - extract VOIs / MIPs
    # ---------------------------------------------------------

    # Temporary metadata write for pipeline scaffolding
    run_params = {
        "scan_id": scan_id,
        "effective_energy_keV": recon_cfg["effective_energy_keV"],
        "energy_bin_size_keV": energy_bin_size_keV,
        "csI_thickness_mm": csI_thickness_mm,
        "vertical_offset_mm": vertical_offset_mm,
        "new_vx_um": new_vx_um,
        "calc_shape": sim_cfg.get("calc_shape", "sphere"),
        "calc_diameter_mm": calc_diameter_mm,
        "num_calcs": num_calcs,
        "cluster_diameter_mm": cluster_diameter_mm,
        "clusterCenteredFLAG": sim_cfg.get("clusterCenteredFLAG", 0),
        "num_SPvois_perbreast": cfg["voi"]["num_SPvois_perbreast"],
        "num_SAvois_perbreast": cfg["voi"]["num_SAvois_perbreast"],
        "voi_size_mm": cfg["voi"]["voi_size_mm"],
        "voi_size_vx": "TBD",
        "flagHU": recon_cfg["flagHU"],
        "mu_water": recon_cfg["mu_water"],
    }

    recon_alg = recon_cfg["algorithms"][0]
    kernel = recon_cfg["kernels"][0]
    folder_suffix = f"scan{scan_id:04d}"

    output_cfg = cfg.get("output", {})
    patch_subdir = output_cfg.get("patch_subdir", "CalcPatches/Patches")
    patch_root = Path(output_dir) / patch_subdir

    write_metadata(
        output_patch_root=patch_root,
        cfg=cfg,
        run_params=run_params,
        density=density,
        geo=geo,
        ang=ang,
        recon_alg=recon_alg,
        folder_suffix=folder_suffix,
        kernel=kernel,
    )

    print("Pipeline scaffold completed.")
