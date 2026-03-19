from __future__ import annotations

from pathlib import Path
import numpy as np

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
from hybrid_bct.simulation.voi import fxn_getVOIcenters
from hybrid_bct.simulation.blur import mtf_blur

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
    voi_size_mm = cfg["voi"]["voi_size_mm"]

    (
        voi_centers_mm_SP,
        voi_centers_mm_SA,
        num_SPvois_perbreast,
        num_SAvois_perbreast,
    ) = fxn_getVOIcenters(
        scanID=scan_id,
        voi_size_mm=voi_size_mm,
        num_SPvois_perbreast=num_SPvois_perbreast,
        num_SAvois_perbreast=num_SAvois_perbreast,
        input_folder=voi_center_dir,
    )

    print(f"Signal-present VOI centers loaded: {len(voi_centers_mm_SP)}")
    print(f"Signal-absent VOI centers loaded: {len(voi_centers_mm_SA)}")

    # Use the first signal-present VOI center for now
    voi_center_mm = voi_centers_mm_SP[0]

    # Convert VOI center from mm to voxel coordinates in the upsampled cropped volume
    voi_center_vx = np.round(voi_center_mm / voxel_size_mm).astype(int)

    print(f"Selected signal-present VOI center (mm): {voi_center_mm}")
    print(f"Selected signal-present VOI center (voxels): {voi_center_vx}")

    # Prepare calc-cluster insertion inputs
    volume_with_calcs = np.zeros_like(seg_volume_HR, dtype=np.uint8)

    cluster_center = voi_center_vx.astype(int)
    half_dim = calc.shape[0] // 2

    logical_sphere_center = np.array(
        [cluster_center[0], cluster_center[1], cluster_center[2]],
        dtype=int,
    )

    # For now, use all nonzero voxels in the calc as candidate offsets
    true_z, true_y, true_x = np.where(calc > 0)

    print(f"Cluster center (voxels): {cluster_center}")
    print(f"Calc half dimension: {half_dim}")
    print(f"Number of candidate calc voxels: {len(true_x)}")

    volume_with_calcs = fxn_insert_calc_cluster_new(
        volume_with_calcs=volume_with_calcs,
        seg_volume_HR=seg_volume_HR,
        true_x=true_x,
        true_y=true_y,
        true_z=true_z,
        logical_sphere_center=logical_sphere_center,
        cluster_center=cluster_center,
        num_calcs=num_calcs,
        calc=calc,
        half_dim=half_dim,
        labels=labels,
    )

    print(f"Inserted calc volume nonzero voxels: {np.count_nonzero(volume_with_calcs)}")
    print(f"Inserted calc unique values: {np.unique(volume_with_calcs)}")

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

    detector_mtf_path = _resolve_cfg_path(cfg, files_cfg["detector_mtf"])
    mtf_data = np.loadtxt(detector_mtf_path, delimiter=",")
    f_MTF = mtf_data[:, 0]
    mtf_MTF = mtf_data[:, 1]

    print(f"Loaded detector MTF points: {len(f_MTF)}")

    # Temporary placeholder for calc-only projection stack until ray tracing is wired in
    prjstack_calcs = np.zeros_like(prjstack, dtype=np.float32)

    prjstack_calcs_blurred = mtf_blur(
        prjstack_calcs=prjstack_calcs,
        dexel_mm=dexel_mm,
        f_MTF=f_MTF,
        mtf_MTF=mtf_MTF,
    )

    print(f"Blurred calc projection stack shape: {prjstack_calcs_blurred.shape}")

    print("Configuration validated.")
    print(f"Loaded scanlog row index: {iscan}")
    print(f"Segmentation shape: {seg_volume.shape}")
    print(f"Projection stack shape: {prjstack.shape}")
    print(f"Spectrum bins: {len(energy_keV)}")
    print(f"Output directory: {output_dir}")
    print(f"VOI center directory: {voi_center_dir}")

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
        "num_SPvois_perbreast": num_SPvois_perbreast,
        "num_SAvois_perbreast": num_SAvois_perbreast,
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
