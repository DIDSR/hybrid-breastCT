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
from hybrid_bct.io.materials import read_material_file, interpolate_material_attenuation
from hybrid_bct.metadata import write_metadata
from hybrid_bct.simulation.calc_models import fxn_generate_calc
from hybrid_bct.simulation.insertion import fxn_insert_calc_cluster_new
from hybrid_bct.simulation.voi import fxn_getVOIcenters, fxn_extract_and_save_vois
from hybrid_bct.simulation.volume import (
    fxn_crop_volume,
    fxn_upsample_volume_in_sections,
)
from hybrid_bct.simulation.projection import (
    build_cluster_candidate_offsets,
    compute_cor_offset_mm,
    generate_hybrid_projection_stack,
)
from hybrid_bct.reconstruction import reconstruct_hybrid_volume

print(">>> ENTERED run_hybrid_simulation", flush=True)

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
    num_SPvois_perbreast = cfg["voi"]["num_SPvois_perbreast"]
    num_SAvois_perbreast = cfg["voi"]["num_SAvois_perbreast"]

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

    material_files = files_cfg.get("material_files", {})

    calc_energy_keV, calc_mu_raw = read_material_file(
        _resolve_cfg_path(cfg, material_files["calc"])
    )
    adipose_energy_keV, adipose_mu_raw = read_material_file(
        _resolve_cfg_path(cfg, material_files["adipose"])
    )
    gland_energy_keV, gland_mu_raw = read_material_file(
        _resolve_cfg_path(cfg, material_files["glandular"])
    )
    csi_energy_keV, csi_mu_raw = read_material_file(
        _resolve_cfg_path(cfg, material_files["csI"])
    )

    mu_cm = {
        "calc": interpolate_material_attenuation(
            calc_energy_keV,
            calc_mu_raw,
            energy_keV,
            density=density["calc"],
        ),
        "adipose": interpolate_material_attenuation(
            adipose_energy_keV,
            adipose_mu_raw,
            energy_keV,
            density=density["adipose"],
        ),
        "glandular": interpolate_material_attenuation(
            gland_energy_keV,
            gland_mu_raw,
            energy_keV,
            density=density["glandular"],
        ),
        "csI": interpolate_material_attenuation(
            csi_energy_keV,
            csi_mu_raw,
            energy_keV,
            density=density["csI"],
        ),
    }

    energy_levels_eV = energy_keV * 1000.0
    photon_fluence = fluence
    QDE_CsI = 1.0 - np.exp(-0.1 * mu_cm["csI"] * csI_thickness_mm)

    detector_mtf_path = _resolve_cfg_path(cfg, files_cfg["detector_mtf"])
    mtf_data = np.loadtxt(detector_mtf_path, delimiter=",")
    f_MTF = mtf_data[:, 0]
    mtf_MTF = mtf_data[:, 1]

    print(f"Loaded detector MTF points: {len(f_MTF)}")
    
    cluster_info = build_cluster_candidate_offsets(
        cluster_diameter_mm=cluster_diameter_mm,
        calc=calc,
        voi_size_mm=voi_size_mm,
        new_vx_um=new_vx_um,
    )

    CORoffset_mm_3D, zstart_mm = compute_cor_offset_mm(
        iscan=0,
        VcropLocs=VcropLocs_scan,
        Nxyz=Nxyz,
        original_vx_um=original_vx_um,
        FLAGchestwall=FLAGchestwall,
    )

    hybrid_prjstack = generate_hybrid_projection_stack(
        volume_with_calcs=volume_with_calcs,
        geo=geo,
        ang=ang,
        photon_fluence=photon_fluence,
        QDE_CsI=QDE_CsI,
        energy_levels_eV=energy_levels_eV,
        mu_cm=mu_cm,
        dexel_mm=dexel_mm,
        f_MTF=f_MTF,
        mtf_MTF=mtf_MTF,
        prjstack=prjstack,
        new_vx_um=new_vx_um,
        vertical_offset_mm=vertical_offset_mm,
        zstart_mm=zstart_mm,
        CORoffset_mm_3D=CORoffset_mm_3D,
    )

    print(f"Hybrid projection stack shape: {hybrid_prjstack.shape}")

    recon_alg = recon_cfg["algorithms"][0]
    kernel = recon_cfg["kernels"][0]
    iterations = recon_cfg.get("iterations", 30)
    folder_suffix = f"scan{scan_id:04d}"

    hybrid_ct_volume = reconstruct_hybrid_volume(
        hybrid_prjstack=hybrid_prjstack,
        geo=geo,
        ang=ang,
        recon_alg=recon_alg,
        kernel=kernel,
        iterations=iterations,
    )

    output_cfg = cfg.get("output", {})
    mip_subdir = output_cfg.get("mip_subdir", "CalcPatches/MIPjpgs")
    patch_subdir = output_cfg.get("patch_subdir", "CalcPatches/Patches")

    loc_save_MIPjpgs = Path(output_dir) / mip_subdir
    loc_save_patches = Path(output_dir) / patch_subdir

    recon_size_mm = geo.dVoxel[0]
    voi_size_vx = round(voi_size_mm / recon_size_mm)
    voi_halfdim_vx = voi_size_vx // 2

    fxn_extract_and_save_vois(
        rec=hybrid_ct_volume,
        recon_alg=recon_alg,
        folder_suffix=folder_suffix,
        scanID=scan_id,
        loc_save_patches=loc_save_patches,
        loc_save_MIPjpgs=loc_save_MIPjpgs,
        calc_diameter_mm=calc_diameter_mm,
        cluster_diameter_mm=cluster_diameter_mm,
        num_calcs=num_calcs,
        density=density,
        num_SPvois_perbreast=num_SPvois_perbreast,
        voi_centers_mm_SP=voi_centers_mm_SP,
        voi_centers_mm_SA=voi_centers_mm_SA,
        recon_size_mm=recon_size_mm,
        voi_halfdim_vx=voi_halfdim_vx,
        flagHU=recon_cfg["flagHU"],
    )

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
        "voi_size_vx": voi_size_vx,
        "flagHU": recon_cfg["flagHU"],
        "mu_water": recon_cfg["mu_water"],
    }

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
