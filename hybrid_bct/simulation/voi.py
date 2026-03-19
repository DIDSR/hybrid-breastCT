from __future__ import annotations

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt


def fxn_getVOIcenters(
    scanID,
    voi_size_mm,
    num_SPvois_perbreast,
    num_SAvois_perbreast,
    input_folder,
):
    txt_filepath = (
        Path(input_folder)
        / f"{voi_size_mm:02d}mm_voi_dim"
        / f"{scanID:04d}_voicenters_mm.txt"
    )

    if txt_filepath.exists():
        voicenters_mm_uncropped = np.loadtxt(txt_filepath)

        # Handle case where only one center is present
        voicenters_mm_uncropped = np.atleast_2d(voicenters_mm_uncropped)

        # Change order to match TIGRE [z, y, x]
        voicenters_mm_uncropped = voicenters_mm_uncropped[:, [2, 1, 0]]

        num_viable_voicenters = voicenters_mm_uncropped.shape[0]

        if num_viable_voicenters < (num_SPvois_perbreast + num_SAvois_perbreast):
            OG_num_SPvois_perbreast = num_SPvois_perbreast
            OG_num_SAvois_perbreast = num_SAvois_perbreast
            num_SPvois_perbreast = num_viable_voicenters // 2
            num_SAvois_perbreast = num_viable_voicenters // 2

            print(
                f"\nCannot generate {OG_num_SPvois_perbreast} signal present and "
                f"{OG_num_SAvois_perbreast} signal absent VOIs because only "
                f"{num_viable_voicenters} viable VOI centers.",
                flush=True,
            )
            print(
                f"New # signal present VOIs: {num_SPvois_perbreast}    "
                f"New # signal absent VOIs: {num_SAvois_perbreast}",
                flush=True,
            )

        np.random.shuffle(voicenters_mm_uncropped)

        voi_centers_mm_SP = voicenters_mm_uncropped[:num_SPvois_perbreast]
        voi_centers_mm_SA = voicenters_mm_uncropped[
            num_SPvois_perbreast:num_SPvois_perbreast + num_SAvois_perbreast
        ]

        return (
            voi_centers_mm_SP,
            voi_centers_mm_SA,
            num_SPvois_perbreast,
            num_SAvois_perbreast,
        )

    raise FileNotFoundError(
        f"File {txt_filepath} does not exist in the input folder."
    )


def fxn_extract_and_save_vois(
    rec,
    recon_alg,
    folder_suffix,
    scanID,
    loc_save_patches,
    loc_save_MIPjpgs,
    calc_diameter_mm,
    cluster_diameter_mm,
    num_calcs,
    density,
    num_SPvois_perbreast,
    voi_centers_mm_SP,
    voi_centers_mm_SA,
    recon_size_mm,
    voi_halfdim_vx,
    flagHU,
):
    num_SAvois_perbreast = num_SPvois_perbreast

    loc_save_patches = Path(loc_save_patches)
    loc_save_MIPjpgs = Path(loc_save_MIPjpgs)
    loc_save_patches.mkdir(parents=True, exist_ok=True)
    loc_save_MIPjpgs.mkdir(parents=True, exist_ok=True)

    # Create directories for this set of parameters
    loc_parameters_MIP = loc_save_patches / (
        f"{calc_diameter_mm:.2f}mmCalc_"
        f"{cluster_diameter_mm:.1f}mmCluster_"
        f"{num_calcs:02d}Calcs_"
        f"CaOx_{density['calc']:1.2f}_{recon_alg}_{folder_suffix}_MIP"
    )
    loc_parameters_VOI = loc_save_patches / (
        f"{calc_diameter_mm:.2f}mmCalc_"
        f"{cluster_diameter_mm:.1f}mmCluster_"
        f"{num_calcs:02d}Calcs_"
        f"CaOx_{density['calc']:1.2f}_{recon_alg}_{folder_suffix}_VOI"
    )

    loc_parameters_MIP.mkdir(parents=True, exist_ok=True)
    loc_parameters_VOI.mkdir(parents=True, exist_ok=True)

    print("Extracting VOIs... ", flush=True)

    # Extract signal-present VOIs
    for ivoi in range(num_SPvois_perbreast):
        voi_local_mm = voi_centers_mm_SP[ivoi]
        voi_zloc_mm = voi_local_mm[0]
        recon_voi_center = np.round(voi_centers_mm_SP[ivoi] / recon_size_mm).astype(int)

        voi_recon = rec[
            recon_voi_center[0] - voi_halfdim_vx : recon_voi_center[0] + voi_halfdim_vx,
            recon_voi_center[1] - voi_halfdim_vx : recon_voi_center[1] + voi_halfdim_vx,
            recon_voi_center[2] - voi_halfdim_vx : recon_voi_center[2] + voi_halfdim_vx,
        ]
        mip = np.max(voi_recon, axis=0)

        filename = loc_parameters_VOI / (
            f"{calc_diameter_mm:.2f}mmCalc_"
            f"{cluster_diameter_mm:.1f}mmCluster_"
            f"{num_calcs:02d}Calcs_"
            f"CaOx_{density['calc']:1.2f}_{recon_alg}_{folder_suffix}_VOI_"
            f"{voi_zloc_mm:06.1f}_{scanID:04d}_{ivoi:03d}_SP.npy"
        )
        np.save(filename, voi_recon)

        filename = loc_parameters_MIP / (
            f"{calc_diameter_mm:.2f}mmCalc_"
            f"{cluster_diameter_mm:.1f}mmCluster_"
            f"{num_calcs:02d}Calcs_"
            f"CaOx_{density['calc']:1.2f}_{recon_alg}_{folder_suffix}_MIP_"
            f"{voi_zloc_mm:06.1f}_{scanID:04d}_{ivoi:03d}_SP.npy"
        )
        np.save(filename, mip)

        if ivoi == 0:
            plt.figure()
            plt.imshow(mip, cmap="gray")
            plt.axis("off")
            plt.colorbar()
            if flagHU == 1:
                plt.title(f"{calc_diameter_mm:.2f} mm calcs | peak: {np.max(mip)} HU")
            else:
                plt.title(f"{calc_diameter_mm:.2f} mm calcs | peak: {np.max(mip):.2f} cm^-1")
            filename = loc_save_MIPjpgs / (
                f"{calc_diameter_mm:.2f}mmCalc_"
                f"{cluster_diameter_mm:.1f}mmCluster_"
                f"{num_calcs:02d}Calcs_"
                f"CaOx_{density['calc']:1.2f}_{recon_alg}_{folder_suffix}_MIP_"
                f"{voi_zloc_mm:06.1f}_{scanID:04d}_{ivoi:03d}_SP.jpg"
            )
            plt.savefig(filename)
            plt.close()

    # Extract signal-absent VOIs
    for ivoi in range(num_SAvois_perbreast):
        voi_local_mm = voi_centers_mm_SA[ivoi]
        voi_zloc_mm = voi_local_mm[0]
        recon_voi_center = np.round(voi_centers_mm_SA[ivoi] / recon_size_mm).astype(int)

        voi_recon = rec[
            recon_voi_center[0] - voi_halfdim_vx : recon_voi_center[0] + voi_halfdim_vx,
            recon_voi_center[1] - voi_halfdim_vx : recon_voi_center[1] + voi_halfdim_vx,
            recon_voi_center[2] - voi_halfdim_vx : recon_voi_center[2] + voi_halfdim_vx,
        ]
        mip = np.max(voi_recon, axis=0)

        filename = loc_parameters_VOI / (
            f"{calc_diameter_mm:.2f}mmCalc_"
            f"{cluster_diameter_mm:.1f}mmCluster_"
            f"{num_calcs:02d}Calcs_"
            f"CaOx_{density['calc']:1.2f}_{recon_alg}_{folder_suffix}_VOI_"
            f"{voi_zloc_mm:06.1f}_{scanID:04d}_{ivoi:03d}_SA.npy"
        )
        np.save(filename, voi_recon)

        filename = loc_parameters_MIP / (
            f"{calc_diameter_mm:.2f}mmCalc_"
            f"{cluster_diameter_mm:.1f}mmCluster_"
            f"{num_calcs:02d}Calcs_"
            f"CaOx_{density['calc']:1.2f}_{recon_alg}_{folder_suffix}_MIP_"
            f"{voi_zloc_mm:06.1f}_{scanID:04d}_{ivoi:03d}_SA.npy"
        )
        np.save(filename, mip)

        if ivoi == 0:
            plt.figure()
            plt.imshow(mip, cmap="gray")
            plt.axis("off")
            plt.colorbar()
            if flagHU == 1:
                plt.title(f"{calc_diameter_mm:.2f} mm calcs | peak: {np.max(mip)} HU")
            else:
                plt.title(f"{calc_diameter_mm:.2f} mm calcs | peak: {np.max(mip):.2f} cm^-1")
            filename = loc_save_MIPjpgs / (
                f"{calc_diameter_mm:.2f}mmCalc_"
                f"{cluster_diameter_mm:.1f}mmCluster_"
                f"{num_calcs:02d}Calcs_"
                f"CaOx_{density['calc']:1.2f}_{recon_alg}_{folder_suffix}_MIP_"
                f"{scanID:04d}_{ivoi:03d}_SA.jpg"
            )
            plt.savefig(filename)
            plt.close()
