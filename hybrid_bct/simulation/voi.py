from __future__ import annotations

from pathlib import Path

import numpy as np


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

        # Randomly shuffle viable VOI centers
        np.random.shuffle(voicenters_mm_uncropped)

        # Assign first num_SPvois_perbreast as signal-present VOI centers
        voi_centers_mm_SP = voicenters_mm_uncropped[:num_SPvois_perbreast]

        # Assign next num_SAvois_perbreast as signal-absent VOI centers
        voi_centers_mm_SA = voicenters_mm_uncropped[
            num_SPvois_perbreast:num_SPvois_perbreast + num_SAvois_perbreast
        ]

        return (
            voi_centers_mm_SP,
            voi_centers_mm_SA,
            num_SPvois_perbreast,
            num_SAvois_perbreast,
        )

    else:
        raise FileNotFoundError(
            f"File {txt_filepath} does not exist in the input folder."
        )
