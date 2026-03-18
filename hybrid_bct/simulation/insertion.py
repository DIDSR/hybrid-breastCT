from __future__ import annotations

import numpy as np


def fxn_insert_calc_cluster_new(
    volume_with_calcs,
    seg_volume_HR,
    true_x,
    true_y,
    true_z,
    logical_sphere_center,
    cluster_center,
    num_calcs,
    calc,
    half_dim,
    labels,
):
    calcs_added = 0
    while calcs_added < num_calcs:
        idx = np.random.randint(len(true_x))
        random_distance_from_cluster_center = (
            np.array([true_z[idx], true_y[idx], true_x[idx]]) - logical_sphere_center
        )
        calc_center = cluster_center - random_distance_from_cluster_center

        z_start = calc_center[0] - half_dim
        z_end = calc_center[0] + (half_dim if calc.shape[0] % 2 == 0 else half_dim + 1)
        y_start = calc_center[1] - half_dim
        y_end = calc_center[1] + (half_dim if calc.shape[1] % 2 == 0 else half_dim + 1)
        x_start = calc_center[2] - half_dim
        x_end = calc_center[2] + (half_dim if calc.shape[2] % 2 == 0 else half_dim + 1)

        # Check if the indices are within the valid range
        if (
            z_start < 0
            or z_end > volume_with_calcs.shape[0]
            or y_start < 0
            or y_end > volume_with_calcs.shape[1]
            or x_start < 0
            or x_end > volume_with_calcs.shape[2]
        ):
            print(f"Skipping calc at index {idx} due to out of bounds indices", flush=True)
            continue

        if volume_with_calcs[z_start:z_end, y_start:y_end, x_start:x_end].shape != calc.shape:
            print("Shape mismatch", flush=True)

        volume_with_calcs[z_start:z_end, y_start:y_end, x_start:x_end] += calc
        calcs_added += 1

    volume_with_calcs[(seg_volume_HR == labels["glandular"]) & (volume_with_calcs != 0)] = 101
    volume_with_calcs[(seg_volume_HR == labels["adipose"]) & (volume_with_calcs != 0)] = 100

    return volume_with_calcs
