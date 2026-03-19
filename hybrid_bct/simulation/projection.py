from __future__ import annotations

import numpy as np
import tigre

from hybrid_bct.simulation.blur import mtf_blur
from hybrid_bct.reconstruction import fxn_alter_geometry


def build_cluster_candidate_offsets(
    cluster_diameter_mm,
    calc,
    voi_size_mm,
    new_vx_um,
):
    cluster_radius_vx = round(cluster_diameter_mm / (2 * (0.001 * new_vx_um)))
    voi_size_vx = round(voi_size_mm / (0.001 * new_vx_um))
    half_dim = calc.shape[0] // 2
    margin_vx = round(np.ceil(cluster_radius_vx) + 1)
    possible_distance_from_voi_center = round((voi_size_vx - 2 * margin_vx) / 2)

    cluster_diameter_vx = int(np.ceil(cluster_diameter_mm / (0.001 * new_vx_um)))
    logical_sphere_center = np.array(
        [cluster_radius_vx, cluster_radius_vx, cluster_radius_vx]
    )

    x, y, z = np.mgrid[
        1 : cluster_diameter_vx + 1,
        1 : cluster_diameter_vx + 1,
        1 : cluster_diameter_vx + 1,
    ]

    distances_from_center = np.sqrt(
        (x - logical_sphere_center[0]) ** 2
        + (y - logical_sphere_center[1]) ** 2
        + (z - logical_sphere_center[2]) ** 2
    )
    logical_sphere = distances_from_center <= (cluster_radius_vx - half_dim)

    true_x, true_y, true_z = np.where(logical_sphere)

    return {
        "cluster_radius_vx": cluster_radius_vx,
        "voi_size_vx": voi_size_vx,
        "half_dim": half_dim,
        "margin_vx": margin_vx,
        "possible_distance_from_voi_center": possible_distance_from_voi_center,
        "cluster_diameter_vx": cluster_diameter_vx,
        "logical_sphere_center": logical_sphere_center,
        "true_x": true_x,
        "true_y": true_y,
        "true_z": true_z,
    }


def compute_cor_offset_mm(
    iscan,
    VcropLocs,
    Nxyz,
    original_vx_um,
    FLAGchestwall,
):
    rowstart = int(VcropLocs[iscan, 0])
    rowend = int(VcropLocs[iscan, 1])
    colstart = int(VcropLocs[iscan, 2])
    colend = int(VcropLocs[iscan, 3])
    zstart = int(VcropLocs[iscan, 4])
    zend = int(VcropLocs[iscan, 5])

    zstart_mm = zstart * (0.001 * original_vx_um)

    COR_mm = np.array(
        [
            Nxyz[1] / 2 * (0.001 * original_vx_um),
            Nxyz[2] / 2 * (0.001 * original_vx_um),
        ]
    )
    newCOR_mm = np.array(
        [
            (rowstart + rowend) / 2 * (0.001 * original_vx_um),
            (colstart + colend) / 2 * (0.001 * original_vx_um),
        ]
    )

    if FLAGchestwall == 1:
        CORoffset_mm_3D = -np.array(
            [0, newCOR_mm[0] - COR_mm[0], newCOR_mm[1] - COR_mm[1]]
        )
    else:
        COR_mm_z = np.squeeze(np.array([Nxyz[2] / 2 * (0.001 * original_vx_um)]))
        newCOR_mm_z = np.squeeze(
            np.array([(zstart + zend) / 2 * (0.001 * original_vx_um)])
        )
        CORoffset_mm_3D = -np.array(
            [
                newCOR_mm_z - COR_mm_z,
                newCOR_mm[0] - COR_mm[0],
                newCOR_mm[1] - COR_mm[1],
            ]
        )

    return CORoffset_mm_3D, zstart_mm


def generate_hybrid_projection_stack(
    volume_with_calcs,
    geo,
    ang,
    photon_fluence,
    QDE_CsI,
    energy_levels_eV,
    mu_cm,
    dexel_mm,
    f_MTF,
    mtf_MTF,
    prjstack,
    new_vx_um,
    vertical_offset_mm,
    zstart_mm,
    CORoffset_mm_3D,
):
    geo_raytracing = fxn_alter_geometry(
        geo,
        volume_with_calcs,
        CORoffset_mm_3D,
        new_vx_um,
        vertical_offset_mm,
        zstart_mm,
    )

    volume_with_calcs_flipped = np.rot90(volume_with_calcs, k=2, axes=(1, 2))

    sum_prjstack_calcs = np.zeros(
        (len(ang), geo_raytracing.nDetector[0], geo_raytracing.nDetector[1]),
        dtype=np.float32,
    )
    total_I_o = 0.0

    for iE, energy_eV in enumerate(energy_levels_eV):
        volume_with_calcs_mu_mm = np.zeros_like(volume_with_calcs, dtype=np.float32)
        volume_with_calcs_mu_mm[volume_with_calcs_flipped == 101] = 0.1 * abs(
            mu_cm["calc"][iE] - mu_cm["glandular"][iE]
        )
        volume_with_calcs_mu_mm[volume_with_calcs_flipped == 100] = 0.1 * abs(
            mu_cm["calc"][iE] - mu_cm["adipose"][iE]
        )

        prjstack_calcs = tigre.Ax(
            np.float32(volume_with_calcs_mu_mm),
            geo_raytracing,
            ang,
            method="interpolated",
        )

        sum_prjstack_calcs += (
            photon_fluence[iE] * QDE_CsI[iE] * np.exp(-prjstack_calcs)
        )
        total_I_o += photon_fluence[iE] * QDE_CsI[iE]

    blurred_prjstack_calcs = mtf_blur(
        prjstack_calcs=sum_prjstack_calcs,
        dexel_mm=dexel_mm,
        f_MTF=f_MTF,
        mtf_MTF=mtf_MTF,
    )

    eps = np.finfo(np.float32).eps
    prjstack_calcs_blur_QDE_lognorm = -np.log(
        np.clip(blurred_prjstack_calcs / total_I_o, eps, None)
    )
    hybrid_prjstack = (prjstack + prjstack_calcs_blur_QDE_lognorm).astype(np.float32)

    return hybrid_prjstack
