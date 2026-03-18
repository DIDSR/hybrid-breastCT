from __future__ import annotations

import numpy as np
from scipy.interpolate import RegularGridInterpolator


def fxn_crop_volume(iscan, volume, original_vx_um, Nxyz, VcropLocs, FLAGchestwall):
    rowstart = int(VcropLocs[iscan, 0])
    rowend = int(VcropLocs[iscan, 1])

    colstart = int(VcropLocs[iscan, 2])
    colend = int(VcropLocs[iscan, 3])

    zstart = int(VcropLocs[iscan, 4])
    zend = int(VcropLocs[iscan, 5])

    if FLAGchestwall == 1:
        volume_cropped = volume[:zend, rowstart:rowend+1, colstart:colend+1]  # [z, y, x]
        shift_mm_3D = [
            0,
            rowstart * (0.001 * original_vx_um),
            colstart * (0.001 * original_vx_um),
        ]  # [Z, Y, X] for tigre
    else:  # cropped a smaller volume that does not include the chest wall
        volume_cropped = volume[zstart:zend+1, rowstart:rowend+1, colstart:colend+1]  # [z, y, x]
        shift_mm_3D = [
            zstart * (0.001 * original_vx_um),
            rowstart * (0.001 * original_vx_um),
            colstart * (0.001 * original_vx_um),
        ]  # [Z, Y, X] for tigre

    Nxyz_cropped = volume_cropped.shape
    print(f"Cropped volume dimensions: {Nxyz_cropped[0]} x {Nxyz_cropped[1]} x {Nxyz_cropped[2]}", flush=True)
    return volume_cropped, shift_mm_3D


def fxn_upsample_volume(seg_volume, original_vx_um, new_vx_um):
    Nxyz = seg_volume.shape
    upsample_rate = original_vx_um / new_vx_um
    print(f"\nOriginal voxel size = {round(original_vx_um)} um.", flush=True)
    print(f"New voxel size      = {round(new_vx_um)} um.", flush=True)

    Nxyz_HR = [
        round(upsample_rate * Nxyz[0]),
        round(upsample_rate * Nxyz[1]),
        round(upsample_rate * Nxyz[2]),
    ]
    print(
        f"Forming upsampled cropped volume with dimensions: {Nxyz_HR[0]} x {Nxyz_HR[1]} x {Nxyz_HR[2]} ... ",
        flush=True,
    )

    # Define the original coordinates
    z_orig = np.linspace(1, Nxyz[0], Nxyz[0], dtype=np.float32)
    y_orig = np.linspace(1, Nxyz[1], Nxyz[1], dtype=np.float32)
    x_orig = np.linspace(1, Nxyz[2], Nxyz[2], dtype=np.float32)

    # Define the new coordinates for the upsampled volume
    z_new = np.linspace(1, Nxyz[0], Nxyz_HR[0], dtype=np.float32)
    y_new = np.linspace(1, Nxyz[1], Nxyz_HR[1], dtype=np.float32)
    x_new = np.linspace(1, Nxyz[2], Nxyz_HR[2], dtype=np.float32)

    # Create a meshgrid for the new coordinates
    z_mesh, y_mesh, x_mesh = np.meshgrid(z_new, y_new, x_new, indexing="ij")

    # Perform the interpolation
    interp = RegularGridInterpolator((z_orig, y_orig, x_orig), seg_volume, method="nearest")
    seg_volume_HR = interp(
        np.array([z_mesh.flatten(), y_mesh.flatten(), x_mesh.flatten()]).T
    ).reshape(z_mesh.shape)

    return seg_volume_HR


def fxn_upsample_volume_in_sections(seg_volume, original_vx_um, new_vx_um):
    Nxyz = seg_volume.shape
    upsample_rate = original_vx_um / new_vx_um
    print(f"\nOriginal voxel size = {round(original_vx_um)} um.", flush=True)
    print(f"New voxel size      = {round(new_vx_um)} um.", flush=True)

    Nxyz_HR = [
        round(upsample_rate * Nxyz[0]),
        round(upsample_rate * Nxyz[1]),
        round(upsample_rate * Nxyz[2]),
    ]
    print(
        f"Forming upsampled cropped volume with dimensions: {Nxyz_HR[0]} x {Nxyz_HR[1]} x {Nxyz_HR[2]} ... ",
        flush=True,
    )

    slices = 50  # Adjust this value based on your available RAM
    sections = list(range(0, Nxyz[0], slices)) + [Nxyz[0]]
    num_sections = len(sections) - 1

    Nz = Nxyz[0]
    Ny = Nxyz[1]
    Nx = Nxyz[2]
    Nz_HR_original = round(upsample_rate * slices)
    Ny_HR = round(upsample_rate * Ny)
    Nx_HR = round(upsample_rate * Nx)
    seg_volume_HR = np.zeros(Nxyz_HR, dtype=np.uint8)

    print("Section ", end="", flush=True)
    for section in range(num_sections):
        print(f"{section+1}/{num_sections}, ", end="", flush=True)
        starting_slice = sections[section]
        ending_slice = sections[section+1]

        subvolume = seg_volume[starting_slice:ending_slice, :, :]

        # Define the original coordinates
        z_orig = np.linspace(1, subvolume.shape[0], subvolume.shape[0], dtype=np.float32)
        y_orig = np.linspace(1, Ny, Ny, dtype=np.float32)
        x_orig = np.linspace(1, Nx, Nx, dtype=np.float32)

        # Define the new coordinates for the upsampled volume
        Nz_HR = round(upsample_rate * subvolume.shape[0])
        z_new = np.linspace(1, subvolume.shape[0], Nz_HR, dtype=np.float32)
        y_new = np.linspace(1, Ny, Ny_HR, dtype=np.float32)
        x_new = np.linspace(1, Nx, Nx_HR, dtype=np.float32)

        # Create a meshgrid for the new coordinates
        z_mesh, y_mesh, x_mesh = np.meshgrid(z_new, y_new, x_new, indexing="ij")

        # Perform the interpolation
        interp = RegularGridInterpolator((z_orig, y_orig, x_orig), subvolume, method="nearest")
        upsampled_volume = interp(
            np.array([z_mesh.flatten(), y_mesh.flatten(), x_mesh.flatten()]).T
        ).reshape(z_mesh.shape)

        if section < num_sections - 1:
            seg_volume_HR[Nz_HR_original * section:Nz_HR_original * section + Nz_HR, :, :] = upsampled_volume
        else:
            # Last section: need to correct for rounding errors
            diff = (
                upsampled_volume.shape[0]
                - seg_volume_HR[
                    Nz_HR_original * section:Nz_HR_original * section + Nz_HR, :, :
                ].shape[0]
            )
            upsampled_volume = upsampled_volume[: (upsampled_volume.shape[0] - diff), :, :]
            seg_volume_HR[(Nz_HR_original * section):((Nz_HR_original * section) + Nz_HR), :, :] = upsampled_volume

    print("done.", flush=True)
    return seg_volume_HR
