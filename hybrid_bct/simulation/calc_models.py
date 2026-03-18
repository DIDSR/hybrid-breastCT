from __future__ import annotations

from pathlib import Path
import math
from typing import Optional, Tuple

import numpy as np
from scipy import ndimage


def rotate_calc(
    calc: np.ndarray,
    angle_deg: float,
    axes: Tuple[int, int] = (1, 2),
    order: int = 1,
) -> np.ndarray:
    """
    Rotate a 3D calcification volume without changing output shape.

    Parameters
    ----------
    calc : np.ndarray
        3D calcification volume.
    angle_deg : float
        Rotation angle in degrees.
    axes : tuple[int, int], optional
        Plane of rotation. Default rotates in the y-x plane.
    order : int, optional
        Interpolation order passed to scipy.ndimage.rotate.

    Returns
    -------
    np.ndarray
        Rotated 3D calcification volume.
    """
    if calc.ndim != 3:
        raise ValueError("rotate_calc expects a 3D array.")

    rotated = ndimage.rotate(
        calc,
        angle=angle_deg,
        axes=axes,
        reshape=False,
        order=order,
        mode="constant",
        cval=0.0,
    )

    return rotated.astype(calc.dtype, copy=False)


def fxn_generate_spherical_calc(
    voxel_size_mm: float,
    default_size_calc_mm: float,
    calc_value: float = 1.0,
) -> np.ndarray:
    """
    Generate a spherical calcification volume.

    Parameters
    ----------
    voxel_size_mm : float
        Isotropic voxel size in mm.
    default_size_calc_mm : float
        Calcification diameter in mm.
    calc_value : float, optional
        Value assigned inside the calcification.

    Returns
    -------
    np.ndarray
        3D volume containing a spherical calcification.
    """
    if voxel_size_mm <= 0:
        raise ValueError("voxel_size_mm must be positive.")
    if default_size_calc_mm <= 0:
        raise ValueError("default_size_calc_mm must be positive.")

    radius_mm = default_size_calc_mm / 2.0
    radius_vox = radius_mm / voxel_size_mm

    margin = 2
    half_width = int(math.ceil(radius_vox)) + margin
    size = 2 * half_width + 1

    z, y, x = np.meshgrid(
        np.arange(size) - half_width,
        np.arange(size) - half_width,
        np.arange(size) - half_width,
        indexing="ij",
    )

    rr = np.sqrt(x**2 + y**2 + z**2)

    calc = np.zeros((size, size, size), dtype=np.float32)
    calc[rr <= radius_vox] = calc_value
    return calc


def _save_calc_volume(
    calc: np.ndarray,
    voxel_size_mm: float,
    default_size_calc_mm: float,
    shape: str,
    save_dir: str | Path,
) -> tuple[Path, Path]:
    """
    Save calc volume and sidecar metadata.

    Returns
    -------
    tuple[Path, Path]
        (raw_file_path, txt_file_path)
    """
    save_dir = Path(save_dir).expanduser().resolve()
    save_dir.mkdir(parents=True, exist_ok=True)

    stem = (
        f"calc_{shape}_"
        f"{default_size_calc_mm:.3f}mm_"
        f"{calc.shape[0]}x{calc.shape[1]}x{calc.shape[2]}"
    )

    raw_path = save_dir / f"{stem}.raw"
    txt_path = save_dir / f"{stem}.txt"

    calc.astype(np.float32).tofile(raw_path)

    txt_path.write_text(
        "\n".join(
            [
                f"raw_file: {raw_path.name}",
                f"shape: {shape}",
                f"calc_diameter_mm: {default_size_calc_mm}",
                f"voxel_size_mm: {voxel_size_mm}",
                f"volume_shape: {calc.shape}",
            ]
        )
        + "\n"
    )

    return raw_path, txt_path


def fxn_generate_calc(
    voxel_size_mm: float,
    default_size_calc_mm: float,
    shape: str = "sphere",
    rotate_angle_deg: Optional[float] = None,
    saveFLAG: int = 0,
    save_dir: Optional[str | Path] = None,
) -> np.ndarray:
    """
    Generate a calcification model.

    This keeps the original function-style naming for easier migration from the
    older codebase, but removes hard-coded internal output paths.

    Parameters
    ----------
    voxel_size_mm : float
        Isotropic voxel size in mm.
    default_size_calc_mm : float
        Calcification diameter in mm.
    shape : str, optional
        Calcification shape. Currently only 'sphere' is supported.
    rotate_angle_deg : float | None, optional
        Optional rotation angle in degrees.
    saveFLAG : int, optional
        If 1, save the generated calcification volume to disk.
    save_dir : str | Path | None, optional
        Directory where the calcification volume should be saved if saveFLAG=1.

    Returns
    -------
    np.ndarray
        3D calcification volume.
    """
    shape = shape.lower()

    if shape == "sphere":
        calc = fxn_generate_spherical_calc(
            voxel_size_mm=voxel_size_mm,
            default_size_calc_mm=default_size_calc_mm,
        )
    else:
        raise ValueError(
            f"Unsupported calcification shape '{shape}'. "
            "Currently supported: ['sphere']"
        )

    if rotate_angle_deg is not None:
        calc = rotate_calc(calc, rotate_angle_deg)

    if saveFLAG == 1:
        if save_dir is None:
            raise ValueError("save_dir must be provided when saveFLAG == 1.")
        _save_calc_volume(
            calc=calc,
            voxel_size_mm=voxel_size_mm,
            default_size_calc_mm=default_size_calc_mm,
            shape=shape,
            save_dir=save_dir,
        )

    return calc
