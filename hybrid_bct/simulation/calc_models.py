from __future__ import annotations

from pathlib import Path
import math
import numpy as np
from scipy import ndimage


def rotate_calc(volume: np.ndarray, angle_deg: float, axes: tuple[int, int] = (0, 1)) -> np.ndarray:
    """
    Rotate a 3D calcification volume without changing output shape.

    Parameters
    ----------
    volume : np.ndarray
        3D binary or numeric volume.
    angle_deg : float
        Rotation angle in degrees.
    axes : tuple[int, int]
        Axes defining the plane of rotation.

    Returns
    -------
    np.ndarray
        Rotated volume with the same shape as the input.
    """
    if volume.ndim != 3:
        raise ValueError("rotate_calc expects a 3D volume.")

    rotated = ndimage.rotate(
        volume,
        angle=angle_deg,
        axes=axes,
        reshape=False,
        order=1,
        mode="constant",
        cval=0.0,
    )
    return rotated


def generate_spherical_calc(
    voxel_size_mm: float,
    calc_diameter_mm: float,
    value: float = 1.0,
) -> np.ndarray:
    """
    Generate a spherical calcification as a 3D array.

    Parameters
    ----------
    voxel_size_mm : float
        Isotropic voxel size in mm.
    calc_diameter_mm : float
        Diameter of spherical calcification in mm.
    value : float
        Value to assign inside the calcification.

    Returns
    -------
    np.ndarray
        3D array containing a spherical calcification.
    """
    if voxel_size_mm <= 0:
        raise ValueError("voxel_size_mm must be positive.")
    if calc_diameter_mm <= 0:
        raise ValueError("calc_diameter_mm must be positive.")

    radius_mm = calc_diameter_mm / 2.0
    radius_vox = radius_mm / voxel_size_mm

    margin = 2
    half_width = int(math.ceil(radius_vox)) + margin
    size = 2 * half_width + 1

    zz, yy, xx = np.meshgrid(
        np.arange(size) - half_width,
        np.arange(size) - half_width,
        np.arange(size) - half_width,
        indexing="ij",
    )

    rr = np.sqrt(xx**2 + yy**2 + zz**2)
    vol = np.zeros((size, size, size), dtype=np.float32)
    vol[rr <= radius_vox] = value
    return vol


def generate_calc(
    voxel_size_mm: float,
    calc_diameter_mm: float,
    shape: str = "sphere",
    rotate_angle_deg: float | None = None,
    save: bool = False,
    save_dir: str | Path | None = None,
) -> np.ndarray:
    """
    Generate a calcification model.

    Currently supported:
    - sphere

    Parameters
    ----------
    voxel_size_mm : float
        Isotropic voxel size in mm.
    calc_diameter_mm : float
        Calcification diameter in mm.
    shape : str
        Calcification shape.
    rotate_angle_deg : float | None
        Optional rotation angle to apply.
    save : bool
        Whether to save the volume to disk.
    save_dir : str | Path | None
        Directory where files are saved if save=True.

    Returns
    -------
    np.ndarray
        3D calcification volume.
    """
    shape = shape.lower()

    if shape == "sphere":
        calc = generate_spherical_calc(
            voxel_size_mm=voxel_size_mm,
            calc_diameter_mm=calc_diameter_mm,
        )
    else:
        raise ValueError(f"Unsupported calcification shape: {shape}")

    if rotate_angle_deg is not None:
        calc = rotate_calc(calc, rotate_angle_deg)

    if save:
        if save_dir is None:
            raise ValueError("save_dir must be provided when save=True")
        save_calc_volume(
            calc=calc,
            voxel_size_mm=voxel_size_mm,
            calc_diameter_mm=calc_diameter_mm,
            save_dir=save_dir,
            shape=shape,
        )

    return calc


def save_calc_volume(
    calc: np.ndarray,
    voxel_size_mm: float,
    calc_diameter_mm: float,
    save_dir: str | Path,
    shape: str = "sphere",
) -> tuple[Path, Path]:
    """
    Save calcification volume and sidecar metadata.

    Returns
    -------
    tuple[Path, Path]
        Paths to the raw volume file and metadata text file.
    """
    save_path = Path(save_dir).expanduser().resolve()
    save_path.mkdir(parents=True, exist_ok=True)

    stem = f"calc_{shape}_{calc_diameter_mm:.3f}mm_{calc.shape[0]}x{calc.shape[1]}x{calc.shape[2]}"
    raw_path = save_path / f"{stem}.raw"
    txt_path = save_path / f"{stem}.txt"

    calc.astype(np.float32).tofile(raw_path)

    txt_path.write_text(
        "\n".join(
            [
                f"raw_file: {raw_path.name}",
                f"shape: {shape}",
                f"calc_diameter_mm: {calc_diameter_mm}",
                f"voxel_size_mm: {voxel_size_mm}",
                f"volume_shape: {calc.shape}",
            ]
        )
        + "\n"
    )

    return raw_path, txt_path
