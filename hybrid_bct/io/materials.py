from __future__ import annotations

from pathlib import Path
import re
import numpy as np


def _extract_numeric_rows(lines: list[str]) -> np.ndarray:
    """
    Extract rows containing at least two numeric values from a text file.
    Returns a 2D float array with shape (N, >=2).
    """
    rows: list[list[float]] = []

    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("#"):
            continue

        parts = re.split(r"\s+", stripped)
        numeric_vals: list[float] = []

        for part in parts:
            try:
                numeric_vals.append(float(part))
            except ValueError:
                continue

        if len(numeric_vals) >= 2:
            rows.append(numeric_vals)

    if not rows:
        raise ValueError("No numeric attenuation data found in material file.")

    min_len = min(len(row) for row in rows)
    trimmed = np.array([row[:min_len] for row in rows], dtype=float)
    return trimmed


def read_material_file(filepath: str | Path) -> tuple[np.ndarray, np.ndarray]:
    """
    Read a material attenuation file and return energy (keV) and attenuation values.

    This function is intentionally permissive to support simple text-based attenuation
    files used in public examples. It looks for rows with at least two numeric columns
    and interprets:
      - column 1 as energy
      - column 2 as attenuation

    Returns
    -------
    energy_keV : np.ndarray
        1D array of energy values in keV.
    attenuation : np.ndarray
        1D array of attenuation values corresponding to each energy.

    Raises
    ------
    FileNotFoundError
        If the file does not exist.
    ValueError
        If numeric data cannot be extracted.
    """
    path = Path(filepath).expanduser().resolve()
    if not path.exists():
        raise FileNotFoundError(f"Material file not found: {path}")

    lines = path.read_text().splitlines()
    data = _extract_numeric_rows(lines)

    energy_keV = np.asarray(data[:, 0], dtype=float)
    attenuation = np.asarray(data[:, 1], dtype=float)

    if energy_keV.size == 0:
        raise ValueError(f"Material file contains no usable data: {path}")

    if np.any(np.diff(energy_keV) < 0):
        order = np.argsort(energy_keV)
        energy_keV = energy_keV[order]
        attenuation = attenuation[order]

    return energy_keV, attenuation

def interpolate_material_attenuation(
    material_energy_keV: np.ndarray,
    material_attenuation: np.ndarray,
    target_energy_keV: np.ndarray,
    density: float | None = None,
) -> np.ndarray:
    """
    Interpolate attenuation values onto a target energy grid.

    Parameters
    ----------
    material_energy_keV : np.ndarray
        Energy grid from the material file.
    material_attenuation : np.ndarray
        Attenuation values from the material file.
    target_energy_keV : np.ndarray
        Energy grid to interpolate onto.
    density : float | None, optional
        Optional density scaling factor. If provided, attenuation values are
        multiplied by density after interpolation.

    Returns
    -------
    np.ndarray
        Interpolated attenuation values on target_energy_keV.
    """
    interp_vals = np.interp(
        target_energy_keV,
        material_energy_keV,
        material_attenuation,
    )

    if density is not None:
        interp_vals = interp_vals * density

    return interp_vals
