from __future__ import annotations

from pathlib import Path
import numpy as np


def read_energy_spectrum(filepath: str | Path) -> tuple[np.ndarray, np.ndarray]:
    """
    Read an energy spectrum file and return energy (keV) and relative fluence.

    Supported formats:
    - Two-column text files: energy_keV, fluence
    - '.spc' files stored as whitespace-delimited two-column text

    Returns
    -------
    energy_keV : np.ndarray
        1D array of energy bin centers in keV.
    fluence : np.ndarray
        1D array of spectral fluence values.

    Raises
    ------
    FileNotFoundError
        If the spectrum file does not exist.
    ValueError
        If the file cannot be parsed into at least two columns.
    """
    path = Path(filepath).expanduser().resolve()
    if not path.exists():
        raise FileNotFoundError(f"Spectrum file not found: {path}")

    try:
        data = np.loadtxt(path, comments="#")
    except Exception as exc:
        raise ValueError(f"Failed to read spectrum file: {path}") from exc

    if data.ndim == 1:
        if data.size < 2:
            raise ValueError(f"Spectrum file must contain at least two values: {path}")
        data = data.reshape(1, -1)

    if data.shape[1] < 2:
        raise ValueError(f"Spectrum file must have at least two columns: {path}")

    energy_keV = np.asarray(data[:, 0], dtype=float)
    fluence = np.asarray(data[:, 1], dtype=float)

    if energy_keV.size == 0:
        raise ValueError(f"Spectrum file is empty: {path}")

    if np.any(np.diff(energy_keV) < 0):
        order = np.argsort(energy_keV)
        energy_keV = energy_keV[order]
        fluence = fluence[order]

    return energy_keV, fluence
