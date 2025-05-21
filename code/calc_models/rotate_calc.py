import numpy as np
from scipy.ndimage import rotate
import random

def rotate_calc(calc):
    # Rotate around z-axis
    calc = rotate(calc, random.randint(0, 360), axes=(1, 0), reshape=True)
    # Rotate around y-axis
    calc = rotate(calc, random.randint(0, 360), axes=(2, 0), reshape=True)
    # Rotate around x-axis
    calc = rotate(calc, random.randint(0, 360), axes=(2, 1), reshape=True)

    new_size = max(calc.shape)
    new_start = np.floor((new_size - np.array(calc.shape)) / 2).astype(int)

    rotated_calc = np.zeros((new_size, new_size, new_size), dtype=np.uint8)
    rotated_calc[
        new_start[0]:new_start[0]+calc.shape[0],
        new_start[1]:new_start[1]+calc.shape[1],
        new_start[2]:new_start[2]+calc.shape[2]
    ] = calc

    return rotated_calc

