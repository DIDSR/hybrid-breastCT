import numpy as np
from calc_models.rotate_calc import rotate_calc
from calc_models.fxn_generate_spherical_calc import fxn_generate_spherical_calc

def fxn_generate_calc(voxel_size_mm, default_size_calc_mm, shape, saveFLAG):
    if shape == 'sphere':
        calc_diameter_mm = default_size_calc_mm
        calc = fxn_generate_spherical_calc(voxel_size_mm, calc_diameter_mm)
    
    calc = rotate_calc(calc)

    ind_x, ind_y, ind_z = np.nonzero(calc)
    min_x, max_x = ind_x.min(), ind_x.max()
    min_y, max_y = ind_y.min(), ind_y.max()
    min_z, max_z = ind_z.min(), ind_z.max()

    max_range = max([max_x - min_x, max_y - min_y, max_z - min_z]) + 1

    center_x = (max_x + min_x) / 2
    center_y = (max_y + min_y) / 2
    center_z = (max_z + min_z) / 2
    new_min_x = round(center_x - max_range / 2)
    new_max_x = new_min_x + max_range - 1
    new_min_y = round(center_y - max_range / 2)
    new_max_y = new_min_y + max_range - 1
    new_min_z = round(center_z - max_range / 2)
    new_max_z = new_min_z + max_range - 1

    new_min_x = max(new_min_x, 1)
    new_max_x = new_max_x #min(new_max_x, calc.shape[0])
    new_min_y = max(new_min_y, 1)
    new_max_y = new_max_y #min(new_max_y, calc.shape[1])
    new_min_z = max(new_min_z, 1)
    new_max_z = new_max_z #min(new_max_z, calc.shape[2])

    calc = calc[new_min_x:new_max_x+1, new_min_y:new_max_y+1, new_min_z:new_max_z+1]

    if saveFLAG == 1:
        filename = f'/mnt/wd_disk/sunny_work/masses/calc_{shape}_{default_size_calc_mm*1000:04d}_{calc.shape[2]:03d}.raw'
        with open(filename, 'wb') as f:
            f.write(calc.tobytes())

        txtfilename = f'/mnt/wd_disk/sunny_work/masses/calc_{shape}_{default_size_calc_mm*1000:04d}_{calc.shape[2]:03d}.txt'
        with open(txtfilename, 'w') as f:
            f.write(f'{filename}\n')
            f.write(f'Calc size (mm): {default_size_calc_mm:.2f}\n')
            f.write(f'Voxel size (mm) = {voxel_size_mm:.3f}\n')

        print('Generated calc.')
    
    return calc
