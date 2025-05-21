import numpy as np

def fxn_generate_spherical_calc(voxel_size_mm, calc_diameter_mm):
    # Calculate the size of the 3D matrix
    matrix_size = round(1.5 * (calc_diameter_mm / voxel_size_mm))
    
    # Create a 3D matrix with all zeros
    calc = np.zeros((matrix_size, matrix_size, matrix_size), dtype=np.uint8)
    
    # Calculate the center of the matrix
    center = round(matrix_size / 2)
    
    # Create a meshgrid for the matrix
    x, y, z = np.meshgrid(np.arange(1, matrix_size+1),
                          np.arange(1, matrix_size+1),
                          np.arange(1, matrix_size+1),
                          indexing='ij')
    
    # Calculate the distance from each point to the center
    distance = np.sqrt((x - center)**2 + (y - center)**2 + (z - center)**2)
    
    # Set the values inside the sphere to 1
    calc[distance <= (calc_diameter_mm / 2) / voxel_size_mm] = 1
    
    return calc
