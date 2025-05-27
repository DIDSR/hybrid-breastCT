import os
import argparse
import pandas as pd
import numpy as np
import time
from datetime import datetime
import matplotlib.pyplot as plt
from glob import glob
import math
import tigre
import tigre.algorithms as algs
import tigre.utilities.gpu as gpu

print(datetime.now(),flush=True)
# USER INPUT  ##################################################################################################

# Define paths and directories
code_dir = '/path/to/code/dir/'                                       # for code and relevant files
data_dir = '/path/to/data/dir/'                                       # larger storage containing patient data. hybrid images will output here
voi_center_dir  = '/path/to/voi/center/dir'                           # pre-generated viable VOI centers for each scan

# Path for energy_spectra and material_files
loc_spectra = os.path.join(code_dir, 'system_specific', 'energy_spectra', 'W60kVp_0.2mmGd.spc')
loc_material_files = {
    'calc': os.path.join(code_dir, 'system_specific', 'material_files', 'CalciumOxalate__5-120keV.mcgpu'),
    #'calc': os.path.join(code_dir, 'system_specific', 'material_files', 'Hydroxyapatite__5-120keV.mcgpu'),
    'adipose': os.path.join(code_dir, 'system_specific', 'material_files', 'adipose__5-120keV.mcgpu'),
    'glandular': os.path.join(code_dir, 'system_specific', 'material_files', 'glandular__5-120keV.mcgpu'),
    'csI': os.path.join(code_dir, 'system_specific', 'material_files', 'CesiumIodide__5-120keV.txt')
}

# Define path for output data
loc_save_MIPjpgs = os.path.join(data_dir, 'CalcPatches', 'MIPjpgs/')   # jpgs of MIP images
loc_save_patches = os.path.join(data_dir, 'CalcPatches', 'Patches/')   # patches = final images displayed in MIP or VOI, in raw form

# Desired voxel size of object (microns)
new_vx_um            = 30                  # voxel size of ray-tracing object

# Ray tracing projection parameters
energy_bin_size_keV  = 1                   # polyenergetic simulation occurs at every ___ keV intervals
csI_thickness_mm     = 0.45                # thickness of CsI scintillator layer
xsiz, ysiz           = 2048, 1536          # detector pixels in x direction, y direction
dexel_mm             = 0.150               # detector pixel size
vertical_offset_mm   = -30                 # vertical offset in system geometry

# VOI parameters
num_SPvois_perbreast = 50                  # max # signal present VOIs to extract per breast. depending on breast size and voi size, actual # may be smaller than this.
num_SAvois_perbreast = 50                  # max # signal absent VOIs to extract per breast. depending on breast size and voi size, actual # may be smaller than this.
voi_size_mm = 8                            # size of output VOI

# Set material densities
density = {
    'calc': 2.12 * 0.84,                   # Pure calcium oxalate * 0.84 (Warren 2013)
    'adipose': 0.920,
    'glandular': 1.035,
    'csI': 4.510
}

# Calc cluster parameters
calc_shape           = 'sphere'
clusterCenteredFLAG  = 0                   # [0 or 1] 0: Cluster not centered within VOI   1: Cluster centered within VOI

# Reconstruction parameters
kernels              = ['hann','shepp_logan', 'cosine', 'ram_lak']  # Refer to TIGRE documentation
recon_algorithms     = ['FDK']             #['FDK', 'SART', 'CGLS', 'MLEM'], refer to TIGRE documentation
iterations           = 30                  # called if using iterative recon algorithm
flagHU               = 0  # [0:mu  1:HU]   # output CT voxel values in mu or HU form
effective_energy_keV = 36.2                # effective energy of current x-ray spectrum
mu_water             = 0.298               # mu of water at effective energy of 36.2keV, used for conversion from mu to HU

# Flag for saving outputs into data_dir
savepatchesFLAG      = 1                

####################################################################################################################

# Read command line arguments
parser = argparse.ArgumentParser(description="Process input parameters for simulation.")
parser.add_argument('--cluster_diameter_mm', type=float, required=True, help="Diameter of the calcification cluster in mm.")
parser.add_argument('--num_calcs', type=int, required=True, help="Number of calcifications.")
parser.add_argument('--calc_diameter_mm', type=float, required=True, help="Diameter of each calcification in mm.")
parser.add_argument('--scanID', type=int, required=True, help="Patient scanID.")

args = parser.parse_args()

if args.cluster_diameter_mm <= 0:
    raise ValueError("Cluster diameter must be a positive number.")
if args.num_calcs <= 0:
    raise ValueError("Number of calcifications must be a positive integer.")
if args.calc_diameter_mm <= 0:
    raise ValueError("Calcification diameter must be a positive number.")
if not args.scanID:
    raise ValueError("scanID must be provided.")

cluster_diameter_mm = args.cluster_diameter_mm
num_calcs = args.num_calcs
calc_diameter_mm = args.calc_diameter_mm
scanID = args.scanID

# Load Doheny patient scan log
file_path = os.path.join(code_dir, 'system_specific', 'scanlog.xlsx')
df = pd.read_excel(file_path, sheet_name='log')

# Extract columns
VscanID = df.iloc[:, 1].values
VcropLocs = df.iloc[:, 2:8].values  # This is in voxels, based on the original CT recon (1 vx = 0.198 mm)
try:
    iscan = np.where(VscanID == scanID)[0][0]
except IndexError:
    raise ValueError(f"scanID {scanID} not found in the scan log.")

from utils import fxn_read_material_file, fxn_read_energy_spectrum

# Load energy spectrum, sampled at bin size "energy_bin_size_keV"
energy_levels_eV, photon_fluence = fxn_read_energy_spectrum(loc_spectra, energy_bin_size_keV)

# Using a dictionary to store mu values and energy levels for different materials
energy_levels = {}
mu_cm = {}

# Load mu values for Cesium Iodide
energy_levels['csI'], mu_cm['csI'] = fxn_read_material_file(loc_material_files['csI'], density['csI'], energy_levels_eV)
energy_levels['calc'], mu_cm['calc'] = fxn_read_material_file(loc_material_files['calc'], density['calc'], energy_levels_eV)
energy_levels['adipose'], mu_cm['adipose'] = fxn_read_material_file(loc_material_files['adipose'], density['adipose'], energy_levels_eV)
energy_levels['glandular'], mu_cm['glandular'] = fxn_read_material_file(loc_material_files['glandular'], density['glandular'], energy_levels_eV)

# Compute quantum detection efficiency (QDE) for Cesium Iodide
QDE_CsI = 1 - np.exp(-mu_cm['csI'] * (csI_thickness_mm / 10))

# Load Doheny detector MTF
datatable = pd.read_csv(os.path.join(code_dir,  'system_specific', 'Doheny_DetectorMTF_2x2_0.4mm_focalspotblur.csv'), header=None)
f_MTF = np.abs(datatable.iloc[:, 0].values)
mtf_MTF = datatable.iloc[:, 1].values

from calc_models.fxn_generate_calc import fxn_generate_calc

calc = fxn_generate_calc(new_vx_um*0.001,calc_diameter_mm,calc_shape,0)
half_dim = calc.shape[0]//2
print(f"Calc shape: {calc.shape} with voxel size {new_vx_um*0.001} mm.\n",flush=True)

# Convert dimensions from mm to voxels
cluster_radius_vx = round(cluster_diameter_mm / (2 * (0.001 * new_vx_um)))
voi_size_vx = round(voi_size_mm / (0.001 * new_vx_um))
margin_vx = round(np.ceil(cluster_radius_vx) + 1)
possible_distance_from_voi_center = round((voi_size_vx - 2 * margin_vx) / 2)  # set to [0,0,0] for cluster to be centered in VOI

cluster_diameter_vx = np.ceil(cluster_diameter_mm / (0.001 * new_vx_um))
logical_sphere_center = np.array([cluster_radius_vx, cluster_radius_vx, cluster_radius_vx])

# Create a grid of coordinates (x, y, z)
x, y, z = np.mgrid[1:cluster_diameter_vx+1, 1:cluster_diameter_vx+1, 1:cluster_diameter_vx+1]

# Calculate the distances from the center
distances_from_center = np.sqrt((x - logical_sphere_center[0])**2 + (y - logical_sphere_center[1])**2 + (z - logical_sphere_center[2])**2)
logical_sphere = distances_from_center <= (cluster_radius_vx - half_dim)

# Find indices of all possible calc centers
true_x, true_y, true_z = np.where(logical_sphere)

from utils import fxn_load_seg

start_time = time.time()
# Load patient segmentation volume
seg_volume, Nxyz, original_vx_um, labels = fxn_load_seg(scanID, data_dir)
print(f"Elapsed time: {time.time() - start_time:.1f} seconds",flush=True)

from utils import fxn_load_projections_and_geometry

start_time = time.time()
#Load patient projection images and define geometry for reconstruction
prjstack, geo, ang = fxn_load_projections_and_geometry(scanID, data_dir, xsiz, ysiz, dexel_mm, vertical_offset_mm);
print(f"Elapsed time: {time.time() - start_time:.1f} seconds",flush=True)

from utils import fxn_crop_volume

FLAGchestwall = 0
seg_volume_cropped, shift_mm_3D = fxn_crop_volume(iscan, seg_volume, original_vx_um, Nxyz, VcropLocs, FLAGchestwall)

from utils import fxn_upsample_volume_in_sections

start_time = time.time()
seg_volume_HR = fxn_upsample_volume_in_sections(seg_volume_cropped, original_vx_um, new_vx_um)
print(f"Elapsed time: {time.time() - start_time:.1f} seconds",flush=True)

from utils import fxn_getVOIcenters

voi_centers_mm_SP, voi_centers_mm_SA, num_SPvois_perbreast, num_SAvois_perbreast = fxn_getVOIcenters(scanID, voi_size_mm, num_SPvois_perbreast, num_SAvois_perbreast, voi_center_dir)

from utils import fxn_insert_calc_cluster_new

volume_with_calcs = np.zeros(seg_volume_HR.shape, dtype=np.uint8)

print(f'Inserting {num_SPvois_perbreast} calc clusters: {cluster_diameter_mm} mm cluster diameter, {calc_diameter_mm:.2f} mm calc diameter, {num_calcs} #calcs...',flush=True)
start_time = time.time()

for ivoi in range(num_SPvois_perbreast):
    if ivoi % 10 == 0:
        print(f'{ivoi}, ', end='',flush=True)

    voi_center_mm_original = voi_centers_mm_SP[ivoi]
    # Adjust VOI centers to be relative to cropped breast volume
    voi_center = np.round((voi_center_mm_original - shift_mm_3D) / (0.001 * new_vx_um)).astype(int)

    if clusterCenteredFLAG == 0:
        # Randomly determine a cluster center within the VOI
        cluster_center = voi_center + np.random.randint(-possible_distance_from_voi_center, possible_distance_from_voi_center + 1, size=(3,))
    else:
        cluster_center = voi_center

    volume_with_calcs = fxn_insert_calc_cluster_new(volume_with_calcs, seg_volume_HR, true_x, true_y, true_z, logical_sphere_center, cluster_center, num_calcs, calc, half_dim, labels)

print('done.',flush=True)
print(f"Elapsed time: {time.time() - start_time:.1f} seconds",flush=True)

# Check that calc was inserted into breast parenchyma only
if np.any(volume_with_calcs == 1) or np.any(volume_with_calcs == 2):
    print('ERROR: calc was added into air or skin.',flush=True)

# Calculate COR offset due to cropping volume
rowstart = int(VcropLocs[iscan, 0])
rowend = int(VcropLocs[iscan, 1])
colstart = int(VcropLocs[iscan, 2])
colend = int(VcropLocs[iscan, 3])
zstart = int(VcropLocs[iscan, 4])
zstart_mm = zstart*(0.001 * original_vx_um)
zend = int(VcropLocs[iscan, 5])

COR_mm = np.array([Nxyz[1] / 2 * (0.001 * original_vx_um), Nxyz[2] / 2 * (0.001 * original_vx_um)])
newCOR_mm = np.array([(rowstart + rowend) / 2 * (0.001 * original_vx_um), (colstart + colend) / 2 * (0.001 * original_vx_um)])

# Center of rotation offset
if FLAGchestwall==1:
    CORoffset_mm_3D = -np.array([0, newCOR_mm[0] - COR_mm[0], newCOR_mm[1] - COR_mm[1]])
else: #smaller cropped volume to save memory
    COR_mm_z = np.squeeze(np.array([Nxyz[2] / 2 * (0.001 * original_vx_um)]))
    newCOR_mm_z = np.squeeze(np.array([(zstart + zend) / 2 * (0.001 * original_vx_um)]))
    CORoffset_mm_3D = -np.array([newCOR_mm_z - COR_mm_z, newCOR_mm[0] - COR_mm[0], newCOR_mm[1] - COR_mm[1]])

from utils import fxn_alter_geometry

geo_raytracing = fxn_alter_geometry(geo, volume_with_calcs, CORoffset_mm_3D, new_vx_um, vertical_offset_mm, zstart_mm)

from utils import fxn_mtf_blur

# Re-orient object so that projection aligns with patient projection
volume_with_calcs_flipped = np.rot90(volume_with_calcs, k=2, axes = (1,2))

# Multi-energy ray tracing simulation
print(f'Starting ray tracing for {len(energy_levels_eV)} energy levels (every {energy_bin_size_keV} keV between {energy_levels_eV[0]/1000:.2f}-{energy_levels_eV[-1]/1000:.2f} keV): ',flush=True)
start_time = time.time()

sum_prjstack_calcs = np.zeros((len(ang), geo_raytracing.nDetector[0], geo_raytracing.nDetector[1]))
total_I_o = 0

# Iterate through energy bins in polyenergetic spectrum
for iE, energy_eV in enumerate(energy_levels_eV):
    energy_eV = energy_levels_eV[iE]
    print(f'{energy_eV/1000:.2f}, ', end='',flush=True)

    # New ray tracing object at every energy
    volume_with_calcs_mu_mm = np.zeros_like(volume_with_calcs,dtype='float32')
    volume_with_calcs_mu_mm[volume_with_calcs_flipped == 101] = 0.1 * abs(mu_cm['calc'][iE] - mu_cm['glandular'][iE])  # calc minus glandular
    volume_with_calcs_mu_mm[volume_with_calcs_flipped == 100] = 0.1 * abs(mu_cm['calc'][iE] - mu_cm['adipose'][iE])  # calc minus adipose

    # Perform ray tracing projection (assume single precision and interpolated method)
    prjstack_calcs = tigre.Ax(np.float32(volume_with_calcs_mu_mm), geo_raytracing, ang, method='interpolated')
    sum_prjstack_calcs += photon_fluence[iE] * QDE_CsI[iE] * np.exp(-prjstack_calcs)  # I_o * QDE * exponent of ray tracing output

    total_I_o += photon_fluence[iE] * QDE_CsI[iE]

print('done. \n',flush=True)
print(f"Elapsed time: {time.time() - start_time:.1f} seconds",flush=True)

# Blur calc projections with Doheny detector MTF
print('\nBlurring all ray-tracing projections with detector MTF...',flush=True)
start_time = time.time()
blurred_prjstack_calcs = fxn_mtf_blur(sum_prjstack_calcs, dexel_mm, f_MTF, mtf_MTF)
print('done. \n',flush=True)
print(f"Elapsed time: {time.time() - start_time:.1f} seconds",flush=True)

# Log normalization: results in line integral of attenuation coeffs over all energies with MTF and QDE incorporated
prjstack_calcs_blur_QDE_lognorm = -np.log(blurred_prjstack_calcs / total_I_o)

# Add blurred calc projections to patient projection images
hybrid_prjstack = (prjstack + prjstack_calcs_blur_QDE_lognorm).astype(np.float32)

del prjstack_calcs_blur_QDE_lognorm, volume_with_calcs_mu_mm, sum_prjstack_calcs, blurred_prjstack_calcs, prjstack_calcs

from utils import fxn_extract_and_save_vois, fxn_write_metadata

# RECONSTRUCTION AND VOI EXTRACTION
recon_size_mm = geo.dVoxel[0]  # Isotropic voxels
voi_size_vx = round(voi_size_mm / (recon_size_mm))
voi_halfdim_vx = round(voi_size_mm / (2 * recon_size_mm))

# DEFAULT UNLESS REDEFINED
kernel = 'nokernel'
folder_suffix = kernel

# Loop through each reconstruction algorithm
for recon_alg in recon_algorithms:
    if recon_alg == 'FDK':
        for kernel in kernels:
            start_time = time.time()
            print(f"Starting {recon_alg} with {kernel} kernel... ", flush=True)
            rec = 10 * algs.fdk(hybrid_prjstack, geo, ang, filter=kernel)
            if flagHU == 1:
                rec = 1000 * (rec - mu_water) / mu_water
            rec = np.rot90(rec, k=2, axes=(1, 2))
            if savepatchesFLAG == 1:
                if kernel == 'ram_lak': #for ease of parsing through folders later
                    kernel = 'ramlak'
                if kernel == 'shepp_logan':
                    kernel = 'shepplogan'
                folder_suffix = kernel
                fxn_extract_and_save_vois(rec, recon_alg, folder_suffix, scanID, loc_save_patches, loc_save_MIPjpgs, calc_diameter_mm, cluster_diameter_mm, num_calcs, density, num_SPvois_perbreast, voi_centers_mm_SP, voi_centers_mm_SA, recon_size_mm, voi_halfdim_vx, flagHU)
                fxn_write_metadata(loc_save_patches, loc_spectra, effective_energy_keV, loc_material_files, density, energy_bin_size_keV, csI_thickness_mm, vertical_offset_mm, new_vx_um, calc_shape, calc_diameter_mm, num_calcs, cluster_diameter_mm, clusterCenteredFLAG,num_SPvois_perbreast, num_SAvois_perbreast, voi_size_mm, voi_size_vx, recon_alg, folder_suffix, kernel, flagHU, mu_water, code_dir, geo, ang)
            elapsed_time = time.time() - start_time
            print(f"Kernel: {kernel} reconstruction completed in {elapsed_time:.1f} seconds", flush=True)
    elif recon_alg == 'SART':
        start_time = time.time()
        print(f"Starting {recon_alg} reconstruction for {iterations} iterations... ", flush=True)
        rec = 10 * algs.ossart(hybrid_prjstack, geo, ang, iterations)
        if flagHU == 1:
            rec = 1000 * (rec - mu_water) / mu_water
        rec = np.rot90(rec, k=2, axes=(1, 2))
        if savepatchesFLAG == 1:
           fxn_extract_and_save_vois(rec, recon_alg, folder_suffix, scanID, loc_save_patches, loc_save_MIPjpgs, calc_diameter_mm, cluster_diameter_mm, num_calcs, density, num_SPvois_perbreast, voi_centers_mm_SP, voi_centers_mm_SA, recon_size_mm, voi_halfdim_vx, flagHU)
           fxn_write_metadata(loc_save_patches, loc_spectra, effective_energy_keV, loc_material_files, density, energy_bin_size_keV, csI_thickness_mm, vertical_offset_mm, new_vx_um, calc_shape, calc_diameter_mm, num_calcs, cluster_diameter_mm, clusterCenteredFLAG,num_SPvois_perbreast, num_SAvois_perbreast, voi_size_mm, voi_size_vx, recon_alg, folder_suffix, kernel, flagHU, mu_water, code_dir, geo, ang)
        elapsed_time = time.time() - start_time
        print(f"{recon_alg} reconstruction completed in {elapsed_time:.1f} seconds", flush=True)
    elif recon_alg == 'MLEM':
        start_time = time.time()
        print(f"Starting {recon_alg} reconstruction for {iterations} iterations... ", flush=True)
        rec = 10 * algs.mlem(hybrid_prjstack, geo, ang, iterations)
        if flagHU == 1:
            rec = 1000 * (rec - mu_water) / mu_water
        rec = np.rot90(rec, k=2, axes=(1, 2))
        if savepatchesFLAG == 1:
           fxn_extract_and_save_vois(rec, recon_alg, folder_suffix, scanID, loc_save_patches, loc_save_MIPjpgs, calc_diameter_mm, cluster_diameter_mm, num_calcs, density, num_SPvois_perbreast, voi_centers_mm_SP, voi_centers_mm_SA, recon_size_mm, voi_halfdim_vx, flagHU)
           fxn_write_metadata(loc_save_patches, loc_spectra, effective_energy_keV, loc_material_files, density, energy_bin_size_keV, csI_thickness_mm, vertical_offset_mm, new_vx_um, calc_shape, calc_diameter_mm, num_calcs, cluster_diameter_mm, clusterCenteredFLAG,num_SPvois_perbreast, num_SAvois_perbreast, voi_size_mm, voi_size_vx, recon_alg, folder_suffix, kernel, flagHU, mu_water, code_dir, geo, ang)
        elapsed_time = time.time() - start_time
        print(f"{recon_alg} reconstruction completed in {elapsed_time:.1f} seconds", flush=True)
    elif recon_alg == 'CGLS':
        start_time = time.time()
        print(f"Starting {recon_alg} reconstruction for {iterations} iterations... ", flush=True)
        rec = 10 * algs.cgls(hybrid_prjstack, geo, ang, iterations)
        if flagHU == 1:
            rec = 1000 * (rec - mu_water) / mu_water
        rec = np.rot90(rec, k=2, axes=(1, 2))
        if savepatchesFLAG == 1:
           fxn_extract_and_save_vois(rec, recon_alg, folder_suffix, scanID, loc_save_patches, loc_save_MIPjpgs, calc_diameter_mm, cluster_diameter_mm, num_calcs, density, num_SPvois_perbreast, voi_centers_mm_SP, voi_centers_mm_SA, recon_size_mm, voi_halfdim_vx, flagHU)
           fxn_write_metadata(loc_save_patches, loc_spectra, effective_energy_keV, loc_material_files, density, energy_bin_size_keV, csI_thickness_mm, vertical_offset_mm, new_vx_um, calc_shape, calc_diameter_mm, num_calcs, cluster_diameter_mm, clusterCenteredFLAG,num_SPvois_perbreast, num_SAvois_perbreast, voi_size_mm, voi_size_vx, recon_alg, folder_suffix, kernel, flagHU, mu_water, code_dir, geo, ang)
        elapsed_time = time.time() - start_time
        print(f"{recon_alg} reconstruction completed in {elapsed_time:.1f} seconds", flush=True)

print('done.\n\n',flush=True)
print(datetime.now(),flush=True)
