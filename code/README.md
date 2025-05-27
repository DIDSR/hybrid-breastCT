## Setup
1. To clone this repository to your local machine, open a terminal and run:
```
git clone https://github.com/DIDSR/hybrid-breastCT.git
```

2. Activate a new Conda environment:
```
conda create --name bct_env python=3.11
conda activate bct_env
```
3. Install [TIGRE](https://github.com/CERN/TIGRE/blob/master/Frontispiece/python_installation.md) for Python
4. Install required Python packages using pip: 
```
pip install -r requirements.txt
```
5. Customize ```hybrid-simulation.py```, ```utils.py```, and directory `/system-specific` for your breast CT system. More info below:

## Requirements For Implementation With Other Breast CT Systems

This simulation framework is designed to be adaptable across different breast CT platforms. While the original implementation uses data from the [Doheny Breast CT system](https://pmc.ncbi.nlm.nih.gov/articles/PMC4376760/), users can integrate their own patient datasets. The following data/information is required for implementation with your system:

#### REQUIREMENTS 
- Patient breast CT projection images (post-flatfielding and other corrections).
- Patient breast CT segmentations: uint8 image volume with each voxel labeled as "air", "adipose tissue", "fibroglandular tissue" and "skin". For reference, the labels used in this work are:
labels = {
    'air': 0,
    'adipose': 1,
    'glandular': 2,
    'skin': 5
}
- X-ray spectrum for your breast CT system
- Information about system geometry and acquisition parameters, for example:
  - `DSD` (Source-to-detector distance)
  - `DSO` (Source-to-object distance)
  - Detector pixel size
  - Number of projection angles in the acquisition
  - Angles of acquisition
- 1D or 2D detector modulation transfer function (MTF) for your system
    note: measure at face of detector

#### OPTIONAL
- CT reconstruction algorithm
- Attenuation coefficient files for materials not found in `/system_specific/material_files`

## Steps: Uploading Data and Customizing Code For Other Breast CT Systems
Below are the directories, parameters, and functions you’ll need to modify or replace for compatibility with your system:

### 1. Directories and Data
- In `hybrid-simulation.py`, define code directory (clone this repository as explained in `setup`):

```
code_dir = '/path/to/code/dir/'                                       # for code and relevant files
```

- In `hybrid-simulation.py`, define or make data directory:

```
data_dir = '/path/to/data/dir/'                                       # larger storage containing patient data. hybrid images will output here
```

- Upload patient projection images, patient segmentations, and any patient metadata  to `data_dir`.
- Adapt `fxn_load_projections_and_geometry(...)` and `fxn_load_seg(...)` in `utils.py` to read patient data and reflect the geometry of your CT system.
- Adapt `fxn_alter_geometry(...)` in `utils.py` if reconstruction geometry differs from acquisition geometry (e.g. binning, recon voxel size).
- Refer to [TIGRE documentation](https://github.com/CERN/TIGRE/blob/master/Python/demos/d01_CreateGeometry.py) for additional help setting geometry 

### 2. X-ray Energy Spectrum and Material Files
- Replace `system_specific/energy_spectra/W60kVp_0.2mmGd.spc` with your system’s X-ray spectrum in `.txt` or `.spc` format.
- In `system_specific/material_files`, use or provide your own attenuation coefficient files for:
  - Calcification material (e.g., calcium oxalate or hydroxyapatite)
  - Adipose and glandular tissue
  - Detector material (e.g., CsI)
- Instructions for generating your own material files and energy spectra are found in `system_specific/material_files/README.md` and `system_specific/energy_spectra/README.md`
- Adapt `fxn_read_material_file(...)` and `fxn_read_energy_spectrum(...)` in `utils.py` to read material files and energy spectra.
 
### 3. Edit parameters in `hybrid-simulation.py`
- In `hybrid-simulation.py`, edit relevant parameters in section
```
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
effective_energy_keV = 36.2                # effective energy of current x-ray spectrum, used for conversion from mu to HU.
mu_water             = 0.298               # mu of water at effective energy of 36.2keV, used for conversion from mu to HU

# Flag for saving outputs into data_dir
savepatchesFLAG      = 1                

####################################################################################################################
```
### 4. Scan Log, VOI Centers, MTF
- If patient database contains more than one patient scan, a scan log may be useful. `system_specific/scanlog.xlsx` contains sample scan log used for Doheny patient database. Replace with relevant files for your system.
- Update `fxn_getVOIcenters(...)` in `utils.py` to generate or supply viable VOI center coordinates for each scan. VOIs surrounding VOI centers should be contained within the breast.
- If you wish to simulate realistic detector blur:
  - Replace `system_specific/Doheny_DetectorMTF_2x2_0.4mm_focalspotblur.csv` with your system’s measured detector MTF.
  - Format the file as `[frequency (lp/mm), MTF value]` and update the path used in the script.
  - NOTE: `fxn_mtf_blur(prjstack_calcs, dexel_mm, f_MTF, mtf_MTF)` in `utils.py` assumes the measured MTF is 1D, and implements MTF blurring first in the X direction, then in the Y direction. If your measured MTF is in 2D, then `fxn_mtf_blur` should be modified to implement 2D MTF blur.
    
#### NOTE: on ray-tracing OBJECT
- To ensure optimal resolution of simulated ray-tracing object, a small object voxel size should be defined in `hybrid-simulation.py`.
```
# Desired voxel size of ray-tracing object (microns)
new_vx_um            = 30                
```
- To achieve desired voxel size, the original segmentation volume is upsampled in `hybrid-simulation.py`
```
seg_volume_HR = fxn_upsample_volume_in_sections(seg_volume_cropped, original_vx_um, new_vx_um)
```
- Because the upsampled segmentation volume can be exceedingly large, it is advised to **crop the original segmentation volume and ray-tracing object prior to upsampling**. The cropped volume should contain the majority of the breast but minimize background (air) and chest wall. Simulated lesions are not to be inserted in the chest wall, so these areas are not within the FOV during ray tracing.
- Cropping occurs in  `hybrid-simulation.py`:
```
seg_volume_cropped, shift_mm_3D = fxn_crop_volume(iscan, seg_volume, original_vx_um, Nxyz, VcropLocs, FLAGchestwall)
```
based on manually-identified cropping locations defined in `system_specific/scanlog.xlsx`.
  - Columns C-H: first and last rows, cols, slices in X, Y, Z respectively were manually identified from breast CT reconstructed volumes for each scanID
  - In `hybrid-simulation.py`, columns C:H are read in as vector: VcropLocs.
- Then, the cropped volume `seg_volume_cropped` can be upsampled, and undergo ray-tracing simulation.
  
### 5. Projection Preprocessing
Check if your projection data is already normalized:
- If raw, you may need to apply flat-field correction (e.g., using I₀ or flood-field data).
- If preprocessed, ensure the log-normalization step in `fxn_load_projections_and_geometry(...)` is either updated or skipped.


### 6. Optional: Reconstruction
If your system does not use TIGRE for CT reconstruction:
- You can still use the calc simulation + patient projection pipeline in `hybrid-simulation.py`.
- Export `hybrid_prjstack` and reconstruct externally using your own algorithms.

## Usage
```hybrid-simulation.py``` can be run in command line. First, activate Python virtual environment. Then, run the following:
```
python hybrid-simulation.py --cluster_diameter_mm <CLUSTER_DIAMETER_MM> --num_calcs <NUM_CALCS> --calc_diameter_mm <CALC_DIAMETER_MM> --scanID <SCAN_ID>
```

where

```<CLUSTER_DIAMETER_MM>``` = diameter of spherical cluster in mm

```<NUM_CALCS>``` = # calcs in cluster

```<CALC_DIAMETER_MM>``` = calc diameter in mm

```<SCAN_ID>``` = scanID of patient images. Calc clusters will be inserted into these images.

## Example

```
python hybrid-simulation.py --cluster_diameter_mm 4.0 --num_calcs 5 --calc_diameter_mm 1.0 --scanID 2878
```
