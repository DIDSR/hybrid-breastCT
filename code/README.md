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
5. Customize ```hybrid-simulation.py``` and ```utils.py``` for your breast CT system. More info below:

## Getting Started with Other Breast CT Systems

This simulation framework is designed to be adaptable across different breast CT platforms. While the original implementation uses data from the [Doheny Breast CT system](https://pmc.ncbi.nlm.nih.gov/articles/PMC4376760/), users can integrate their own patient datasets by updating a few key components of the pipeline.

Below are the parameters and functions you’ll need to modify or replace in ```hybrid-simulation.py``` or ```utils.py``` for compatibility with your system:


### 1. File Paths and Data Formats
Update file I/O functions to match your system’s file structure and formats:
- Replace **segmentation loading** logic in `fxn_load_seg(...)` if your data is in DICOM, NIfTI, or other formats.
- Adapt `fxn_load_projections_and_geometry(...)` to read your projection files, angle data, and calibration metadata.
- Update `fxn_getVOIcenters(...)` to generate or supply VOI center coordinates (e.g., from external software or manually). VOIs should be contained within the breast.


### 2. Geometry and Detector Configuration
Your CT system’s geometric setup must be reflected in the code:
- Update parameters such as:
  - `DSD` (Source-to-detector distance)
  - `DSO` (Source-to-object distance)
  - Detector pixel size and resolution
  - Voxel size of reconstructed volume
- Modify how these are initialized in `fxn_load_projections_and_geometry(...)` and `fxn_alter_geometry(...)`.


### 3. X-ray Spectrum and Material Properties
The simulation depends on accurate modeling of energy-dependent attenuation:
- Replace `W60kVp_0.2mmGd.spc` with your system’s X-ray spectrum in `.txt` or `.spc` format.
- Provide your own attenuation coefficient files for:
  - Calcification material (e.g., calcium oxalate or hydroxyapatite)
  - Adipose and glandular tissue
  - Detector material (e.g., CsI)
- Update these inputs in the `material_files/` directory and ensure compatibility with `fxn_read_material_file(...)`.


### 4. Tissue Segmentation Labels
Ensure the code uses the correct label IDs for your dataset:

```python
labels = {
    'air': 0,
    'adipose': 1,
    'glandular': 2,
    'skin': 5
}
```


### 5. Projection Preprocessing
Check if your projection data is already normalized:
- If raw, you may need to apply flat-field correction (e.g., using I₀ or flood-field data).
- If preprocessed, ensure the log-normalization step in `fxn_load_projections_and_geometry(...)` is either updated or skipped.


### 6. Detector MTF
If you wish to simulate realistic detector blur:
- Replace `Doheny_DetectorMTF_2x2_0.4mm_focalspotblur.csv` with your system’s measured MTF curve.
- Format the file as `[frequency (lp/mm), MTF value]` and update the path used in the script.


### 7. Optional Reconstruction
If your system does not use TIGRE for reconstruction:
- You can still use the calc simulation + projection pipeline.
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
