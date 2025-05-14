# hybrid-breastCT
This repository contains code for generating hybrid breast CT images combining simulated pathology with real anatomical background.
![image](https://github.com/user-attachments/assets/539ca4a7-c2c6-46b8-8069-080405025048)

Each of these steps is embedded in ```hybrid-simulation.py```, found in the ```code``` directory.
## Code Features


### Define calc cluster model:
- Generates calcification clusters with flexible configuration of:
  - Number of calcs
  - Cluster diameter
  - Calc diameter
- Identifies background tissue type at insertion site (VOI centers) based on breast segmentation volume. Object values assigned according to density and linear attenuation coefficients of calcification and background tissue.
- OPTIONAL: Object may be cropped (to reduce computational burden) and resampled (to increase object resolution)
  

### Simulate x-ray projections of clusters:
- Supports polyenergetic X-ray spectra with energy binning.
- Applies detector Quantum Detection Efficiency (QDE) using CsI thickness.
- Generates projections of inserted calcs using [TIGRE reconstruction toolbox](https://github.com/CERN/TIGRE/tree/master) GPU-accelerated ray tracing.


### Blurred x-ray projections of clusters: 
- Applies detector blur using Modulation Transfer Function (MTF).


### Patient x-ray projections (anatomical background):
- Loads real segmentation and patient breast CT projection data.
- Extracts scan parameters, system geometry, and viable VOI centers for calc insertion.


### Hybrid projection images:
- Combines simulated calc projections with real patient projections to create hybrid projection images.


### Hybrid breast CT volume:
- Utilizes [TIGRE](https://github.com/CERN/TIGRE/tree/master) reconstruction algorithms to generate hybrid CT volume.
  - Multiple algorithms supported: FDK, SART, CGLS, and MLEM.
  - Apodization filters can be modified for FDK algorithm (e.g., Ram-Lak, Hann).


### Hybrid breast CT ROI/VOI:
- Extracts signal-present and signal-absent VOIs and maximum intensity projections (MIPs) for training, testing, or reader studies.


## Requirements
1. Install [TIGRE](https://github.com/CERN/TIGRE/blob/master/Frontispiece/python_installation.md) for Python
2. Install required packages in a new Python environment
   
```
conda create --name bct_env python=3.11
conda activate bct_env
pip install -r requirements.txt
```
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


## Citation

"_Hybrid simulation of breast CT for assessing microcalcification detectability_". Lyu SH, Makeev A, Li D, Badal A, Hernandez AM, Boone JM, Glick SJ.
Journal of Medical Imaging, 2025.

## More information

[Catalog of Regulatory Science Tools to Help Assess New Medical Devices](https://www.fda.gov/medical-devices/science-and-research-medical-devices/catalog-regulatory-science-tools-help-assess-new-medical-devices)


### About the Catalog of Regulatory Science Tools

The enclosed tool is part of the Catalog of Regulatory Science Tools, which provides a peer-reviewed resource for stakeholders to use where standards and qualified Medical Device Development Tools (MDDTs) do not yet exist. These tools do not replace FDA-recognized standards or MDDTs. This catalog collates a variety of regulatory science tools that the FDA's Center for Devices and Radiological Health's (CDRH) Office of Science and Engineering Labs (OSEL) developed. These tools use the most innovative science to support medical device development and patient access to safe and effective medical devices. If you are considering using a tool from this catalog in your marketing submissions, note that these tools have not been qualified as Medical Device Development Tools and the FDA has not evaluated the suitability of these tools within any specific context of use. You may request feedback or meetings for medical device submissions as part of the Q-Submission Program.

For more information about the Catalog of Regulatory Science Tools, email OSEL_CDRH@fda.hhs.gov. 

## Disclaimer
This software and documentation (the "Software") were developed at the Food and Drug Administration (FDA) by employees of the Federal Government in the course of their official duties. Pursuant to Title 17, Section 105 of the United States Code, this work is not subject to copyright protection and is in the public domain. Permission is hereby granted, free of charge, to any person obtaining a copy of the Software, to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, or sell copies of the Software or derivatives, and to permit persons to whom the Software is furnished to do so. FDA assumes no responsibility whatsoever for use by other parties of the Software, its source code, documentation or compiled executables, and makes no guarantees, expressed or implied, about its quality, reliability, or any other characteristic. Further, use of this code in no way implies endorsement by the FDA or confers any advantage in regulatory decisions. Although this software can be redistributed and/or modified freely, we ask that any derivative works bear some notice that they are derived from it, and any modified versions bear some notice that they have been modified.
