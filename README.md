# hybrid-breastCT

This repository contains a Python-based framework for generating hybrid breast CT images by combining simulated microcalcification clusters with real patient breast CT projection data.

![image](https://github.com/user-attachments/assets/539ca4a7-c2c6-46b8-8069-080405025048)



---

## Overview

The hybrid simulation framework performs the following steps:

#### 1. Define calc cluster model
- Generates calcification clusters with configurable:
  - Number of calcs
  - Cluster diameter
  - Calc diameter
    
#### 2. Simulate X-ray projections of clusters
- Supports polyenergetic spectra with energy binning
- Uses [TIGRE](https://github.com/CERN/TIGRE) for GPU-accelerated ray tracing

#### 3. Apply detector blur
- Uses detector Modulation Transfer Function (MTF)

#### 4. Load patient projection data
- Uses real breast CT projections and segmentation volumes
- Extracts geometry and viable VOI centers

#### 5. Generate hybrid projections
- Combines simulated lesion projections with real projections

#### 6. Reconstruct hybrid CT volume
- Uses TIGRE reconstruction algorithms:
  - FDK, SART, CGLS, MLEM
- Supports multiple reconstruction filters

#### 7. Extract VOIs and MIPs
- Generates signal-present and signal-absent volumes
- Outputs MIPs for analysis or reader studies

---

## Quick Start

#### 1. Clone repository
```bash
git clone https://github.com/DIDSR/hybrid-breastCT.git
cd hybrid-breastCT
```

#### 2. Create environment
```bash
python -m venv bct_env
source bct_env/bin/activate
pip install -r requirements.txt
```

#### 3.  Install [TIGRE](https://github.com/CERN/TIGRE/blob/master/Frontispiece/python_installation.md) for Python (requires GPU)

#### 4. Configuration

All system-specific inputs and paths are defined in: `configs/doheny_example.yaml`

##### Required paths
```yaml
paths:
  patient_data_dir: /path/to/patient_data_root
  output_dir: /path/to/output
  voi_center_dir: /path/to/voi_centers
```

##### Required files
```yaml
files:
  spectrum: ...
  detector_mtf: ...
  scanlog: ...
  material_files:
    calc: ...
    adipose: ...
    glandular: ...
    csI: ...
```

#### 5. Validate Inputs

Check that all required files and paths are accessible:

```bash
python -m hybrid_bct.cli validate --config configs/doheny_example.yaml
```

#### 6. Run Simulation

Run the hybrid simulation pipeline:

```bash
python -m hybrid_bct.cli run \
  --config configs/doheny_example.yaml \
  --scan-id 2878 \
  --cluster-diameter-mm 4.0 \
  --num-calcs 5 \
  --calc-diameter-mm 1.0
```

Outputs from the pipeline include:
- VOI `.npy` volumes
- MIP `.npy` arrays
- MIP `.jpg` images
- Metadata files

## Requirements For Implementation With Other Breast CT Systems

This simulation framework is designed to be adaptable across different breast CT platforms. While the original implementation uses data from the [Doheny Breast CT system](https://pmc.ncbi.nlm.nih.gov/articles/PMC4376760/), users can integrate their own patient datasets. The following data/information is required for implementation with your system:

To run the pipeline with your own system, you must supply:

#### REQUIREMENTS 
- Patient breast CT projection images 
- Patient breast CT segmentation volumes
- X-ray spectrum
- System geometry parameters
- 1D or 2D detector modulation transfer function (MTF) for your system
  - NOTE: measure at face of detector
- VOI center files

#### OPTIONAL
- CT reconstruction algorithm
- Attenuation coefficient files for additional materials

To adapt this framework to a new breast CT system, users must modify **two components**:

#### 1. System-specific loading code

Modify: `hybrid_bct/systems/doheny.py`. 
Specifically, replace these functions with implementations compatible with your system:
- `load_seg_doheny(...)`
- `load_projections_and_geometry_doheny(...)`

These functions are responsible for:
- loading your segmentation volume
- loading your projection data
- defining TIGRE geometry (`geo`)
- defining projection angles (`ang`)

This is the where system-specific data formats and geometry are handled.

#### 2. Configuration file

Update: `configs/doheny_example.yaml` to reflect your system:

- file paths (patient data, VOI centers, outputs)
- detector parameters (pixel size, dimensions, offsets)
- X-ray spectrum
- material files
- reconstruction settings

This file controls all runtime parameters and inputs without modifying code.


### What you do NOT need to modify

You generally do NOT need to change:

- `pipeline.py` (core workflow)
- simulation modules (`simulation/`)
- reconstruction code (`reconstruction/`)
- CLI interface


### Data format expectations

#### Segmentation volume
- 3D numpy array (uint8)
- Labels:
```python
labels = {
    'air': 0,
    'adipose': 1,
    'glandular': 2,
    'skin': 5
}
```
#### Projection data
- preprocessed projection images
- log-normalized or raw (consistent with your loader)
  
#### Geometry
Defined using TIGRE: refer to https://github.com/CERN/TIGRE/blob/master/Python/demos/d01_CreateGeometry.py

### Notes on resolution and performance
- The segmentation volume is **upsampled** before simulation
- Cropping is recommended to reduce memory usage
- GPU acceleration (TIGRE) is required for projection simulation

### Optional customization
Advanced users may modify:
- `simulation/projection.py` --> forward projection model
- `simulation/blur.py` --> detector blur model
- `reconstruction.py` --> reconstruction algorithms
  
## Citation

"_Hybrid simulation of breast CT for assessing microcalcification detectability_". Lyu SH, Makeev A, Li D, Badal A, Hernandez AM, Boone JM, Glick SJ.
Journal of Medical Imaging, 2025. https://doi.org/10.1117/1.JMI.12.S2.S22015

## Regulatory Science Tool Reference

RST Reference Number: RST26MD04.01

Date of Publication: 05/04/2026

Recommended Citation: U.S. Food and Drug Administration. (2026). Hybrid Simulation Framework for Breast CT Virtual Imaging Trials (RST26MD04.01). https://cdrh-rst.fda.gov/hybrid-simulation-framework-breast-ct-virtual-imaging-trials

#### About the Catalog of Regulatory Science Tools

The enclosed tool is part of the[Catalog of Regulatory Science Tools to Help Assess New Medical Devices](https://www.fda.gov/medical-devices/science-and-research-medical-devices/catalog-regulatory-science-tools-help-assess-new-medical-devices)
, which provides a peer-reviewed resource for stakeholders to use where standards and qualified Medical Device Development Tools (MDDTs) do not yet exist. These tools do not replace FDA-recognized standards or MDDTs. This catalog collates a variety of regulatory science tools that the FDA's Center for Devices and Radiological Health's (CDRH) Office of Science and Engineering Labs (OSEL) developed. These tools use the most innovative science to support medical device development and patient access to safe and effective medical devices. If you are considering using a tool from this catalog in your marketing submissions, note that these tools have not been qualified as Medical Device Development Tools and the FDA has not evaluated the suitability of these tools within any specific context of use. You may request feedback or meetings for medical device submissions as part of the Q-Submission Program.

For more information about the Catalog of Regulatory Science Tools, email RST_CDRH@fda.hhs.gov.

## Disclaimer
This software and documentation (the "Software") were developed at the Food and Drug Administration (FDA) by employees of the Federal Government in the course of their official duties. Pursuant to Title 17, Section 105 of the United States Code, this work is not subject to copyright protection and is in the public domain. Permission is hereby granted, free of charge, to any person obtaining a copy of the Software, to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, or sell copies of the Software or derivatives, and to permit persons to whom the Software is furnished to do so. FDA assumes no responsibility whatsoever for use by other parties of the Software, its source code, documentation or compiled executables, and makes no guarantees, expressed or implied, about its quality, reliability, or any other characteristic. Further, use of this code in no way implies endorsement by the FDA or confers any advantage in regulatory decisions. Although this software can be redistributed and/or modified freely, we ask that any derivative works bear some notice that they are derived from it, and any modified versions bear some notice that they have been modified.
