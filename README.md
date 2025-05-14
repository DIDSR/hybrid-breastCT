# hybrid-breastCT
This repository contains code for generating hybrid breast CT images combining simulated pathology with real anatomical background.
![image](https://github.com/user-attachments/assets/539ca4a7-c2c6-46b8-8069-080405025048)

Each of these steps is embedded in ```hybrid-simulation.py```, found in the ```code``` directory.
#### Define calc cluster model:
- Generates calcification clusters with flexible configuration of:
  - Number of calcs, Cluster diameter, Calc diameter
- Identifies background tissue type at insertion site based on breast segmentation volume. Assigns values according to density and linear attenuation coefficients of calcification and background tissue.
  

#### Simulate x-ray projections of clusters:
- Supports polyenergetic X-ray spectra with energy binning.
- Generates projections of inserted calcs using [TIGRE reconstruction toolbox](https://github.com/CERN/TIGRE/tree/master) GPU-accelerated ray tracing.


#### Blurred x-ray projections of clusters: 
- Applies detector blur using Modulation Transfer Function (MTF).


#### Patient x-ray projections (anatomical background):
- Loads real segmentation and patient breast CT projection data.
- Extracts scan parameters, system geometry, and viable VOI centers for calc insertion.


#### Hybrid projection images:
- Combines simulated calc projections with real patient projections to create hybrid projection images.


#### Hybrid breast CT volume:
- Utilizes [TIGRE](https://github.com/CERN/TIGRE/tree/master) reconstruction algorithms to generate hybrid CT volume.
  - Multiple algorithms supported: FDK, SART, CGLS, and MLEM.
  - Apodization filters can be modified for FDK algorithm (e.g., Ram-Lak, Hann).


#### Hybrid breast CT ROI/VOI:
- Extracts signal-present and signal-absent VOIs and maximum intensity projections (MIPs) for training, testing, or reader studies.

## Usage
Instructions for setup, installation, and customizing code can be found in ```code```


## Citation

"_Hybrid simulation of breast CT for assessing microcalcification detectability_". Lyu SH, Makeev A, Li D, Badal A, Hernandez AM, Boone JM, Glick SJ.
Journal of Medical Imaging, 2025.

## More information

[Catalog of Regulatory Science Tools to Help Assess New Medical Devices](https://www.fda.gov/medical-devices/science-and-research-medical-devices/catalog-regulatory-science-tools-help-assess-new-medical-devices)


#### About the Catalog of Regulatory Science Tools

The enclosed tool is part of the Catalog of Regulatory Science Tools, which provides a peer-reviewed resource for stakeholders to use where standards and qualified Medical Device Development Tools (MDDTs) do not yet exist. These tools do not replace FDA-recognized standards or MDDTs. This catalog collates a variety of regulatory science tools that the FDA's Center for Devices and Radiological Health's (CDRH) Office of Science and Engineering Labs (OSEL) developed. These tools use the most innovative science to support medical device development and patient access to safe and effective medical devices. If you are considering using a tool from this catalog in your marketing submissions, note that these tools have not been qualified as Medical Device Development Tools and the FDA has not evaluated the suitability of these tools within any specific context of use. You may request feedback or meetings for medical device submissions as part of the Q-Submission Program.

For more information about the Catalog of Regulatory Science Tools, email OSEL_CDRH@fda.hhs.gov. 

## Disclaimer
This software and documentation (the "Software") were developed at the Food and Drug Administration (FDA) by employees of the Federal Government in the course of their official duties. Pursuant to Title 17, Section 105 of the United States Code, this work is not subject to copyright protection and is in the public domain. Permission is hereby granted, free of charge, to any person obtaining a copy of the Software, to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, or sell copies of the Software or derivatives, and to permit persons to whom the Software is furnished to do so. FDA assumes no responsibility whatsoever for use by other parties of the Software, its source code, documentation or compiled executables, and makes no guarantees, expressed or implied, about its quality, reliability, or any other characteristic. Further, use of this code in no way implies endorsement by the FDA or confers any advantage in regulatory decisions. Although this software can be redistributed and/or modified freely, we ask that any derivative works bear some notice that they are derived from it, and any modified versions bear some notice that they have been modified.
