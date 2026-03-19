# hybrid-breastCT

This repository contains a Python-based framework for generating **hybrid breast CT images** by combining simulated microcalcification clusters with real patient breast CT projection data.

![image](https://github.com/user-attachments/assets/539ca4a7-c2c6-46b8-8069-080405025048)

The method preserves realistic anatomical background while ensuring physics-based image formation through projection-domain lesion insertion.

---

## Workflow Overview

The hybrid simulation pipeline performs the following steps:

### 1. Define calc cluster model
- Generates calcification clusters with configurable:
  - Number of calcs
  - Cluster diameter
  - Calc diameter
- Assigns material properties based on segmentation labels

### 2. Simulate X-ray projections of clusters
- Supports polyenergetic spectra
- Uses [TIGRE](https://github.com/CERN/TIGRE) for GPU-accelerated ray tracing

### 3. Apply detector blur
- Uses detector Modulation Transfer Function (MTF)

### 4. Load patient projection data
- Uses real breast CT projections and segmentation volumes
- Extracts geometry and viable VOI centers

### 5. Generate hybrid projections
- Combines simulated lesion projections with real projections

### 6. Reconstruct hybrid CT volume
- Uses TIGRE reconstruction algorithms:
  - FDK, SART, CGLS, MLEM
- Supports multiple reconstruction filters

### 7. Extract VOIs and MIPs
- Generates signal-present and signal-absent volumes
- Outputs MIPs for analysis or reader studies

---

## Quick Start

### 1. Clone repository
```bash
git clone https://github.com/DIDSR/hybrid-breastCT.git
cd hybrid-breastCT
