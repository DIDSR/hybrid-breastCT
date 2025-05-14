### Setup
1. Install [TIGRE](https://github.com/CERN/TIGRE/blob/master/Frontispiece/python_installation.md) for Python
2. Setup and activate Python environment:
   
```
conda create --name bct_env python=3.11
conda activate bct_env
pip install -r requirements.txt
```
3. Adapt ```hybrid-simulation.py``` and ```utils/utils.py``` for your breast CT system. More info on homepage. 

### Usage
```hybrid-simulation.py``` can be run in command line. First, activate Python virtual environment. Then, run the following:
```
python hybrid-simulation.py --cluster_diameter_mm <CLUSTER_DIAMETER_MM> --num_calcs <NUM_CALCS> --calc_diameter_mm <CALC_DIAMETER_MM> --scanID <SCAN_ID>
```

where

```<CLUSTER_DIAMETER_MM>``` = diameter of spherical cluster in mm

```<NUM_CALCS>``` = # calcs in cluster

```<CALC_DIAMETER_MM>``` = calc diameter in mm

```<SCAN_ID>``` = scanID of patient images. Calc clusters will be inserted into these images.

### Example

```
python hybrid-simulation.py --cluster_diameter_mm 4.0 --num_calcs 5 --calc_diameter_mm 1.0 --scanID 2878
```
