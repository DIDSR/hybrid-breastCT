## RUNNING THE SOFTWARE
hybrid-simulation.py can be run in command line. 

First, activate the Python virtual environment.

Then, run the following with the command line arguments:
```
python hybrid-simulation.py --cluster_diameter_mm <CLUSTER_DIAMETER_MM> --num_calcs <NUM_CALCS> --calc_diameter_mm <CALC_DIAMETER_MM> --scanID <SCANID>
```

<CLUSTER_DIAMETER_MM> = diameter of spherical cluster 
<NUM_CALCS> = # calcs in cluster
<CALC_DIAMETER_MM> = calc diameter 
<SCANID> = scanID of patient images. Calc clusters will be inserted into these images.

## EXAMPLE

```
python hybrid-simulation.py --cluster_diameter_mm 4.0 --num_calcs 5 --calc_diameter_mm 1.0 --scanID 2878
```
