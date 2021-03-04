# PIP WEIGHT TWOPCF for OnePercent
### Files needed for the calculation of the two-point correlation function using the PIP weight scheme

This repository contains the files and scripts necessary to generate and calculate the Two-Point Correlation Function using the PIP weight scheme. The entire procedure runs at NERSC. The steps to generate our results are as follows:

#### First Step

Run python `python joining_tiles.py` to read all available tiles for each catalog: "/global/cfs/cdirs/desi/datachallenge/onepercent/catalogs/{cat}/" to build a join summary table. These files will be needed for fiberassign runs.
