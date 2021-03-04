# PIP WEIGHT TWOPCF for OnePercent
### Files needed for the calculation of the two-point correlation function using the PIP weight scheme

This repository contains the files and scripts necessary to generate and calculate the Two-Point Correlation Function using the PIP weight scheme. The entire procedure runs at NERSC. The steps to generate our results are as follows:

#### First Step

Run `python joining_tiles.py` to read all available tiles for each catalog: "/global/cfs/cdirs/desi/datachallenge/onepercent/catalogs/{cat}/" (with cat={'bright', 'dark' or 'gray'}) to build a join summary table. These files will be needed for fiberassign runs.

#### Second Step

We now need to run the DESI fiberassign code for each catalog. To do this, we will first generate a set of N(=128 in this case) mtl-files for each catalog that differ only in the "subpriority" column using a random uniform distribution. Run `python generate_multipass.py` to generate these files and a new folder "fiberassign_pip/" will be created.

#### Third step

