# PIP WEIGHT TWOPCF for OnePercent
### Files needed for the calculation of the two-point correlation function using the PIP weight scheme

This repository contains the files and scripts necessary to generate and calculate the Two-Point Correlation Function using the PIP weight scheme. The entire procedure runs at NERSC. The steps to generate our results are as follows:

#### Step 1

Run `python joining_tiles.py` to read all available tiles for each catalog: "/global/cfs/cdirs/desi/datachallenge/onepercent/catalogs/{cat}/" (with cat={'bright', 'dark' or 'gray'}) to build a join summary table. These files will be needed for fiberassign runs.

#### Step 2

We now need to run the DESI fiberassign code for each catalog. To do this, we will first generate a set of N(=128 in this case) mtl-files for each catalog that differ only in the "subpriority" column using a random uniform distribution. Run `python generate_multipass.py` to generate these files and a new folder "fiberassign_pip/" will be created.

#### Step 3

Run `python run_pip_fba.py` to read all mtl-files for each catalog and run fiberassign code for each sample. A folder "fiberassign_pip/{cat}/fba-{time}-{cat}" (with time={'bright', 'dark' or 'gray'} and cat={'BGS', 'LRG' and 'ELG'}) will be created for each sample and then all the mtl-files will be deleted to free up disk space.

#### Step 4

Run `python generate_pipfiles.py` to generate the input files needed to calculate the two-point correlation function using the pip weight scheme. For each catalog the three files will be created in a new directory "clus_files_pip/parent_{cat}.hdf5", "clus_files_pip/randoms_{cat}.hdf5" and "clus_files_pip/targeted_{cat}.hdf5". These files contain RA, DEC, Z and a PIP weight vector (BITWEIGHT0 and BITWEIGHT1). The PIP weights are stored as a 64-bit signed integer (for 128=2^7 fibre assignment realizations). It's possible to have longer weight vectors by adding more datasets (BITWEIGHT2, BITWEIGHT3, etc). In the parent and random files, the PIP weights are all set to 1 (which corresponds to the integer -1), and the file of targeted objects contains the actual PIP weight vectors from fibre assignment.

##### TWOPCF code

There is a code that can calculate the correlation function with the PIP corrections on github here
https://github.com/lstothert/two_pcf

This code runs on NERSC. I have put a `makefile` file that has been modified so that the executable can be built in the NERSC system.

It should be possible to get it to compile using `make TWOPCF`.

Then to run the code with a parameter file, do `./TWOPCF param.txt`.

In the makefile, you can turn on/off the PIP weights with the option D_USE_INV_WEIGHT. Also, the option D_N_MASK_INTS sets the number of columns of bitwise weights in the hdf5 files (so in this example, with 1 64-bit integer, it is set to 1, in this case it is set to 2 corresponding with 2 64-bit integer).

To apply the PIP weights, you first need to calculate the angular pair counts for the parent sample. If you compile (with the D_USE_INV_WEIGHT option turned on), then run TWOPCF with param_parent.txt, it will create the files angDD.hdf5 and angDR.hdf5 with these angular pair counts. It will also output the monopole for the parent sample.

Next, run TWOPCF with param_pip.txt. It will read in the angular pair counts, and also calculate the angular pair counts that use the PIP weights. It will then use these to output the monopole of the targeted sample with the PIP correction applied.

To compare this to the monopole without the PIP correction, you can re-compile (without D_USE_INV_WEIGHT), then run TWOPCF using param_targeted.txt. The python script `plot_xi.ipynb` will make a plot to check that the PIP correction is working for each sample.

For this trge samples (BGS, LRG and ELG) we obtain the folloging plots:
![Texto alternativo](correlfunctPIPBGS.pdf)
