import numpy as np
import fitsio as F
from astropy.table import Table
import glob, sys, os

def joinning_tiles(cat='bright'):
    # 'bright': BGS tiles, 'dark': LRG tiles, 'gray': ELG tiles
    path_to_tiles_e2e="/global/cfs/cdirs/desi/datachallenge/onepercent/catalogs/{}/".format(cat)
    tiles_files=glob.glob(path_to_tiles_e2e+"e2etiles_run*")
    alltiles=Table(np.concatenate([F.read(f) for f in tiles_files]))
    print("We have {} tiles for {} time".format(len(tiles_files),cat))
    print("writing tiles for {} time".format(cat))
    alltiles.write("alle2etiles_{}.fits".format(cat),overwrite=True)
    print("DOne..!")
    return

print("Writing tiles for catalogs...")
joinning_tiles('bright')
joinning_tiles('dark')
joinning_tiles('gray')
print("DOne..!")
