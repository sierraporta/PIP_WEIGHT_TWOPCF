from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import fitsio as F
from fitsio import FITS,FITSHDR
import desimodel.io
import desitarget.mtl
import desisim.quickcat
from astropy.io import fits
from astropy.table import Table, Column, vstack
import json
import shutil
import healpy
from desitarget.targetmask import desi_mask, bgs_mask, obsconditions
from collections import Counter
import subprocess

N_real=128

def generate_data(cat):
    path="/global/cfs/cdirs/desi/datachallenge/onepercent/targets/"
    columns = ['TARGETID', 'DESI_TARGET', 'MWS_TARGET', 'BGS_TARGET', 'SUBPRIORITY', 'NUMOBS_INIT', 'PRIORITY_INIT', 'RA', 'DEC', 'HPXPIXEL', 'BRICKNAME', 'FLUX_R', 'MW_TRANSMISSION_R','OBSCONDITIONS']
    if cat=='BGS':
        filename="mtl-bright.fits"
        targ=F.read(path+filename)
        sel = (targ["DESI_TARGET"] & desi_mask["BGS_ANY"]) > 0
        targ = targ[sel]
    if cat=='LRG':
        filename="mtl-dark.fits"
        targ=F.read(path+filename)
        sel = (targ["DESI_TARGET"] & desi_mask["LRG"]) > 0
        targ = targ[sel]
    if cat=='QSO':
        filename="mtl-dark.fits"
        targ=F.read(path+filename)
        sel = (targ["DESI_TARGET"] & desi_mask["QSO"]) > 0
        targ = targ[sel]
    if cat=='ELG':
        filename="mtl-dark.fits"
        targ=F.read(path+filename)
        sel = (targ["DESI_TARGET"] & desi_mask["ELG"]) > 0
        targ = targ[sel]

    data2=[]

    outdir="fiberassign_pip/"
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    print(len(targ))
    for i in range(N_real):
        subp=np.random.uniform(0,1,len(targ))
        targ["SUBPRIORITY"]=subp
        F.write(outdir+filename.split(".")[0]+str("_")+cat+str(i)+str(".fits"),targ,clobber=True)
        data2.append(subp)
    names=["SUBP_"+str(i) for i in range(N_real)]
    F.write(outdir+"subpriorities_targets"+str("_")+cat+str(N_real)+".fits",data2,clobber=True,names=names)
    return

generate_data('BGS')
generate_data('LRG')
generate_data('QSO')
generate_data('ELG')

