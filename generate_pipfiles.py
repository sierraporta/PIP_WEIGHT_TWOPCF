#!/usr/bin/env python
# coding: utf-8

import fitsio as F
import numpy as np
import h5py
import matplotlib.pyplot as plt
from astropy.table import Table, join, unique
from desitarget.targetmask import obsconditions, bgs_mask, desi_mask
import random
from random import randint
from random import seed
random.seed(30)
import glob, os, sys
import fitsio
from random import random
from tqdm import tqdm 

def select_cat_dat(cat):
    path="/global/cfs/cdirs/desi/datachallenge/onepercent/targets/"
    if cat=='BGS':
        targ=F.read(path+"mtl-bright.fits",columns=["TARGETID","RA","DEC","FLUX_G","FLUX_R","FLUX_Z","NOBS_G","NOBS_R","NOBS_Z","DESI_TARGET","BGS_TARGET","MWS_TARGET"])
        truth=F.read(path+"truth-bright.fits",ext=1,columns=["TARGETID","TRUEZ","TRUESPECTYPE","TEMPLATETYPE"])
        truth=truth[truth["TEMPLATETYPE"]=='BGS       ']
        sel = (targ["DESI_TARGET"] & desi_mask["BGS_ANY"]) > 0
        targ = targ[sel]
        #targ=targ[targ["DESI_TARGET"]==desi_mask.mask("BGS_ANY")]
    if cat=='LRG':
        targ=F.read(path+"mtl-dark.fits",columns=["TARGETID","RA","DEC","FLUX_G","FLUX_R","FLUX_Z","NOBS_G","NOBS_R","NOBS_Z","DESI_TARGET","BGS_TARGET","MWS_TARGET"])
        truth=F.read(path+"truth-dark.fits",ext=1,columns=["TARGETID","TRUEZ","TRUESPECTYPE","TEMPLATETYPE"])
        truth=truth[truth["TEMPLATETYPE"]=='LRG       ']
        sel = (targ["DESI_TARGET"] & desi_mask["LRG"]) > 0
        targ = targ[sel]
        #targ=targ[targ["DESI_TARGET"]==desi_mask.mask("LRG")]
    if cat=='QSO':
        targ=F.read(path+"mtl-dark.fits",columns=["TARGETID","RA","DEC","FLUX_G","FLUX_R","FLUX_Z","NOBS_G","NOBS_R","NOBS_Z","DESI_TARGET","BGS_TARGET","MWS_TARGET"])
        truth=F.read(path+"truth-dark.fits",ext=1,columns=["TARGETID","TRUEZ","TRUESPECTYPE","TEMPLATETYPE"])
        truth=truth[truth["TEMPLATETYPE"]=='QSO       ']
        sel = (targ["DESI_TARGET"] & desi_mask["QSO"]) > 0
        targ = targ[sel]
        #targ=targ[targ["DESI_TARGET"]==desi_mask.mask("QSO")]
    if cat=='ELG':
        targ=F.read(path+"mtl-dark.fits",columns=["TARGETID","RA","DEC","FLUX_G","FLUX_R","FLUX_Z","NOBS_G","NOBS_R","NOBS_Z","DESI_TARGET","BGS_TARGET","MWS_TARGET"])
        truth=F.read(path+"truth-dark.fits",ext=1,columns=["TARGETID","TRUEZ","TRUESPECTYPE","TEMPLATETYPE"])
        truth=truth[truth["TEMPLATETYPE"]=='ELG       ']
        sel = (targ["DESI_TARGET"] & desi_mask["ELG"]) > 0
        targ = targ[sel]
        #targ=targ[targ["DESI_TARGET"]==desi_mask.mask("ELG")]
    targ=Table(targ)
    truth=Table(truth)
    all_tar=join(targ,truth)
    # Keep cols: ["RA","DEC","TARGETID","TRUEZ"]
    #all_tar=all_tar["RA","DEC","TARGETID","TRUEZ"]
    all_tar.rename_column('TRUEZ', 'Z')
    rmag=np.where(all_tar["FLUX_R"]>0,22.-2.5*np.log10(all_tar["FLUX_R"]),0)
    all_tar["RMAG"]=rmag
    return all_tar

def select_cat_ran(cat):
    path="/global/cfs/cdirs/desi/datachallenge/onepercent/catalogs/"
    columns=['TARGETID','RA','DEC','BRICKNAME','NOBS_G','NOBS_R','NOBS_Z','HPXPIXEL','DESI_TARGET','NUMOBS_INIT','PRIORITY','OBSCONDITIONS']
    if cat=='BGS':
        randoms=Table(F.read(path+"bright/randoms_mtl_cuttod.fits",columns=columns))
    if cat=='LRG':
        randoms=Table(F.read(path+"dark/randoms_mtl_cuttod.fits",columns=columns))
    if cat=='QSO':
        randoms=Table(F.read(path+"dark/randoms_mtl_cuttod.fits",columns=columns))    
    if cat=='ELG':
        randoms=Table(F.read(path+"gray/randoms_mtl_cuttod.fits",columns=columns))
    return randoms

def get_bitweights(array):
    """
    Creates an array of bitwise weights stored as 64-bit signed integers
    
    Input: a 2D boolean array of shape (Ngal, Nreal), where Ngal is the total number 
           of target galaxies, and Nreal is the number of fibre assignment realizations.
           
    Output: returns a 2D array of 64-bit signed integers. 
    """

    Nbits=64
    dtype=np.int64

    Ngal, Nreal = array.shape           # total number of realizations and number of target galaxies
    Nout = (Nreal + Nbits - 1) // Nbits # number of output columns

    # intermediate arrays
    bitw8 = np.zeros((Ngal, 8), dtype="i")   # array of individual bits of 8 realizations
    bitweights = np.zeros(Ngal, dtype=dtype) # array of 64-bit integers

    # array to store final output
    output_array = np.zeros((Ngal, Nout), dtype=dtype)

    idx_out = 0 # initial column in output_array

    # loop through realizations to build bitwise weights
    for i in range(Nreal):
        bitw8[array[:,i], i%8] = 1
        arr = np.array(np.packbits(bitw8[:,::-1]), dtype=dtype)
        bitweights = np.bitwise_or(bitweights, np.left_shift(arr, 8*((i%Nbits)//8)))

        if (i+1)%Nbits == 0 or i+1 == Nreal:
            output_array[:,idx_out] = bitweights
            bitweights[:] = 0
            idx_out += 1

        if (i+1)%8 == 0:
            bitw8[:] = 0

    return output_array

def generate_pip_files(cat,realizations=2):
    data=select_cat_dat(cat)
    randoms=select_cat_ran(cat)
    if cat=='BGS':
        path='fiberassign_pip/bgs/fba_bright_targets_BGS'
        #pathtofba_ran="fiberassign_cat/fba_bright_randoms/fba*"
        #fba_files_ran=glob.glob(pathtofba_ran)
    if cat=='LRG':
        path='fiberassign_pip/lrg/fba_dark_targets_LRG'
        #pathtofba_ran="fiberassign_cat/fba_dark_randoms/fba*"
        #fba_files_ran=glob.glob(pathtofba_ran)
    if cat=='QSO':
        path='fiberassign_pip/qso/fba_dark_targets_QSO'
        #pathtofba_ran="fiberassign_cat/fba_dark_randoms/fba*"
        #fba_files_ran=glob.glob(pathtofba_ran)
    if cat=='ELG':
        path='fiberassign_pip/elg/fba_dark_targets_ELG'
        #pathtofba_ran="fiberassign_cat/fba_gray_randoms/fba*"
        #fba_files_ran=glob.glob(pathtofba_ran)
    for i in range(realizations):
        fba_files = glob.glob(path+str(i)+"/fba*.fits")
        allavail_tar=np.concatenate([F.read(f,ext="FAVAIL") for f in fba_files])
        allassign_tar=np.concatenate([F.read(f,ext="FASSIGN",columns=["FIBER","TARGETID","LOCATION","FIBERSTATUS"]) for f in fba_files])
        isas=np.zeros(len(data),dtype=bool)
        idas=np.isin(data["TARGETID"], allassign_tar["TARGETID"])
        isas[idas] = True
        data["ISASS_"+str(i)]=isas
    isav=np.zeros(len(data),dtype=bool)
    idav=np.isin(data["TARGETID"], allavail_tar["TARGETID"])
    isav[idav] = True
    data["ISAVAIL"]=isav
    
    parent=data[data["ISAVAIL"]==True]
    targeted=parent[parent["ISASS_0"]==True]
    col=["ISASS_"+str(i) for i in range(realizations)]
    targeted["wi"]=np.sum(targeted["ISASS_"+str(i)] for i in range(realizations))
    parent["wi"]=np.sum(parent["ISASS_"+str(i)] for i in range(realizations))
    
    vec1,vec2=[],[]
    for i in targeted[col]:
    #    print(np.packbits(list(i)).view(np.int)[0], np.packbits(list(i)).view(np.int)[1])
        vec1.append(np.packbits(list(i)).view(np.int)[0])
        vec2.append(np.packbits(list(i)).view(np.int)[1])
    vec1,vec2=np.array(vec1),np.array(vec2)

    targeted["BITWEIGHT0"]=vec1
    targeted["BITWEIGHT1"]=vec2
    parent["BITWEIGHT0"]=np.ones(len(parent),dtype=int)*-1
    parent["BITWEIGHT1"]=np.ones(len(parent),dtype=int)*-1

    nd=len(parent)
    randoms["Z"]=np.zeros(len(randoms))

    for i in range(0,len(randoms)):
        ind=int(random()*nd)
        randoms["Z"][i]=parent["Z"][ind]

    randoms["BITWEIGHT0"]=np.ones(len(randoms),dtype=int)*-1
    randoms["BITWEIGHT1"]=np.ones(len(randoms),dtype=int)*-1

    cols1=["RA","DEC","Z","BITWEIGHT0","BITWEIGHT1"]
    cols2=["RA","DEC","Z","BITWEIGHT0","BITWEIGHT1","RMAG"]
    parent0=parent[cols1]
    targeted0=targeted[cols2]
    randoms0=randoms[cols1]
    #targeted1=targeted0[targeted0["BITWEIGHT0"]%2==1]
    return parent0, targeted0, randoms0

def writefiles(parent,targeted,randoms,cat,rmag=25):
    rmag_cut=targeted["RMAG"]<=rmag
    targeted=targeted[rmag_cut]
    outdir="clus_files_pip/"
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    with h5py.File(outdir+'parent_'+cat+'.hdf5', 'w') as f:
        f.create_dataset("RA", data=parent["RA"])
        f.create_dataset("DEC", data=parent["DEC"])
        f.create_dataset("Z", data=parent["Z"])
        f.create_dataset("BITWEIGHT0", data=parent["BITWEIGHT0"])
        f.create_dataset("BITWEIGHT1", data=parent["BITWEIGHT1"])
    with h5py.File(outdir+'randoms_'+cat+'.hdf5', 'w') as f:
        f.create_dataset("RA", data=randoms["RA"])
        f.create_dataset("DEC", data=randoms["DEC"])
        f.create_dataset("Z", data=randoms["Z"])
        f.create_dataset("BITWEIGHT0", data=randoms["BITWEIGHT0"])
        f.create_dataset("BITWEIGHT1", data=randoms["BITWEIGHT1"])
    with h5py.File(outdir+'targeted_'+cat+'.hdf5', 'w') as f:
        f.create_dataset("RA", data=targeted["RA"])
        f.create_dataset("DEC", data=targeted["DEC"])
        f.create_dataset("Z", data=targeted["Z"])
        f.create_dataset("BITWEIGHT0", data=targeted["BITWEIGHT0"])
        f.create_dataset("BITWEIGHT1", data=targeted["BITWEIGHT1"])

bgs=generate_pip_files('BGS',128)
parent_bgs=bgs[0]
targeted_bgs=bgs[1]
randoms_bgs=bgs[2]
writefiles(parent_bgs,targeted_bgs,randoms_bgs,cat='bgs')
print(len(parent_bgs),len(targeted_bgs),len(randoms_bgs))

lrg=generate_pip_files('LRG',128)
parent_lrg=lrg[0]
targeted_lrg=lrg[1]
randoms_lrg=lrg[2]
writefiles(parent_lrg,targeted_lrg,randoms_lrg,cat='lrg')
print(len(parent_lrg),len(targeted_lrg),len(randoms_lrg))

qso=generate_pip_files('QSO',128)
parent_qso=qso[0]
targeted_qso=qso[1]
randoms_qso=qso[2]
writefiles(parent_qso,targeted_qso,randoms_qso,cat='qso')
print(len(parent_qso),len(targeted_qso),len(randoms_qso))

elg=generate_pip_files('ELG',128)
parent_elg=elg[0]
targeted_elg=elg[1]
randoms_elg=elg[2]
writefiles(parent_elg,targeted_elg,randoms_elg,cat='elg')
print(len(parent_elg),len(targeted_elg),len(randoms_elg))

print("DOne..")
