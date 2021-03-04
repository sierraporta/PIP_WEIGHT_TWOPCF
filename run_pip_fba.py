import sys, os
# import subprocess library
import subprocess
os.system('source /global/common/software/desi/desi_environment.sh master')

path="fiberassign_pip/"

N_real=128

for i in range(N_real):
    os.system('nohup fba_run --targets fiberassign_pip/mtl-bright_BGS{}.fits --sky /global/cfs/cdirs/desi/datachallenge/onepercent/targets/skies_fixedid.fits --footprint alle2etiles_bright.fits --rundate 2020-01-01T00:00:00 --dir fiberassign_pip/bgs/fba_bright_targets_BGS{} --overwrite'.format(i,i))
os.system('rm fiberassign_pip/mtl-bright_BGS*')

for i in range(N_real):
    os.system('nohup fba_run --targets fiberassign_pip/mtl-dark_LRG{}.fits --sky /global/cfs/cdirs/desi/datachallenge/onepercent/targets/skies_fixedid.fits --footprint alle2etiles_dark.fits --rundate 2020-01-01T00:00:00 --dir fiberassign_pip/lrg/fba_dark_targets_LRG{} --overwrite'.format(i,i))
os.system('rm fiberassign_pip/mtl-dark_LRG*')

for i in range(N_real):
    os.system('nohup fba_run --targets fiberassign_pip/mtl-dark_QSO{}.fits --sky /global/cfs/cdirs/desi/datachallenge/onepercent/targets/skies_fixedid.fits --footprint alle2etiles_dark.fits --rundate 2020-01-01T00:00:00 --dir fiberassign_pip/qso/fba_dark_targets_QSO{} --overwrite'.format(i,i))
os.system('rm fiberassign_pip/mtl-dark_QSO*')

for i in range(N_real):
    os.system('nohup fba_run --targets fiberassign_pip/mtl-dark_ELG{}.fits --sky /global/cfs/cdirs/desi/datachallenge/onepercent/targets/skies_fixedid.fits --footprint alle2etiles_gray.fits --rundate 2020-01-01T00:00:00 --dir fiberassign_pip/elg/fba_dark_targets_ELG{} --overwrite'.format(i,i))
os.system('rm fiberassign_pip/mtl-dark_ELG*')


