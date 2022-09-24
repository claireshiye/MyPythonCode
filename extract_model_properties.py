import numpy as np
import pandas as pd
import re
import os,sys
import subprocess
import gzip
from glob import glob
import matplotlib.pyplot as plt
import random
import scipy.optimize as opt
import unit_convert as uc
import bh_wd_hdf5 as bw_hdf
import ns_tidalcapture_hdf5 as ntc_hdf
import ns_hdf5 as ns_hdf
import dynamics as dyn
import ns

sys.path.insert(1, '/projects/b1095/syr904/MyCodes/cmctoolkit')
import cmctoolkit as cmct
sys.path.append('/projects/b1095/syr904/MyCodes/PythonCodes/makeSBP')
import make_SBP_and_vel_profile_one_hdf5 as mone_hdf


modelpath = '/projects/b1095/syr904/cmc/CMC-COSMIC/wdmerger_update/CMC-COSMIC_modified/CMC/runs/n8-rv0.5-rg8-z0.0002_tcoff/'
modelN = 800000
modelZ = 0.0002
modelrv = 0.5
modelrg = 8.
modeldist  = 4.0
lowtime = 9000.; hightime = 13100.
deltas = -2

def read_keys(thekey):
    return re.findall(r'\d+\.\d+|\d+', thekey)


def make_snap_key_file(sourcepath):
    snap_h5 = 'initial.snapshots.h5'
    with pd.HDFStore(sourcepath+snap_h5) as snap_hdf:
        snap_keys = np.sort(snap_hdf.keys())

    #print(snap_keys)    
    #max_snapno = len(snap_keys)    
    
    #t_conv = dyn.conv('t', sourcepath+'initial.conv.sh')


    snapno = []; snaptime = []
    for ii in range(len(snap_keys)):
        theno = read_keys(snap_keys[ii])[0]; thetime = read_keys(snap_keys[ii])[1]
        snapno.append(int(theno)); snaptime.append(thetime)
    
    snapno_sort, snaptime_sort = (np.array(t) for t in zip(*sorted(zip(snapno, snaptime))))
    #print(snapno_sort)
    np.savetxt(sourcepath+'snap_keys.txt', np.c_[snapno_sort, snaptime_sort], fmt = '%s %s', header = '1.snap_no 2.snap_time', comments = '#')

#make_snap_key_file(modelpath)
#print('make snap_keys.txt')

#os.system('rm '+modelpath+'initial.ns.dat')
#ns.print_Nns_psrfile([modelpath], 0, 1, 0)

mone_hdf.main(modelpath, modelN, modelZ, modelrv, modelrg, modeldist, lowtime, hightime, deltas)
print('make projection files, surface brightness, velocity dispersions, and parameter files')

#datasnapkey = np.genfromtxt(modelpath+'snap_keys.txt')
#t_conv = dyn.conv('t', modelpath+'initial.conv.sh')
#snapnos = datasnapkey[:,0]; snaptimes = datasnapkey[:,1]*t_conv
#for xx in range(len(snaptimes)-1, 0, deltas):
#    if snaptimes[xx] > hightime: 
#        continue
#    if snaptimes[xx] < lowtime:
#        print('stop!')
#        break
#    
#    ns_hdf.get_allpsr_atsnap(modelpath, 'MSP', int(snapnos[xx]), modeldist, modelZ)
#    ns_hdf.get_allpsr_atsnap(modelpath, 'PSR', int(snapnos[xx]), modeldist, modelZ)
#     
#    ntc_hdf.find_NS_XX_atsnap(modelpath, 0, 1, int(snapnos[xx]), 'NSMS', modeldist, modelZ)
#    ntc_hdf.find_NS_XX_atsnap(modelpath, 10, 12, int(snapnos[xx]), 'NSWD', modeldist, modelZ)
#    ntc_hdf.find_NS_XX_atsnap(modelpath, 2, 9, int(snapnos[xx]), 'NSG', modeldist, modelZ)
#
#    bw_hdf.find_BH_XX_atsnap_hdf5(modelpath, 0, 1, int(snapnos[xx]), 'BHMS', modeldist, modelZ)
#    bw_hdf.find_BH_XX_atsnap_hdf5(modelpath, 10, 12, int(snapnos[xx]), 'BHWD', modeldist, modelZ)
#    bw_hdf.find_BH_XX_atsnap_hdf5(modelpath, 2, 9, int(snapnos[xx]), 'BHG', modeldist, modelZ)
#
#    bw_hdf.find_WD_XX_atsnap_hdf5(modelpath, 0, 1, int(snapnos[xx]), 'WDMS', modeldist, modelZ)
