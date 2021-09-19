import numpy as np
import pandas as pd
import matplotlib
print(matplotlib.__version__)
matplotlib.use('PDF')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import matplotlib.lines as mlines
from glob import glob
from collections import Counter
import ns
import history_cmc as hic
import math
import scipy
from scipy.interpolate import interp1d
from scipy import stats
import matplotlib.cm as cm
import matplotlib as mpl
import random
from random import shuffle
import gzip
import sys
import astropy
from astropy import units

import ecc_calc as gwcalc
import unit_convert as uc
import merger_rate_calculator as mr
import ns_tidalcapture as tc
import conversions
import dynamics as dyn

#%matplotlib inline

matplotlib.rcParams.update({'font.size': 24})
      
twopi=2.*np.pi
yearsc=3.1557*10**7
Kconst=9.87*10**-48 ##yr/G^2
Gconst=6.674*10**-8 ##cm3*g-1*s-2
Gconst_sun = 4.30091*10**-3 ##pc*M_sun**-1*(km/s)2
clight=3*10**10 ##cm/s
Msun=2*10**33 ##gram
AU=1.496*10**13  ##cm
AU_Rsun=214.93946938362 ##AU to R_sun
PC=3.086*10**18  ##cm
pC_Rsun = 44334448.0068964 ##pc to R_sun


##Finding eclipsing binaries in the cluster
#path = '/projects/b1095/syr904/cmc/47Tuc/rundir/47Tuc/best_fits/MOCHA47Tuc_elson_rv4_3e6_tcon/'
path = '/projects/b1095/syr904/cmc/47Tuc/rundir/47Tuc/best_fits/MOCHA47Tuc_elson_rv4_3e6_tcon_fixbin_ic/'
sys.path.insert(1, '/projects/b1095/syr904/MyCodes/cmctoolkit_oldsnap')

#path = '/projects/b1095/syr904/cmc/CMC-COSMIC/master_tc_test/ver_0601/MOCHA47Tuc_elson_rv4_3e6_tcon_fb10/'
#sys.path.insert(1, '/projects/b1095/syr904/MyCodes/cmctoolkit')

import cmctoolkit as cmct

t_obs = 8.29; t_det = 2*t_obs  ##in days
def eclipsing_binary(a_au, ecc, r0, r1, m0, m1, k0, k1):
    ecl_flag = []
    period = uc.au_to_period(a_au, m0, m1)
    
    theta = np.random.uniform(low=0., high=2.*np.pi, size=len(r0))
            
    a_rsun = a_au*AU_Rsun*(1-ecc)  ##use the minimum separation
    
    for xx in range(len(period)):
        if theta[xx] < (1./2.)*np.pi:
            beta = (1./2.)*np.pi - theta[xx]
        elif (1./2.)*np.pi <= theta[xx] < np.pi:
            beta = theta[xx] - np.pi/2.
        elif np.pi <= theta[xx] < (3./2.)*np.pi:
            beta = (3./2.)*np.pi - theta[xx]
        else:
            beta = theta[xx] - (3./2.)*np.pi
        
        if period[xx] > t_det or k0[xx] >=10 or k1[xx]>=10:
            ecl_flag.append(0)
        
        else:
            if r0[xx] > r1[xx]:
                l0 = a_rsun[xx] - r1[xx]/np.sin(beta)
                r_cut = l0 * np.sin(beta)
                
                if r_cut <= 0:
                    r_cov = 0
                elif r_cut >= r0[xx]:
                    r_cov = 0
                elif r0[xx]-r_cut > 2*r1[xx]:
                    r_cov = 2*r1[xx]
                else:
                    r_cov = r0[xx]-r_cut

                cov_ratio = r_cov/r0[xx]
                    
            else:
                l1 = a_rsun[xx] - r0[xx]/np.sin(beta)
                r_cut = l1 * np.sin(beta)
                
                if r_cut <= 0:
                    r_cov = 0
                elif r_cut >= r1[xx]:
                    r_cov = 0
                elif r1[xx]-r_cut > 2*r0[xx]:
                    r_cov = 2*r0[xx]
                else:
                    r_cov = r1[xx]-r_cut

                cov_ratio = r_cov/r1[xx]
            
            if cov_ratio >= 0.1:
                ecl_flag.append(1)
            else:
                ecl_flag.append(0)
                
    #print(period)
                
    return period, ecl_flag
    
    
def dect_prob(r0, r1, m0, m1, a_au):
    zeta = 0.9
    period = uc.au_to_period(a_au, m0, m1)
    
    mtot = m0+m1
    
    dp = []
    for xx in range(len(r0)):
        if r0[xx] > r1[xx]:
            dp.append((zeta*r0[xx]+r1[xx])/(4.21*mtot[xx]**(1/3)*period[xx]**(2/3))*max(0., min(1., 2*t_obs/period[xx]-1.)))
        else:
            dp.append((zeta*r1[xx]+r0[xx])/(4.21*mtot[xx]**(1/3)*period[xx]**(2/3))*max(0., min(1., 2*t_obs/period[xx]-1.)))

    return dp
    

#id0_bin = []; id1_bin = []; k0_bin = []; k1_bin = []
#m0_bin = []; m1_bin = []; r0_bin = []; r1_bin = []
#sma_bin = []; ecc_bin = []
#
#with gzip.open(path+'initial.snap0507.dat.gz', 'r') as fsnap:
#    next(fsnap); next(fsnap)
#    for line in fsnap:
#        data = line.split()
#        if int(data[7])==1:
#            id0_bin.append(int(data[10])); id1_bin.append(int(data[11]))
#            k0_bin.append(int(data[17])); k1_bin.append(int(data[18]))
#            m0_bin.append(float(data[8])); m1_bin.append(float(data[9]))
#            r0_bin.append(float(data[21])); r1_bin.append(float(data[22]))
#            sma_bin.append(float(data[12])); ecc_bin.append(float(data[13]))
#            

###################CMC-COSMIC################################
#snap = cmct.Snapshot(fname=path+'initial.snapshots.h5', snapshot_name='/446(t=0.58030745)', 
#                         conv=path+'initial.conv.sh', 
#                         dist=4.52, # distance to cluster in kpc
#                         z=0.0038)


#binflag = snap.data['binflag']
#m0 = snap.data['m0_MSUN']; m1 = snap.data['m1_MSUN']
#r0 = snap.data['bin_star_radius0_RSUN']; r1 = snap.data['bin_star_radius1_RSUN']
#k0 = snap.data['bin_startype0']; k1 = snap.data['bin_startype1']
#sma = snap.data['a_AU']; ecc = snap.data['e']
#id0 = snap.data['id0']; id1 = snap.data['id1']
#id_star = snap.data['id']

##################cmc-3.2####################################
snap = cmct.Snapshot(fname=path+'initial.snap0651.dat.gz',
                         conv=path+'initial.conv.sh', 
                         dist=4.52, # distance to cluster in kpc
                         z=0.0038)

binflag = snap.data['binflag']
m0 = snap.data['m0[MSUN]']; m1 = snap.data['m1[MSUN]']
r0 = snap.data['bin_star_radius0[RSUN]']; r1 = snap.data['bin_star_radius1[RSUN]']
k0 = snap.data['bin_startype0']; k1 = snap.data['bin_startype1']
sma = snap.data['a[AU]']; ecc = snap.data['e']
id0 = snap.data['id0']; id1 = snap.data['id1']
id_star = snap.data['id']
print('read data')

snap.add_photometry('/projects/b1095/syr904/MyCodes/cmctoolkit_oldsnap/filt_index.txt')
print('add photometry')

tot_obsmag = list(snap.data['tot_obsMag_V'])

############################################################
sin_obsmag = list(snap.data['obsMag_V'][binflag != 1])
bin_obsmag0 = list(snap.data['bin_obsMag0_V'][binflag == 1])
bin_obsmag1 = list(snap.data['bin_obsMag1_V'][binflag == 1])
#tot_obsmag = list(snap.data['tot_obsMag_V'])

id0_bin = list(id0[binflag == 1]); id1_bin = list(id1[binflag == 1])
sma_bin = list(sma[binflag == 1]); ecc_bin = list(ecc[binflag == 1])
m0_bin = list(m0[binflag == 1]); m1_bin = list(m1[binflag == 1])
k0_bin = list(k0[binflag == 1]); k1_bin = list(k1[binflag == 1])
r0_bin = list(r0[binflag == 1]); r1_bin = list(r1[binflag == 1])
id_sin = list(id_star[binflag != 1])

#projfile = np.sort(glob(path+'SNAP446_*.2Dproj.dat.gz'))
projfile = np.sort(glob(path+'SNAP651_*.2Dproj.dat.gz'))
#print(projfile)

period_bin, eclipsing = eclipsing_binary(np.array(sma_bin), np.array(ecc_bin), np.array(r0_bin), np.array(r1_bin), np.array(m0_bin), np.array(m1_bin), np.array(k0_bin), np.array(k1_bin))
detection_prob = dect_prob(np.array(r0_bin), np.array(r1_bin), np.array(m0_bin), np.array(m1_bin),
                              np.array(sma_bin))
print(np.sum(eclipsing), len(eclipsing))
    
#eclips_bin = [[],[],[],[],[],[]]
n_eclips = []; n_all = []
n_all_brt = []; n_bin_brt = []; n_eclips_brt = []
for ii in range(len(projfile)):
    print(ii)
    m0_eclips = []; m1_eclips = []
    n_eclips.append(0); n_all.append(0)
    n_all_brt.append(0); n_bin_brt.append(0); n_eclips_brt.append(0)
    with gzip.open(projfile[ii], 'r') as fproj:
        next(fproj); next(fproj)
        for line in fproj:
            data = line.split()
            r2d = uc.pc2arcsec(4.52, float(data[0]))
            #print(r2d, float(data[0]))
            if r2d > 90: break
                
            #id_index = id0_bin.index(int(float(data[13])))
            #eclips_flag = eclipsing[id_index]
            #detec_flag = detection_prob[id_index]
            #Pday = period_bin[id_index]
            #tot_obsmag_index = tot_obsmag[id_index]
            #bin_obsmag0_index = bin_obsmag0[id_index]; bin_obsmag1_index = bin_obsmag1[id_index]
            #sin_obsmag_index = sin_obsmag[id_index]
            
            if int(data[2]) == 1:
                #print(r2d)
                id_index = id0_bin.index(int(float(data[13])))
                eclips_flag = eclipsing[id_index]
                detec_flag = detection_prob[id_index]
                Pday = period_bin[id_index]
                Ecc = ecc_bin[id_index]
                #tot_obsmag_index = tot_obsmag[id_index]
                bin_obsmag0_index = bin_obsmag0[id_index]; bin_obsmag1_index = bin_obsmag1[id_index]
                
                if float(data[10])>float(data[11]):
                    q_frac = float(data[11])/float(data[10])
                else:
                    q_frac = float(data[10])/float(data[11])
                #if tot_obsmag_index <= 25:
                #    n_all_brt[ii]+=1
                #    n_bin_brt[ii]+=1
                #print(data)
                if Pday < 4. or q_frac < 0.95 or Ecc > 0.5:  #detec_flag < 0.1 or
                    continue
                if eclips_flag == 1:
                    n_eclips[ii]+=1
                    #eclips_bin[0].append(int(data[13])); eclips_bin[1].append(int(data[14]))
                    #eclips_bin[2].append(int(data[5])); eclips_bin[3].append(int(data[6]))
                    #eclips_bin[4].append(float(data[10])); eclips_bin[5].append(float(data[11]))
                    
                    if 17 <= bin_obsmag0_index <= 25 and 17 <= bin_obsmag1_index <= 25:
                        n_eclips_brt[ii]+=1
                        print(Ecc)
                        #m0_eclips.append(float(data[10])); m1_eclips.append(float(data[11]))
            #else:
            #    n_all[ii]+=1
            #    id_index = id_sin.index(int(float(data[12]))) 
            #    sin_obsmag_index = sin_obsmag[id_index]                     
              
            #    if sin_obsmag_index >= 25:
            #        n_all_brt[ii]+=1

    print(n_eclips)
    print(n_eclips_brt)
    #print(m0_eclips, m1_eclips)                  
    #print(eclips_bin)
                                      
    #print(ii)
#tot_obsmag = np.array(tot_obsmag)
#print(len(tot_obsmag[17 <= tot_obsmag <= 25][0]))                    
print(np.mean(n_eclips))
print(np.mean(n_eclips_brt))

