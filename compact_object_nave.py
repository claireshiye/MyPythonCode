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

twopi=2.*np.pi
yearsc=3.1557*10**7
Kconst=9.87*10**-48 ##yr/G^2
Gconst=6.674*10**-8 ##cm3*g-1*s-2
Gconst_sun = 4.30091*10**-3 ##pc*M_sun**-1*(km/s)^2
clight=3*10**10 ##cm/s
Msun=2*10**33 ##gram
AU=1.496*10**13  ##cm
AU_Rsun=214.93946938362 ##AU to R_sun
PC=3.086*10**18  ##cm
PC_Rsun = 44334448.0068964 ##pc to R_sun

##Calculate average properties and standard deviations from best-fit and near-fit models

sb_const = 5.670374419*10**(-8)  ##W*m^−2*K^−4, Stefan-Boltzmann Constant
Lsun = 3.828*10**26  ##Watts
Rsun = 6.957*10**8 ##meters

import scripts, scripts1, scripts2, scripts3

def LtoT(lumi, radi):
    r_meter = radi*Rsun
    Area = 2*twopi*r_meter**2
    return pow((lumi*Lsun)/sb_const/Area, 1./4.)



def med3(array):
    if len(array)>4: return np.median(array)
    else: return 0
    

def get_turnoff(snapshot,Nbins=200):
    '''Purpose: Find the MS turnoff, defined as the luminosity (and temperature)
    of the highest-temperature upper-MS single that is NOT a blue straggler.
    Inputs: snapshot array and bin resolution for both luminosity and temperature axes of an HR diagram.
    Outputs: turnoff point (luminosity,temperature). Units: (Lsun,Kelvin).'''
    
    T0, L0 = [], []
    L_bin_edges = np.logspace(-0.5,0.5,int(Nbins+1))
    #print(L_bin_edges)
    L_bins = [[] for i in range(int(Nbins))]
    T_bins = [[] for i in range(int(Nbins))]
    
    ##For CMC-COSMIC ver 1.0.0
    binflag = snapshot.data['binflag']
    L = snapshot.data['luminosity_LSUN'][binflag != 1]
    ktype = snapshot.data['startype'][binflag != 1]
    R = snapshot.data['radius_RSUN'][binflag != 1]
    
    L0 = L[ktype == 1]; R0 = R[ktype == 1]
    T0 = LtoT(L0, R0)
    L0 = np.array(L0); T0 = np.array(T0)
    #print(L0, T0)
       
    #print(L0, T0)
    for j,Lj in enumerate(L0):
        for b in range(int(Nbins)):
            if L_bin_edges[b] <= Lj < L_bin_edges[b+1]: T_bins[b].append(T0[j]); L_bins[b].append(L0[j])

    T_meds = [med3(bin) for bin in T_bins]                  # List the median temperature of each bin
    turnoff_bin_index = T_meds.index(np.max(T_meds))        # Find the bin with the highest median temperature
    
    return (np.median(L_bins[turnoff_bin_index]), T_meds[turnoff_bin_index])

sys.path.insert(1, '/projects/b1095/syr904/MyCodes/cmctoolkit')
import cmctoolkit as cmct


paths = ['/projects/b1095/syr904/cmc/CMC-COSMIC/master_tc_test/ver_0601/MOCHA47Tuc_elson_rv4_3e6_tcon/',
         '/projects/b1095/syr904/cmc/CMC-COSMIC/master_tc_test/ver_0601/MOCHA47Tuc_elson_rv3.5_3e6_tcon_fb10/',
         '/projects/b1095/syr904/cmc/CMC-COSMIC/master_tc_test/ver_0601/MOCHA47Tuc_elson_rv4_3e6_tcon_rt171/']

#t_span_up = [1207., 13600., 13100.]; t_span_down = [8980., 9950., 8990.]
t_span_up = [13550., 13600., 13100.]; t_span_down = [8300., 9950., 8990.]

r_c = []; r_h = []; mass = []
r_cobs = []; r_hl = []
Nbh_inter = []; Nbhbh_inter = []; Nbhbin_inter = []
Nns_inter = []; Npsr_inter = []; Nmsp_inter = []; Nnsns_inter = []
n_cv = []
n_tot = [[],[],[],[],[],[]]
n_bss = []; bss_id = []
nsms_tc = []
for xx in range(3):
    print(paths[xx])
    l_conv = dyn.conv('l', paths[xx]+'initial.conv.sh')
    m_conv = dyn.conv('m', paths[xx]+'initial.conv.sh')
    t_conv = dyn.conv('t', paths[xx]+'initial.conv.sh')
    
    data_key = np.genfromtxt(paths[xx]+'snap_keys.txt', dtype = 'str')
    snapno = data_key[:,0]; snaptime = (data_key[:,1]).astype(np.float)*t_conv
    snapcodet = data_key[:,1]
    ##print(snaptime[10])

    ###################################################
    ##Average rc and rh
    with open(paths[xx]+'initial.dyn.dat', 'r') as fdyn:
        next(fdyn); next(fdyn)
        for line in fdyn:
            datadyn = line.split()
            if t_span_down[xx]/t_conv<=float(datadyn[0])<=t_span_up[xx]/t_conv:
                r_c.append(float(datadyn[7])*l_conv)
                r_h.append(float(datadyn[20])*l_conv)
                mass.append(float(datadyn[4])*m_conv)

    params_files = np.sort(glob(paths[xx]+'initial.snap*.cluster_params.dat'))
    if xx == 0:
        params_files = np.sort(glob(paths[xx]+'observed_profiles/initial.snap*.cluster_params.dat'))
    print(len(params_files))
    for vv in range(len(params_files)):
        datapara = np.genfromtxt(params_files[vv])
        if t_span_down[xx] <= datapara[0,0] <= t_span_up[xx]:
            r_cobs.append(datapara[0,9]); r_hl.append(datapara[0,10])
            
    print('rc, rh done')   
        
    ###################################################
    ##Average number of BHs
    bhs = np.genfromtxt(paths[xx]+'initial.bh.dat')
    bh_time = bhs[:,1]; Nbh_tot = bhs[:,2]; Nbhbh = bhs[:,5]; Nbhbin = bhs[:,4]

    for vv in range(len(bh_time)):
        thetime = bh_time[vv]*t_conv
        if t_span_down[xx] <= thetime <= t_span_up[xx]:
            Nbh_inter.append(Nbh_tot[vv])
            Nbhbh_inter.append(Nbhbh[vv])
            Nbhbin_inter.append(Nbhbin[vv])
            
    print('bh done')
            
    ###################################################
    ##Average number of NSs, young pulsars and MSPs
    nss = np.genfromtxt(paths[xx]+'initial.ns.dat')
    ns_time = nss[:,0]; Nns_tot = nss[:,1]; Npsr = nss[:,5]; Nmsp = nss[:,6]; Nnsns = nss[:,7]

    for vv in range(len(ns_time)):
        thetime = ns_time[vv]
        if t_span_down[xx] <= thetime <= t_span_up[xx]:
            Nns_inter.append(Nns_tot[vv])
            Npsr_inter.append(Npsr[vv])
            Nmsp_inter.append(Nmsp[vv])
            Nnsns_inter.append(Nnsns[vv])
        
    print('ns done')
                 
    ###################################################
    ##Average number of CVs
    wdms_files = np.sort(glob(paths[xx]+'WDMS*.dat'))
    print(len(wdms_files))
    for vv in range(len(wdms_files)):
        wdms_temp = wdms_files[vv].replace(paths[xx]+'WDMS', '')
        wdms_snapno = wdms_temp.split('.')[0]
        #print(wdms_snapno)
        t_wdms = snaptime[int(wdms_snapno)]
        
        Ncv = 0
        if t_span_down[xx] <= t_wdms <= t_span_up[xx]:
            #print(wdms_files[vv])
            wd_file = np.genfromtxt(wdms_files[vv])
            m1 = wd_file[:,3]; rad1 = wd_file[:,9]
            id0 = wd_file[:,0]; id1 = wd_file[:,1]
    
            for kk in range(len(rad1)):
                if rad1[kk]>=1.:
                    Ncv+=1
                
            n_cv.append(Ncv)

    print('cv done')
    
    ###################################################
    ##Average number of LMXBs
    nswd_files = np.sort(glob(paths[xx]+'NSWD*.dat'))
    nsg_files = np.sort(glob(paths[xx]+'NSGiant*.dat'))
    nsms_files = np.sort(glob(paths[xx]+'NSMS*.dat'))
    bhwd_files = np.sort(glob(paths[xx]+'BHWD*.dat'))
    bhg_files = np.sort(glob(paths[xx]+'BHGiant*.dat'))
    bhms_files = np.sort(glob(paths[xx]+'BHMS*.dat'))
    
    msp_files = np.sort(glob(paths[xx]+'MSP*.dat'))
    
    print(len(nswd_files), len(nsg_files), len(nsms_files))
    for vv in range(len(nswd_files)):
        nswd_temp = nswd_files[vv].replace(paths[xx]+'NSWD', '')
        nswd_snapno = nswd_temp.split('.')[0]
        #print(wdms_snapno)
        t_nswd = snaptime[int(nswd_snapno)]

        if t_span_down[xx] > t_nswd or t_nswd > t_span_up[xx]:
            continue
        
        print(nswd_files[vv])
        nswd_data = np.genfromtxt(nswd_files[vv])
        nsg_data = np.genfromtxt(nsg_files[vv])
        nsms_data = np.genfromtxt(nsms_files[vv])
        
        bhwd_data = np.genfromtxt(bhwd_files[vv])
        bhg_data = np.genfromtxt(bhg_files[vv])
        bhms_data = np.genfromtxt(bhms_files[vv])
        
        msp_data = np.genfromtxt(msp_files[vv])
        id0_msp = msp_data[:,10]; id1_msp = msp_data[:,11]
        
        Ntot = [0, 0, 0, 0 ,0 ,0]
        n_nsms_tc = 0

        all_files = [nswd_data, nsg_data, nsms_data, bhwd_data, bhg_data, bhms_data]
        for kk in range(len(all_files)):
            if len(all_files[kk])==0:
                n_tot[kk].append(Ntot[kk])
                continue
            
            if isinstance(all_files[kk][0], float):
                m1 = all_files[kk][3]; rad1= all_files[kk][9]
                id0 = all_files[kk][0]; id1 = all_files[kk][1]; tcflag = all_files[kk][12]
                
                if rad1>=1. and id0 not in id0_msp:
                    Ntot[kk]+=1
                    if kk==2 and tcflag==91:
                        Ntot[kk]-=1
                        
                if kk==2 and tcflag==91 and id0 not in id0_msp: 
                    n_nsms_tc+=1
                        
            else:
                m1 = all_files[kk][:,3]; rad1= all_files[kk][:,9]
                id0 = all_files[kk][:,0]; id1 = all_files[kk][:,1]; tcflag = all_files[kk][:,12]
                
                for ii in range(len(rad1)):  
                    if rad1[ii]>=1. and id0[ii] not in id0_msp:
                        Ntot[kk]+=1
                        if kk==2 and tcflag[ii]==91:
                            Ntot[kk]-=1
                            
                    if kk==2 and tcflag[ii]==91 and id0[ii] not in id0_msp: 
                        n_nsms_tc+=1
                            
            n_tot[kk].append(Ntot[kk])
        nsms_tc.append(n_nsms_tc)

    print('lmxb done, tidal capture nsms done')

    ###################################################
    ####Average number of BSSs 
    for ii in range(0, len(snaptime), 3):
        if t_span_down[xx] <= snaptime[ii] <= t_span_up[xx]:
            thekey = '/'+str(int(snapno[ii]))+'(t='+snapcodet[int(snapno[ii])]+')'
            snap_h5 = cmct.Snapshot(fname=paths[xx]+'initial.snapshots.h5', snapshot_name=thekey, 
                         conv=paths[xx]+'initial.conv.sh', 
                         dist=4.52, # distance to cluster in kpc
                         z=0.0038)

            Lto, Tto = get_turnoff(snap_h5,400)
            #print(Lto, Tto)
            print(snapno[ii])
            
            L_bss = []; T_bss = []; k0_bss = []; k1_bss = []
            Temp = []; Ltot = []
    
    
            ##For CMC-COSMIC ver 1.0.0
            binflag = snap_h5.data['binflag']
            rgc = snap_h5.data['r']*l_conv
            rgcsin = np.array(rgc[binflag != 1]); rgcbin = np.array(rgc[binflag == 1])
            
            Lsin = snap_h5.data['luminosity_LSUN'][binflag != 1]
            Lbin0 = np.array(snap_h5.data['bin_star_lum0_LSUN'][binflag == 1])
            Lbin1 = np.array(snap_h5.data['bin_star_lum1_LSUN'][binflag == 1])
            ktype = snap_h5.data['startype'][binflag != 1]
            kbin0 = np.array(snap_h5.data['bin_startype0'][binflag == 1])
            kbin1 = np.array(snap_h5.data['bin_startype1'][binflag == 1])
            Rsin = snap_h5.data['radius_RSUN'][binflag != 1]
            Rbin0 = np.array(snap_h5.data['bin_star_radius0_RSUN'][binflag == 1])
            Rbin1 = np.array(snap_h5.data['bin_star_radius1_RSUN'][binflag == 1])
            idsin = snap_h5.data['id'][binflag != 1]
            idbin0 = np.array(snap_h5.data['id0'][binflag == 1]); idbin1 = np.array(snap_h5.data['id1'][binflag == 1])
            
            ###For single stars
            Lsin_nobh = np.array(Lsin[ktype != 14]); Rsin_nobh = np.array(Rsin[ktype != 14])
            Tsin_nobh = LtoT(Lsin_nobh, Rsin_nobh)
            ktype_nobh = ktype[ktype != 14]
            rgcsin_nobh = rgcsin[ktype != 14]
            idsin_nobh = idsin[ktype != 14]
            
            ##Selecting BSS
            Ltemp1 = Lsin_nobh[Lsin_nobh >= 2*Lto]; Ttemp1 = Tsin_nobh[Lsin_nobh >= 2*Lto]
            ktype_temp1 = ktype_nobh[Lsin_nobh >= 2*Lto]
            Ltemp2 = Ltemp1[Ttemp1 >= Tto]; Ttemp2 = Ttemp1[Ttemp1 >= Tto]
            ktype_temp2 = ktype_temp1[Ttemp1 >= Tto]
            idsin_temp1 = idsin_nobh[Lsin_nobh >= 2*Lto]
            idsin_temp2 = idsin_temp1[Ttemp1 >= Tto]
            r_temp1 = rgcsin_nobh[Lsin_nobh >= 2*Lto]; r_temp2 = r_temp1[Ttemp1 >= Tto]
            
            r_bss = np.concatenate((r_temp2[ktype_temp2 == 0], r_temp2[ktype_temp2 == 1]), axis=None)
            L_bss = L_bss + list(Ltemp2[ktype_temp2 == 0]) + list(Ltemp2[ktype_temp2 == 1])
            T_bss = T_bss + list(Ttemp2[ktype_temp2 == 0]) + list(Ttemp2[ktype_temp2 == 1])
            bss_id = bss_id + list(idsin_temp2[ktype_temp2 == 0]) + list(idsin_temp2[ktype_temp2 == 1])
            bss_id = [int(i) for i in bss_id]
            
            L_bss = np.array(L_bss); T_bss = np.array(T_bss)
            L_bss = list(L_bss); T_bss = list(T_bss)
            
            ###For binary stars
            for kk in range(len(Lbin0)):
                if kbin0[kk] != 14 or kbin1[kk] != 14:
                    temperature0 = LtoT(Lbin0[kk], Rbin0[kk])
                    temperature1 = LtoT(Lbin1[kk], Rbin1[kk])
                    temp_eff = (temperature0*Lbin0[kk]+temperature1*Lbin1[kk])/(Lbin0[kk] + Lbin1[kk])
                    Temp.append(temp_eff)
                    Ltot.append(Lbin0[kk] + Lbin1[kk])
                        
                    if Lbin0[kk] > 2*Lto and temperature0 > Tto and (kbin0[kk] == 1 or kbin0[kk] == 0):
                        L_bss.append(Lbin0[kk]); T_bss.append(temperature0)
                        bss_id.append(int(idbin0[kk]))
                        
                    if Lbin1[kk] > 2*Lto and temperature1 > Tto and (kbin1[kk] == 1 or kbin1[kk] == 0):
                        L_bss.append(Lbin1[kk]); T_bss.append(temperature1)
                        bss_id.append(int(idbin1[kk]))
    
            n_bss.append(len(L_bss))
    
    print('bss done')
    
    
print(np.mean(mass), np.mean(r_c), np.mean(r_h))
print(np.std(mass), np.std(r_c), np.std(r_h))
print(np.mean(r_cobs), np.mean(r_hl))
print(uc.pc2arcsec(4.52, np.mean(r_cobs))/60., uc.pc2arcsec(4.52, np.mean(r_hl))/60.)
print(uc.pc2arcsec(4.52, np.std(r_cobs))/60., uc.pc2arcsec(4.52, np.std(r_hl))/60.)
r_cobs_arcmin = uc.pc2arcsec(4.52, np.array(r_cobs))/60.
r_hl_arcmin = uc.pc2arcsec(4.52, np.array(r_hl))/60.
print(np.mean(r_cobs_arcmin), np.mean(r_hl_arcmin))
print(np.std(r_cobs_arcmin), np.std(r_hl_arcmin))

print(np.mean(Nbh_inter), np.mean(Nbhbh_inter), np.mean(Nbhbin_inter))
print(np.std(Nbh_inter), np.std(Nbhbh_inter), np.std(Nbhbin_inter))

#print(len(Nmsp))
print(Nns_inter)
print(np.mean(Nns_inter), np.mean(Npsr_inter),  np.mean(Nmsp_inter), np.mean(Nnsns_inter))
print(np.std(Nns_inter), np.std(Npsr_inter),  np.std(Nmsp_inter), np.std(Nnsns_inter))
print(np.mean(np.array(Npsr_inter)-np.array(Nmsp_inter)))
print(np.std(np.array(Npsr_inter)-np.array(Nmsp_inter)))

print(np.mean(n_cv), np.std(n_cv))

print(np.sum(n_tot, axis = 0))
print(np.mean(np.sum(n_tot, axis=0)), np.std(np.sum(n_tot, axis=0)))

print('n_tc_nsms', nsms_tc, np.mean(nsms_tc), np.std(nsms_tc))

print(n_bss)
print(np.mean(n_bss),np.std(n_bss))
