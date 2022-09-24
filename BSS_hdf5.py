import numpy as np
from glob import glob
import collections
from collections import Counter
import os, sys
import re, gzip
import scipy.stats as ss
import scripts
import ns
import dynamics as dyn
import scripts, scripts1, scripts2, scripts3

sys.path.insert(1, '/projects/b1095/syr904/MyCodes/cmctoolkit')
import cmctoolkit as cmct


sb_const = 5.670374419*10**(-8)  ##W*m^−2*K^−4, Stefan-Boltzmann Constant
Lsun = 3.828*10**26  ##Watts
Rsun = 6.957*10**8 ##meters

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


path = '/projects/b1095/syr904/cmc/CMC-COSMIC/master_tc_test/ver_0601/MOCHA47Tuc_elson_rv4_3e6_tcon/'
l_conv = dyn.conv('l', path+'initial.conv.sh')
t_conv = dyn.conv('t', path+'initial.conv.sh')
all_keys = np.genfromtxt(path+'snap_keys.txt', dtype = 'str')
all_snapno = all_keys[:,0]; all_snaptime = all_keys[:,1]


##BSS at the present-day
snap_h5_726 = cmct.Snapshot(fname=path+'initial.snapshots.h5', snapshot_name='/733(t=0.47867046)', 
                     conv=path+'initial.conv.sh', 
                     dist=4.52, # distance to cluster in kpc
                     z=0.0038)
Lto_726, Tto_726 = get_turnoff(snap_h5_726,400)

bss_726 = np.genfromtxt(path+'observed_profiles/allBSS_733.dat')
id0_726 = bss_726[:,0]; id1_726 = bss_726[:,1]
k0_726 = bss_726[:,2]; k1_726 = bss_726[:,3]
T_726 = bss_726[:,4]; L_726 = bss_726[:,5]
tform_726 = np.zeros(len(id0_726))


nbss_snap_726 = []; nbss_snap = []; tsnap = []; Tto_snap = []; Lto_snap = []
allbss_id = []
for ii in range(785):
    print(ii)
    thekey = '/'+str(ii)+'(t='+all_snaptime[ii]+')'
    snap_h5 = cmct.Snapshot(fname=path+'initial.snapshots.h5', snapshot_name=thekey, 
                     conv=path+'initial.conv.sh', 
                     dist=4.52, # distance to cluster in kpc
                     z=0.0038)

    Lto, Tto = get_turnoff(snap_h5,400)

    L_bss = []; T_bss = []; k0_bss = []; k1_bss = []; id_bss = []
    Temp = []; Ltot = []
    L_bss_726 = []; T_bss_726 = []
    
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
    id_bss = id_bss + list(idsin_temp2[ktype_temp2 == 0]) + list(idsin_temp2[ktype_temp2 == 1])
    id_bss = [int(i) for i in id_bss]

    allbss_id = allbss_id + list(idsin_temp2[ktype_temp2 == 0]) + list(idsin_temp2[ktype_temp2 == 1])
    allbss_id = [int(i) for i in allbss_id]


    ##Selecting BSS using the turnoff values at snapshot 726
    Ltemp1_726 = Lsin_nobh[Lsin_nobh >= 2*Lto_726]; Ttemp1_726 = Tsin_nobh[Lsin_nobh >= 2*Lto_726]
    ktype_temp1_726 = ktype_nobh[Lsin_nobh >= 2*Lto_726]
    Ltemp2_726 = Ltemp1_726[Ttemp1_726 >= Tto_726]; Ttemp2_726 = Ttemp1_726[Ttemp1_726 >= Tto_726]
    ktype_temp2_726 = ktype_temp1_726[Ttemp1_726 >= Tto_726]
    L_bss_726 = L_bss_726 + list(Ltemp2_726[ktype_temp2_726 == 0]) + list(Ltemp2_726[ktype_temp2_726 == 1])
    T_bss_726 = T_bss_726 + list(Ttemp2_726[ktype_temp2_726 == 0]) + list(Ttemp2_726[ktype_temp2_726 == 1])

    #L_bss = np.array(L_bss); T_bss = np.array(T_bss)
    #L_bss = list(L_bss); T_bss = list(T_bss)
    
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
                id_bss.append(int(idbin0[kk]))
                allbss_id.append(int(idbin0[kk]))

               
            if Lbin1[kk] > 2*Lto and temperature1 > Tto and (kbin1[kk] == 1 or kbin1[kk] == 0):
                L_bss.append(Lbin1[kk]); T_bss.append(temperature1)
                id_bss.append(int(idbin1[kk]))
                allbss_id.append(int(idbin1[kk]))


            ##Selecting BSS using the turnoff values at snapshot 726
            if Lbin0[kk] > 2*Lto_726 and temperature0 > Tto_726 and (kbin0[kk] == 1 or kbin0[kk] == 0):
                L_bss_726.append(Lbin0[kk]); T_bss_726.append(temperature0)
                
            if Lbin1[kk] > 2*Lto_726 and temperature1 > Tto_726 and (kbin1[kk] == 1 or kbin1[kk] == 0):
                L_bss_726.append(Lbin1[kk]); T_bss_726.append(temperature1)


    nbss_snap.append(len(L_bss)); nbss_snap_726.append(len(L_bss_726))
    tsnap.append(t_conv*float(all_snaptime[ii]))
    Tto_snap.append(Tto); Lto_snap.append(Lto)


    ##Found formation time of BSS_726
    for xx in range(len(id0_726)):
        for yy in range(len(id_bss)):
            if int(id_bss[yy]) == int(id0_726[xx]) and tform_726[xx] == 0.:
                tform_726[xx] = t_conv*float(all_snaptime[ii])


np.savetxt(path+'observed_profiles/allBSS_733.dat', np.c_[id0_726, id1_726, k0_726, k1_726, T_726, L_726, tform_726],
          fmt = '%d %d %d %d %f %f %f', delimiter = ' ',
          header = '1.id0 2.id1 3.k0 4.k1 5.T(K) 6.L(Lsun) 7.Tform(Myr)', comments = '#')
#np.savetxt(path+'nbss_time.dat', np.c_[tsnap, nbss_snap, nbss_snap_726, Tto_snap, Lto_snap],
#    fmt = '%f %d %d %f %f', delimiter = ' ', header = '1.Time(Myr) 2.Nbss 3.Nbss_726(use turnoffs at snapshot726) 4.Tto(K) 5.Lto(Lsun)', comments = '#')


#print(bss_id)
bss_id_unique = Counter(bss_id).keys()
collfile = scripts1.readcollfile(path+'initial.collision.log')
n_bss_coll = 0
bss_coll = []
for ii in range(len(collfile)):
    line = collfile[ii].split()
    if int(line[3]) in bss_id_unique:
        n_bss_coll+=1
        bss_coll.append(int(line[3]))
#bss_coll =  list(Counter(bss_coll).keys())
    
bss_id_unique = list(bss_id_unique)
for kk in range(len(bss_id_unique)):
    if bss_id_unique[kk] in bss_coll:
        continue
    else:
        print(bss_id_unique[kk])
print(np.sort(bss_id_unique))
print(np.sort(bss_coll))
                
    
print(n_bss_coll, len(bss_coll))
print(len(bss_id_unique))
