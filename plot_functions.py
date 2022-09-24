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
from scipy import stats
import matplotlib.cm as cm
import matplotlib as mpl
import random
from random import shuffle
import json
import scripts
import scripts1
import scripts2
import scripts3
import dynamics as dyn

import ecc_calc as gwcalc
import unit_convert as uc
import merger_rate_calculator as mr
import ns_tidalcapture as tc
import psr_catalog as pc

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
PC_Rsun = 44334448.0068964 ##pc to R_sun
PC_AU = 206265 ##pc to AU


###Get model data
##Pulsars at the last timestep
datamsp=np.genfromtxt('/projects/b1095/syr904/cmc/cmc-mpi-tidalcapture/rvgrid/standard_models_tcon/msp_last_tcon.dat')
datapsr=np.genfromtxt('/projects/b1095/syr904/cmc/cmc-mpi-tidalcapture/rvgrid/standard_models_tcon/normalpsr_last_tcon.dat')
m0_msp=datamsp[:,10]; m1_msp=datamsp[:,11]; sma_msp=datamsp[:,16]; ecc_msp=datamsp[:,17]
k0_msp=datamsp[:,14]; k1_msp=datamsp[:,15]; B_msp=datamsp[:,4]; P_msp=datamsp[:,5]
id0_msp=datamsp[:,12]; id1_msp=datamsp[:,13]; model_msp=datamsp[:,0]

m0_psr=datapsr[:,10]; m1_psr=datapsr[:,11]; sma_psr=datapsr[:,16]; ecc_psr=datapsr[:,17]
k0_psr=datapsr[:,14]; k1_psr=datapsr[:,15]; B_psr=datapsr[:,4]; P_psr=datapsr[:,5]
id0_psr=datapsr[:,12]; id1_psr=datapsr[:,13]; model_psr=datapsr[:,0]

Pdot_psr=Kconst*yearsc*B_psr*B_psr/P_psr
Pdot_msp=Kconst*yearsc*B_msp*B_msp/P_msp
Pdot_msp=Pdot_msp *10**19*1000; Pdot_psr=Pdot_psr*10**19*1000
P_msp = P_msp*1000; P_psr = P_psr*1000


##Create distinc ids for different models
id0_new_psr = []
for x in range(len(model_psr)):
    id0_new_psr.append(str(int(model_psr[x]))+str(int(id0_psr[x])))
    
id0_new_msp = []
for y in range(len(model_msp)):
    id0_new_msp.append(str(int(model_msp[y]))+str(int(id0_msp[y])))



##NSWD and NSMS binaries that are not MSPs or young pulsars at the last timestep
data13ms=np.genfromtxt('/projects/b1095/syr904/cmc/cmc-mpi-tidalcapture/rvgrid/standard_models_tcon/NSMS_last_rvgrid_tcon.dat')
data13wd=np.genfromtxt('/projects/b1095/syr904/cmc/cmc-mpi-tidalcapture/rvgrid/standard_models_tcon/NSWD_last_rvgrid_tcon.dat')

id0_13ms=data13ms[:,2]; id1_13ms=data13ms[:,3]; m0_13ms=data13ms[:,4]; m1_13ms=data13ms[:,5]
sma_13ms=data13ms[:,8]; ecc_13ms=data13ms[:,9]; tcflag_13ms=data13ms[:,14]
model_13ms=data13ms[:,0]
k0_13ms=data13ms[:,6]; k1_13ms=data13ms[:,7]

id0_13wd=data13wd[:,2]; id1_13wd=data13wd[:,3]; m0_13wd=data13wd[:,4]; m1_13wd=data13wd[:,5]
sma_13wd=data13wd[:,8]; ecc_13wd=data13wd[:,9]; tcflag_13wd=data13wd[:,14]
model_13wd=data13wd[:,0]
k0_13wd=data13wd[:,6]; k1_13wd=data13wd[:,7]

##Create distinct ids for different models
id0_new_13ms = []
for x in range(len(model_13ms)):
    id0_new_13ms.append(str(int(model_13ms[x]))+str(int(id0_13ms[x])))
    
id0_new_13wd = []
for y in range(len(model_13wd)):
    id0_new_13wd.append(str(int(model_13wd[y]))+str(int(id0_13wd[y])))

                                        
###Get observed data
P, Pdot, Binflag, Namespin, Period, Ecc, Mc, Names, Pall, Bfall, Nameall = pc.readdata_freire()
print(len(Pall), len(Bfall))
print(Ecc)
datarb=np.genfromtxt('/projects/b1095/syr904/projects/PULSAR/data_observed/GCredback.dat', dtype=str)
databw=np.genfromtxt('/projects/b1095/syr904/projects/PULSAR/data_observed/GCblackwidow.dat', dtype=str)
namerb=datarb[:,0]; namebw=databw[:,0]


###Plotting
plt.rcParams['figure.figsize'] = [8,16]
rdot=mlines.Line2D([], [],  linestyle = 'None', color='r', marker='o',
                  markersize=7, label='Redbacks')
kdot=mlines.Line2D([], [],  linestyle = 'None', color='k', marker='o',
                  markersize=7, label='Black widows')
bdot=mlines.Line2D([], [],  linestyle = 'None', color='b', marker='o',
                  markersize=7, label='Others')
gtri=mlines.Line2D([], [],  linestyle = 'None', color='g', marker='^',
                  markersize=7, label='Others')
star=mlines.Line2D([], [],  linestyle = 'None', color='g', marker='*',
                  markersize=7, label='Tidal capture binaries')
square=mlines.Line2D([], [],  linestyle = 'None', color='g', marker='s',
                  markersize=7, label='Giant collision binaries')
tri=mlines.Line2D([], [],  linestyle = 'None', color='g', marker='^',
                  markersize=7, label='Others')


fig, axs = plt.subplots(nrows=2, sharey=False)
fig.subplots_adjust(wspace=0.05)

##Model data
for ii in range(len(id0_new_13ms)):
    pb_13ms = uc.au_to_period(sma_13ms[ii], m0_13ms[ii], m1_13ms[ii])
    if m1_13ms[ii]<0.04:
        print(m1_13ms[ii], pb_13ms, id0_new_13ms[ii], id1_13ms[ii], k1_13ms[ii], tcflag_13ms[ii])
    if tcflag_13ms[ii]==91:
        if id0_new_13ms[ii] in id0_new_msp:# or id0_new_13ms[ii] in id0_new_psr:
            axs[0].scatter(m1_13ms[ii], pb_13ms, color='red', marker = '*', s=50, lw = 1)
            #axs[0].scatter(m1_13ms[ii], pb_13ms, color='red', marker = 'o', 
            #               lw = 1, facecolor = 'none')
            axs[1].scatter(ecc_13ms[ii], pb_13ms,  color='red', marker = '*', s=50, lw = 1)
        else:
            axs[0].scatter(m1_13ms[ii], pb_13ms, color='red', marker = '*', s=50, 
                           lw = 1, facecolors='none')
            axs[1].scatter(ecc_13ms[ii], pb_13ms,  color='red', marker = '*', s=50, 
                           lw = 1, facecolors='none')
    elif tcflag_13ms[ii]==81:
        if id0_new_13ms[ii] in id0_new_msp:# or id0_new_13ms[ii] in id0_new_psr:
            axs[0].scatter(m1_13ms[ii], pb_13ms, color='red', marker = 's', lw = 1)
            axs[1].scatter(ecc_13ms[ii], pb_13ms,  color='red', marker = 's', lw = 1)
        else:
            axs[0].scatter(m1_13ms[ii], pb_13ms, color='red', marker = 's', 
                           lw = 1, facecolors='none')
            axs[1].scatter(ecc_13ms[ii], pb_13ms,  color='red', marker = 's', 
                           lw = 1, facecolors='none')
        
    else:
        if id0_new_13ms[ii] in id0_new_msp:# or id0_new_13ms[ii] in id0_new_psr:
            print(id0_new_13ms[ii])
            axs[0].scatter(m1_13ms[ii], pb_13ms, color='red', marker = '^', 
                           lw = 1, alpha=0.6)
            axs[1].scatter(ecc_13ms[ii], pb_13ms,  color='red', marker = '^', lw = 1)
        #else:
        #    axs[0].scatter(m1_13ms[ii], pb_13ms, color='red', marker = '^', lw = 1, facecolors='none')
        #    axs[1].scatter(ecc_13ms[ii], pb_13ms,  color='red', marker = '^', lw = 1, facecolors='none')

            
for hh in range(len(id0_new_13wd)):
    pb_13wd = uc.au_to_period(sma_13wd[hh], m0_13wd[hh], m1_13wd[hh])
    if m1_13wd[hh]<0.04:
        print(m1_13wd[hh], pb_13wd, id0_new_13wd[hh], id1_13wd[hh], k1_13wd[hh], tcflag_13wd[hh])
    if ecc_13wd[hh]>1:
        print(id0_13wd[hh], id1_13wd[hh])
    if tcflag_13wd[hh]==91:
        if k1_13wd[hh]==10:
            if id0_new_13wd[hh] in id0_new_msp:# or id0_new_13wd[hh] in id0_new_psr:
                axs[0].scatter(m1_13wd[hh], pb_13wd, color='b', marker = '*', s=50, lw = 1)
                axs[0].scatter(m1_13wd[hh], pb_13wd, color='b', marker = 'o', 
                               lw = 1, facecolor = 'none', s=70)
                axs[1].scatter(ecc_13wd[hh], pb_13wd,  color='b', marker = '*', s=50, lw = 1) 
            else:
                axs[0].scatter(m1_13wd[hh], pb_13wd, color='b', marker = '*', s=50, 
                               lw = 1, facecolors='none')
                axs[0].scatter(m1_13wd[hh], pb_13wd, color='b', marker = 'o', 
                               lw = 1, facecolor = 'none', s=70)
                axs[1].scatter(ecc_13wd[hh], pb_13wd,  color='b', marker = '*', s=50, 
                               lw = 1, facecolors='none') 
    
        elif k1_13wd[hh]>10:
            if id0_new_13wd[hh] in id0_new_msp:# or id0_new_13wd[hh] in id0_new_psr:
                axs[0].scatter(m1_13wd[hh], pb_13wd, color='b', marker = '*', s=50, lw = 1)
                axs[1].scatter(ecc_13wd[hh], pb_13wd,  color='b', marker = '*', s=50, lw = 1) 
            else:
                axs[0].scatter(m1_13wd[hh], pb_13wd, color='b', marker = '*', s=50,
                               lw = 1, facecolors='none')
                axs[1].scatter(ecc_13wd[hh], pb_13wd,  color='b', marker = '*', s=50,
                               lw = 1, facecolors='none') 
    
    elif tcflag_13wd[hh]==81:
        if k1_13wd[hh]==10:
            if id0_new_13wd[hh] in id0_new_msp:# or id0_new_13wd[hh] in id0_new_psr:
                axs[0].scatter(m1_13wd[hh], pb_13wd, color='b', marker = 's', lw = 1)
                axs[0].scatter(m1_13wd[hh], pb_13wd, color='b', marker = 'o', 
                               lw = 1, facecolor = 'none', s=70)
                axs[1].scatter(ecc_13wd[hh], pb_13wd,  color='b', marker = 's', lw = 1) 
            else:
                axs[0].scatter(m1_13wd[hh], pb_13wd, color='b', marker = 's', 
                               lw = 1, facecolors='none')
                axs[0].scatter(m1_13wd[hh], pb_13wd, color='b', marker = 'o', 
                               lw = 1, facecolor = 'none', s=70)
                axs[1].scatter(ecc_13wd[hh], pb_13wd,  color='b', marker = 's', s=50, 
                               lw = 1, facecolors='none') 
    
        elif k1_13wd[hh]>10:
            if id0_new_13wd[hh] in id0_new_msp:# or id0_new_13wd[hh] in id0_new_psr:
                axs[0].scatter(m1_13wd[hh], pb_13wd, color='b', marker = 's', lw = 1)
                axs[1].scatter(ecc_13wd[hh], pb_13wd,  color='b', marker = 's', lw = 1) 
            else:
                axs[0].scatter(m1_13wd[hh], pb_13wd, color='b', marker = 's', 
                               lw = 1, facecolors='none')
                axs[1].scatter(ecc_13wd[hh], pb_13wd,  color='b', marker = 's', 
                               lw = 1, facecolors='none') 
    else: 
        if k1_13wd[hh]==10:
            if id0_new_13wd[hh] in id0_new_msp:# or id0_new_13wd[hh] in id0_new_psr:
                axs[0].scatter(m1_13wd[hh], pb_13wd, color='b', marker = '^', 
                               lw = 1, alpha=0.6)
                axs[0].scatter(m1_13wd[hh], pb_13wd, color='b', marker = 'o', 
                               lw = 1, facecolor = 'none', s=70)
                axs[1].scatter(ecc_13wd[hh], pb_13wd,  color='b', marker = '^', lw = 1)  
        elif k1_13wd[hh]>10:
            if id0_new_13wd[hh] in id0_new_msp:# or id0_new_13wd[hh] in id0_new_psr:
                axs[0].scatter(m1_13wd[hh], pb_13wd, color='b', marker = '^', 
                               lw = 1, alpha=0.6)
                axs[1].scatter(ecc_13wd[hh], pb_13wd,  color='b', marker = '^', lw = 1)  

            
##Observed data
for i in range(len(Mc)):
    if Mc[i][0] != '>':
        if Names[i] in namerb:
            axs[0].scatter(float(Mc[i]), float(Period[i]), color='red', s = 10, 
                           alpha=0.7, zorder=2)
            if Ecc[i][0]!='<' and Ecc[i][0]!='-':
                axs[1].scatter(float(Ecc[i]), float(Period[i]),  color='red', alpha=0.7, zorder=2)
            elif Ecc[i][0]!='-':
                axs[1].errorbar(float(Ecc[i].split('<')[1]), float(Period[i]), xerr=0.02, 
                                xuplims=True,  marker='o', color='red', linestyle='none')
                
            
        elif Names[i] in namebw:
            axs[0].scatter(float(Mc[i]), float(Period[i]), color='k', s = 10, 
                           alpha=0.7, zorder=2)
            if Ecc[i][0]!='<':
                axs[1].scatter(float(Ecc[i]), float(Period[i]), color='k', alpha=0.7, zorder=2)
            else:
                axs[1].errorbar(float(Ecc[i].split('<')[1]), float(Period[i]), xerr=0.02, 
                                xuplims=True,  marker='o', color='k', linestyle='none')
                
            
        else:
            axs[0].scatter(float(Mc[i]), float(Period[i]), color='b', s = 10, 
                           alpha=0.7, zorder=1)
            if Ecc[i][0]!='<':
                axs[1].scatter(float(Ecc[i]), float(Period[i]), color='b', alpha=0.7, zorder=1)
            else:
                #print(Ecc[i])
                axs[1].errorbar(float(Ecc[i].split('<')[1]), float(Period[i]), xerr=0.02, 
                                xuplims=True,  marker='o', color='b', linestyle='none')
                
    
    else:
        axs[0].errorbar(float(Mc[i].split('>')[1]), float(Period[i].split('>')[1]), xerr=0.002, yerr=0.2, lolims=True, 
                        xlolims=True, marker='o', color='b', markersize = 4, linestyle='none')
        axs[1].errorbar(float(Ecc[i].split('>')[1]), float(Period[i].split('>')[1]), xerr=0.02, yerr=0.2, lolims=True, 
                        xlolims=True, marker='o', markerfacecolor='none', color='b', linestyle='none')
        
        
            
        
axs[0].set_yscale('log')
axs[0].set_xscale('log')
axs[0].set_ylim(0.01, 3000)
axs[0].set_xlim(0.001, 20)
axs[0].set_xlabel(r'$M_c\ (M_\odot)$')
axs[0].set_ylabel(r'$P_{orb}\ (d)$')
axs[0].legend(handles=[rdot, kdot, bdot, star, square, tri], loc='upper left', prop={'size': 12}, numpoints=1, frameon=True)

axs[1].set_yscale('log')
#axs[1].set_xscale('symlog')
axs[1].set_ylim(0.01, 3000)
#axs[1].set_xlim(-0.1, 1)
axs[1].set_xlabel(r'$eccentricity$')
#axs[1].set_ylabel(r'$P_{orb}\ (d)$')
