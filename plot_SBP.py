import numpy as np
import matplotlib.pyplot as plt
import gzip
import scripts
import constants
import blackhole
import conversions
import math
from scipy import integrate
import extract_observed_prop

cluster = 'NGC3201'

if cluster == 'm10':
	d_helio = 4.4
	rc_obs = 0.99 # pc
	rh_obs = 2.50
if cluster == 'm22':
	d_helio = 3.2
        rc_obs = 1.24 # pc
        rh_obs = 3.12
if cluster == 'NGC3201':
        #d_helio = 4.9
        d_helio = 5.7
	rc_obs = 1.85 # pc
        rh_obs = 4.42

#simulation = 'm10_400000e5_5.0_1.0_0.05_FULL'
#path = '/projects/b1011/sourav/new_runs/kick_grid/rv2/kickscale_0.1/'
path = '/projects/b1011/sourav/new_runs/kick_grid/rv1/kickscale_0.05/'
snapno = '0894'
#snapno = '0731'
#path = '/projects/b1011/sourav/new_runs/kick_grid/rv1/kickscale_0.1/'
#snapno = '1056'

########  MAKE SBP ##############

data = np.genfromtxt(path+'initial.snap'+snapno+'.2D_SBPLcut15.dat')
data2 = np.genfromtxt(cluster+'_data_2.dat')
arcsec = conversions.pc_to_arcsec(data[:,1],d_helio)
SB = conversions.SB_converter_tot(data[:,3])
SB_err = data[:,6]/data[:,5]*SB

fig = plt.figure()
ax = fig.add_subplot(1,1,1)

rc_obs = conversions.pc_to_arcsec(rc_obs,d_helio) 
rh_obs = conversions.pc_to_arcsec(rh_obs,d_helio)

##### get observed properties #####
props = np.genfromtxt(path+'initial.snap'+snapno+'.obs_params.dat')
rc_model = conversions.pc_to_arcsec(props[0,7],d_helio)
rh_model = conversions.pc_to_arcsec(props[0,8],d_helio)
M_total = np.str(props[0,12])
time = np.str(props[0,10])
N_BH = np.str(props[0,13])
######################

plt.scatter(10**data2[:,2],data2[:,3], marker='o', c='red',edgecolor='red',label=r'$\rm{Observed}$')
#plt.scatter(10**data3[:,2],data3[:,3], marker='o', c='yellow',edgecolor='yellow',label='Noyola and Gebhardt 2006')
plt.scatter(arcsec, SB, color='blue',label=r'$\rm{Model}$')
#plt.errorbar(arcsec,SB,yerr=SB_err, fmt='o',c='blue')
#plt.plot(np.unique(arcsec), np.poly1d(np.polyfit(arcsec,SB, 8))(np.unique(arcsec)),color='black',lw=1)
#plt.plot([rc_obs,rc_obs],[-10,100],lw=2.0,color='black')#,label='$r_c$')
#plt.plot([rc_model,rc_model],[-10,100],lw=2.0,color='blue')

#plt.plot([rh_obs,rh_obs],[-10,100],lw=2.0,ls='--',color='black')#,label='$r_h$')
#plt.plot([rh_model,rh_model],[-10,100],lw=2.0,ls='--',color='blue')

plt.xlim(0.3,3000)
plt.ylim(30,15)
ax.set_xscale('log')

plt.xlabel(r'$r\,(\rm{arcsec})$',fontsize=20)
plt.ylabel(r'$\mu_v$',fontsize=20)

plt.legend(loc=1,frameon='False')

plt.show()


