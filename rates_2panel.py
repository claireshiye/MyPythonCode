import numpy as np
import os,sys
import subprocess
import gzip
import scripts
import matplotlib.pyplot as plt
import history_cmc
import LISA_calculations as LISA
from scipy import stats
from scipy import interpolate
from scipy.signal import savgol_filter
from scipy import ndimage

data = np.genfromtxt('BHMScollisions_more_2.dat')
data2 = np.genfromtxt('NSMScollisions_2.dat')
data3 = np.genfromtxt('collapsars_2.dat')


data_obs = np.genfromtxt('observed_GRB.dat')
density = np.genfromtxt('observed_GRBdensity.dat')
density_L14 = np.genfromtxt('observed_GRB_L14.dat')

t_array = np.linspace(350,13700,400) #350

N = len(t_array)-2

##### initialize BH-MS arrays

z_array = []
R_array = []
R_8_array = []
scriptR_array = []
R_cum_array = []
R_8_cum_array = []
scriptR_8_array = []

R_array_upper = []
R_8_array_upper = []
scriptR_array_upper = []
R_cum_array_upper = []
R_8_cum_array_upper = []
scriptR_8_array_upper = []

R_array_lower = []
R_8_array_lower = []
scriptR_array_lower = []
R_cum_array_lower = []
R_8_cum_array_lower = []
scriptR_8_array_lower = []

R_cum = 0

##### initialize NS-MS arrays

R_array2 = []
scriptR_array2 = []
R_cum_array2 = []

R_array_upper2 = []
scriptR_array_upper2 = []
R_cum_array_upper2 = []

R_array_lower2 = []
scriptR_array_lower2 = []
R_cum_array_lower2 = []

##### initialize collapsar arrays

R_array3 = []
scriptR_array3 = []
R_cum_array3 = []

R_array_upper3 = []
scriptR_array_upper3 = []
R_cum_array_upper3 = []

R_array_lower3 = []
scriptR_array_lower3 = []
R_cum_array_lower3 = []


for i in range(1,len(t_array)-1):
	t_lower = t_array[i]
	t_upper = t_array[i+1]
	t_mid = 0.5*(t_lower+t_upper)
	z_upper = LISA.zAtLookbackTime((13700.-t_lower)/1000.)
	z_lower = LISA.zAtLookbackTime((13700.-t_upper)/1000.)
	z_mid = LISA.zAtLookbackTime((13700.-t_mid)/1000.)
	D_upper = LISA.comovingDistance(z_upper)
	D_lower = LISA.comovingDistance(z_lower)

	#### BH-MS
	count_MS = 0
	count_MS_8 = 0
	count_G = 0
	
	count_MS_upper = 0
        count_MS_8_upper = 0
        count_G_upper = 0

	count_MS_lower = 0
        count_MS_8_lower = 0
        count_G_lower = 0
        #### NS-MS
        count_MS2 = 0
	count_G2 = 0

        count_MS_upper2 = 0
	count_G_upper2 = 0	

        count_MS_lower2 = 0
        count_G_lower2 = 0
        #### collapsar
        count_MS3 = 0
        count_G3 = 0

        count_MS_upper3 = 0
        count_G_upper3 = 0

        count_MS_lower3 = 0
        count_G_lower3 = 0


	##### LOOK THROUGH BH-MS COLLISION #######
	for j in range(len(data)):
		t = data[j,0]
		if t > t_lower and t < t_upper:
			#2190.86035452 122.887084122 0.203042 52.5302 0.0 14.0 11.0
			k_star = data[j,4]
			m_BH = data[j,3]
			m_star = data[j,2]
			model = data[j,6]
			if k_star <= 1:
				count_MS += 1
				if m_star >= 5.0:
					count_MS_8 += 1
			else:
				count_G += 1
			if model <= 13:
				if k_star <= 1:
					count_MS_upper += 1
					if m_star >= 5.0:
						count_MS_8_upper += 1
				else:
					count_G_upper += 1
			if model >= 17:
                                if k_star <= 1:
                                        count_MS_lower += 1
                                        if m_star >= 5.0:
                                                count_MS_8_lower += 1
                                else:
                                        count_G_lower += 1
        ##### LOOK THROUGH NS-MS COLLISION #######
        for j in range(len(data2)):
                t = data2[j,0]
                if t > t_lower and t < t_upper:
                        #2190.86035452 122.887084122 0.203042 52.5302 0.0 14.0 11.0
                        k_star = data2[j,4]
                        m_BH = data2[j,3]
                        m_star = data2[j,2]
                        model = data2[j,6]
                        if k_star <= 1:
                                count_MS2 += 1
                        else:
                                count_G2 += 1
                        if model <= 13:
                                if k_star <= 1:
                                        count_MS_upper2 += 1
                                else:
                                        count_G_upper2 += 1
                        if model >= 17:
                                if k_star <= 1:
                                        count_MS_lower2 += 1
                                else:
                                        count_G_lower2 += 1	
        ##### LOOK THROUGH COLLAPSARS #######
        for j in range(len(data3)):
                t = data3[j,0]
                if t > t_lower and t < t_upper:
			model = data3[j,11]
			#9926.28262997 4.36546830797 0.501201 0.0 943982.0 0.118514 30.7634 27.6871 0.0 1.0 -0.0 21.0
                        count_MS3 += 1
                        if model <= 13:
                        	count_MS_upper3 += 1
                        if model >= 17:
				count_MS_lower3 += 1
	
	delta_t = (t_upper - t_lower)*1.e6  # in years

	weight = 1./(1000.*11.)	# number of time draws (100) times number of models (11)
	
	weight_lower = 1./(1000.*5.)
	weight_upper = 1./(1000.*3.)


########## BH-MS collisions ###############
	###### fidicial
	R = 1. * count_MS/delta_t * weight * 4./3.*np.pi * 0.77 * (D_upper**3. - D_lower**3.) / (1+z_mid)	
	scriptR = 1. * count_MS/delta_t * weight * 0.77 * (1.e3)**3.

	R_8 = 1. * count_MS_8/delta_t * weight * 4./3.*np.pi * 0.77 * (D_upper**3. - D_lower**3.) / (1+z_mid)
	scriptR_8 = 1. * count_MS_8/delta_t * weight * 0.77 * (1.e3)**3.

	##### only rv < 0.7
        R_upper = 1. * count_MS_upper/delta_t * weight_upper * 4./3.*np.pi * 2.31 * (D_upper**3. - D_lower**3.) / (1+z_mid)
        scriptR_upper = 1. * count_MS_upper/delta_t * weight_upper * 2.31 * (1.e3)**3.

        R_8_upper = 1. * count_MS_8_upper/delta_t * weight_upper * 4./3.*np.pi * 2.31 * (D_upper**3. - D_lower**3.) / (1+z_mid)
        scriptR_8_upper = 1. * count_MS_8_upper/delta_t * weight_upper * 2.31 * (1.e3)**3.
	
	##### only rv > 2
        R_lower = 1. * count_MS_lower/delta_t * weight_lower * 4./3.*np.pi * 0.32 * (D_upper**3. - D_lower**3.) / (1+z_mid)
        scriptR_lower = 1. * count_MS_lower/delta_t * weight_lower * 0.32 * (1.e3)**3.

        R_8_lower = 1. * count_MS_8_lower/delta_t * weight_lower * 4./3.*np.pi * 0.32 * (D_upper**3. - D_lower**3.) / (1+z_mid)
        scriptR_8_lower = 1. * count_MS_8_lower/delta_t * weight_lower * 0.32 * (1.e3)**3

	### now append the arrays

        z_array.append(z_mid)
        R_array.append(R)
        scriptR_array.append(scriptR)
        R_8_array.append(R_8)
        scriptR_8_array.append(scriptR_8)

        R_array_upper.append(R_upper)
        scriptR_array_upper.append(scriptR_upper)
        R_8_array_upper.append(R_8_upper)
        scriptR_8_array_upper.append(scriptR_8_upper)

        R_array_lower.append(R_lower)
        scriptR_array_lower.append(scriptR_lower)
        R_8_array_lower.append(R_8_lower)
        scriptR_8_array_lower.append(scriptR_8_lower)


######## NS-MS collision ##################

        ###### fidicial
        R2 = 1. * count_MS2/delta_t * weight * 4./3.*np.pi * 0.77 * (D_upper**3. - D_lower**3.) / (1+z_mid)
        scriptR2 = 1. * count_MS2/delta_t * weight * 0.77 * (1.e3)**3.

        ##### only rv < 0.7
        R_upper2 = 1. * count_MS_upper2/delta_t * weight_upper * 4./3.*np.pi * 2.31 * (D_upper**3. - D_lower**3.) / (1+z_mid)
        scriptR_upper2 = 1. * count_MS_upper2/delta_t * weight_upper * 2.31 * (1.e3)**3.

        ##### only rv > 2
        R_lower2 = 1. * count_MS_lower2/delta_t * weight_lower * 4./3.*np.pi * 0.32 * (D_upper**3. - D_lower**3.) / (1+z_mid)
        scriptR_lower2 = 1. * count_MS_lower2/delta_t * weight_lower * 0.32 * (1.e3)**3.

        ### now append the arrays

        R_array2.append(R2)
        scriptR_array2.append(scriptR2)

        R_array_upper2.append(R_upper2)
        scriptR_array_upper2.append(scriptR_upper2)

        R_array_lower2.append(R_lower2)
        scriptR_array_lower2.append(scriptR_lower2)

######## collapsar ##################

        ###### fidicial
        R3 = 1. * count_MS3/delta_t * weight * 4./3.*np.pi * 0.77 * (D_upper**3. - D_lower**3.) / (1+z_mid)
        scriptR3 = 1. * count_MS3/delta_t * weight * 0.77 * (1.e3)**3.
        
        ##### only rv < 0.7
        R_upper3 = 1. * count_MS_upper3/delta_t * weight_upper * 4./3.*np.pi * 2.31 * (D_upper**3. - D_lower**3.) / (1+z_mid)
        scriptR_upper3 = 1. * count_MS_upper3/delta_t * weight_upper * 2.31 * (1.e3)**3.

        ##### only rv > 2
        R_lower3 = 1. * count_MS_lower3/delta_t * weight_lower * 4./3.*np.pi * 0.32 * (D_upper**3. - D_lower**3.) / (1+z_mid)
        scriptR_lower3 = 1. * count_MS_lower3/delta_t * weight_lower * 0.32 * (1.e3)**3.

        ### now append the arrays

        R_array3.append(R3)
        scriptR_array3.append(scriptR3)

        R_array_upper3.append(R_upper3)
        scriptR_array_upper3.append(scriptR_upper3)

        R_array_lower3.append(R_lower3)
        scriptR_array_lower3.append(scriptR_lower3)

##########################################

	print 't =', t_lower, t_upper, 'z =', z_lower, z_upper, 'D =',  D_lower, D_upper, R, scriptR, R_8, scriptR_8


# initialize cumulative BH-MS arrays
R_cum = 0
R_cum_upper = 0
R_cum_lower = 0
R_cum_array_opt = []
R_cum_array_pess = []
R_8_cum = 0
R_8_cum_upper = 0
R_8_cum_lower = 0
R_8_cum_array_opt = []
R_8_cum_array_pess = []
# initialize cumulative NS-MS arrays
R_cum2 = 0
R_cum_upper2 = 0
R_cum_lower2 = 0
R_cum_array_opt2 = []
R_cum_array_pess2 = []
# initialize cumulative collapsar arrays
R_cum3 = 0
R_cum_upper3 = 0
R_cum_lower3 = 0
R_cum_array_opt3 = []
R_cum_array_pess3 = []

for i in range(len(z_array)-1,-1,-1):
	R_cum += R_array[i]
	R_cum_array.append(R_cum)
	R_cum_array_opt.append(R_cum*2.31/0.77)
	R_cum_array_pess.append(R_cum*0.32/0.77)
        R_8_cum += R_8_array[i]
        R_8_cum_array.append(R_8_cum)
        R_8_cum_array_opt.append(R_8_cum*2.31/0.77)
        R_8_cum_array_pess.append(R_8_cum*0.32/0.77)

	R_cum_upper += R_array_upper[i]
        R_cum_array_upper.append(R_cum_upper)
	R_8_cum_upper += R_8_array_upper[i]
        R_8_cum_array_upper.append(R_8_cum_upper)

	R_cum_lower += R_array_lower[i]
        R_cum_array_lower.append(R_cum_lower)
        R_8_cum_lower += R_8_array_lower[i]
        R_8_cum_array_lower.append(R_8_cum_lower)

	print z_array[i], R_array[i], R_cum, R_8_cum

	## do the same for NS-MS
        R_cum2 += R_array2[i]
        R_cum_array2.append(R_cum2)
        R_cum_array_opt2.append(R_cum2*2.31/0.77)
        R_cum_array_pess2.append(R_cum2*0.32/0.77)

        R_cum_upper2 += R_array_upper2[i]
        R_cum_array_upper2.append(R_cum_upper2)

        R_cum_lower2 += R_array_lower2[i]
        R_cum_array_lower2.append(R_cum_lower2)

	## do the same for collapsars
        R_cum3 += R_array3[i]
        R_cum_array3.append(R_cum3)
        R_cum_array_opt3.append(R_cum3*2.31/0.77)
        R_cum_array_pess3.append(R_cum3*0.32/0.77)

        R_cum_upper3 += R_array_upper3[i]
        R_cum_array_upper3.append(R_cum_upper3)

        R_cum_lower3 += R_array_lower3[i]
        R_cum_array_lower3.append(R_cum_lower3)


## now reverse arrays for plotting purposes
R_cum_array = R_cum_array[::-1]
R_cum_array_upper = R_cum_array_upper[::-1]
R_cum_array_lower = R_cum_array_lower[::-1]
R_cum_array_opt = R_cum_array_opt[::-1]
R_cum_array_pess = R_cum_array_pess[::-1]
R_8_cum_array = R_8_cum_array[::-1]
R_8_cum_array_upper = R_8_cum_array_upper[::-1]
R_8_cum_array_lower = R_8_cum_array_lower[::-1]
R_8_cum_array_opt = R_8_cum_array_opt[::-1]
R_8_cum_array_pess = R_8_cum_array_pess[::-1]
#### NS_MS
R_cum_array2 = R_cum_array2[::-1]
R_cum_array_upper2 = R_cum_array_upper2[::-1]
R_cum_array_lower2 = R_cum_array_lower2[::-1]
R_cum_array_opt2 = R_cum_array_opt2[::-1]
R_cum_array_pess2 = R_cum_array_pess2[::-1]
#### collapsar
R_cum_array3 = R_cum_array3[::-1]
R_cum_array_upper3 = R_cum_array_upper3[::-1]
R_cum_array_lower3 = R_cum_array_lower3[::-1]
R_cum_array_opt3 = R_cum_array_opt3[::-1]
R_cum_array_pess3 = R_cum_array_pess3[::-1]


##### now make the plots... 8 panel


### BH-MS collisions
ax1 = plt.subplot(211)
plt.plot(z_array,R_cum_array,lw=2,color='black',ls='--',label=r'$\rm{All\,\,collisions}$')
#plt.plot(z_array,R_cum_array_opt,lw=2,color='black',ls='--')
#plt.plot(z_array,R_cum_array_pess,lw=2,color='black',ls=':')
plt.fill_between(z_array,R_cum_array_pess, R_cum_array_opt,color='black',alpha=0.2)
plt.fill_between(z_array,R_cum_array_lower, R_cum_array_upper,color='black',alpha=0.1)

#plt.title(r'$\rm{BH-MS\,\,collisions}$',fontsize=20)

plt.plot(z_array,R_8_cum_array,lw=2,color='blue',ls='--',label=r'$M_{\star}>5\,M_{\odot}$')
#plt.plot(z_array,R_8_cum_array_opt,lw=2,color='blue',ls='--')
#plt.plot(z_array,R_8_cum_array_pess,lw=2,color='blue',ls=':')

plt.fill_between(z_array,R_8_cum_array_pess, R_8_cum_array_opt,color='blue',alpha=0.2)
plt.fill_between(z_array,R_8_cum_array_lower, R_8_cum_array_upper,color='blue',alpha=0.1)

#plt.scatter(1e7,1e7,color='red',facecolor='none',s=30,label=r'$\it{Swift}-\rm{WP10}$')

#plt.xlabel(r'$\rm{Redshift}\,(\it{z})$',fontsize=20)
plt.ylabel(r'$\rm{Cumulative\,\,rate}\,(\rm{yr}^{-1})$',fontsize=16)
plt.yscale('log')
plt.ylim(.01,1e5)
plt.xlim(0.05,8.3)

plt.setp(ax1.get_xticklabels(), visible=False)

plt.legend(loc=4,scatterpoints=1,prop={'size': 12})

ax2 = plt.subplot(212)

yhat = savgol_filter(scriptR_array, 51, 3)
yhat_upper = yhat*3.
yhat_lower = yhat/2.
yhat_upper_upper = savgol_filter(scriptR_array_upper, 51, 3)
yhat_lower_lower = savgol_filter(scriptR_array_lower, 51, 3)

yhat_8 = savgol_filter(scriptR_8_array, 51, 3)
yhat_8_upper = yhat_8*3.
yhat_8_lower = yhat_8/2.
yhat_8_upper_upper = savgol_filter(scriptR_8_array_upper, 51, 3)
yhat_8_lower_lower = savgol_filter(scriptR_8_array_lower, 51, 3)


#plt.plot(z_array,scriptR_array,lw=2,color='black')
plt.plot(z_array,yhat,lw=2,color='black',ls='--')
plt.fill_between(z_array,yhat_lower, yhat_upper,color='black',alpha=0.2)
plt.fill_between(z_array,yhat_lower_lower, yhat_upper_upper,color='black',alpha=0.1)

plt.plot(z_array,yhat_8,lw=2,color='blue',ls='--')
plt.fill_between(z_array,yhat_8_lower, yhat_8_upper,color='blue',alpha=0.2)
plt.fill_between(z_array,yhat_8_lower_lower, yhat_8_upper_upper,color='blue',alpha=0.1)

#plt.plot(z_array,scriptR_8_array,lw=2,color='gray')

y_lower = density[:,2] - density[:,1]
y_upper = density[:,3] - density[:,2]
#plt.errorbar(density[:,0], density[:,2], yerr = [y_lower, y_upper],ecolor='red',fmt='none',zorder=0)
#plt.scatter(density[:,0], density[:,2],marker='o',facecolor='white',color='red',s=30,zorder=1)

#plt.plot(density_L14[:,0], density_L14[:,1], color='darkgreen',lw=2)

plt.xlabel(r'$\rm{Redshift}\,(\it{z})$',fontsize=18)
#plt.ylabel(r'$R\,(\rm{yr}^{-1})$',fontsize=20)
plt.ylabel(r'$\rm{Comoving\,\,rate}\,(\rm{Gpc}^{-3}\,\rm{yr}^{-1})$',fontsize=16)
plt.xlim(0.05,8.3)
plt.yscale('log',nonposy='clip')
plt.ylim(.01,400)

ax1.xaxis.set_tick_params(width=1,length=6,which='major')
ax1.xaxis.set_tick_params(width=1,length=3.5,which='minor')
ax1.yaxis.set_tick_params(width=1,length=6, which='major')
ax1.yaxis.set_tick_params(width=1,length=3.5, which='minor')

for tick in ax1.xaxis.get_major_ticks():
        tick.label.set_fontsize(14)
for tick in ax1.yaxis.get_major_ticks():
        tick.label.set_fontsize(14)

ax2.xaxis.set_tick_params(width=1,length=6,which='major')
ax2.xaxis.set_tick_params(width=1,length=3.5,which='minor')
ax2.yaxis.set_tick_params(width=1,length=6, which='major')
ax2.yaxis.set_tick_params(width=1,length=3.5, which='minor')

for tick in ax2.xaxis.get_major_ticks():
        tick.label.set_fontsize(14)
for tick in ax2.yaxis.get_major_ticks():
        tick.label.set_fontsize(14)

plt.subplots_adjust(wspace=0, hspace=0)
plt.show()
