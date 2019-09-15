import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import conversions
import finder
import os,sys
import subprocess
import gzip
import scripts
import history_cmc

# 108 total bins (3*3*3*4)

f = open('model_weights.dat','w')

mass_array = [0,5e4,1e5,2e5,7e5]
rcrh_array = [0,0.4,0.8,1.3]
rgc_array = [0,5,14,100]
z_array = [0,0.0005,0.005,0.03]
#rgc_array = [0,100]
#z_array = [0,0.03]

data = np.genfromtxt('harris.dat')
data2 = np.genfromtxt('cluster_params_12gyr.dat')
obs_total = 60.

obs_count_total = 0
weight_total = 0
for i in range(len(mass_array)-1):  # Total mass bins
	m_low = mass_array[i]
	m_high = mass_array[i+1]

	for j in range(len(rcrh_array)-1): # rc/rh bins
		rcrh_low = rcrh_array[j]
		rcrh_high = rcrh_array[j+1]

		for k in range(len(rgc_array)-1): # Rgc bins
			rgc_low = rgc_array[k]
			rgc_high = rgc_array[k+1]			

			for l in range(len(z_array)-1): # Metallicity bins
				z_low = z_array[l]
				z_high = z_array[l+1]

				obs_count = 0
				for x in range(len(data)):
				#0)n 1)ID 2)Rsun(kpc) 3)Rgc(kpc) 4)MVt 5)Fe/H 6)c 7)Rc(arcmin) 8)Rh(arcmin) 9)muV(mag/arcsec^2) #10)core-collapsed?
					R_sun_obs = data[x,2]
					M_bolo = data[x,4]-0.107
					mass_obs = 2.0*10**(0.4*(4.74-M_bolo))
					rc_obs = data[x,7]*(R_sun_obs*1000)/3437.75
					rh_obs = data[x,8]*(R_sun_obs*1000)/3437.75
					rcrh_obs = rc_obs/rh_obs
					rgc_obs = data[x,3]
					z_obs = conversions.metallicity(data[x,5],'fe/htoz')
					if m_low <= mass_obs <= m_high and rcrh_low <= rcrh_obs <= rcrh_high and rgc_low <= rgc_obs <= rgc_high and z_low <= z_obs <= z_high:
						obs_count += 1
						#obs_count_total += 1

				model_count = 0
                for y in range(len(data2)):
				#25 11923.8755321 2.0 0.5 20.0 0.0002 23240.0406083 0.383022426497 2.08378323199
                    mass_mod = data2[y,6]
                    rc_mod = data2[y,7]
                    rh_mod = data2[y,8]
                    rcrh_mod = rc_mod/rh_mod
                    rgc_mod = data2[y,4]
                    z_mod = data2[y,5]
                    if m_low <= mass_mod <= m_high and rcrh_low <= rcrh_mod <= rcrh_high and rgc_low <= rgc_mod <= rgc_high and z_low <= z_mod <= z_high:
                            model_count += 1		

				if model_count > 0:	
					weight = np.float(obs_count)/np.float(model_count)/obs_total
					obs_count_total += obs_count
					#weight_total += (weight*model_count)
				else:
					weight = 'nan'
				print i,j,k,l, 'obs:', obs_count, 'model:', model_count, 'weight:', weight


				### Now make the new cluster params file with weight added
				for y in range(len(data2)):
                #25 11923.8755321 2.0 0.5 20.0 0.0002 23240.0406083 0.383022426497 2.08378323199
                    mass_mod = data2[y,6]
					rv_mod = data2[y,3]
					N_mod = data2[y,2]
                    rc_mod = data2[y,7]
                    rh_mod = data2[y,8]
                    rcrh_mod = rc_mod/rh_mod
                    rgc_mod = data2[y,4]
                    z_mod = data2[y,5]
                                        if m_low <= mass_mod <= m_high and rcrh_low <= rcrh_mod <= rcrh_high and rgc_low <= rgc_mod <= rgc_high and z_low <= z_mod <= z_high:
                                                print>>f, data2[y,0], N_mod, rv_mod, rgc_mod, z_mod, weight, mass_mod, rc_mod, rh_mod, rgc_mod, z_mod
						weight_total += weight

### Now re-sort the weights file and also assign weights for dissolved clusters to be zero
f.close()
data = np.genfromtxt('model_weights.dat')
f = open('model_weights.dat','w')
print>>f,'#0)model 1)weight 2)mass 3)rc 4)rh 5)rgc 6)Z'
weight_total = 0
for i in range(1,145):
	flag = 0
	for j in range(len(data)):
		if data[j,0] == i:
			print>>f, i, data[j,5], data[j,6], data[j,7], data[j,8], data[j,9], data[j,10], data[j,1], data[j,2], data[j,3], data[j,4]
			weight_total += data[j,5]
			flag = 1
	if flag == 0:
		print>>f, i, 0.0, 0, 0, 0, 0, 0 ,0 ,0 ,0 ,0

print obs_count_total		
print weight_total
