import numpy as np
import conversions
import finder
import os,sys
import subprocess
import gzip
import scripts
import history_cmc

# Define bins: 108 total bins (3*3*3*4)
mass_array = [0,5e4,1e5,2e5,7e5]
rcrh_array = [0,0.4,0.8,1.3]
rgc_array = [0,5,14,100]
z_array = [0,0.0005,0.005,0.03]
#rgc_array = [0,100]
#z_array = [0,0.03]

##Read Data
data_harris = np.genfromtxt('/projects/b1095/kylekremer/python_code/CMC_Grid_March2019/harris.dat')
R_sun_obs = data_harris[:,2]; M_bolo = data_harris[:,4]-0.107
rc_obs = data_harris[:,7]*(R_sun_obs*1000)/3437.75; rh_obs = data_harris[:,8]*(R_sun_obs*1000)/3437.75
rgc_obs = data_harris[:,3]; z_obs = conversions.metallicity(data_harris[:,5],'fe/htoz')
mass_obs = 2.0*10**(0.4*(4.74-M_bolo)); rcrh_obs = rc_obs/rh_obs

data_model = np.genfromtxt('/projects/b1095/kylekremer/python_code/CMC_Grid_March2019/cluster_params_12gyr.dat')
mass_mod = data_model[:,6]; rc_mod = data_model[:,7]; rh_mod = data_model[:,8]
rcrh_mod = rc_mod/rh_mod
rgc_mod = data_model[:,4]; z_mod = data_model[:,5]


obs_total = 60.   ##Why 60?

obs_bin_count=[]
mod_bin_count=[]
for i in range(len(mass_array)-1):
    for j in range(len(rcrh_array)-1):
        for k in range(len(rgc_array)-1):
            for l in range(len(z_array)-1):
                obscount=0; modcount=0
                for x in range(len(mass_obs)):
                    if mass_array[i]<=mass_obs[x]<mass_array[i+1] and rcrh_array[j]<=rcrh_obs[x]<rcrh_array[j+1] and rgc_array[k]<=rgc_obs<rgc_array[k+1] and z_array[l]<=z_obs[x]<z_array[l+1]:
                        obscount+=1

                obs_bin_count.append(obscount)


