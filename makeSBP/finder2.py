import numpy as np
import constants
import math
import conversions
from scipy import integrate
import matplotlib.pyplot as plt
import os
import sys
import glob
import gzip
import scripts

def find_num_binaries_half_light(filestring, r_h, snapno = '0010', num = 587943):
	"""finds number of binaries outside of half-light radius. Inputs are simulation folder, r_h (pc)"""
	data = np.loadtxt('/projects/b1011/kyle/cmc/branches/cmc-mpi/m22_project/rundir/'+filestring+'/initial.snap'+snapno+'.2Dproj.dat') 
	count = 0.0
	count_bin = 0.0
	for i in range(0,num):
		if data[i,0] > r_h:
			count = count + 1  # Counts number of objects outside r_h
			if data[i,2] == 1: # Asks if it is a binary
				count_bin = count_bin + 1	# Counts num of binaries outside r_H
	bin_frac = count_bin/count
	return count_bin, count, bin_frac

def find_rh(cluster, filestring, snapno):
	"""finds half-light radius"""
	data = np.loadtxt('/projects/b1011/kyle/cmc/branches/cmc-mpi/m22_project/rundir/'+cluster+'/'+filestring+'/initial.snap'+snapno+'.2Dproj.dat')	
	L = 0
	L_half = 0
	for i in range(0,len(data)):
    		L = L + data[i,1]
	for i in range(0,len(data)):
		L_half = L_half + data[i,1]
		if L_half > L/2:
			return data[i-1,0], L_half, L
			break


def find_bin_half_light_ALT(filestring, r_h, snapno = '0010', num = 587943):
	"""finds number of binaries outside r_h by finding al binaries with q>0.5"""
	data = np.loadtxt('/projects/b1011/kyle/cmc/branches/cmc-mpi/m22_project/rundir/'+filestring+'/initial.snap'+snapno+'.2Dproj.dat')
	count = 0.0
	count_bin = 0.0
	for i in range(0,num):
		if data[i,0] > r_h:
			if data[i,2] == 0:	# if single star
				if data[i,3] < 10:	# if non-compact object
					count = count + 1 # counts number of non-compact object single stars outside r_h
			if data[i,2] == 1:
				if data[i,5] < 10 and data[i,6] < 10:  # asks if both binary components are non-compact objects
					if data[i,10] > data[i,11]:
						q = data[i,11]/data[i,10]
					else:
						q = data[i,10]/data[i,11]
					if q > 0.5:
						count_bin = count_bin + 1
	bin_frac = count_bin/count*2
	return count_bin, count, bin_frac

def find_MS_turnoff(t):
	"""given the time in Myr it finds the MS turn-off mass in Solar masses.  Very simple now.  Need to make the MS lifetime formula better. """
	t_yr = t*10**6
	lm = (9.921 - np.log10(t_yr))/3.6648
	m = 10**lm
	return(m)	

def find_t_ms(z, m):
	eta = np.log10(z/0.02)
	a1 = 1.593890e3+2.053038e3*eta+1.231226e3*eta**2.+2.327785e2*eta**3.
	a2 = 2.706708e3+ 1.483131e3*eta+ 5.772723e2*eta**2.+ 7.411230e1*eta**3.
	a3 = 1.466143e2 - 1.048442e2*eta - 6.795374e1*eta**2. - 1.391127e1*eta**3.
	a4 = 4.141960e-2 + 4.564888e-2*eta + 2.958542e-2*eta**2 + 5.571483e-3*eta**3.
	a5 = 3.426349e-1
	a6 = 1.949814e1 + 1.758178*eta - 6.008212*eta**2. - 4.470533*eta**3.
	a7 = 4.903830
	a8 = 5.212154e-2 + 3.166411e-2*eta - 2.750074e-3*eta**2. - 2.271549e-3*eta**3.
	a9 = 1.312179 - 3.294936e-1*eta + 9.231860e-2*eta**2. + 2.610989e-2*eta**3.
	a10 = 8.073972e-1


	m_hook = 1.0185 + 0.16015*eta + 0.0892*eta**2.
	m_HeF = 1.995 + 0.25*eta + 0.087*eta**2.
	m_FGB = 13.048*(z/0.02)**0.06/(1+0.0012*(0.02/z)**1.27)

	t_BGB = (a1+a2*m**4.+a3*m**5.5+m**7.)/(a4*m**2.+a5*m**7.)
	x = max([0.95,min([0.95-0.03*(eta+0.30103)]),0.99])
	mu = max(0.5, 1.0-0.01* max(a6/(m**a7) , a8+a9/m**a10))
	t_hook = mu*t_BGB

	t_MS = max(t_hook, x*t_BGB)

	return (t_MS)

def find_MS_TO(t, z, mguess):
	tguess = find_t_ms(z, mguess)
	#print mguess, tguess, (t-tguess)/t

	while (t-tguess)/t > 0.00005:
		mguess -= 0.00001
		tguess = find_t_ms(z, mguess)
		#print mguess, tguess, (t-tguess)/t
	mto = mguess
	return mto

def find_num_BS(cluster,filestring, snapno,  m_to = 0.828):
	"""find number of blue stragglers. Outputs total number, total singles, total in binaries"""
	data = np.genfromtxt('/projects/b1011/kyle/cmc/branches/cmc-mpi/m22_project/rundir/'+cluster+'/'+filestring+'/initial.snap'+snapno+'.dat.gz')
	count = 0
	count_bin = 0
	for i in range(0,len(data)):
		if data[i,7] == 0:	# If single star
			if data[i,14] == 0 or data[i,14] == 1:	# k-type
				if data[i,1] > 1.05*m_to:
					count = count + 1
		else:			# If binary star
			if data[i,17] == 0 or data[i,17] == 1:
				if data[i,8] > 1.05*m_to:
					count_bin = count_bin + 1
			if data[i,18] == 0 or data[i,18] == 1:
				if data[i,9] > 1.05*m_to:
                                        count_bin = count_bin + 1
	total = count + count_bin
	return total, count, count_bin

def find_BHs(cluster,filestring,snapno):
	"""finds all BH (single and in binaries) retained and ejected. Returns num single BHs retained, num BBHs retained, num single BBHS retained, num single BHs ejected, BBHs ejected, single BBHs ejected. num is number of lines in snap file. num2 is number of lines in .esc file."""
	f = open('BH_'+filestring+'.dat','a')
	data = np.genfromtxt('/projects/b1011/kyle/cmc/branches/cmc-mpi/m22_project/rundir/'+cluster+'/'+filestring+'/initial.snap'+snapno+'.dat.gz')
	data2 = np.genfromtxt('/projects/b1011/kyle/cmc/branches/cmc-mpi/m22_project/rundir/'+cluster+'/'+filestring+'/initial.esc.dat')
	count = 0
	count_bin_single = 0
	count_bin_double = 0
        count2 = 0
        count2_bin_single = 0
        count2_bin_double = 0
	xray_cand = 0
	for i in range(0,len(data)):
		if data[i,7] == 0:
			if data[i,14] == 14:
				count = count + 1
		if data[i,7] == 1:
			if data[i,17] == 14 and data[i,18] == 14:
				count_bin_double = count_bin_double + 1
			if data[i,17] == 14 or data[i,18] == 14:
				count_bin_single = count_bin_single + 1
				if data[i,17] < 10 or data[i,18] < 10:
					xray_cand = xray_cand + 1
	count_bin_single = count_bin_single - count_bin_double
	for i in range(0,len(data2)):
		if data2[i,14] == 0:
			if data2[i,21] == 14:
				count2 = count2 + 1
		if data2[i,14] == 1:
			if data2[i,22] == 14 and data2[i,23] == 14:
				count2_bin_double = count2_bin_double + 1
			if data2[i,22] == 14 or data2[i,23] == 14:
				count2_bin_single = count2_bin_single + 1
	count2_bin_single = count2_bin_single - count2_bin_double
	print>>f, snapno, count, count_bin_double, count_bin_single, xray_cand
	#return count, count_bin_double, count_bin_single, count2, count2_bin_double, count2_bin_single, "num single BHs retained, num double BH binaries retained, num single BH binaries retained, num single BHs ejected, num double BH binaries ejected, num single BH binaries ejected"

def find_half_mass_radius(cluster, filestring, snapno):
	data = np.loadtxt('/projects/b1011/kyle/cmc/branches/cmc-mpi/m22_project/rundir/'+cluster+'/'+filestring+'/initial.snap'+snapno+'.2Dproj.dat')
	total_mass = 0
	total_mass_temp = 0
	for i in range(0,len(data)):
		total_mass = total_mass + data[i,9]
	for j in range(0,len(data)):
		total_mass_temp = total_mass_temp + data[j,9]
		if total_mass_temp > 0.5*total_mass:
			r_h = data[j-1,0]
			break
	return total_mass, r_h

def find_binfrac(cluster, filestring, snapno, rc, rh, flag):
	if flag == 0: # Want to calculate for MS stars
		data = np.loadtxt('/projects/b1011/kyle/cmc/branches/cmc-mpi/m22_project/rundir/'+cluster+'/'+filestring+'/initial.snap'+snapno+'.2Dproj_I.dat')
		Mag_lower = 20.3	#set lower and upper bounds on magnitude cut
		Mag_upper = 23.0
	if flag == 1: # Want to calculate for RGB stars
		data = np.loadtxt('/projects/b1011/kyle/cmc/branches/cmc-mpi/m22_project/rundir/'+cluster+'/'+filestring+'/initial.snap'+snapno+'.2Dproj.dat')
		Mag_lower = 12.017
                Mag_upper = 15.017
	total_rc = 0.0
	total_int = 0.0
	total_rh = 0.0
	bin_total_rc = 0.0
	bin_total_int = 0.0
	bin_total_rh = 0.0
	count = 0.0
	primary = 0.0

	for i in range(0,len(data)):
		if data[i,2] == 1:
			count = count + 1.0
			if data[i,10] > data[i,11]:          #this determines mass ratio if object is a binary
				q = data[i,11]/data[i,10]
				primary = data[i,10]
			if data[i,11] >= data[1,10]:
				q = data[i,10]/data[i,11]
				primary = data[i,11]
		Mag_I = -2.5*np.log10(data[i,1]) + 13.217     #converts to apparent I-band magnitude
		#if Mag_I >= 0 and Mag_I <= 100:
		if Mag_I >= Mag_lower and Mag_I <= Mag_upper: # Determine if in I-band range
			#print primary, data[i,3], data[i,5], data[i,6]
			if data[i,0] < rc:
				if data[i,2] == 0:
					total_rc = total_rc + 1
				if data[i,2] == 1:
					if q >= 0:
						total_rc = total_rc + 1
						bin_total_rc = bin_total_rc + 1	
			if data[i,0] < rh and data[i,0] > rc:
				if data[i,2] == 0:
                                        total_int = total_int + 1
				if data[i,2] == 1:
					if q >= 0:
						total_int = total_int + 1
						bin_total_int = bin_total_int + 1
			if data[i,0] > rh:
				if data[i,2] == 0:
					total_rh = total_rh + 1
				if data[i,2] == 1:
					if q >= 0:
						total_rh = total_rh + 1
						bin_total_rh = bin_total_rh + 1
	binfrac_rc = bin_total_rc/total_rc
	binfrac_int = bin_total_int/total_int
	binfrac_rh = bin_total_rh/total_rh
	total = total_rc+total_int+total_rh
	return total, binfrac_rc, binfrac_int, binfrac_rh, count		

def convert_to_3d(r, vr, vt):
	#costheta = np.random.uniform(-1, 1)
	#sintheta = (1-costheta**2.)**0.5

	if np.shape(r)==():
		#print 'came here'
		sintheta = np.random.uniform(low=-1., high=1.)
		phi = np.random.uniform(low=0., high=2.*np.pi)
		anglev = np.random.uniform(low=0., high=2.*np.pi)
	else:
		#print 'came here too'
		r = np.array(r)
		vr = np.array(vr)
		vt = np.array(vt)
	
		sintheta = np.random.uniform(low=-1., high=1., size=len(r))
		phi = np.random.uniform(low=0., high=2.*np.pi, size=len(r))
		anglev = np.random.uniform(low=0., high=2.*np.pi, size=len(r))
		
	costheta = (1-sintheta**2.)**0.5
	#costheta = (1-sintheta**2.)**0.5
	#phi = np.random.uniform(low=0., high=2.*np.pi, size=len(r))
	
	rz = r*sintheta
	rx = r*costheta*np.cos(phi)
	ry = r*costheta*np.sin(phi)
	
	#anglev = np.random.uniform(low=0., high=2.*np.pi, size=len(r))
	magv = (vr*vr + vt*vt)**0.5	
	thetadot = np.cos(anglev) * vt/r
	phidot = np.sin(anglev)*vt/(r*costheta)
	
	#vx = vr * np.sin(np.arccos(costheta)) * np.cos(phi) + r * thetadot * costheta * np.cos(phi) - r * phidot * np.sin(np.arccos(costheta)) * np.sin(phi)
	vx = vr * costheta * np.cos(phi) - r * phidot * costheta * np.sin(phi) - r * thetadot * sintheta * np.cos(phi)
	#vy = vr * np.sin(np.arccos(costheta)) * np.sin(phi) + r * thetadot * costheta * np.sin(phi) + r * phidot * np.sin(np.arccos(costheta)) * np.cos(phi)
	vy = vr * costheta * np.sin(phi) + r * phidot * costheta * np.cos(phi) - r * thetadot * sintheta * np.sin(phi)
	#vz = vr * costheta - r * thetadot * np.sin(np.arccos(costheta))
	vz = vr * sintheta + r * thetadot * costheta


	r3d = np.array([rx, ry, rz])
	v3d = np.array([vx, vy, vz])

	return r3d, v3d

def velocity_dispersion(path,string,snapno,ALL=1,Bin_no=25):
	import scripts
	import random
	units=scripts.read_units(path+string)
	km = units[0]['l_cgs']*1.0e-5
	time = units[0]['nbt_cgs']
	f = open('RV_model.dat','w')
	f1 = open('RV_model_sigma.dat','a')
	if ALL == 1:
		f4 = open(path+string+'.snap'+snapno+'.vel_dispersion_tight.dat','w')
	if ALL == 0:
		f4 = open(path+string+'.snap'+snapno+'.vel_dispersion_giants_25.dat','w')
	f2 = open('RV_obs_sigma.dat','w')
	f3 = open('giants.dat','w')
	print>>f4,'#0)r2D(pc) 1)sigma_v(1D; km/s) 2) sigma_err (km/s)'
####################
	f55 = gzip.open(path+string+'.snap'+snapno+'.dat.gz','r')
	lines55 = f55.readlines()
	#data = np.genfromtxt(path+string+'.snap'+snapno+'.dat.gz')
#####################
###################
	Vr = []; Vpm_r = []; Vpm_t = []; Vpm = []
	R = []
	bin_count = 0
	count = 0
	bin_array = []
	for i in range(2,len(lines55)):
		line55 = lines55[i]
		data = line55.split(' ')
		for k in range(0,21):
			data[k] = np.float(data[k])			
		bin = 0
		if ALL == 1:
			if data[1] <= 1e6:	# Looks at all stars, not a cut.
			#if data[i,3] <= 1.8 and data[i,3] >= -1.2:  # asks if in V-band mag range given by Caretta et al. 2009.
				r_xyz = [data[15], data[16], data[17]]  #x,y,z components of r in pc
				v_xyz = [data[18], data[19], data[20]]  #x,y,z components of velocity in km/s
				count = count + 1
				if data[2] == 1:
					bin_count = bin_count + 1
					bin = 1
				bin_array.append(bin)
				r = np.sqrt(r_xyz[0]**2. + r_xyz[1]**2.)
				#R_temp = conversions.pc_to_arcsec(r,d_heliocentric)  #Converts radius in pc to arcsec
				R.append(r)
				Vr.append(data[20])		# Just use the z component as your LOS direction
		if ALL == 0:
			v_r = data[3]*km/time # converts from code units to km/s
			v_t = data[4]*km/time
			r_km = data[2]*km  #convert r from code units to km
			## EXCLUDE BINARIES HERE ######
			#if data[i,7] == 1:
			#	if data[i,17] >= 2 and data[i,17] <= 9:
			#		for k in range(0,70):
			#			r,v = convert_to_3d(r_km, v_r, v_t)
			#			Vr.append(v[0])    # Use this if you want 1d vel dispersion
			#			r = np.sqrt(r[1]**2. + r[2]**2.)/3.086e13  # convert r from km to pc
			#			R.append(r)
			#			count = count + 1
				
			#	if data[i,18] >= 2 and data[i,18] <= 9:
                        #                for k in range(0,70):
                        #                       r,v = convert_to_3d(r_km, v_r, v_t)
                        #                        Vr.append(v[0])    # Use this if you want 1d vel dispersion
                        #                        r = np.sqrt(r[1]**2. + r[2]**2.)/3.086e13  # convert r from km to pc
                        #                        R.append(r)
                        #                        count = count + 1			
			if data[7] != 1:
				if data[14] >= 2 and data[14] <= 9:
                    for k in range(0,1):
					#for k in range(0,25):
                        r,v = convert_to_3d(r_km, v_r, v_t)
                        Vr.append(v[0])    # Use this if you want 1d vel dispersion
                        Vpm_r.append(v[1]); Vpm_t.append(v[2])
                        Vpm.append(np.sqrt(v[1]**2 + v[2]**2))
                        r = np.sqrt(r[1]**2. + r[2]**2.)/3.086e13  # convert r from km to pc
                        R.append(r)
                        count = count + 1
################################################
	mean_model = np.mean(Vr[:]) #Find global mean of all model RVs
	mean_model_pmr = np.mean(Vpm_r[:])
	mean_model_pmt = np.mean(Vpm_t[:])
	mean_model_pm = np.mean(Vpm[:])
	print 'mean=',mean_model
	array = np.zeros((len(Vr),5))
	for i in range(0,len(Vr)):
		array[i,0] = R[i]
		array[i,1] = Vr[i]
		array[i,2] = Vpm_r[i]
		array[i,3] = Vpm_t[i]
		array[i,4] = Vpm[i]
	array = array[array[:,0].argsort()]  # sorts each measurement in order of radial position
	for i in range(0,len(Vr)):
		print>>f, array[i,0], array[i,1], array[i,2], array[i,3], array[i,4] #Print each radius/LOS velocity pair in order of radial position	

#####################################	
	sigma_vel_array = []
	sigma_pmr_array = []
	sigma_pmt_array = []
	sigma_pm_array = []
	R_array = []
	#mean_array = []
	#flag = 0
	sum_vel = 0
	sum_pmr = 0
	sum_pmt = 0
	sum_pm = 0
	print>>f1, snapno,
        ## Define each bin as having 2000 stars. Every 2000 stars, start new bin

	vel_array = []  # Makes an array with velocites of all stars within each bin
	pmr_array = []; pmt_array = []; pm_array = []
	r_array = []
	#bin_count = 0
	total_count = 0
	for j in range(0,len(array)):   
		if total_count <= Bin_no:
		#if total_count <= 500:
			vel_array.append(array[j,1])
			pmr_array.append(array[j,2])
			pmt_array.append(array[j,3])
			pm_array.append(array[j,4])
			r_array.append(array[j,0])
			total_count = total_count + 1
		else:
			count = 0.
			for k in range(0,len(vel_array)):
				r = np.mean(r_array)
				sum_vel = sum_vel + (vel_array[k]-mean_model)**2.
				sum_pmr = sum_pmr + (pmr_array[k]-mean_model_pmr)**2.
				sum_pmt = sum_pmt + (pmt_array[k]-mean_model_pmt)**2.
				sum_pm = sum_pm + (pm_array[k]-mean_model_pm)**2.
				count = count + 1.
			sigma_vel = np.sqrt(sum_vel/count)
			error_vel = np.sqrt(sigma_vel**2.0/(2.*count))
            sigma_pmr = np.sqrt(sum_pmr/count)
			error_pmr = np.sqrt(sigma_pmr**2.0/(2.*count))
			sigma_pmt = np.sqrt(sum_pmt/count)
			error_pmt = np.sqrt(sigma_pmt**2.0/(2.*count))
			sigma_pm = np.sqrt(sum_pm/count)
			error_pm = np.sqrt(sigma_pm**2.0/(2.*count))

			print>>f1, r, sigma, error, sigma_pmr, error_pmr, sigma_pmt, error_pmt, sigma_pm, error_pm
			print>>f4, r, sigma, error, sigma_pmr, error_pmr, sigma_pmt, error_pmt, sigma_pm, error_pm
			#print r, 'sigma=',sigma, 'error=',error, 'count=',count, 'N_binaries=',bin_count,'N_stars=',total_count,'binary fraction=',float(bin_count)/float(total_count)
			sigma_vel_array.append(sigma_vel)
			sigma_pmr_array.append(sigma_pmr)
			sigma_pmt_array.append(sigma_pmt)
			sigma_pm_array.append(sigma_pm)
			R_array.append(np.mean(r_array))
			sum_vel = 0; sum_pmr = 0; sum_pmt = 0; sum_pm = 0
			#flag = 0
			#bin_count = 0
			total_count = 0
			vel_array = []; pmr_array = []; pmt_array = []; pm_array = []
			r_array = []
	print>>f1,' '

def velocity_dispersion_type(path,string, snapno,ALL=0,type=14):
        import scripts
        import random
        units=scripts.read_units(path+string)
        km = units[0]['l_cgs']*1.0e-5
        time = units[0]['nbt_cgs']
        f = open('RV_model.dat','w')
        f1 = open('RV_model_sigma.dat','a')
	if type == 14:
		type_label = 'BH'
	if type == 13:
                type_label = 'NS'
        if ALL == 1:
                f4 = open(path+string+'.snap'+snapno+'.vel_dispersion_tight.dat','w')
        if ALL == 0:
                f4 = open(path+string+'.snap'+snapno+'.vel_dispersion_'+type_label+'.dat','w')
        f2 = open('RV_obs_sigma.dat','w')
        f3 = open('giants.dat','w')
        print>>f4,'#0)r2D(pc) 1)sigma_v(1D; km/s) 2) sigma_err (km/s)'
####################
        f55 = gzip.open(path+string+'.snap'+snapno+'.dat.gz','r')
        lines55 = f55.readlines()
        #data = np.genfromtxt(path+string+'.snap'+snapno+'.dat.gz')
#####################
###################
        Vr = []
        R = []
        bin_count = 0
        count = 0
        bin_array = []
        for i in range(2,len(lines55)):
                line55 = lines55[i]
                data = line55.split(' ')
                for k in range(0,21):
                        data[k] = np.float(data[k])
                bin = 0
                if ALL == 1:
			v_r = data[3]*km/time # converts from code units to km/s
                        v_t = data[4]*km/time
                        r_km = data[2]*km  #convert r from code units to km

                        if data[7] != 1:
				for k in range(0,1):
					r,v = convert_to_3d(r_km, v_r, v_t)
					Vr.append(v[0])    # Use this if you want 1d vel dispersion
					r = np.sqrt(r[1]**2. + r[2]**2.)/3.086e13  # convert r from km to pc
					R.append(r)
					count = count + 1

                if ALL == 0:
                        v_r = data[3]*km/time # converts from code units to km/s
                        v_t = data[4]*km/time
                        r_km = data[2]*km  #convert r from code units to km

                        if data[7] != 1:
                                if data[14] == type:
                                        for k in range(0,1):
                                        #for k in range(0,25):
                                                r,v = convert_to_3d(r_km, v_r, v_t)
                                                Vr.append(v[0])    # Use this if you want 1d vel dispersion
                                                r = np.sqrt(r[1]**2. + r[2]**2.)/3.086e13  # convert r from km to pc
                                                R.append(r)
                                                count = count + 1
################################################
        mean_model = np.mean(Vr[:]) #Find global mean of all model RVs
        print 'mean=',mean_model, len(Vr)
        array = np.zeros((len(Vr),2))
        for i in range(0,len(Vr)):
                array[i,0] = R[i]
                array[i,1] = Vr[i]
        array = array[array[:,0].argsort()]  # sorts each measurement in order of radial position
        for i in range(0,len(Vr)):
                print>>f, array[i,0], array[i,1] #Print each radius/LOS velocity pair in order of radial position       

#####################################   
        sigma_array = []
        R_array = []
        R_array2 = []
        mean_array = []
        mean_array2 = []
        sigma_array2 = []
        flag = 0
        sum = 0
        print>>f1, snapno,
        ## Define each bin as having 2000 stars. Every 2000 stars, start new bin

        vel_array = []  # Makes an array with velocites of all stars within each bin
        r_array = []
        bin_count = 0
        total_count = 0
        for j in range(0,len(array)):
		if ALL == 0:
			count_max = 70
		if ALL == 1:
                        count_max = 1000
                if total_count <= count_max:
                        vel_array.append(array[j,1])
                        r_array.append(array[j,0])
                        total_count = total_count + 1
                else:
                        count = 0.
                        for k in range(0,len(vel_array)):
                                r = np.mean(r_array)
                                sum = sum + (vel_array[k]-mean_model)**2.
                                count = count + 1.
                        sigma = np.sqrt(sum/count)
                        error = np.sqrt(sigma**2.0/(2.*count))
                        print>>f1, r, sigma, error,
                        print>>f4, r, sigma, error
                        #print r, 'sigma=',sigma, 'error=',error, 'count=',count, 'N_binaries=',bin_count,'N_stars=',total_count,'binary fraction=',float(bin_count)/float(total_count)
                        sigma_array.append(sigma)
                        R_array.append(np.mean(r_array))
                        sum = 0
                        flag = 0
                        bin_count = 0
                        total_count = 0
                        vel_array = []
                        r_array = []
        print>>f1,' '

	 	
def angular_distance(RA_1,RA_2,RA_3,Dec_1,Dec_2,Dec_3):
	RA_m10 = 16.0+57.0/60.0+9.05/3600.0
	Dec_m10 = -4.0 + 6.0/60.0 + 1.1/3600.0
	RA_star = RA_1+RA_2/60.0+RA_3/3600.0
	Dec_star = Dec_1 + Dec_2/60.0+Dec_3/3600.0
 	#print RA_m10, Dec_m10, RA_star, Dec_star	
	
	RA_m10 = np.pi/12.0*RA_m10  #converts hours to radians
	RA_star = np.pi/12.0*RA_star
	Dec_m10 = np.pi/180.0*Dec_m10   #converts degrees to radians
	Dec_star = np.pi/180.0*Dec_star

	cos_angdistance = np.cos(np.pi/2.0 - Dec_m10)*np.cos(np.pi/2.0 - Dec_star) + np.sin(np.pi/2.0 - Dec_m10)*np.sin(np.pi/2.0 - Dec_star)*np.cos(RA_m10-RA_star)
	#print cos_angdistance
	angdistance = 180.0/np.pi*np.arccos(cos_angdistance)
	arcsec = angdistance*3600.0
	#print arcsec
	return arcsec	
	
def hard_soft(cluster,filestring,snapno):	
        import scripts
        units=scripts.read_units('/projects/b1011/kyle/cmc/branches/cmc-mpi/m22_project/rundir/'+cluster+'/'+filestring+'/initial')
        data2 = np.loadtxt('RV_model_sigma_'+filestring+'.dat')
	km = units[0]['l_cgs']*1.0e-5
	pc = units[0]['l_pc']
        time = units[0]['nbt_cgs']

	########## make array with velocity dispersions for each bin for the particular snapshot you're interested in	
	sigma = []
	s = int(snapno)
	for i in range(0,len(data2)):
		if data2[i,0] == s:
			for j in range(2,19,3):
				sigma.append(data2[i,j])
	#######################
	
	data = np.loadtxt('/projects/b1011/kyle/cmc/branches/cmc-mpi/m22_project/rundir/'+cluster+'/'+filestring+'/initial.snap'+snapno+'.hrdiag_modified.dat')
	data2 = np.loadtxt('/projects/b1011/kyle/cmc/branches/cmc-mpi/m22_project/rundir/'+cluster+'/'+filestring+'/initial.snap'+snapno+'.dat.gz')
	data[:,0] = data[:,0]*pc   # converts from code units to pc
	data[:,0] = conversions.pc_to_arcsec(data[:,0],4.4)   # converts radius from pc to arcsec

	M_G_array = []
        M_MS_array = []
        m_array = []
	count_G = 0
	count_MS = 0
	count = 0	
	for j in range(0,6):
		M_G = 0
		M_MS = 0
		m = 0
		count_G = 0
		count_MS = 0
		count = 0
		if j == 0:
			r_min = 26.67
			r_max = 79.85
		if j == 1:
			r_min = 80.173
			r_max = 123.844
                if j == 2:
                        r_min = 127.816
                        r_max = 168.111
                if j == 3:
                        r_min = 169.443
                        r_max = 234.509
                if j == 4:
                        r_min = 237.539
                        r_max = 326.939
                if j == 5:
                        r_min = 335.386
                        r_max = 611.081
		for i in range(0,len(data)):
			if data[i,0] > r_min and data[i,0] <= r_max:
				m = m + data[i,8] 	# If any object (single or binary)
				count = count + 1	# counts all objects
				if data[i,1] == 0: 	# if a single star
					if data[i,2] >= 3 and data[i,2] <= 6:	# If a giant star						
						M_G = M_G + data[i,8]
						count_G = count_G + 1
					if data[i,2] <= 1:	# If a MS star
						M_MS = M_MS + data[i,8]
						count_MS = count_MS + 1
				else:			# If a binary
					if data[i,4] >= 3 and data[i,4] <= 6:	# If component 0 is a giant
						M_G = M_G + data[i,9]
						count_G = count_G + 1
					if data[i,5] >= 3 and data[i,5] <= 6:	# If component 1 is a giant
                                                M_G = M_G + data[i,10]
                                                count_G = count_G + 1
                                        if data[i,4] <= 1:			# If component 0 is a MS
                                                M_MS = M_MS + data[i,9]
                                                count_MS = count_MS + 1
                                        if data[i,5] <= 1:			# If component 1 is a MS
                                                M_MS = M_MS + data[i,10]
                                                count_MS = count_MS + 1

		M_G_array.append(M_G/count_G)
		M_MS_array.append(M_MS/count_MS)
		m_array.append(m/count)
		print M_G/count_G, M_MS/count_MS, m/count, count_G, count_MS, count
	
	for i in range(0,len(m_array)):
		# Convert the arrays from Msun to kg 
		M_G_array[i] = M_G_array[i] * 1.99e30
		M_MS_array[i] = M_MS_array[i] * 1.99e30
		m_array[i] = m_array[i] * 1.99e30
		sigma[i] = sigma[i] * 1000.	# convert sigma to m/s
	G = 6.67e-11
	a = []
	for i in range(0,len(m_array)):
		a.append(G*M_G_array[i]*M_MS_array[i]/(2.0*m_array[i]*sigma[i]**2.)/1.496e11)
	print a
	print sigma

	# find mean giant radius for entire cluster
	R = 0
	count = 0
	for i in range(0,len(data2)):
		if data2[i,7] == 0:
			if data2[i,14] >= 3 and data2[i,14] <= 6:
				R = R + data2[i,16]
				print data2[i,16], data2[i,14], data2[i,1]
				count = count + 1
		else:
			if data2[i,17] >= 3 and data2[i,17] <= 6:
				R = R + data2[i,21]
				count = count + 1
			if data2[i,18] >= 3 and data2[i,18] <= 6:
                                R = R + data2[i,22]
                                count = count + 1
	print R/count*1.0, count

def giant_radius(cluster, filestring, snapno):
	data = np.loadtxt('/projects/b1011/kyle/cmc/branches/cmc-mpi/m22_project/rundir/'+cluster+'/'+filestring+'/initial.snap'+snapno+'.dat.gz')
	array = []
	count = 0
	for i in range(0,len(data)):
		if data[i,14] == 4:
			print data[i,14], data[i,1], data[i,16]
			array.append(data[i,16])
			count = count + 1
	print np.mean(array), np.min(array), np.max(array), count	

def core_radius_obs(filestring,string,snapno):
        import blackhole

        filestring = filestring+string
        filename2dproj = filestring+'.snap'+snapno+'.2Dproj.dat'


        if not os.path.exists(filename2dproj):
                print "Making 2D proj file..."
                blackhole.make_2D_projection(filestring, snapno, seedy=100, proj=(0,1))

        hllow, hlhigh, hllowcut, hlhighcut, t_myr = blackhole.get_half_light_obs(filename2dproj, Lcut=20.)

        print hllow, hlhigh, hllowcut, hlhighcut, t_myr

        rhl = (hllow+hlhigh)/2.

        r2d, Lcum, p_opt, p_cov, t_myr = blackhole.get_kingfit_params(filename2dproj, rhl, p0guess=[1e5, 1.])

        rc = 1.169*np.abs(p_opt[1])
        return rc
        #print 'rc =', rc, 'time =', t_myr

def find_max_snap_time(path,string):
        snapstring = path+string+'.snap*.dat.gz'
        snapfiles = np.sort(glob.glob(snapstring))

        snapno_max = len(snapfiles) - 1
        if snapno_max < 10:
                snapno_max_str = '000'+str(snapno_max)
        if snapno_max < 100 and snapno_max >= 10:
                snapno_max_str = '00'+str(snapno_max)
        if snapno_max < 1000 and snapno_max >= 100:
                snapno_max_str = '0'+str(snapno_max)
        if snapno_max < 10000 and snapno_max >= 1000:
                snapno_max_str = str(snapno_max)
        if snapno_max >= 10000:
                snapno_max_str = str(snapno_max)
        return snapno_max_str

def find_snap_time_array(path,string):
        snapno_max = int(find_max_snap_time(path,string))
        units=scripts.read_units(path+string)
        km = units[0]['l_cgs']*1.0e-5
        time_units = units[0]['t_myr']

        time_array = []
        snap_array = []
        for j in range(0,snapno_max+1):
                if j < 10:
                        snapno = '000'+str(j)
                if j < 100 and j >= 10:
                        snapno = '00'+str(j)
                if j < 1000 and j >= 100:
                        snapno = '0'+str(j)
                if j < 10000 and j >= 1000:
                        snapno = str(j)
                if j >= 10000:
                        snapno = str(j)

                f2 = gzip.open(path+string+'.snap'+snapno+'.dat.gz','r')
                line = f2.readline()
                line = line.split(' ')
                line = line[1]
                parsed = line.split('=')
                time = float(parsed[1])*time_units      # DEFINE THE TIME OF THE SNAPFILE
                time_array.append(time)
                snap_array.append(snapno)
        return time_array, snap_array
	
