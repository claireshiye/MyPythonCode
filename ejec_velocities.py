import numpy as np
import os,sys
import subprocess
import gzip
import scripts
import history_cmc
import finder
import random
import gzip
import matplotlib.pyplot as plt
import glob

f = open('path.dat')
lines = f.readlines()

dataBH = np.genfromtxt('Accretion_BH.dat')
dataNS = np.genfromtxt('Accretion_NS.dat')

ID_array = []

f = open('Accretion_BH2.dat','w')
#print>>f,'# 1) ID 2) t_eject 3) t_MT_start 4) t_MT_finish 5) M1 6) M2 7) k1 8) k2 9) a 10) e 11) v_ejection'

data = dataBH

for i in range(len(data)):
	flag = 0
	#for j in range(len(ID_array)):
	#	if ID_array[j] == dataNS[i,1]:
	#		flag = 1
	#		break

	if flag == 0:
		if data[i,3] != data[i,4]:
			ID_array.append(data[i,1])
			modelID = int(data[i,22])
			line = lines[modelID]
			line = line.split('\n')
			line = line[0]
			line = line.split(' ')
			path = line[0]
			string = line[1]
			print i, modelID, path, string,
			
			# IMPORT THE CONVERSION FROM CODE UNITS TO CGS FROM INITIAL.CONV.SH
			units=scripts.read_units(path+string)
			m = units[0]['l_cgs']*1.0e-2
			kg = units[0]['m_cgs']*1.0e-3
			time = units[0]['t_cgs']	
			nbt = units[0]['nbt_cgs']	

			### LOOK HERE CLAIRE!!!
			E = data[i,20]*m**2./nbt**2.
			t_ejection = data[i,2]
			R_ejection = data[i,12]		# in code units
			R_tidal = data[i,17]
			vr = data[i,13]*m/nbt
			vt = data[i,14]*m/nbt
			M1 = data[i,5]*1.99e30
			M2 = data[i,6]*1.99e30
			a = data[i,10]*1.5e11
			KE_i = 0.5*(M1+M2)*(vr**2.+vt**2.)
			E_int = -6.67e-11*M1*M2/(2*a)
			phi_tidal = data[i,18]*m**2./nbt**2.
			phi_zero = data[i,19]*m**2./nbt**2.*(M1+M2)

                        dyn = open(path+string+'.dyn.dat','r')
			lines10 = dyn.readlines()
			for x in range(2,len(lines10)):
				line = lines10[x]
				line = line.split(' ')
				t_temp = float(line[0])*time/3.15e7/1.e6
				if t_temp >= t_ejection:
					N = float(line[3])
					break
			
			#snap_time_array,snap_array = finder.find_snap_time_array(path,string)
			#for j in range(len(snap_time_array)):
			#	if t_ejection < snap_time_array[j]:
			#		snapno_eject = snap_array[j-1]
                        #                time_before = snap_time_array[j-1]
			#		time_after = snap_time_array[j]
			#		snapno_eject_2 = snap_array[j]
			#		break
			#print j, time_before, time_after, snapno_eject, snapno_eject_2

			#f5 = gzip.open(path+string+'.snap'+snapno_eject+'.dat.gz','r')
			#lines5 = f5.readlines()
			#Menc = 0
			#for j in range(2,len(lines5)):
			#	line5 = lines5[j]
			#	line5 = line5.split(' ')
			#	R = float(line5[2])
			#	M = float(line5[1])
			#	if R <= R_ejection:			
			#		Menc = Menc + M
			#	else:
			#		break
			#print Menc,
			#phi_i = -6.67e-11*Menc*1.99e30*(M1+M2)/(R_ejection*m) + phi_zero
			#E_i = KE_i + phi_i

			#################
			#f6 = gzip.open(path+string+'.snap'+snapno_eject_2+'.dat.gz','r')
			#lines6 = f6.readlines()
			#Menc_2 = 0
                        #for j in range(2,len(lines6)):
                        #        line6 = lines6[j]
                        #        line6 = line6.split(' ')
                        #        R = float(line6[2])
                        #        M = float(line6[1])
                        #        if R <= R_ejection:
                        #                Menc_2 = Menc_2 + M
			#	else:
                        #                break
                        #print Menc_2,
                        #phi_i_2 = -6.67e-11*Menc_2*1.99e30*(M1+M2)/(R_ejection*m) + phi_zero
                        #E_i_2 = KE_i + phi_i_2
			##################
			
			#### INTERPOLATE THE POTENTIAL FROM THE TWO SNAPSHOTS  #####
			#times = [time_before,time_after]
			#potentials = [phi_i,phi_i_2]
			#phi_new = np.interp(t_ejection,times,potentials)
			#E_i = KE_i + phi_new
			#######################################	
			
			Nstar = N
			GAMMA = 0.01
			gierszalpha = 1.5 - 3.0 * (np.log(GAMMA * Nstar) / Nstar)**0.25
			print 'Nstar=', Nstar, 'gierszalpha=',gierszalpha, t_ejection
			#v_eject = np.sqrt(2*(E_i/(M1+M2) - gierszalpha*phi_tidal))

                        ################
                        PE_tidal = phi_tidal*gierszalpha
                        KE = 2*(E-PE_tidal)
			if KE >= 0.:
				v_alt = np.sqrt(2*(E-PE_tidal))
                        else:
				v_alt = -1
			#################	
			if v_alt >= 0:
				#print 'KE=',KE_i, 'PE=',phi_new, 'E_i=',E_i/(M1+M2), 'phi_tidal=',phi_tidal*gierszalpha, 'v_eject=',v_eject, '///////', 'E=',E, 'PE_rtidal=',PE_tidal, 'v_alt=',v_alt
        			print>>f, data[i,22], data[i,2], data[i,3],data[i,4],data[i,5],data[i,6],data[i,7],data[i,8],data[i,10], data[i,11], v_alt, data[i,16], data[i,17]

	
