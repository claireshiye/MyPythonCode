import numpy as np
import finder
import finder2
import finder3
import numpy as np
import os,sys
import subprocess
import gzip
import scripts
import matplotlib.pyplot as plt
import history_cmc
import blackhole
import extract_observed_prop as OBS


path = '/projects/b1095/syr904/cmc/47Tuc/rundir/47Tuc_size/MOCHA47Tuc_elson_rv4_3e6/'
N=3000000
Z=0.0038
rv=4
rg=7.4

string = 'initial'
units=scripts.read_units(path+string)
km = units[0]['l_cgs']*1.0e-5
time_units = units[0]['t_myr']
m_cgs = units[0]['l_cgs']*1.0e-2
kg = units[0]['m_cgs']*1.0e-3
time_cgs = units[0]['t_cgs']
nbt = units[0]['nbt_cgs']
M_total_units = units[0]['m_msun']
pc = units[0]['l_pc']

time_array, snap_array = finder.find_snap_time_array(path,string)
#print snap_array
snapno_max_str = snap_array[-1]
snapno_max = int(snapno_max_str)
Delta = -5  #default -5000 for only making the last snapshot

for k in range(len(time_array)-1,-1,Delta):
	time = time_array[k]
	snapno = snap_array[k]
	print 'time=', time, 'snapno=', snapno,
	if time < 10000:
		print 'stop!'
		break
	try:
		f = open(path+'initial.snap'+snapno+'.vel_dispersion_giants_25.dat','r')
		print snapno, 'is done'
		continue
	except:
		try:
			OBS.make_2D_projection(path+'initial', snapno, units)
			os.system('gzip '+path+'initial.snap'+snapno+'.2Dproj.dat')
                        print 'made 2D projection'
			###### make cluster params file
			f2 = open(path+'initial.snap'+snapno+'.cluster_params.dat','w')
			f5 = gzip.open(path+'initial.snap'+snapno+'.dat.gz','r')
			lines5 = f5.readlines()
			props = OBS.get_obs_props(path+'initial', snapno, FAC=1.)
                        print 'props=', props
			rc = props['rc']
			print 'rc=', rc

			##Initialization
			count_obj=0; count=0; BHnonBH=0; BBH=0; NSnonNS=0; BNS=0	
			P=0; MS=0; G=0; WD=0; NS=0; BH=0

			for j in range(2,len(lines5)):
				line5 = lines5[j]
				data = line5.split(' ')
				#print data
				for j in range(0,len(data)-2):
					if data[j] != 'na':
						data[j] = np.float(data[j])   # Convert strings to floats
				r = data[2]*pc
				if r <= rc:
					count_obj += 1
				if data[7] == 1:
					if r <= rc:
						count += 2
					k1 = data[17]
					k2 = data[18]
					M1 = data[8]
					M2 = data[9]
					ID1 = data[10]
					ID2 = data[11]
					a = data[12]
					e = data[13]
					if k1 == 14 and k2 != 14:
						BHnonBH += 1
					if k2 == 14 and k1 != 14:
						BHnonBH += 1
					if k1 == 14 and k2 == 14:
						BBH += 1
					if k1 == 13 and k2 < 13:
						NSnonNS += 1
					if k2 == 13 and k1 < 13:
						NSnonNS += 1
					if k1 == 13 and k2 == 13:
						BNS += 1
				else:
					if r <= rc:
						count += 1
					M = data[1]
					k = data[14]
					if k == 0 and M == 0.001:
						P += 1
					if k <= 1 and M > 0.01:
						MS += 1
					if k>= 10 and k <= 12:
						WD += 1
					if k == 13:
						NS += 1
					if k == 14:
						BH += 1
					if k >= 2 and k <= 9:
						G += 1
			print 'end of loop'
			number_density = count/(rc**3.)
			number_density2 = count_obj/(rc**3.)
			#print 'number_density', number_density

			print>>f2, "#time, number_density, props['Ltot'], props['M/L'], props['Mave'], props['drc'], props['drhl'], props['dsigmac'], props['dvsigmac_rv'], props['rc'], props['rhl'], props['sigmac'],  props['vsigmac_rv'], number_density2"
			print>>f2, time, number_density, props['Ltot'], props['M/L'], props['Mave'], props['drc'], props['drhl'], props['dsigmac'], props['dvsigmac_rv'], props['rc'], props['rhl'], props['sigmac'],  props['vsigmac_rv'], number_density2, MS, G, BH+BHnonBH+2.*BBH, BHnonBH, BBH, NS+NSnonNS+2.*BNS, NSnonNS, BNS, N, Z, rv, rg
			print>>f2, time, number_density, props['Ltot'], props['M/L'], props['Mave'], props['drc'], props['drhl'], props['dsigmac'], props['dvsigmac_rv'], props['rc'], props['rhl'], props['sigmac'],  props['vsigmac_rv'], number_density2, MS, G, BH+BHnonBH+2.*BBH, BHnonBH, BBH, NS+NSnonNS+2.*BNS, NSnonNS, BNS, N, Z, rv, rg
			f2.close()
			f5.close()
			print 'cluster_params done'	
	
			###############
			print 'made params file'
			blackhole.get_sbp_from_2D_projection(path+string, snapno)
			print snapno, 'made SBP'
			finder2.velocity_dispersion(path,string, snapno, ALL=0, Bin_no=100)
			print 'Made vel dispersion for',snapno
		except:
			print snapno, 'failed'




