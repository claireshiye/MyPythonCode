import numpy as np
import conversions
import sys
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt

cluster = np.str(sys.argv[1])

bf_number = np.float(sys.argv[2])

data = np.genfromtxt('cluster_params_all_final2.dat')
f4 = open('cluster_params_all_final2.dat','r')
lines4 = f4.readlines()

f = open('harris.dat','r')
lines = f.readlines()

data_h = np.genfromtxt('harris.dat')

f2 = open('trager.dat','r')
lines2 = f2.readlines()

f5 = open('vel_dispersion_profiles.dat','r')
lines5 = f5.readlines()

f3 = open('path.dat','r')
lines3 = f3.readlines()

flag = 0
#### find observed SBP from trager file
arcsec_t = []
SB_t = []
for k in range(len(lines2)):
	line2 = lines2[k]
	line2 = line2.split('\n')
	line2 = line2[0]
	line2 = line2.split('\t')
	cluster_tregar = line2[1]
	if cluster_tregar == cluster:
		arcsec_t.append(np.float(line2[2]))
		SB_t.append(np.float(line2[3]))
		flag = 1
arcsec_t = np.array(arcsec_t[:])
SB_t = np.array(SB_t[:])
if flag == 0:
        print 'Trager doesnt have the SBP or sigma v profile'
        sys.exit()

flag = 0
#### find observed sigma_v profile
data5 = np.genfromtxt('vel_dispersion_profiles.dat')
R_obs = [] # in arcsec
sigma_obs = []
sigma_err_obs_up = []
sigma_err_obs_down = []
for k in range(len(lines5)):
        line5 = lines5[k]
	line5 = line5.split('\n')
	line5 = line5[0]
        line5 = line5.split(' ')
        cluster_sigma = line5[0]
        if cluster_sigma == cluster:
		R_obs.append(data5[k,1])
		sigma_obs.append(data5[k,2])
		sigma_err_obs_up.append(data5[k,3])
		sigma_err_obs_down.append(data5[k,4])
                flag = 1
R_obs = np.array(R_obs[:])
sigma_obs = np.array(sigma_obs[:])
sigma_err_obs_up = np.array(sigma_err_obs_up[:])
sigma_err_obs_down = np.array(sigma_err_obs_down[:])
############
if flag == 0:
	print 'No observed sigmav profile'
        sys.exit()

#0)n 1)ID 2)Rsun(kpc) 3)Rgc(kpc) 4)MVt 5)Fe/H 6)c 7)Rc(arcmin) 8)Rh(arcmin) 9)muV(mag/arcsec^2)
for i in range(0,len(lines)-1):
	line = lines[i+1]
	line = line.split('\n')
	line = line[0]
	data2 = line.split('\t')
	ID_h = data2[1]
	Rgc_h = data_h[i,3]
	R_sun_obs = data_h[i,2]
	Z_h = conversions.metallicity(data_h[i,5],'fe/htoz')  
	rc_h = np.float(data_h[i,7])*(R_sun_obs*1000)/3437.75
	rh_h = np.float(data_h[i,8])*(R_sun_obs*1000)/3437.75
	M_bolo = np.float(data_h[i,4])#-0.107       #BC = -0.107 for sun
	Mass_h = 2.0*10**(0.4*(4.74-M_bolo))      

	if cluster == ID_h:
		if Z_h < 0.0006:
			Z_bin = 0.0002
		if Z_h >= 0.0006 and Z_h < 0.006:
			Z_bin = 0.002
		if Z_h >= 0.006:
			Z_bin = 0.02

		if Rgc_h < 4:
			Rgc_bin = 2
		if Rgc_h >= 4 and Rgc_h < 15:
			Rgc_bin = 8
		if Rgc_h >= 15:
			Rgc_bin = 20

		Mass_arr = []
		t_arr = []
		rc_arr = []
		rh_arr = []
		Nbh_arr = []
		Nns_arr = []
		chi_sq_arr = []
		rv_arr = []
		N_arr = []
		t_arr = []
		snapno_arr = []
		path_arr = []
		print Z_bin, Rgc_bin
		for j in range(0,len(data)):
			#0)i, 1)time, 2)number_density, 3)props['Ltot'], 4)props['M/L'], 5)props['Mave'], 6)props['drc'], 7)props['drhl'], 8)props['dsigmac'], 9)props['dvsigmac_rv'], 10)props['rc'], 11)props['rhl'], 12)props['sigmac'],  13)props['vsigmac_rv'], 14)number_density2, 15)MS, 16)G, 17)Nbh_tot, 18)BHnonBH, 19)BBH, 20)Nns_tot, 21)NSnonNS, 22)BNS, 23)N, 24)Z, 25)rv, 26)rg
			line4 = lines4[j]
			line4 = line4.split('\n')
			line4 = line4[0]
			line4 = line4.split(' ')
			snapno = line4[27]
			t = data[j,1]
			Z = data[j,24]
			Rgc = data[j,26]
			Mass = data[j,3]*data[j,4]
			rc = data[j,10]
			rh = data[j,11]
			Nbh = data[j,17]
			Nns = data[j,20]
			RV = data[j,25]		
			Number = data[j,23]
			if Z == Z_bin and Rgc == Rgc_bin and t > 10000.0:
			#if Z == Z_bin and t > 10000.0 and RV == 0.5:
			#if t > 10000.0 and Number == 32 and Z == 0.0002:
				Mass_arr.append(Mass)
				t_arr.append(t)
				rc_arr.append(rc)
				rh_arr.append(rh)
				Nbh_arr.append(Nbh)
				Nns_arr.append(Nns)
				rv_arr.append(data[j,25])
				N_arr.append(data[j,23])
				t_arr.append(t)
				snapno_arr.append(snapno)
				path_num = np.int(data[j,0])
				line3 = lines3[path_num]
				line3 = line3.split(' ')
				path = line3[0]
				path_arr.append(path)
				#print path, snapno, Z, Rgc, Z_bin, Rgc_bin,	
				####### FIND SBP CHI_SQUARED ###############    
                                #data = np.genfromtxt(path+'initial.snap'+snapno+'.2D_SBPLcut15_vband.dat')
				try:
					data5 = np.genfromtxt(path+'initial.snap'+snapno+'.2D_SBPLcut15.dat')
					arcsec = conversions.pc_to_arcsec(data5[:,1],R_sun_obs)
					SB = conversions.SB_converter(data5[:,3])
					SBerr = data5[:,6]/data5[:,5]*SB

					arcsec_cut = []
					SB_cut = []
					SBerr_cut = []
					for k in range(len(SB)):
						if arcsec[k] < 10 and SB[k] > 20:
							Nothing = 0
						else:
							arcsec_cut.append(arcsec[k])
							SB_cut.append(SB[k])
							#SBerr_cut.append(SBerr[k])
							SBerr_cut.append(1)
					total_chi_squared_SBP = 0
					for k in range(len(arcsec_cut)):
						Obs_mag = np.interp(arcsec_cut[k], 10**arcsec_t[:], SB_t[:])

						chi_squared = pow((SB_cut[k] - Obs_mag),2.0)/pow(SBerr_cut[k],2.0)
						total_chi_squared_SBP = total_chi_squared_SBP + chi_squared
				except:
					total_chi_squared_SBP = 1000.0
				try:
					####### FIND SIGMA_V CHI_SQUARED ################
					dataG = np.genfromtxt(path+'initial.snap'+snapno+'.vel_dispersion_giants_25.dat')
					R_model = conversions.pc_to_arcsec(dataG[:,0],R_sun_obs)
					sigma_model = dataG[:,1]
					sigma_err_model = dataG[:,2]

					total_chi_squared_sigma = 0
					count = 0
					sigma_sq_array = []
					mean_array = []
					var_array = []
					for k in range(len(R_obs)):
						## Find model data point closest to each obs data point
						diff = 1000.
						loc = 0
						for n in range(len(R_model)):
							diff_temp = np.abs(R_obs[k]-R_model[n])
							if diff_temp < diff:
								diff = diff_temp
								loc = n
						chi_squared = pow(sigma_model[loc] - sigma_obs[k],2.0)/(sigma_err_model[loc]**2. + sigma_err_obs_up[k]**2.)
						total_chi_squared_sigma = total_chi_squared_sigma + chi_squared
				except:
					total_chi_squared_sigma=1000.0
					
				total_chi_squared = total_chi_squared_SBP + total_chi_squared_sigma
				#print total_chi_squared_SBP, total_chi_squared_sigma, total_chi_squared
                                #############################################
				chi_sq_arr.append(total_chi_squared)
		chi_sq_arr = np.array(chi_sq_arr[:])
		#loc = np.argmin(chi_sq_arr)
		ind = np.argpartition(chi_sq_arr, bf_number)[:bf_number]
                #print ind
		N_bh_array = []
                #print 'check'
                fig, (ax1, ax2) = plt.subplots(nrows = 2)
		#ax1 = plt.subplot(211)
		#ax2 = plt.subplot(212)
                #print 'check'
                for n in range(len(ind)):
                        loc = ind[n]
			#print ID_h, Mass_h, Mass[loc], rc_h, rc[loc], rh_h, rh[loc], np.len(chi_sq_arr)
                        finfo = open('/projects/b1095/syr904/projects/massive_clusters/matching_SBP/'+cluster+'_info.txt','w')
			print '####################'
			print >> finfo, ID_h, i
			print >> finfo, 'Mass:', Mass_h, Mass_arr[loc]
			print >> finfo, 'rc:', rc_h, rc_arr[loc]
			print >> finfo, 'rh:', rh_h, rh_arr[loc]
			print >> finfo, 'Nbh:', Nbh_arr[loc]
			print >> finfo, 'Nns:', Nns_arr[loc]
			print >> finfo, 'Number of fits:', len(chi_sq_arr), 'Goodness of fit:', chi_sq_arr[loc]
			print >> finfo, 'rv:', rv_arr[loc], 'N:',N_arr[loc], 'time:', t_arr[loc], 'Z:', Z_bin, 'Rgc:', Rgc_bin
			print >> finfo, 'path:', path_arr[loc], 'snapno:', snapno_arr[loc]

			## Get SBP for best-fit model
			data_bf = np.genfromtxt(path_arr[loc]+'initial.snap'+snapno_arr[loc]+'.2D_SBPLcut15.dat')
			arcsec = conversions.pc_to_arcsec(data_bf[:,1],R_sun_obs)
			SB = conversions.SB_converter(data_bf[:,3])
			arcsec_cut = []
			SB_cut = []
			for k in range(len(SB)):
				if arcsec[k] < 10 and SB[k] > 20:
					Nothing = 0
				else:
					arcsec_cut.append(arcsec[k])
					SB_cut.append(SB[k])
			ax1.plot(arcsec_cut, SB_cut, color='black',lw=2)
			## Get sigma_v profile for best-fit model
			try:
				data_bf2 = np.genfromtxt(path_arr[loc]+'initial.snap'+snapno_arr[loc]+'.vel_dispersion_giants_25.dat')
				arcsec2 = conversions.pc_to_arcsec(data_bf2[:,0],R_sun_obs)
				sigma_model = data_bf2[:,1]
				sigma_model_err = data_bf2[:,2]
				ax2.scatter(arcsec2, sigma_model,color='black',s=15,zorder=2,alpha=0.5)
				ax2.errorbar(arcsec2,sigma_model,yerr=2*sigma_model_err, fmt='o',markersize=0.01,color='black',zorder=1,alpha=0.5)
                        except:
                                Nothing = 0			
			ax2.scatter(R_obs, sigma_obs,c='gold',s=30, edgecolor='black')
			ax2.errorbar(R_obs,sigma_obs,yerr=2*sigma_err_obs_up, fmt='o',c='gold', lw=2.0)
			#############################
			N_bh_array.append(Nbh_arr[loc])

		print >> finfo, 'BH mean+range:',np.mean(N_bh_array), np.max(N_bh_array), np.min(N_bh_array)
	        finfo.close()	
                #ax1.title(cluster,fontsize=24)
		ax1.text(300,15,ID_h,fontsize=20)
		ax1.plot([10000,10000],[-5,-5],lw=2,color='black',label=r'$\rm{Best-fit\,models}$')
		ax1.scatter(10**arcsec_t[:], SB_t[:],facecolor='gold',edgecolor='black',s=20,label=r'$\rm{Trager\,et\,al.\,1995}$')
		ax1.set_xscale('log')
		ax1.set_xlim(0.06,5000)
		ax1.set_ylim(30,11)
		ax1.set_xlabel(r'$r\,(\rm{arcsec})$',fontsize=20)
		ax1.set_ylabel(r'$\mu_v$',fontsize=24)
		ax1.legend(loc=3,scatterpoints=1)

		ax2.scatter([10000,10000],[-5,-5],lw=2,color='black',label=r'$\rm{Best-fit\,models}$')
		ax2.scatter([10000,10000],[-5,-5],c='gold',s=30, edgecolor='black',label=r'$\rm{Baumgardt\,&\,Hilker\,(2018)}$')
		ax2.set_xscale('log')
		ax2.set_xlim(0.06,5000)
	        #ax2.set_ylim(0,10)
                ax2.set_xlabel(r'$r\,(\rm{arcsec})$',fontsize=20)
                ax2.set_ylabel(r'$\sigma_v\,(\rm{km\,s^{-1}})$',fontsize=24)
		ax2.legend(loc=3,scatterpoints=1)
		plt.savefig('/projects/b1095/syr904/projects/massive_clusters/matching_SBP/'+cluster+'.pdf', dpi=300)
		break
