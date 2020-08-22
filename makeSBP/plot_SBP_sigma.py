import numpy as np
import conversions
import sys
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
matplotlib.use('PDF')

path = np.str(sys.argv[1])

snapno = np.str(sys.argv[2])

cluster = np.str(sys.argv[3])

bf_number = 1

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
        print('Trager doesnt have the SBP or sigma v profile')
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

		#data = np.genfromtxt(path+'initial.snap'+snapno+'.cluster_params.dat')
		#t = data[0,1]
		#Z = data[0,24]
		#Rgc = data[0,26]
		#Mass = data[0,3]*data[0,4]
		#rc = data[0,10]
		#rh = data[0,11]
		#Nbh = data[0,17]
		#Nns = data[0,20]
		#RV = data[0,25]
		#Number = data[0,23]		

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

		dataG = np.genfromtxt(path+'initial.snap'+snapno+'.vel_dispersion_giants_25.dat')
		R_model = conversions.pc_to_arcsec(dataG[:,0],R_sun_obs)
		sigma_model = dataG[:,1]
		sigma_err_model = dataG[:,2]

		fig, (ax1, ax2)=plt.subplots(2, 1)
		#ax1 = plt.subplot(211)
		#ax2 = plt.subplot(212)
		#print ID_h, Mass_h, Mass[loc], rc_h, rc[loc], rh_h, rh[loc], np.len(chi_sq_arr)
		#print '####################'
		#print ID_h, i
		#print 'Mass:', Mass_h, Mass
		#print 'rc:', rc_h, rc
		#print 'rh:', rh_h, rh
		#print 'Nbh:', Nbh
		#print 'Nns:', Nns
		#print 'rv:', data[0,25], 'N:',data[0,23], 'time:', t, 'Z:', Z_bin, 'Rgc:', Rgc_bin

		ax1.plot(arcsec_cut, SB_cut, color='black',lw=2)
		ax2.scatter(R_model, sigma_model,color='black',s=15,zorder=2,alpha=0.5)
		ax2.errorbar(R_model,sigma_model,yerr=2*sigma_err_model, fmt='o',markersize=0.01,color='black',zorder=1,alpha=0.5)
		ax2.scatter(R_obs, sigma_obs,c='gold',s=30, edgecolor='black')
		ax2.errorbar(R_obs,sigma_obs,yerr=2*sigma_err_obs_up, fmt='o',c='gold', lw=2.0)
		#############################

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
		#plt.show()
		plt.savefig(path+'vel_sbp_'+snapno+'.pdf', dpi=300)	
		break
