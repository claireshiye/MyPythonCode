import numpy as np
import conversions
import sys
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
matplotlib.use('PDF')

#path = np.str(sys.argv[1])

#snapno = np.str(sys.argv[2])

cluster = np.str(sys.argv[1])

ft = open('trager.dat','r')
linest = ft.readlines()

fv = open('vel_dispersion_profiles.dat','r')
linesv = fv.readlines()

R_sun_obs=4.3


flag = 0
#### find observed SBP from trager file
arcsec_t = []
SB_t = []
for k in range(len(linest)):
        linet = linest[k]
        linet = linet.split('\n')
        linet = linet[0]
        linet = linet.split('\t')
        cluster_tregar = linet[1]
        if cluster_tregar == cluster:
                arcsec_t.append(np.float(linet[2]))
                SB_t.append(np.float(linet[3]))
                flag = 1
arcsec_t = np.array(arcsec_t[:])
SB_t = np.array(SB_t[:])
if flag == 0:
        print('Trager doesnt have the SBP or sigma v profile')
        sys.exit()

flag = 0
#### find observed sigma_v profile
datav = np.genfromtxt('vel_dispersion_profiles.dat')
R_obs = [] # in arcsec
sigma_obs = []
sigma_err_obs_up = []
sigma_err_obs_down = []
for k in range(len(linesv)):
        linev = linesv[k]
        linev = linev.split('\n')
        linev = linev[0]
        linev = linev.split(' ')
        cluster_sigma = linev[0]
        if cluster_sigma == cluster:
                R_obs.append(datav[k,1])
                sigma_obs.append(datav[k,2])
                sigma_err_obs_up.append(datav[k,3])
                sigma_err_obs_down.append(datav[k,4])
                flag = 1
R_obs = np.array(R_obs[:])
sigma_obs = np.array(sigma_obs[:])
sigma_err_obs_up = np.array(sigma_err_obs_up[:])
sigma_err_obs_down = np.array(sigma_err_obs_down[:])
############
if flag == 0:
        print 'No observed sigmav profile'
        sys.exit()

#pathfile=np.genfromtxt('/projects/b1095/syr904/cmc/47Tuc/rundir/lowIMF/path_lowimf.dat', dtype=str)
#snap2d=pathfile[:,1]; paths=pathfile[:,0]; mass=pathfile[:,-1]
paths=['/projects/b1095/syr904/cmc/47Tuc/rundir/rv1.0_rg7.4_z0.0038_N7.5e5_alpha32.3_fb0.05/', '/projects/b1095/syr904/cmc/47Tuc/rundir/rv1.0_rg7.4_z0.0038_N7.5e5_alpha32.3_fb0.05/', '/projects/b1095/syr904/cmc/47Tuc/rundir/rv1.0_rg7.4_z0.0038_N7.5e5_alpha32.3_fb0.05/', '/projects/b1095/syr904/cmc/47Tuc/rundir/elson_profile/rv1.0_N7.5e5_elson/', '/projects/b1095/syr904/cmc/47Tuc/rundir/elson_profile/rv1.0_N7.5e5_elson/', '/projects/b1095/syr904/cmc/47Tuc/rundir/elson_profile/rv1.0_N7.5e5_elson/']; snap2d=['0018','0123', '0398', '0050', '0160', '0385']
#rv=[0.75, 0.75, 1.0, 1.0, 1.0, 1.0, 3.0, 3.0, 3.0, 3.0]
#alpha1=[0.1, 0.4, 0.7, 0.99, 1.3, 1.6]
#cs=['orange', 'blue', 'red', 'black', 'gold', 'cyan']
#lss=['-', '-', '-', '-', '-', '-']
#rv=[1.84, 1.84]; rg=[5.47, 7.4]
prof=['King-100Myr', 'King-1Gyr', 'King-14Gyr', 'Elson-100Myr', 'Elson-1Gyr', 'Elson-14Gyr']
cs=['blue','blue','blue','orange','orange','orange']
lss=['-.','--', '-', '-.', '--', '-']

fig, (ax1, ax2)=plt.subplots(2, 1, sharex=True, figsize=(10,16))
ax2.scatter(R_obs, sigma_obs,c='gold',s=30, edgecolor='black')
ax2.errorbar(R_obs,sigma_obs,yerr=2*sigma_err_obs_up, fmt='o',c='gold', lw=2.0)
ax1.scatter(10**arcsec_t[:], SB_t[:],facecolor='gold',edgecolor='black',s=20,label=r'$\rm{Trager\,et\,al.\,1995}$')
ax2.scatter([10000,10000],[-5,-5],c='gold',s=30, edgecolor='black',label=r'$\rm{Baumgardt\,&\,Hilker\,(2018)}$')

for i in range(0, len(paths)):
    path=paths[i]
    snapno=snap2d[i]
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

    ax1.plot(arcsec_cut, SB_cut, lw=2, color=cs[i], ls=lss[i])
    ax2.scatter(R_model, sigma_model,s=15,zorder=2,alpha=0.5, color=cs[i])
    ax2.errorbar(R_model,sigma_model,yerr=2*sigma_err_model, fmt='o',markersize=0.01,zorder=1,alpha=0.5, color=cs[i])
    #ax2.scatter(R_obs, sigma_obs,c='gold',s=30, edgecolor='black')
    #ax2.errorbar(R_obs,sigma_obs,yerr=2*sigma_err_obs_up, fmt='o',c='gold', lw=2.0)
    #############################

    #ax1.title(cluster,fontsize=24)
    #ax1.text(300,15,cluster,fontsize=20)
    ax1.plot([10000,10000],[-5,-5], lw=2, ls=lss[i], label='profile:'+prof[i], color=cs[i])
    #ax1.scatter(10**arcsec_t[:], SB_t[:],facecolor='gold',edgecolor='black',s=20,label=r'$\rm{Trager\,et\,al.\,1995}$')
    ax1.set_xscale('log')
    ax1.set_xlim(0.03,5000)
    ax1.set_ylim(30,6)
    #ax1.set_xlabel(r'$r\,(\rm{arcsec})$',fontsize=20)
    ax1.set_ylabel(r'$\mu_v$',fontsize=24)
    ax1.legend(loc='best',scatterpoints=1, ncol=2)

    ax2.scatter([10000,10000],[-5,-5],lw=2,label='profile:'+prof[i], color=cs[i])
    #ax2.scatter([10000,10000],[-5,-5],c='gold',s=30, edgecolor='black',label=r'$\rm{Baumgardt\,&\,Hilker\,(2018)}$')
    ax2.set_xscale('log')
    ax2.set_xlim(0.03,5000)
    #ax2.set_ylim(0,10)
    ax2.set_xlabel(r'$r\,(\rm{arcsec})$',fontsize=20)
    ax2.set_ylabel(r'$\sigma_v\,(\rm{km\,s^{-1}})$',fontsize=24)
    ax2.legend(loc=3,scatterpoints=1, ncol=2)
    #plt.show()
    plt.savefig('/projects/b1095/syr904/cmc/47Tuc/rundir/elson_profile/vel_sbp_elson_multipletime.pdf', dpi=300)
