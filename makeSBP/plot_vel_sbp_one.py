import numpy as np
import conversions
import sys
from glob import glob
import gzip
import matplotlib
import matplotlib.cm as cm
matplotlib.use('agg')
import matplotlib.pyplot as plt
matplotlib.use('PDF')


path = np.str(sys.argv[1])

#snapno = np.str(sys.argv[2])

cluster = np.str(sys.argv[2])

ft = open('trager.dat','r')
linest = ft.readlines()

fv = open('vel_dispersion_profiles.dat','r')
linesv = fv.readlines()

R_sun_obs=4.5


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


snap2D = np.sort(glob(path+'initial.snap*.2Dproj.dat.gz'))
print snap2D
snap2D_L15 = np.sort(glob(path+'initial.snap*.2D_SBPLcut15.dat'))
snap_giant = np.sort(glob(path+'initial.snap*.vel_dispersion_giants_25.dat'))

cs = cm.Blues(np.linspace(0.5, 1, len(snap2D)))
#cmap=matplotlib.colors.ListedColormap([c[0], c[1], c[2], c[3]])


fig, (ax1, ax2)=plt.subplots(2, 1, sharex=True, figsize=(10,16))
ax2.scatter(R_obs, sigma_obs,c='gold',s=30, edgecolor='black')
ax2.errorbar(R_obs,sigma_obs,yerr=2*sigma_err_obs_up, fmt='o',c='gold', lw=2.0)
ax1.scatter(10**arcsec_t[:], SB_t[:],facecolor='gold',edgecolor='black',s=20,label=r'$\rm{Trager\,et\,al.\,1995}$')
ax2.scatter([10000,10000],[-5,-5],c='gold',s=30, edgecolor='black',label=r'$\rm{Baumgardt\,&\,Hilker\,(2018)}$')

for i in range(0, len(snap2D)):
    with gzip.open(snap2D[i], 'r') as f2D:
        first_line=f2D.readline()
            
    t_gyr = first_line.strip().split('=')[-1]
    print t_gyr 
    if t_gyr>12000.: continue

    data5 = np.genfromtxt(snap2D_L15[i])
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

    dataG = np.genfromtxt(snap_giant[i])
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

    ax1.plot(arcsec_cut, SB_cut, lw=2, color=cs[i])
    ax2.scatter(R_model, sigma_model,s=15,zorder=2,alpha=0.5, color=cs[i])
    ax2.errorbar(R_model,sigma_model,yerr=2*sigma_err_model, fmt='o',markersize=0.01,zorder=1,alpha=0.5, color=cs[i])
    #ax2.scatter(R_obs, sigma_obs,c='gold',s=30, edgecolor='black')
    #ax2.errorbar(R_obs,sigma_obs,yerr=2*sigma_err_obs_up, fmt='o',c='gold', lw=2.0)
    #############################

    #ax1.title(cluster,fontsize=24)
    #ax1.text(300,15,cluster,fontsize=20)
    ax1.plot([10000,10000],[-5,-5], lw=2, label='time = '+str(t_gyr), color=cs[i])
    #ax1.scatter(10**arcsec_t[:], SB_t[:],facecolor='gold',edgecolor='black',s=20,label=r'$\rm{Trager\,et\,al.\,1995}$')
    ax1.set_xscale('log')
    ax1.set_xlim(0.03,5000)
    ax1.set_ylim(30,6)
    #ax1.set_xlabel(r'$r\,(\rm{arcsec})$',fontsize=20)
    ax1.set_ylabel(r'$\mu_v$',fontsize=24)
    ax1.legend(loc='best',scatterpoints=1, ncol=2)

    ax2.scatter([10000,10000],[-5,-5],lw=2,label='time = '+str(t_gyr), color=cs[i])
    #ax2.scatter([10000,10000],[-5,-5],c='gold',s=30, edgecolor='black',label=r'$\rm{Baumgardt\,&\,Hilker\,(2018)}$')
    ax2.set_xscale('log')
    ax2.set_xlim(0.03,5000)
    #ax2.set_ylim(0,10)
    ax2.set_xlabel(r'$r\,(\rm{arcsec})$',fontsize=20)
    ax2.set_ylabel(r'$\sigma_v\,(\rm{km\,s^{-1}})$',fontsize=24)
    ax2.legend(loc=3,scatterpoints=1, ncol=2)
    #plt.show()
    plt.savefig(path+'sbp_vel_timeserie.pdf', dpi=300)
