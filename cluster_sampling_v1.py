import numpy as np
import scipy
from scipy.special import gamma, gammainc
import scipy.integrate as integrate
import math
import pandas as pd
import matplotlib.pyplot as plt
import random
from random import choices

import dynamics as dyn

twopi = 2*np.pi
Gconst = 4.30091*10**(-3)   #pc*Msun^-1*(km/s)^2
H = 69.6 #km/s*Mpc^-1


##Stellar mass density profile without dark matter halo
def smdf(r, Mtot = 5.*1e10, Re = 4., ns = 2.2):  ##density in Msun/kpc^3
    ##stellar density
    g2n = gamma(2*ns)
    p = 1.0-0.6097/ns+0.05563/ns**2
    b = 2*ns-1./3.+0.009876/ns
    gn3p = gamma(ns*(3-p))
    L = twopi*Re**2*ns*b**(-2*ns)*g2n

    rho_0 = (Mtot/L)*b**(ns*(1-p))*g2n/(2*Re*gn3p)
    exponent = -b*(r/Re)**(1/ns)
    rho_star = rho_0*(r/Re)**(-p)*np.exp(exponent)

    return rho_star

##Galaxy density profile with dark matter halo
def gdf(r, Mtot = 5.*1e10, Re = 4., ns = 2.2, Mdm=1.e12, rs=20.):  ##density in Msun/kpc^3
    ##stellar density
    g2n = gamma(2*ns)
    p = 1.0-0.6097/ns+0.05563/ns**2
    b = 2*ns-1./3.+0.009876/ns
    gn3p = gamma(ns*(3-p))
    L = twopi*Re**2*ns*b**(-2*ns)*g2n

    rho_0 = (Mtot/L)*b**(ns*(1-p))*g2n/(2*Re*gn3p)
    exponent = -b*(r/Re)**(1/ns)
    rho_star = rho_0*(r/Re)**(-p)*np.exp(exponent)

    ##dark matter halo
    rho_crit = 3*(H**2/10**6)/(4*twopi*Gconst*0.001)
    r200 = (3*Mdm/(200*2*twopi*rho_crit))**(1./3.)
    conc = r200/rs
    delta_c = (200./3.)*conc**3/(np.log(1+conc)-conc/(1+conc))

    rho_0_dm = rho_crit*delta_c
    rho_dm = rho_0_dm/(r/rs)/(1+r/rs)**2

    ##Total
    rho = rho_star+rho_dm

    return rho

##Generate random initial model galactocentric distance
def gdf_frommodel():  ##density in Msun/kpc^3
    init_rg = [2, 8, 20]
    prp_rho = [gdf(x) for x in init_rg]
    sum_weight = sum(prp_rho)
    weights = [y/sum_weight for y in prp_rho]
    #print(weights)
    return choices(init_rg, weights)

#def NFW_DM(r, Mdm=1.e12, rs=20.):    ##density in Msun/kpc^3
#    rho_crit = 3*(H**2/10**6)/(4*twopi*Gconst*0.001)
#    r200 = (3*Mdm/(200*2*twopi*rho_crit))**(1./3.)
#    conc = r200/rs
#    delta_c = (200./3.)*conc**3/(np.log(1+conc)-conc/(1+conc))
#
#    rho_0_dm = rho_crit*delta_c
#    rho_dm = rho_0_dm/(r/rs)/(1+r/rs)**2
#
#    return rho_dm


##Calculate circular velocity
def circ_vel(ra, rho_func):  ##v_circ in km/s; ra in kpc
    Mr = 2*twopi*integrate.quad(lambda x: rho_func(x)*x*x, 0, ra)[0]

    v_circ = np.sqrt(Gconst*Mr*0.001/ra)

    #print(Mr)
    return v_circ


##Calculate galaxy mass density at r_g
def galaxy_rho_r(r_g, vel_func):   ##r_g in kpc; return value is Msun/pc^3
    rg_pc = r_g*1000.
    #print(vel_func(r_g, gdf)**2/(twopi*Gconst*rg_pc**2))
    return vel_func(r_g, gdf)**2/(twopi*Gconst*rg_pc**2)



##Cluster mass function
def gc_mf(mgc, beta=2): ##mmin = 1e4, mmax = 1e7
    return mgc**(-beta)


##Generate random initial model mass
def gc_mf_frommodel():
    init_n = ['5e4', '1e5', '2e5', '4e5', '8e5', '1.6e6', '3.2e6', '8.6e6']
    init_mass = [29655.3, 59638.8, 1.197630e+05, 2.423500e+05, 4.848440e+05, 9.703820e+05, 1.93472e+06, 5.32797e+06]
    inverse_m2 = [1/x**2 for x in init_mass]
    sum_weight = sum(inverse_m2)
    weights = [y/sum_weight for y in inverse_m2]
    #print(weights)
    return choices(init_n, weights)


##Dynamical friction
def df(dt, r_old, Mgc):  ##t_df in Gyr; r_old in kpc
    t_df = 0.23*circ_vel(r_old, gdf)*r_old**2*(1e5/Mgc)  ##default 0.23, fast: 0.1755
    dr = -r_old*dt/(2*t_df)

    r_new = r_old+dr

    #print(r_new)
    return r_new


##Galaxy tidal field cluster disruption
def t_tidal(Mgc, r_gd, vel_func):  ##r_gd in kpc, vel in km/s
    Pr = 41.4*r_gd/vel_func(r_gd, gdf)
    t_tid = 10*(Mgc/1e5)**(2/3)*Pr
    return t_tid


##Mass loss from stellar winds, dynamical ejection of stars through two-body relaxation and stripping of stars by the Galactic tidal field.
def dmdt(dt, mold, rgd):
    tidaltime = t_tidal(mold, rgd, circ_vel)
    dm = -mold*dt/tidaltime
    mnew = dm+mold

    return mnew


##calculate the initial tidal radius of a cluster
def r_tidal_initial(m0_cluster, rg0):  ##m0 in Msun and rg0 in kpc
    r_tidal = (Gconst*m0_cluster/(2.*220.**2.))**(1./3.)*(rg0*1000)**(2./3.)
    return r_tidal  ##in pc



##Generate galaxy density and velocity profile
def galaxy_profile():
    xga = np.logspace(np.log10(0.0001), np.log10(100), 5000)
    print(xga)
    gdf_array = []; vel_array = []; ga_rho_array = []
    for ii in range(len(xga)):
        gdf_array.append(gdf(xga[ii]))
        vel_array.append(circ_vel(xga[ii], gdf))
        ga_rho_array.append(galaxy_rho_r(xga[ii], circ_vel))

        print(ii)

    #print(gdf_array, vel_array, ga_rho_array)

    np.savetxt('/projects/b1095/syr904/projects/GCE/catalog/galaxy_profile.txt', np.c_[xga, gdf_array, vel_array, ga_rho_array], fmt = '%.8f %f %f %f', header = '1.rg(kpc) 2.rho(Msun/kpc^3) 3.vel(km/s) 4.rho(Msun/pc^3)', comments = '#')




##Random sampling
def random_sampling(function, nums, batch=100):
    samples = []
    #samples_x = []; samples_y = []
    #rej_samples_x = []; rej_samples_y = []

    while len(samples) < nums:  
        samples += function()

        #samples_x+= x[y < function(x)].tolist(); samples_y+= y[y < function(x)].tolist()
        #rej_samples_x+= x[y >= function(x)].tolist(); rej_samples_y+= y[y >= function(x)].tolist()

    #print(len(samples_x), len(rej_samples_x))

    return samples[:nums]
    #return samples_x, samples_y, rej_samples_x, rej_samples_y



##Von Neumann rejection sampling
def VN_sampling(function, xmin, xmax, fmax, nums, fmin=0, batch=1000):
    samples = []
    #samples_x = []; samples_y = []
    #rej_samples_x = []; rej_samples_y = []

    while len(samples) < nums:
        x = np.random.uniform(low=xmin, high=xmax, size=batch)
        y = np.random.uniform(low=fmin, high=fmax, size=batch)

        samples += x[y < function(x)].tolist()

        #samples_x+= x[y < function(x)].tolist(); samples_y+= y[y < function(x)].tolist()
        #rej_samples_x+= x[y >= function(x)].tolist(); rej_samples_y+= y[y >= function(x)].tolist()

    #print(len(samples_x), len(rej_samples_x))

    return samples[:nums]
    #return samples_x, samples_y, rej_samples_x, rej_samples_y


##MCMC
#def MCMC():


def plot_raw_distr(function, xmin, xmax, thelabel):
    xs = np.linspace(xmin, xmax, 200)
    ys = function(xs)

    plt.plot(xs, ys, label=thelabel) 
    plt.fill_between(xs, ys, 0, alpha=0.2)
    plt.xlim(xmin/2, xmax*2)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel("r(kpc)")
    plt.ylabel(r'$mass\ density\ (M_{\odot}/kpc^3)$')
    plt.legend()
    plt.show()


def plot_VNsampling_distr(function, xmin_sample, xmax_sample, fmax_sample, num_sample):
    samps = VN_sampling(function, xmin_sample, xmax_sample, fmax_sample, num_sample)

    xs = np.linspace(xmin_sample, xmax_sample, 200)
    ys = function(xs)

    plt.plot(xs, ys, label="Function")
    #plt.hist(samps, bins=np.logspace(np.log10(xmin_sample), np.log10(xmax_sample), 30),
    #    density=False, alpha=0.2, label="Sample distribution")
    plt.hist(samps, bins=50,
        density=False, alpha=0.2, label="Sample distribution")
    #bins = np.logspace(np.log10(xmin_sample), np.log10(xmax_sample), 10),
    plt.xlim(xmin_sample/2, xmax_sample*1.1)
    #plt.xscale('log')
    plt.yscale('log')
    plt.xlabel("x")
    plt.ylabel("f(x)") 
    plt.legend()
    plt.show()


def plot_sampling_distr(function_ran, function, num_sample, xmin_sample, xmax_sample):
    samps = random_sampling(function_ran, num_sample)
    #print(samps)

    xs = np.linspace(xmin_sample, xmax_sample, 1000)
    ys = function(xs)

    plt.plot(xs, ys, label="Function")
    #plt.hist(samps, bins=np.logspace(np.log10(xmin_sample), np.log10(xmax_sample), 100),
    #    density=True, alpha=0.2, label="Sample distribution")
    plt.hist(samps, density=False, alpha=0.2, label="Sample distribution")
    #bins = np.logspace(np.log10(xmin_sample), np.log10(xmax_sample), 10),
    #plt.xlim(xmin_sample/2, xmax_sample*1.2)
    #plt.xscale('log')
    plt.yscale('log')
    plt.xlabel("x")
    plt.ylabel("f(x)") 
    plt.legend()
    plt.show()


##Read rho_rh at initial time
def read_rho_rh(modelpath):
    mstar_conv = dyn.conv('mstar', modelpath+'initial.conv.sh')
    l_conv = dyn.conv('l', modelpath+'initial.conv.sh')

    #datarho = np.genfromtxt(modelpath+'initial.rho_lagrad.dat')
    #rho_rh0 = datarho[:,20][0]*mstar_conv/l_conv**3
    filerho = modelpath+'initial.rho_lagrad.dat'
    with open(filerho, 'r') as frho:
        for i, line in enumerate(frho):
            if i == 2:   # 3th line
                data = line.split()
                rho_rh0 = float(data[20])*mstar_conv/l_conv**3
                break

    return rho_rh0


##Find number of MSP, gamma-ray luminosity and cluster mass at a time
def read_property_attime(modelpath, time, behemoth_flag):   ##time in Myr
    Num_msp = -100; Mass = -100; Lum_gamma = -100
    Num_msp_old = -100; Mass_old = -100; Lum_gamma_old = -100


    if behemoth_flag == 0:
        t_conv = dyn.conv('t', modelpath+'initial.conv.sh')
        m_conv = dyn.conv('m', modelpath+'initial.conv.sh')
        nsfile = modelpath+'initial.ns.dat'
        dynfile = modelpath+'initial.dyn.dat'

        s=modelpath.split('/')
        n_star=s[-2]
        z=s[-3][1:]
        rg=s[-4][2:]
        rv=s[-5][2:]
        lgamafile = '/projects/b1095/syr904/projects/GCE/catalog/data_lgamma/model_rv'+rv+'_rg'+rg+'_z'+z+'_'+n_star+'.txt'

    else:
        t_conv = dyn.conv('t', modelpath+'behemoth.conv.sh')
        m_conv = dyn.conv('m', modelpath+'behemoth.conv.sh')
        nsfile = modelpath+'behemoth.ns.dat'
        dynfile = modelpath+'behemoth.dyn.dat'

        lgamafile = modelpath+'Lgamma_alltime_behemoth.dat'


    ##Extract number of MSP
    with open(nsfile, 'r') as fns:
        next(fns)
        for i, line in enumerate(fns):
            datans = line.split()
            t_ns = float(datans[0])*t_conv
            if i == 0: n_msp_old = int(datans[6])
            if round(t_ns,6) >= round(time,6):
                print(t_ns, time)
                Num_msp = int(datans[6])
                Num_msp_old = n_msp_old
                break
            n_msp_old = int(datans[6])

    if Num_msp == -100:
        Num_msp = int(datans[6]); Num_msp_old = int(datans[6])


    ##Extract mass of cluster
    #datadyn = np.genfromtxt(modelpath+'initial.dyn.dat')
    with open(dynfile, 'r') as fdyn:
        next(fdyn); next(fdyn)
        for i, line in enumerate(fdyn):
            datadyn = line.split()
            t_dyn = float(datadyn[0])*t_conv
            if i == 0: m_dyn_old = float(datadyn[4])*m_conv
            if round(t_dyn,6) >= round(time,6):
                print(t_dyn, time)
                Mass = float(datadyn[4])*m_conv
                Mass_old = m_dyn_old
                break
            m_dyn_old = float(datadyn[4])*m_conv

    if Mass == -100:
        Mass = float(datadyn[4])*m_conv; Mass_old = float(datadyn[4])*m_conv


    ##Extract Lgamma
    if Num_msp != 0 or Num_msp_old != 0:
        #sourcedir = np.genfromtxt('/projects/b1095/syr904/projects/GCE/catalog/path_allfinished_newruns_maingrid.dat', dtype=str)
        #paths = sourcedir[:,0]
        #model_no = list(paths).index(modelpath)

        datagama = np.genfromtxt(lgamafile)
        model_gama = datagama[:,0]; t_gama = datagama[:,1]; lgama = datagama[:,2]

        #t_model = t_gama[model_gama==model_no]; lum_model = lgama[model_gama==model_no]
        t_model = t_gama; lum_model = lgama
        for xx in range(len(t_model)):
            if round(t_model[xx],6) >= round(time,6):
                Lum_gamma = lum_model[xx]
                Lum_gamma_old = lum_model[xx-1]


    print(Num_msp, Mass, Lum_gamma, Num_msp_old, Mass_old, Lum_gamma_old)
    return Num_msp, Mass, Lum_gamma, Num_msp_old, Mass_old, Lum_gamma_old


####Find number of MSP, gamma-ray luminosity and cluster mass at a time for all models
def read_property_all(start, end, disrupfile, proptyfile, file_ver):
    sample_disrp = np.genfromtxt('/projects/b1095/syr904/projects/GCE/catalog/'+disrupfile+'.dat', dtype = str)
    nstar_disrp = sample_disrp[:,0]; rv_disrp = sample_disrp[:,1]
    rg_model_disrp = sample_disrp[:,2]; z_disrp = sample_disrp[:,3]
    rg_init = sample_disrp[:,5].astype(np.float)
    rg_disrp = sample_disrp[:,6].astype(np.float)
    t_disrp = sample_disrp[:,7].astype(np.float)
    print(np.where(t_disrp == 0))

    prop = [[],[],[]]; prop_old = [[],[],[]]
    for ii in range(start, end):
        if float(nstar_disrp[ii])<=3.2e6:
            thepath = '/projects/b1091/CMC_Grid_March2019/rundir/rv'+rv_disrp[ii]+'/rg'+rg_model_disrp[ii]+'/z'+z_disrp[ii]+'/'+nstar_disrp[ii]+'/'
            bhm_flag = 0
        else:
            thepath = '/projects/b1095/syr904/projects/GCE/behemoth/'
            bhm_flag = 1
        
        if float(nstar_disrp[ii])>1e5:
            Nmsp, Mgc, Lmsp, Nmsp_old, Mgc_old, Lmsp_old = read_property_attime(thepath, t_disrp[ii], bhm_flag)
        else:
            Nmsp = 0.03; Lmsp = 3.5e34
            Nmsp_old = 0.03; Lmsp_old = 3.5e34
            Mass = -100

            t_conv = dyn.conv('t', thepath+'initial.conv.sh')
            m_conv = dyn.conv('m', thepath+'initial.conv.sh')

            ##Extract mass of cluster
            with open(thepath+'initial.dyn.dat', 'r') as fdyn:
                next(fdyn); next(fdyn)
                for i, line in enumerate(fdyn):
                    datadyn = line.split()
                    t_dyn = float(datadyn[0])*t_conv
                    if i == 0: m_dyn_old = float(datadyn[4])*m_conv
                    if round(t_dyn,6) >= round(t_disrp[ii],6):
                        print(t_dyn, t_disrp[ii])
                        Mass = float(datadyn[4])*m_conv
                        Mass_old = m_dyn_old
                        break
                    m_dyn_old = float(datadyn[4])*m_conv

            if Mass == -100:
                Mass = float(datadyn[4])*m_conv; Mass_old = float(datadyn[4])*m_conv

        if Mass == 0 or Mass_old == 0:
            print(modelpath)
        prop[0].append(Nmsp); prop[1].append(Mgc); prop[2].append(Lmsp)
        prop_old[0].append(Nmsp_old); prop_old[1].append(Mgc_old); prop_old[2].append(Lmsp_old)
        print(ii)

    print('finished searching')

    np.savetxt('/projects/b1095/syr904/projects/GCE/catalog/'+proptyfile+file_ver+'.dat', np.c_[nstar_disrp[start:end], rv_disrp[start:end], rg_model_disrp[start:end], z_disrp[start:end], t_disrp[start:end], prop[0], prop[1], prop[2], prop_old[0], prop_old[1], prop_old[2]], fmt = '%s %s %s %s %s %s %s %s %s %s %s', header = '1.N_star 2.RV 3.RG 4.Z 5.T_DISRUP(Myr) 6.Nmsp 7.Mgc(Msun) 8.Lmsp(erg/s) 9.Nmsp_old 10.Mgc_old(Msun) 11.Lmsp_old(erg/s)', comments = '#')


##Extract rho and dyn data
def extract_cluster_data(cluster_path):
    #prefix = 'behemoth'
    prefix = 'initial'
    t_conv = dyn.conv('t', cluster_path+prefix+'.conv.sh')
    l_conv = dyn.conv('l', cluster_path+prefix+'.conv.sh')
    m_conv = dyn.conv('m', cluster_path+prefix+'.conv.sh')
    mstar_conv = dyn.conv('mstar', cluster_path+prefix+'.conv.sh')

    s=cluster_path.split('/')
    n_star=s[-2]
    z=s[-3][1:]
    rg=s[-4][2:]
    rv=s[-5][2:]

    filerho = cluster_path+prefix+'.rho_lagrad.dat'
    filedyn = cluster_path+prefix+'.dyn.dat'

    #with open(filerho, 'r') as fr:
    #    next(fr); next(fr)
    #    for i_rho, l_rho in enumerate(fr):
    #        pass
    #line_count = i_rho+1
    #print(nline_rho)
    #line_cut = int(0.7*line_count)


    data_array = [[],[],[],[]]
    with open(filerho, 'r') as frho, open(filedyn, 'r') as fdyn:
        next(frho); next(frho)
        next(fdyn); next(fdyn)
        for i, lines in enumerate(zip(frho, fdyn)):
            #if i <= line_cut and i % 5 != 0: 
            #    continue
            #if i> line_cut and i % 2 != 0:
            #    continue
            if i % 5 != 0:
                continue

            datarho = lines[0].split()
            datadyn = lines[1].split()

            t_curr = float(datarho[0])
            rho_halfm_curr = float(datarho[20])*(mstar_conv/l_conv**3)

            m_curr = float(datadyn[4])*m_conv

            rh_curr = float(datadyn[20])*l_conv

            data_array[0].append(t_curr)
            data_array[1].append(rho_halfm_curr)
            data_array[2].append(m_curr)
            data_array[3].append(rh_curr)

    
    #np.savetxt('/projects/b1095/syr904/projects/GCE/behemoth/behemoth_rho_dyn.txt', np.c_[data_array[0], data_array[1], data_array[2], data_array[3]], fmt = '%f %f %f %f', header = '1.Time(Code Unit) 2.rho(Msun/pc^3) 3.mass(Msun) 4.rh(pc) 5.t_conv('+str(t_conv)+')', comments = '#')
    np.savetxt('/projects/b1095/syr904/projects/GCE/catalog/data_rho_dyn/model_rv'+rv+'_rg'+rg+'_z'+z+'_'+n_star+'.txt', np.c_[data_array[0], data_array[1], data_array[2], data_array[3]], fmt = '%f %f %f %f', header = '1.Time(Code Unit) 2.rho(Msun/pc^3) 3.mass(Msun) 4.rh(pc) 5.t_conv('+str(t_conv)+')', comments = '#')


##Sampling initial cluster distribution
def cluster_initialization_rg2only(sample_num, dissolflag):
    f_gc = 0.012
    xmin_samps = 0.1; xmax_samps = 5.
    samps_den = VN_sampling(smdf, xmin_samps, xmax_samps, smdf(xmin_samps), sample_num)
    samps_mass = random_sampling(gc_mf_frommodel, sample_num)
    #samps_rg = random_sampling(gdf_frommodel, sample_num)
    samps_rv = np.random.randint(4, size=sample_num)
    samps_z = np.random.randint(3, size=sample_num)

    bin_num = 400
    rg_bin = np.linspace(xmin_samps, xmax_samps, bin_num+1)
    cumu_mass_bin = []
    #rho_bin = []
    #rho_bin = smdf(rg_bin[:-1])
    for kk in range(len(rg_bin)-1):
        rg_med = (rg_bin[kk]+rg_bin[kk+1])/2.
        #rho_bin.append(smdf(rg_med)*f_gc)
        Mr = 2*twopi*integrate.quad(lambda x: smdf(x)*x*x, 0, rg_med)[0]
        cumu_mass_bin.append(Mr*f_gc)

    print('finished sampling')

    modeln = ['2e5', '4e5', '8e5', '1.6e6']
    modelrv = ['0.5', '1', '2', '4']; modelrg = ['2', '8', '20']; modelz = ['0.0002', '0.002', '0.02']
    modelmass = [1.197630e+05, 2.423500e+05, 4.848440e+05, 9.703820e+05]
    rg_range = [[0., 5.], [5., 15.], [15., 100.]]
    print(len(rg_range))

    sourcedir = np.genfromtxt('/projects/b1095/syr904/projects/GCE/catalog/path_allfinished_newruns_maingrid.dat', dtype=str)
    paths = list(sourcedir[:,0]); status = sourcedir[:,1]


    allsamps = [[],[],[],[],[]]
    rg_initial = []
    model_status = []
    time_tidal = []
    model_mass = []
    numgc_tot = 0; massgc_tot = 0; massmodel_tot = 0

    mass_init = np.zeros(bin_num); mass_density_init = np.zeros(bin_num)
    rg_randint = np.ones(bin_num)
    for ii in range(len(samps_mass)):
        check = 0

        N_star = samps_mass[ii]
        rv = modelrv[int(samps_rv[ii])]
        z = modelz[int(samps_z[ii])]

        for xx in range(len(rg_range)):
            #if xx == 0: rg_randint = 1
            #elif xx == 1: rg_randint = 1#np.random.randint(2)
            #elif xx == 2: rg_randint = 1#np.random.randint(17)
            if rg_range[xx][0] <= samps_den[ii] <= rg_range[xx][1]:
                rg = modelrg[xx]

        thepath = '/projects/b1091/CMC_Grid_March2019/rundir/rv'+rv+'/rg'+rg+'/z'+z+'/'+N_star+'/'
        model_no = int(paths.index(thepath))
        thestatus = int(status[model_no])

        if dissolflag == 1 and thestatus !=1: continue


        for xx in range(len(rg_bin)-1):
            if rg_bin[xx] <= samps_den[ii] < rg_bin[xx+1] and rg_randint[xx] == 1:
                check = 1

        if check == 0: continue


        mass = modelmass[modeln.index(N_star)]

        tid_time = t_tidal(mass, samps_den[ii], circ_vel)
        rho_rh_0 = read_rho_rh(thepath)
        rho_ga = galaxy_rho_r(samps_den[ii], circ_vel)

        

        if rho_rh_0 > rho_ga:
            allsamps[0].append(N_star); allsamps[1].append(rv); allsamps[2].append(rg); allsamps[3].append(z); allsamps[4].append(rho_rh_0)
            rg_initial.append(samps_den[ii])
            model_status.append(thestatus)
            time_tidal.append(tid_time)
            model_mass.append(mass)
            numgc_tot += 1
            massgc_tot += mass; massmodel_tot += mass

            for yy in range(len(rg_bin)-1):
                if rg_bin[yy]<=samps_den[ii]<rg_bin[yy+1]:
                    rg_bin_med = (rg_bin[yy+1]+rg_bin[yy])/2.
                    mass_init[yy]+=mass
                    #vol_init = (2.*twopi/3.)*(rg_bin[yy+1]**3-rg_bin[yy]**3)
                    #mass_density_init[yy] = mass_init[yy]/vol_init
    
                    mass_cumu_init = np.cumsum(mass_init)

                    if mass_cumu_init[yy] > cumu_mass_bin[yy]: ##or mass_density_init[yy] > rho_bin[yy]
                        rg_randint[yy] = 0

        if massgc_tot > 5e8: break

        print(ii, thestatus)

    
    print(numgc_tot, massgc_tot, massmodel_tot)
    np.savetxt('/projects/b1095/syr904/projects/GCE/catalog/cluster_sample_initial_M_RG_rg2_nondissol.dat', np.c_[allsamps[0], allsamps[1], allsamps[2], allsamps[3], allsamps[4], rg_initial, model_status, time_tidal, model_mass], fmt = '%s %s %s %s %s %s %s %s %s', header = '1.N_star 2.RV 3.RG 4.Z 5.RHO_RH(Msun/pc^3) 6.RG_INITIAL(kpc) 7.Model_status 8.Ttid(Gyr) 9.Mass(Msun)', comments = '#')

    #return allsamps, rg_initial, model_status


##Sampling initial cluster distribution
def cluster_initialization(sample_num, dissolflag, frac_gc_close, frac_gc_far, xcut):
    #f_gc_close = frac_gc_close; f_gc_far = frac_gc_far ##0.012
    xmin_samps = 0.01; xmax_samps = 100.
    samps_den = VN_sampling(smdf, xmin_samps, xmax_samps, smdf(xmin_samps), sample_num)
    samps_mass = random_sampling(gc_mf_frommodel, sample_num)
    #samps_rg = random_sampling(gdf_frommodel, sample_num)
    samps_rv = np.random.randint(4, size=sample_num)
    samps_z = np.random.randint(3, size=sample_num)

    #samps_den = np.sort(samps_den)

    bin_num = 400
    rg_bin = np.logspace(np.log10(xmin_samps), np.log10(xmax_samps), bin_num+1)
    cumu_mass_bin = []
    #rho_bin = []
    #rho_bin = smdf(rg_bin[:-1])
    check_xcut = 0
    for kk in range(len(rg_bin)-1):
        rg_med = (rg_bin[kk]+rg_bin[kk+1])/2.
        #rho_bin.append(smdf(rg_med)*f_gc)
        if rg_med < xcut:
            Mr = frac_gc_close*2*twopi*integrate.quad(lambda x: smdf(x)*x*x, 0, rg_med)[0]
            rg_med_prev = rg_med
        else:
            if check_xcut == 0:
                xcut = rg_med_prev
                Mr_cut = frac_gc_close*2*twopi*integrate.quad(lambda x: smdf(x)*x*x, 0, xcut)[0]
                check_xcut = 1
            Mr = Mr_cut + frac_gc_far*2*twopi*integrate.quad(lambda x: smdf(x)*x*x, xcut, rg_med)[0]

        cumu_mass_bin.append(Mr)
    
    print('finished sampling')


    modeln = ['5e4', '1e5', '2e5', '4e5', '8e5', '1.6e6']
    modelrv = ['0.5', '1', '2', '4']; modelrg = ['2', '8', '20']; modelz = ['0.0002', '0.002', '0.02']
    modelmass = [29655.3, 59638.8, 1.197630e+05, 2.423500e+05, 4.848440e+05, 9.703820e+05]
    rg_range = [[0., 5.], [5., 15.], [15., 100.]]
    print(len(rg_range))


    sourcedir = np.genfromtxt('/projects/b1095/syr904/projects/GCE/catalog/path_allfinished_newruns_maingrid.dat', dtype=str)
    paths = list(sourcedir[:,0]); status = sourcedir[:,1]


    allsamps = [[],[],[],[],[]]
    rg_initial = []
    model_status = []
    time_tidal = []
    model_mass = []
    numgc_tot = 0; massgc_tot = 0; massmodel_tot = 0

    mass_init = np.zeros(bin_num); mass_density_init = np.zeros(bin_num)
    rg_randint = np.ones(bin_num)
    for ii in range(len(samps_mass)):
        check = 0
        
        if float(samps_mass[ii]) <= 1e5:
            N_star = samps_mass[ii]
            rv = str(int(choices([2,4])[0]))
            z = modelz[int(samps_z[ii])]
            rg_random_select = samps_den[ii]

            for xx in range(len(rg_range)):
                #if xx == 0: rg_randint = 1
                #elif xx == 1: rg_randint = 1#np.random.randint(2)
                #elif xx == 2: rg_randint = 1#np.random.randint(17)
                if rg_range[xx][0] < samps_den[ii] <= rg_range[xx][1]:
                    rg = modelrg[xx]

            if rv == '4' and rg == '20' and z == '0.002' and N_star == '1e5':
                continue

            thepath = '/projects/b1091/CMC_Grid_March2019/rundir/rv'+rv+'/rg'+rg+'/z'+z+'/'+N_star+'/'
            model_no = int(paths.index(thepath))
            thestatus = int(status[model_no])

            mass = modelmass[modeln.index(N_star)]
            rho_rh_0 = read_rho_rh(thepath)


        elif float(samps_mass[ii]) <= 1.6e6:
            N_star = samps_mass[ii]
            rv = modelrv[int(samps_rv[ii])]
            z = modelz[int(samps_z[ii])]
            rg_random_select = samps_den[ii]

            for xx in range(len(rg_range)):
                #if xx == 0: rg_randint = 1
                #elif xx == 1: rg_randint = 1#np.random.randint(2)
                #elif xx == 2: rg_randint = 1#np.random.randint(17)
                if rg_range[xx][0] < samps_den[ii] <= rg_range[xx][1]:
                    rg = modelrg[xx]

            thepath = '/projects/b1091/CMC_Grid_March2019/rundir/rv'+rv+'/rg'+rg+'/z'+z+'/'+N_star+'/'
            model_no = int(paths.index(thepath))
            thestatus = int(status[model_no])

            mass = modelmass[modeln.index(N_star)]
            rho_rh_0 = read_rho_rh(thepath)


        elif float(samps_mass[ii])==3.2e6:
            N_star = samps_mass[ii]
            rv = str(int(choices([1, 2])[0]))
            z = '0.02'
            #rg_random_select = random.uniform(15., 100.)
            rg = choices(['8', '20'])[0]
            if rg == '8': rg_random_select = random.uniform(5., 15.)
            if rg == '20': rg_random_select = random.uniform(15., 100.)

            ##reset rg back
            rg = '20'

            thepath = '/projects/b1091/CMC_Grid_March2019/rundir/rv'+rv+'/rg'+rg+'/z'+z+'/'+N_star+'/'
            model_no = int(paths.index(thepath))
            thestatus = int(status[model_no])

            mass = 1.93472e+06
            rho_rh_0 = read_rho_rh(thepath)

        elif float(samps_mass[ii])==8.6e6:
            N_star = samps_mass[ii]
            rv = '-100'
            z = '0.013'
            rg = choices(['8', '20'])[0]

            if rg == '8': rg_random_select = random.uniform(5., 15.)
            if rg == '20': rg_random_select = random.uniform(15., 100.)

            model_no = int(147)
            thestatus = int(2)

            mass = 5.32797e+06

            filerho = '/projects/b1095/syr904/projects/GCE/behemoth/rho_dyn_behemoth.txt'
            with open(filerho, 'r') as frho:
                for i, line in enumerate(frho):
                    if i == 1:   # 2th line
                        data = line.split()
                        rho_rh_0 = float(data[1])
                        break


        if dissolflag == 1 and thestatus !=1: continue

        for xx in range(len(rg_bin)-1):
            if rg_bin[xx] <= rg_random_select < rg_bin[xx+1] and rg_randint[xx] == 1:
                check = 1

        if check == 0: continue


        tid_time = t_tidal(mass, rg_random_select, circ_vel)
        rho_ga = galaxy_rho_r(rg_random_select, circ_vel)


        if rho_rh_0 > rho_ga:
            for yy in range(len(rg_bin)-1):
                if rg_bin[yy]<=rg_random_select<rg_bin[yy+1]:
                    mass_cumu_init = np.cumsum(mass_init)
                    if mass_cumu_init[yy] > cumu_mass_bin[yy]: ##or mass_density_init[yy] > rho_bin[yy]
                        rg_randint[yy] = 0

                    else:
                        #rg_bin_med = (rg_bin[yy+1]+rg_bin[yy])/2.
                        mass_init[yy]+=mass
                        #vol_init = (2.*twopi/3.)*(rg_bin[yy+1]**3-rg_bin[yy]**3)
                        #mass_density_init[yy] = mass_init[yy]/vol_init
                        mass_cumu_init = np.cumsum(mass_init)
                        if mass_cumu_init[yy] <= cumu_mass_bin[yy]:
                            allsamps[0].append(N_star); allsamps[1].append(rv); allsamps[2].append(rg); allsamps[3].append(z); allsamps[4].append(rho_rh_0)
                            rg_initial.append(rg_random_select)
                            model_status.append(thestatus)
                            time_tidal.append(tid_time)
                            model_mass.append(mass)
                            numgc_tot += 1
                            massgc_tot += mass; massmodel_tot += mass


        if massgc_tot > 2e9: break

        print(ii, thestatus)

    
    print(numgc_tot, massgc_tot, massmodel_tot)
    print(cumu_mass_bin)
    np.savetxt('/projects/b1095/syr904/projects/GCE/catalog/cluster_sample_initial_M_RG_dissol'+str(dissolflag)+'_fcl'+str(frac_gc_close)+'_ffa'+str(frac_gc_far)+'_xcut'+str(round(xcut, 1))+'_xmin'+str(xmin_samps)+'_massive_small.dat', np.c_[allsamps[0], allsamps[1], allsamps[2], allsamps[3], allsamps[4], rg_initial, model_status, time_tidal, model_mass], fmt = '%s %s %s %s %s %s %s %s %s', header = '1.N_star 2.RV 3.RG 4.Z 5.RHO_RH(Msun/pc^3) 6.RG_INITIAL(kpc) 7.Model_status 8.Ttid(Gyr) 9.Mass(Msun)', comments = '#')

    #return allsamps, rg_initial, model_status



##Calculate cluster inspiral
def cluster_inspiral(r_initial, cluster_path, rho_initial, cluster_status, td_time, bhm_flag):
    if bhm_flag == 0: 
        t_conv = dyn.conv('t', cluster_path+'initial.conv.sh')
        l_conv = dyn.conv('l', cluster_path+'initial.conv.sh')
        m_conv = dyn.conv('m', cluster_path+'initial.conv.sh')
        mstar_conv = dyn.conv('mstar', cluster_path+'initial.conv.sh')

        s=cluster_path.split('/')
        n_star=s[-2]
        z=s[-3][1:]
        rg=s[-4][2:]
        rv=s[-5][2:]

        data_rhodyn = np.genfromtxt('/projects/b1095/syr904/projects/GCE/catalog/data_rho_dyn/model_rv'+rv+'_rg'+rg+'_z'+z+'_'+n_star+'.txt')
        t_all = data_rhodyn[:,0]; rho_all = data_rhodyn[:,1]; 
        mass_all = data_rhodyn[:,2]; rh_all = data_rhodyn[:,3]

    else:
        t_conv = dyn.conv('t', '/projects/b1095/syr904/projects/GCE/behemoth/behemoth.conv.sh')
        l_conv = dyn.conv('l', '/projects/b1095/syr904/projects/GCE/behemoth/behemoth.conv.sh')
        m_conv = dyn.conv('m', '/projects/b1095/syr904/projects/GCE/behemoth/behemoth.conv.sh')
        mstar_conv = dyn.conv('mstar', '/projects/b1095/syr904/projects/GCE/behemoth/behemoth.conv.sh')

        data_rhodyn = np.genfromtxt('/projects/b1095/syr904/projects/GCE/behemoth/rho_dyn_behemoth.txt')
        t_all= data_rhodyn[:,0]; rho_all = data_rhodyn[:,1]; 
        mass_all = data_rhodyn[:,2]; rh_all = data_rhodyn[:,3]


    rg_disrupt = -100; t_disrupt = -100; rho_disrupt = -100


    data_gaprof = np.genfromtxt('/projects/b1095/syr904/projects/GCE/catalog/galaxy_profile.txt')
    xga = data_gaprof[:,0]; rho_ga = data_gaprof[:,3]


    ##Initialize the spiral in
    rg_old = r_initial; t_old = 0.; m_old = mass_all[0]; rho_halfm_old = rho_initial
    r0_tidal = r_tidal_initial(m_old, r_initial)
    #print(r0_tidal)

    ##Starting spiral in
    check = 0
    for ii in range(1, len(t_all)):
        t_curr = t_all[ii]; rhm_curr = rh_all[ii]/1000.
        rho_halfm_curr = rho_all[ii]; m_curr = mass_all[ii]

        delta_t = (t_curr-t_old)*t_conv/1000.
        rg_new = df(delta_t, rg_old, m_old)
        #print(rg_new)
        #td_time = t_tidal(m_old, rg_old, circ_vel)
        r_tidal_new = r_tidal_initial(m_curr, rg_new)
        #print(r_tidal_new)

        if rg_new < 0.001:#rhm_curr:
            check = 1
            print('cluster reached center', t_curr*t_conv, rg_new, rho_halfm_curr)
            rg_disrupt = rg_old; t_disrupt = t_old*t_conv; rho_disrupt = rho_halfm_old
            type_disrupt = 1
            break

        if t_curr*t_conv/1000. > td_time:
            check = 1
            print('cluster disrupted by tidal field', t_curr*t_conv, rg_new, rho_halfm_curr)
            rg_disrupt = rg_old; t_disrupt = t_old*t_conv; rho_disrupt = rho_halfm_old
            type_disrupt = 2
            break

        td_time = t_tidal(m_curr, rg_new, circ_vel)

        bin_right = xga[xga>=rg_new][0]
        index_right = np.where(xga == bin_right)[0][0]
        index_left = index_right - 1
        rho_ga_new_right = rho_ga[index_right]
        if rho_halfm_curr <= rho_ga_new_right:
            check = 1
            print('cluster disrupted', t_curr*t_conv, rg_new, rho_halfm_curr)
            rg_disrupt = rg_new; t_disrupt = t_curr*t_conv; rho_disrupt = rho_halfm_curr
            type_disrupt = 3
            break

        if r_tidal_new < 0.2*r0_tidal:
            check = 1
            print('cluster disrupted by tidal field', t_curr*t_conv, rg_new, rho_halfm_curr)
            rg_disrupt = rg_old; t_disrupt = t_old*t_conv; rho_disrupt = rho_halfm_old
            type_disrupt = 4
            break

        t_old = t_curr; m_old = m_curr; rg_old = rg_new; rho_halfm_old = rho_halfm_curr

        #print(i)
        #print(rho_halfm_curr, rho_ga_new, m_curr, rg_new)
    
    if check == 0:
        rg_disrupt = rg_new; t_disrupt = t_curr*t_conv; rho_disrupt = rho_halfm_curr
        type_disrupt = 5

    print(rg_disrupt, t_disrupt, rho_disrupt)

    return rg_disrupt, t_disrupt, rho_disrupt, type_disrupt


##Calculate cluster inspiral
def cluster_inspiral_old(r_initial, cluster_path, rho_initial, cluster_status, td_time):
    #datadyn = np.genfromtxt(cluster_path+'initial.dyn.dat')
    #datarho = np.genfromtxt(cluster_path+'initial.rho_lagrad.dat')
    #m_gc = np.array(datadyn[:,4])*m_conv
    #rho_halfm = np.array(datarho[:,20])*(mstar_conv/l_conv**3)
    #t_rho = np.array(datarho[:,0])*t_conv

    #for kk in range(len(t_rho)-1):
    #    delta_t = (t_rho[kk+1]-t_rho[kk])/1000.
    #    rg_new = df(delta_t, rg_old, m_gc[kk])

    #    rho_ga_new = galaxy_rho_r(rg_new, circ_vel)
    #    if rho_halfm[kk+1] <= rho_ga_new:
    #        print('cluster disrupted', t_rho[kk+1], rg_new, rho_halfm[kk+1])
    #        rg_disrupt = rg_new; t_disrupt = t_rho[kk+1]; rho_disrupt = rho_halfm[kk+1]
    #        break

    #    rg_old = rg_new

    t_conv = dyn.conv('t', cluster_path+'initial.conv.sh')
    l_conv = dyn.conv('l', cluster_path+'initial.conv.sh')
    m_conv = dyn.conv('m', cluster_path+'initial.conv.sh')
    mstar_conv = dyn.conv('mstar', cluster_path+'initial.conv.sh')

    s=cluster_path.split('/')
    n_star=s[-2]
    z=s[-3][1:]
    rg=s[-4][2:]
    rv=s[-5][2:]

    rg_disrupt = -100; t_disrupt = -100; rho_disrupt = -100

    #filerho = cluster_path+'initial.rho_lagrad.dat'
    #filedyn = cluster_path+'initial.dyn.dat'

    #with open(filerho, 'r') as fr:
    #    next(fr); next(fr)
    #    for i_rho, l_rho in enumerate(fr):
    #        pass
    #nline_rho = i_rho+1
    ##print(nline_rho)
    #line_cut = int(0.7*nline_rho)

    #with open(filedyn, 'r') as fd:
    #    next(fd); next(fd)
    #    for l_dyn in fd:
    #        data = l_dyn.split()
    #        m_old = float(data[4])
    #        break

    #print(m_old)

    data_rhodyn = np.genfromtxt('/projects/b1095/syr904/projects/GCE/catalog/data_rho_dyn/model_rv'+rv+'_rg'+rg+'_z'+z+'_'+n_star+'.txt')
    t_all = data_rhodyn[:,0]; rho_all = data_rhodyn[:,1]; 
    mass_all = data_rhodyn[:,2]; rh_all = data_rhodyn[:,3]


    ##Initialize the spiral in
    rg_old = r_initial; t_old = 0.; m_old = mass_all[0]; rho_halfm_old = rho_initial

    ##Starting spiral in
    check = 0
    for ii in range(1, len(t_all)):
        t_curr = t_all[ii]; rhm_curr = rh_all[ii]/1000.
        rho_halfm_curr = rho_all[ii]; m_curr = mass_all[ii]

        delta_t = (t_curr-t_old)*t_conv/1000.
        rg_new = df(delta_t, rg_old, m_old)
        #td_time = t_tidal(m_old, rg_old, circ_vel)

        if rg_new < rhm_curr:
            check = 1
            print('cluster reached center', t_curr*t_conv, rg_new, rho_halfm_curr)
            rg_disrupt = rg_old; t_disrupt = t_old*t_conv; rho_disrupt = rho_halfm_old
            type_disrupt = 1
            break

        if t_curr*t_conv/1000. > td_time:
            check = 1
            print('cluster disrupted by tidal field', t_curr*t_conv, rg_new, rho_halfm_curr)
            rg_disrupt = rg_old; t_disrupt = t_old*t_conv; rho_disrupt = rho_halfm_old
            type_disrupt = 2
            break

        rho_ga_new = galaxy_rho_r(rg_new, circ_vel)
        if rho_halfm_curr <= rho_ga_new:
            check = 1
            print('cluster disrupted', t_curr*t_conv, rg_new, rho_halfm_curr)
            rg_disrupt = rg_new; t_disrupt = t_curr*t_conv; rho_disrupt = rho_halfm_curr
            type_disrupt = 3
            break

        t_old = t_curr; m_old = m_curr; rg_old = rg_new; rho_halfm_old = rho_halfm_curr

        #print(i)
        #print(rho_halfm_curr, rho_ga_new, m_curr, rg_new)
    
    if check == 0:
        rg_disrupt = rg_new; t_disrupt = t_curr*t_conv; rho_disrupt = rho_halfm_curr
        type_disrupt = 4

    print(rg_disrupt, t_disrupt, rho_disrupt)

    return rg_disrupt, t_disrupt, rho_disrupt, type_disrupt


##Main function: Initialization and disruption time and position
def main(sample_num, start, end, initfile, disrupfile, file_version):
    savepath = '/projects/b1095/syr904/projects/GCE/catalog/'
    #gc_samps, gc_rg_init, gc_status = cluster_initialization(sample_num)
    data_samps = np.genfromtxt(savepath+initfile+'.dat', dtype = str)
    gc_samps = [[],[],[],[],[]]
    gc_samps[0] = data_samps[:,0]; gc_samps[1] = data_samps[:,1]; gc_samps[2] = data_samps[:,2]
    gc_samps[3] = data_samps[:,3]; gc_samps[4] = data_samps[:,4]
    gc_rg_init = data_samps[:,5]; gc_status = data_samps[:,6]
    gc_ttid = data_samps[:,7]

    
    gc_disrupt = [[],[],[],[]]
    f_disrt = open(savepath+disrupfile+file_version+'.dat', 'a+')
    f_disrt.write('#1.N_star 2.RV 3.RG 4.Z 5.RHO_RH(Msun/pc^3) 6.RG_INITIAL(kpc) 7.RG_DISRUP(kpc) 8.T_DISRUP(Myr) 9.RHO_DISRUP(Msun/pc^3) 10.TYPE_DISRUP 11.GC_STATUS 12.Ttid(Gyr)\n')

    print('start iteration')
    for jj in range(start, end): #len(gc_rg_init)
        if float(gc_samps[0][jj])<=3.2e6: 
            gc_path = '/projects/b1091/CMC_Grid_March2019/rundir/rv'+gc_samps[1][jj]+'/rg'+gc_samps[2][jj]+'/z'+gc_samps[3][jj]+'/'+gc_samps[0][jj]+'/'
            behemoth_flag = 0
        else:
            gc_path = '/projects/b1095/syr904/projects/GCE/behemoth/'
            behemoth_flag = 1

        rg_drpt, t_drpt, rho_drpt, type_drpt = cluster_inspiral(float(gc_rg_init[jj]), gc_path, float(gc_samps[4][jj]), gc_status[jj], float(gc_ttid[jj]), behemoth_flag)
        gc_disrupt[0].append(rg_drpt); gc_disrupt[1].append(t_drpt); gc_disrupt[2].append(rho_drpt)
        gc_disrupt[3].append(type_drpt)
        #if rg_drpt!=-100: drpt_flag = 1
        #else: drpt_flag = 0
        f_disrt.write('%s %s %s %s %s %s %s %s %s %s %s %s\n'%(gc_samps[0][jj], gc_samps[1][jj], gc_samps[2][jj], gc_samps[3][jj], gc_samps[4][jj], gc_rg_init[jj], gc_disrupt[0][jj-start], gc_disrupt[1][jj-start], gc_disrupt[2][jj-start], gc_disrupt[3][jj-start], gc_status[jj], gc_ttid[jj]))

        print(jj)
        
        #num_msp, mass_gc, lum_gamma = read_property_attime(gc_path, t_drpt, drpt_flag)
        #gc_disrupt[3].append(num_msp); gc_disrupt[4].append(mass_gc); gc_disrupt[5].append(lum_gamma)

    #np.savetxt('/projects/b1095/syr904/projects/GCE/catalog/cluster_sample_disrupt_M_RG.dat', np.c_[gc_samps[0], gc_samps[1], gc_samps[2], gc_samps[3], gc_samps[4], gc_rg_init, gc_disrupt[0], gc_disrupt[1], gc_disrupt[2], gc_disrupt[3], gc_status], fmt = '%s %s %s %s %s %s %s %s %s %s %s', header = '1.N_star 2.RV 3.RG 4.Z 5.RHO_RH(Msun/pc^3) 6.RG_INITIAL(kpc) 7.RG_DISRUP(kpc) 8.T_DISRUP(Myr) 9.RHO_DISRUP(Msun/pc^3) 10.TYPE_DISRUP 11.GC_STATUS', comments = '#')

    f_disrt.close()


##Analytical models following Fragione et al. 2018 but use the number of MSPs from the catalog models
def analytical_main(xmin_samps, xmax_samps, mmin_samps, mmax_samps, frac_gc_close, frac_gc_far, xcut, sample_num):
    ###################################################
    ##Initialization
    samps_dist = VN_sampling(smdf, xmin_samps, xmax_samps, smdf(xmin_samps), sample_num)
    samps_mass = VN_sampling(gc_mf, mmin_samps, mmax_samps, gc_mf(mmin_samps), sample_num)
    samps_rv = np.random.randint(4, size=sample_num)
    samps_z = np.random.randint(3, size=sample_num)

    #samps_dist = np.sort(samps_dist)

    bin_num = 400
    rg_bin = np.logspace(np.log10(xmin_samps), np.log10(xmax_samps), bin_num+1)
    cumu_mass_bin = []
    #rho_bin = []
    #rho_bin = smdf(rg_bin[:-1])
    check_xcut = 0
    for kk in range(len(rg_bin)-1):
        rg_med = (rg_bin[kk]+rg_bin[kk+1])/2.
        #rho_bin.append(smdf(rg_med)*f_gc)
        if rg_med < xcut:
            Mr = frac_gc_close*2*twopi*integrate.quad(lambda x: smdf(x)*x*x, 0, rg_med)[0]
            rg_med_prev = rg_med
        else:
            if check_xcut == 0:
                xcut = rg_med_prev
                Mr_cut = frac_gc_close*2*twopi*integrate.quad(lambda x: smdf(x)*x*x, 0, xcut)[0]
                check_xcut = 1
            Mr = Mr_cut + frac_gc_far*2*twopi*integrate.quad(lambda x: smdf(x)*x*x, xcut, rg_med)[0]

        cumu_mass_bin.append(Mr)
    
    print('finished sampling')

    allsamps = [[],[],[],[]]
    numgc_tot = 0; massgc_tot = 0

    mass_init = np.zeros(bin_num); mass_density_init = np.zeros(bin_num)
    rg_randint = np.ones(bin_num)
    for ii in range(len(samps_mass)):
        check = 0

        for xx in range(len(rg_bin)-1):
            if rg_bin[xx] <= samps_dist[ii]< rg_bin[xx+1] and rg_randint[xx] == 1:
                check = 1

        if check == 0: continue

        rho_ga = galaxy_rho_r(samps_dist[ii], circ_vel)
        if samps_mass[ii] <= 1e5:
            rho_rh_0 = 1000.  ##Msun/pc^3
        elif 1e5 < samps_mass[ii] < 1e6:
            rho_rh_0 = 1000.*(samps_mass[ii]/1e5)**2
        else:
            rho_rh_0 = 100000.

        tid_time = t_tidal(samps_mass[ii], samps_dist[ii], circ_vel)
     
        if rho_rh_0 > rho_ga:
            for yy in range(len(rg_bin)-1):
                if rg_bin[yy]<=samps_dist[ii]<rg_bin[yy+1]:
                    mass_cumu_init = np.cumsum(mass_init)
                    if mass_cumu_init[yy] > cumu_mass_bin[yy]: ##or mass_density_init[yy] > rho_bin[yy]
                        rg_randint[yy] = 0

                    else:
                        #rg_bin_med = (rg_bin[yy+1]+rg_bin[yy])/2.
                        mass_init[yy]+=samps_mass[ii]
                        #vol_init = (2.*twopi/3.)*(rg_bin[yy+1]**3-rg_bin[yy]**3)
                        #mass_density_init[yy] = mass_init[yy]/vol_init
                        mass_cumu_init = np.cumsum(mass_init)
                        if mass_cumu_init[yy] <= cumu_mass_bin[yy]:
                            allsamps[0].append(samps_mass[ii]); allsamps[1].append(samps_dist[ii]); allsamps[2].append(rho_rh_0); allsamps[3].append(tid_time)
                            numgc_tot += 1
                            massgc_tot += samps_mass[ii]


        if massgc_tot > 2e9: break

        print(ii)

    
    print(numgc_tot, massgc_tot)
    print(cumu_mass_bin)

    np.savetxt('/projects/b1095/syr904/projects/GCE/analytical_model/cluster_analytical_initial_M_RG_fcl'+str(frac_gc_close)+'_ffa'+str(frac_gc_far)+'_xcut'+str(round(xcut, 1))+'_xmin'+str(xmin_samps)+'_xmax'+str(xmax_samps)+'_mmin'+str(int(mmin_samps))+'_mmax'+str(int(mmax_samps))+'.dat', np.c_[allsamps[0], allsamps[1], allsamps[2], allsamps[3]], fmt = '%f %f %f %f', header = '1.Mass(Msun) 2.RG_INITIAL(kpc) 3.RHO_RH(Msun/pc^3) 4.Ttid(Gyr)', comments = '#')

    ###################################################
    ##Inspiral
    data_gaprof = np.genfromtxt('/projects/b1095/syr904/projects/GCE/catalog/galaxy_profile.txt')
    xga = data_gaprof[:,0]; rho_ga = data_gaprof[:,3]

    time_steps = np.linspace(0.0, 11.5, 3000) ##Gyr
    gc_disrupt = [[],[],[],[],[],[],[],[],[],[],[],[],[]]
    for kk in range(numgc_tot):
        ##Initialize the spiral in
        rg_old = allsamps[1][kk]; t_old = 0.; m_old = allsamps[0][kk]; rho_halfm_old = allsamps[2][kk]
        mcut = 0.1*m_old
        r0_tidal = r_tidal_initial(m_old, rg_old)
        rt_cut = 0.2*r0_tidal

        ##Starting spiral in
        for zz in range(1, len(time_steps)):
            check = 0
            delta_t = time_steps[zz]-t_old
            rg_new = df(delta_t, rg_old, m_old)
            m_new = dmdt(delta_t, m_old, rg_old)
            rt_new = r_tidal_initial(m_old, rg_old)

            if m_new <= 1e5:
                rho_halfm_new = 1000.  ##Msun/pc^3
            elif 1e5 < m_new < 1e6:
                rho_halfm_new = 1000.*(m_new/1e5)**2
            else:
                rho_halfm_new = 100000.


            if rg_new < 0.001:#rhm_curr:
                check = 1
                print('cluster reached center', time_steps[zz], rg_new, rho_halfm_new)
                rg_disrupt = rg_new; t_disrupt = time_steps[zz]; rho_disrupt = rho_halfm_new; m_disrupt = m_new
                type_disrupt = 1
                break

            if m_new < mcut:
                check = 1
                print('cluster mass negative', time_steps[zz], rg_new, rho_halfm_new)
                rg_disrupt = rg_new; t_disrupt = time_steps[zz]; rho_disrupt = rho_halfm_new;m_disrupt = m_new
                type_disrupt = 2
                break

            if rt_new < rt_cut:
                check = 1
                print('cluster mass negative', time_steps[zz], rg_new, rho_halfm_new)
                rg_disrupt = rg_new; t_disrupt = time_steps[zz]; rho_disrupt = rho_halfm_new;m_disrupt = m_new
                type_disrupt = 3
                break

            #print(rg_new, m_new)
            bin_right = xga[xga>=rg_new][0]
            index_right = np.where(xga == bin_right)[0][0]
            index_left = index_right - 1
            rho_ga_new_right = rho_ga[index_right]
            if rho_halfm_new <= rho_ga_new_right:
                check = 1
                print('cluster disrupted', time_steps[zz], rg_new, rho_halfm_new)
                rg_disrupt = rg_new; t_disrupt = time_steps[zz]; rho_disrupt = rho_halfm_new; m_disrupt = m_new
                type_disrupt = 4
                break

            t_old = time_steps[zz]; m_old = m_new; rg_old = rg_new; rho_halfm_old = rho_halfm_new

        if check == 0:
            rg_disrupt = rg_new; t_disrupt = time_steps[zz]; rho_disrupt = rho_halfm_new
            type_disrupt = 5; m_disrupt = m_new

        gc_disrupt[0].append(rg_disrupt); gc_disrupt[1].append(t_disrupt); gc_disrupt[2].append(rho_disrupt); gc_disrupt[3].append(type_disrupt); gc_disrupt[4].append(m_disrupt)

        #print(rg_disrupt, t_disrupt, rho_disrupt)


        #############################################
        ##Match cluster with model
        modeln = ['5e4', '1e5', '2e5', '4e5', '8e5', '1.6e6', '3.2e6', '8.6e6']
        modelrv = ['0.5', '1', '2', '4']; modelrg = ['2', '8', '20']; modelz = ['0.0002', '0.002', '0.02']
        modelmass = [29655.3, 59638.8, 1.197630e+05, 2.423500e+05, 4.848440e+05, 9.703820e+05, 1.93472e+06, 5.32797e+06]
        rg_range = [[0., 5.], [5., 15.], [15., 100.]]

        
        behemoth_flag = 0; small_mass = 0
        for xx in range(len(rg_range)):
                if rg_range[xx][0] <= allsamps[1][kk] < rg_range[xx][1]:
                    rg = modelrg[xx]

        rv = modelrv[int(choices([0,1,2,3])[0])]
        z = modelz[int(choices([0,1,2])[0])]

        
        if allsamps[0][kk] <= (modelmass[0]+modelmass[1])/2:
            rv = str(int(choices([2,4])[0]))
            n_star=modeln[0]

            small_mass = 1

            #if (rv == '4' and rg == '2' and z == '0.002') or (rv == '2' and rg == '2' and z == '0.02'):
            #    z = choices(modelz)[0]
            #    rv = str(int(choices([2,4])[0]))

        elif allsamps[0][kk] <= (modelmass[1]+modelmass[2])/2:
            rv = str(int(choices([2,4])[0]))
            n_star=modeln[1]

            small_mass = 1
            
            #if rv == '4' and rg == '20' and z == '0.002':
            #    z = choices(modelz)[0]
            #    rv = str(int(choices([2,4])[0]))


        elif allsamps[0][kk] <= (modelmass[2]+modelmass[3])/2:
            n_star=modeln[2]

        elif allsamps[0][kk] <= (modelmass[3]+modelmass[4])/2:
            n_star=modeln[3]

        elif allsamps[0][kk] <= (modelmass[4]+modelmass[5])/2:
            n_star=modeln[4]

        elif allsamps[0][kk] <= (modelmass[5]+modelmass[6])/2:
            n_star=modeln[5]

        elif allsamps[0][kk] <= (modelmass[6]+modelmass[7])/2:
            z = '0.02'
            rg = '20'
            rv = str(int(choices([1,2])[0]))
            n_star=modeln[6]
        else:
            behemoth_flag = 1


        if behemoth_flag == 0:
            if small_mass == 0:
                file_name = 'model_rv'+rv+'_rg'+rg+'_z'+z+'_'+n_star+'.txt'
                print(file_name)
                data_rhodyn = np.genfromtxt('/projects/b1095/syr904/projects/GCE/catalog/data_rho_dyn/'+file_name)

                thepath = '/projects/b1091/CMC_Grid_March2019/rundir/rv'+rv+'/rg'+rg+'/z'+z+'/'+n_star+'/'
                t_conv = dyn.conv('t', thepath+'initial.conv.sh')
                l_conv = dyn.conv('l', thepath+'initial.conv.sh')
                m_conv = dyn.conv('m', thepath+'initial.conv.sh')
                mstar_conv = dyn.conv('mstar', thepath+'initial.conv.sh')

                last_time = 0
                mass_index = np.where(data_rhodyn[:,2]<=m_new)[0]
                print(mass_index)
                if len(mass_index) == 0:
                    last_time = 1
                    t_model_disrupt = data_rhodyn[:,0][-1]*t_conv
                    #lgamma = data_lgamma[:,2][-1]
                else:
                    themass_index = mass_index[0]
                    t_model_disrupt = data_rhodyn[:,0][themass_index-1]*t_conv
                    print(t_model_disrupt)
                    #time_index = np.where(data_lgamma[:,1]>t_model_disrupt)[0]
                    #print(time_index)
                    #lgamma = data_lgamma[:,2][time_index-1]
                    #print(lgamma)
                
                nstime_old = 0; msp_disrupt = 0; nmsp_old = 0
                with open(thepath+'initial.ns.dat', 'r') as fns:
                    next(fns)
                    for line in fns:
                        datans = line.split()
                        if nstime_old>t_model_disrupt:
                            msp_disrupt = nmsp_old
                            break

                        nstime_old = float(datans[0])*t_conv; nmsp_old = int(datans[6])

                ##Extract Lgamma
                if msp_disrupt != 0:
                    data_lgamma = np.genfromtxt('/projects/b1095/syr904/projects/GCE/catalog/data_lgamma/'+file_name)
                    if len(data_lgamma.shape) != 1:
                        if len(mass_index) == 0:
                            lgamma = data_lgamma[:,2][-1]
                        else:
                            time_index = np.where(data_lgamma[:,1]>t_model_disrupt)[0]
                            print(time_index)
                            if len(time_index) == 0:
                                lgamma = data_lgamma[:,2][-1]
                                print(lgamma)
                            else:
                                lgamma = data_lgamma[:,2][time_index[0]-1]
                                print(lgamma)
                    else:
                        lgamma = data_lgamma[2]

                else:
                    lgamma = 0

            else:
                t_model_disrupt = -100
                lgamma = 3.5e34; msp_disrupt = 0.03
                last_time = -100
        
        else:
            data_rhodyn = np.genfromtxt('/projects/b1095/syr904/projects/GCE/behemoth/rho_dyn_behemoth.txt')
            data_lgamma = np.genfromtxt('/projects/b1095/syr904/projects/GCE/behemoth/Lgamma_alltime_behemoth.dat')

            t_conv = dyn.conv('t', '/projects/b1095/syr904/projects/GCE/behemoth/behemoth.conv.sh')
            l_conv = dyn.conv('l', '/projects/b1095/syr904/projects/GCE/behemoth/behemoth.conv.sh')
            m_conv = dyn.conv('m', '/projects/b1095/syr904/projects/GCE/behemoth/behemoth.conv.sh')
            mstar_conv = dyn.conv('mstar', '/projects/b1095/syr904/projects/GCE/behemoth/behemoth.conv.sh')

            last_time = 0
            mass_index = np.where(data_rhodyn[:,2]<=m_disrupt)[0]
            if len(mass_index) == 0:
                last_time = 1
                t_model_disrupt = data_rhodyn[:,0][-1]*t_conv
                lgamma = data_lgamma[:,2][-1]
            else:
                mass_index = mass_index[0]
                t_model_disrupt = data_rhodyn[:,0][mass_index-1]*t_conv
                time_index = np.where(data_lgamma[:,0]>t_model_disrupt)[0]
                lgamma = data_lgamma[:,2][time_index[0]-1]

            nstime_old = 0
            with open('/projects/b1095/syr904/projects/GCE/behemoth/behemoth.ns.dat', 'r') as fns:
                next(fns)
                for line in fns:
                    datans = line.split()
                    if nstime_old>t_model_disrupt:
                        msp_disrupt = nmsp_old
                        break

                    nstime_old = float(datans[0])*t_conv; nmsp_old = int(datans[6])


        gc_disrupt[5].append(lgamma); gc_disrupt[6].append(msp_disrupt); gc_disrupt[7].append(t_model_disrupt); gc_disrupt[8].append(last_time)
        gc_disrupt[9].append(float(n_star)); gc_disrupt[10].append(float(rv)); gc_disrupt[11].append(float(rg)); gc_disrupt[12].append(float(z)) 


    np.savetxt('/projects/b1095/syr904/projects/GCE/analytical_model/semi_analytic_model_fcl'+str(frac_gc_close)+'_ffa'+str(frac_gc_far)+'_xcut'+str(round(xcut, 1))+'_xmin'+str(xmin_samps)+'_xmax'+str(xmax_samps)+'_mmin'+str(int(mmin_samps))+'_mmax'+str(int(mmax_samps))+'.txt', np.c_[allsamps[0], allsamps[1], allsamps[2], allsamps[3], gc_disrupt[9], gc_disrupt[10], gc_disrupt[11], gc_disrupt[12], gc_disrupt[0], gc_disrupt[1], gc_disrupt[2], gc_disrupt[4], gc_disrupt[3], gc_disrupt[5], gc_disrupt[6], gc_disrupt[7], gc_disrupt[8]], fmt = '%f %f %f %f %e %f %f %f %f %f %f %f %d %e %f %f %d', header = '1.M_init(Msun) 2.Rgc_init(kpc) 3.Rho_rh_init(Msun/pc^3) 4.t_disrupt_init(Gyr) 5.N_star 6.RV 7.RG 8.Z 9.Rgc_distupt(kpc) 10.T_disrupt(Gyr) 11.Rho_rh_disrupt(Msun/pc^3), 12.M_disrupt(Msun) 13.Type_disrupt 14.Lgamma 15.MSP 16.T_model_disrupt(Myr), 17.Last_time_flag', comments = '#', delimiter = ' ')




##Running the functions
#read_property_all()
#cluster_initialization(8000)
#main(5000, 0, 200)
#cluster_inspiral(0.1, '/projects/b1091/CMC_Grid_March2019/rundir/rv4/rg20/z0.02/8e5/')

#df(0.2, 0.1, 1e5)
#galaxy_rho_r(0.1, circ_vel)

#plot_raw_distr(gdf, 0.001, 40, 'Galaxy Mass Density Distribution')
#plot_raw_distr(gc_mf, 1e4, 1e7, 'Cluster Mass Function')
#plot_sampling_distr(gc_mf, 1e4, 1e7, 8000, 300)
#plot_sampling_distr(gdf, 0.001, 40, 8*1e11, 300)

##Checking plots
#samps_x, samps_y, rejsamps_x, rejsamps_y = VN_sampling(gc_mf, 1e4, 1e7, 5*1e3, 500)
#xs = np.linspace(1e4, 1e7, 1000)
#ys = gc_mf(xs)
#
#plt.plot(xs, ys, label="Function")
#plt.scatter(samps_x, samps_y, color='green', s=5)
#plt.scatter(rejsamps_x, rejsamps_y, color='red', s=5)
##plt.xscale('log')
##plt.yscale('log')
#plt.show()


##Plot velocity distribution
#xs = np.linspace(0.001, 40, 1000)
#ys = []
#for ii in range(len(xs)):
#    ys.append(circ_vel(xs[ii], gdf, NFW_DM))
#
#plt.plot(xs, ys, label='Velocity Distribution') 
#plt.fill_between(xs, ys, 0, alpha=0.2)
##plt.xlim(0, 30)
#plt.xscale('log')
#plt.yscale('log')
#plt.xlabel("r(kpc)")
#plt.ylabel(r'$circular\ velocity\ (km/s)$')
#plt.legend()
#plt.show()


