import numpy as np
import scipy
from scipy.special import gamma, gammainc
import scipy.integrate as integrate
import math
import pandas as pd
import matplotlib.pyplot as plt
from random import choices

import dynamics as dyn

twopi = 2*np.pi
Gconst = 4.30091*10**(-3)   #pc*Msun^-1*(km/s)^2
H = 69.6 #km/s*Mpc^-1

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
def gc_mf(xm, constn = 8*1e7, beta=2):
    return constn*xm**(-(beta-1))


##Generate random initial model mass
def gc_mf_frommodel():
    init_n = ['2e5', '4e5', '8e5', '1.6e6']
    init_mass = [1.197630e+05, 2.423500e+05, 4.848440e+05, 9.703820e+05]
    inverse_m2 = [1/x**2 for x in init_mass]
    sum_weight = sum(inverse_m2)
    weights = [y/sum_weight for y in inverse_m2]
    #print(weights)
    return choices(init_n, weights)


##Dynamical friction
def df(dt, r_old, Mgc):  ##t_df in Gyr; r_old in kpc
    t_df = 0.23*circ_vel(r_old, gdf)*r_old**2*(1e5/Mgc)
    dr = -r_old*dt/(2*t_df)

    r_new = r_old+dr

    #print(r_new)
    return r_new


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
def VN_sampling(function, xmin, xmax, fmax, nums, fmin=0, batch=100):
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
    xs = np.linspace(xmin, xmax, 50)
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

    xs = np.linspace(xmin_sample, xmax_sample, 50)
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
def read_property_attime(modelpath, time):   ##time in Myr
    Num_msp = -100; Mass = -100; Lum_gamma = -100
    Num_msp_old = -100; Mass_old = -100; Lum_gamma_old = -100

    t_conv = dyn.conv('t', modelpath+'initial.conv.sh')
    m_conv = dyn.conv('m', modelpath+'initial.conv.sh')

    ##Extract number of MSP
    with open(modelpath+'initial.ns.dat', 'r') as fns:
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
    datadyn = np.genfromtxt(modelpath+'initial.dyn.dat')
    with open(modelpath+'initial.dyn.dat', 'r') as fdyn:
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
        sourcedir = np.genfromtxt('/projects/b1095/syr904/projects/GCE/catalog/path_allfinished_newruns_maingrid.dat', dtype=str)
        paths = sourcedir[:,0]
        model_no = list(paths).index(modelpath)

        datagama = np.genfromtxt('/projects/b1095/syr904/projects/GCE/catalog/Lgamma_alltime_maingrid.dat')
        model_gama = datagama[:,0]; t_gama = datagama[:,1]; lgama = datagama[:,2]

        t_model = t_gama[model_gama==model_no]; lum_model = lgama[model_gama==model_no]
        for xx in range(len(t_model)):
            if round(t_model[xx],6) >= round(time,6):
                Lum_gamma = lum_model[xx]
                Lum_gamma_old = lum_model[xx-1]

    #print(Num_msp, Mass, Lum_gamma, Num_msp_old, Mass_old, Lum_gamma_old)
    return Num_msp, Mass, Lum_gamma, Num_msp_old, Mass_old, Lum_gamma_old


####Find number of MSP, gamma-ray luminosity and cluster mass at a time for all models
def read_property_all():
    sample_disrp = np.genfromtxt('/projects/b1095/syr904/projects/GCE/catalog/cluster_sample_disrupt_M_RG.dat', dtype = str)
    nstar_disrp = sample_disrp[:,0]; rv_disrp = sample_disrp[:,1]
    rg_model_disrp = sample_disrp[:,2]; z_disrp = sample_disrp[:,3]
    rg_init = sample_disrp[:,5].astype(np.float)
    rg_disrp = sample_disrp[:,6].astype(np.float)
    t_disrp = sample_disrp[:,7].astype(np.float)
    print(np.where(t_disrp == 0))

    prop = [[],[],[]]; prop_old = [[],[],[]]
    for ii in range(len(nstar_disrp)):
        thepath = '/projects/b1091/CMC_Grid_March2019/rundir/rv'+rv_disrp[ii]+'/rg'+rg_model_disrp[ii]+'/z'+z_disrp[ii]+'/'+nstar_disrp[ii]+'/'
        Nmsp, Mgc, Lmsp, Nmsp_old, Mgc_old, Lmsp_old = read_property_attime(thepath, t_disrp[ii])
        if Mgc == 0 or Mgc_old == 0:
            print(modelpath)
        prop[0].append(Nmsp); prop[1].append(Mgc); prop[2].append(Lmsp)
        prop_old[0].append(Nmsp_old); prop_old[1].append(Mgc_old); prop_old[2].append(Lmsp_old)
        print(ii)

    print('finished searching')

    np.savetxt('/projects/b1095/syr904/projects/GCE/catalog/cluster_sample_property_M_RG.dat', np.c_[nstar_disrp, rv_disrp, rg_model_disrp, z_disrp, t_disrp, prop[0], prop[1], prop[2], prop_old[0], prop_old[1], prop_old[2]], fmt = '%s %s %s %s %s %s %s %s %s %s %s', header = '1.N_star 2.RV 3.RG 4.Z 5.T_DISRUP(Myr) 6.Nmsp 7.Mgc(Msun) 8.Lmsp(erg/s) 9.Nmsp_old 10.Mgc_old(Msun) 11.Lmsp_old(erg/s)', comments = '#')



##Sampling initial cluster distribution
def cluster_initialization(sample_num):
    samps_den = VN_sampling(gdf, 0.001, 40., 8*1e11, sample_num)
    samps_mass = random_sampling(gc_mf_frommodel, sample_num)
    #samps_rg = random_sampling(gdf_frommodel, sample_num)
    samps_rv = np.random.randint(4, size=sample_num)
    samps_z = np.random.randint(3, size=sample_num)

    print('finished sampling')

    modeln = ['2e5', '4e5', '8e5', '1.6e6']
    modelrv = ['0.5', '1', '2', '4']; modelrg = ['2', '8', '20']; modelz = ['0.0002', '0.002', '0.02']
    modelmass = [1.197630e+05, 2.423500e+05, 4.848440e+05, 9.703820e+05]
    rg_range = [[0., 5.], [5., 15.], [15., 40.]]
    print(len(rg_range))

    allsamps = [[],[],[],[],[]]
    rg_initial = []
    model_status = []
    numgc_tot = 0; massgc_tot = 0; massmodel_tot = 0


    for ii in range(len(samps_mass)):
        check = 0
        for xx in range(len(rg_range)):
            if xx == 0: rg_randint = 1
            elif xx == 1: rg_randint = np.random.randint(2)
            elif xx == 2: rg_randint = np.random.randint(5)
            if rg_range[xx][0] <= samps_den[ii] <= rg_range[xx][1] and rg_randint == 1:
                check = 1
                rg = modelrg[xx]

        if check == 0: continue

        N_star = samps_mass[ii]
        rv = modelrv[int(samps_rv[ii])]
        z = modelz[int(samps_z[ii])]


        mass = modelmass[modeln.index(N_star)]


        thepath = '/projects/b1091/CMC_Grid_March2019/rundir/rv'+rv+'/rg'+rg+'/z'+z+'/'+N_star+'/'
        rho_rh_0 = read_rho_rh(thepath)

        rho_ga = galaxy_rho_r(samps_den[ii], circ_vel)

        sourcedir = np.genfromtxt('/projects/b1095/syr904/projects/GCE/catalog/path_allfinished_newruns_maingrid.dat', dtype=str)
        paths = list(sourcedir[:,0]); status = sourcedir[:,1]
        model_no = int(paths.index(thepath))
        thestatus = int(status[model_no])
        #print(thestatus)

        if rho_rh_0 > rho_ga:
            allsamps[0].append(N_star); allsamps[1].append(rv); allsamps[2].append(rg); allsamps[3].append(z); allsamps[4].append(rho_rh_0)
            rg_initial.append(samps_den[ii])
            model_status.append(thestatus)
            numgc_tot += 1
            massgc_tot += mass; massmodel_tot += mass

        if massgc_tot > 3e8: break

        print(ii, thestatus)

    
    print(numgc_tot, massgc_tot, massmodel_tot)
    np.savetxt('/projects/b1095/syr904/projects/GCE/catalog/cluster_sample_initial_M_RG.dat', np.c_[allsamps[0], allsamps[1], allsamps[2], allsamps[3], allsamps[4], rg_initial, model_status], fmt = '%s %s %s %s %s %s %s', header = '1.N_star 2.RV 3.RG 4.Z 5.RHO_RH(Msun/pc^3) 6.RG_INITIAL(kpc) 7.Model_status', comments = '#')

    return allsamps, rg_initial, model_status


##Calculate cluster inspiral
def cluster_inspiral(r_initial, cluster_path, rho_initial, cluster_status):
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

    rg_disrupt = -100; t_disrupt = -100; rho_disrupt = -100

    filerho = cluster_path+'initial.rho_lagrad.dat'
    filedyn = cluster_path+'initial.dyn.dat'

    with open(filerho, 'r') as fr:
        next(fr); next(fr)
        for i_rho, l_rho in enumerate(fr):
            pass
    nline_rho = i_rho+1
    #print(nline_rho)
    line_cut = int(0.7*nline_rho)

    with open(filedyn, 'r') as fd:
        next(fd); next(fd)
        for l_dyn in fd:
            data = l_dyn.split()
            m_old = float(data[4])
            break

    print(m_old)

    ##Initialize the spiral in
    rg_old = r_initial; t_old = 0.; m_old = m_old*m_conv; rho_halfm_old = rho_initial

    ##Starting spiral in
    check = 0
    with open(filerho, 'r') as frho, open(filedyn, 'r') as fdyn:
        next(frho); next(frho)
        next(fdyn); next(fdyn)
        for i, lines in enumerate(zip(frho, fdyn)):
            if i==0: continue
            if i <= line_cut and i % 10 != 0: 
                continue
            datarho = lines[0].split()
            datadyn = lines[1].split()
            t_curr = float(datarho[0])
            rho_halfm_curr = float(datarho[20])*(mstar_conv/l_conv**3)

            m_curr = float(datadyn[4])*m_conv

            delta_t = (t_curr-t_old)*t_conv/1000.
            rg_new = df(delta_t, rg_old, m_old)

            if rg_new <= 0:
                check = 1
                print('cluster reached center', t_curr*t_conv, rg_new, rho_halfm_curr)
                rg_disrupt = rg_old; t_disrupt = t_old*t_conv; rho_disrupt = rho_halfm_old
                type_disrupt = 1
                break

            rho_ga_new = galaxy_rho_r(rg_new, circ_vel)
            if (rho_halfm_curr <= rho_ga_new):
                check = 1
                print('cluster disrupted', t_curr*t_conv, rg_new, rho_halfm_curr)
                rg_disrupt = rg_new; t_disrupt = t_curr*t_conv; rho_disrupt = rho_halfm_curr
                type_disrupt = 2
                break

            t_old = t_curr; m_old = m_curr; rg_old = rg_new; rho_halfm_old = rho_halfm_curr

            #print(i)
            #print(rho_halfm_curr, rho_ga_new, m_curr, rg_new)
    
    if check == 0:
        rg_disrupt = rg_new; t_disrupt = t_curr*t_conv; rho_disrupt = rho_halfm_curr
        type_disrupt = 3

    print(rg_disrupt, t_disrupt, rho_disrupt)

    return rg_disrupt, t_disrupt, rho_disrupt, type_disrupt


##Mian function: Initialization and disruption time and position
def main(sample_numm, start, end):
    #gc_samps, gc_rg_init, gc_status = cluster_initialization(sample_num)
    data_samps = np.genfromtxt('/projects/b1095/syr904/projects/GCE/catalog/cluster_sample_initial_M_RG.dat', dtype = str)
    gc_samps = [[],[],[],[],[]]
    gc_samps[0] = data_samps[:,0]; gc_samps[1] = data_samps[:,1]; gc_samps[2] = data_samps[:,2]
    gc_samps[3] = data_samps[:,3]; gc_samps[4] = data_samps[:,4]
    gc_rg_init = data_samps[:,5]; gc_status = data_samps[:,6]

    
    gc_disrupt = [[],[],[],[]]
    f_disrt = open('/projects/b1095/syr904/projects/GCE/catalog/cluster_sample_disrupt_M_RG_v2.dat', 'a+')
    f_disrt.write('#1.N_star 2.RV 3.RG 4.Z 5.RHO_RH(Msun/pc^3) 6.RG_INITIAL(kpc) 7.RG_DISRUP(kpc) 8.T_DISRUP(Myr) 9.RHO_DISRUP(Msun/pc^3) 10.TYPE_DISRUP 11.GC_STATUS\n')
    for jj in range(start, end): #len(gc_rg_init)
        gc_path = '/projects/b1091/CMC_Grid_March2019/rundir/rv'+gc_samps[1][jj]+'/rg'+gc_samps[2][jj]+'/z'+gc_samps[3][jj]+'/'+gc_samps[0][jj]+'/'
        rg_drpt, t_drpt, rho_drpt, type_drpt = cluster_inspiral(float(gc_rg_init[jj]), gc_path, float(gc_samps[4][jj]), gc_status[jj])
        gc_disrupt[0].append(rg_drpt); gc_disrupt[1].append(t_drpt); gc_disrupt[2].append(rho_drpt)
        gc_disrupt[3].append(type_drpt)
        #if rg_drpt!=-100: drpt_flag = 1
        #else: drpt_flag = 0
        f_disrt.write('%s %s %s %s %s %s %s %s %s %s %s\n'%(gc_samps[0][jj], gc_samps[1][jj], gc_samps[2][jj], gc_samps[3][jj], gc_samps[4][jj], gc_rg_init[jj], gc_disrupt[0][jj-200], gc_disrupt[1][jj-200], gc_disrupt[2][jj-200], gc_disrupt[3][jj-200], gc_status[jj]))

        print(jj)
        
        #num_msp, mass_gc, lum_gamma = read_property_attime(gc_path, t_drpt, drpt_flag)
        #gc_disrupt[3].append(num_msp); gc_disrupt[4].append(mass_gc); gc_disrupt[5].append(lum_gamma)

    #np.savetxt('/projects/b1095/syr904/projects/GCE/catalog/cluster_sample_disrupt_M_RG.dat', np.c_[gc_samps[0], gc_samps[1], gc_samps[2], gc_samps[3], gc_samps[4], gc_rg_init, gc_disrupt[0], gc_disrupt[1], gc_disrupt[2], gc_disrupt[3], gc_status], fmt = '%s %s %s %s %s %s %s %s %s %s %s', header = '1.N_star 2.RV 3.RG 4.Z 5.RHO_RH(Msun/pc^3) 6.RG_INITIAL(kpc) 7.RG_DISRUP(kpc) 8.T_DISRUP(Myr) 9.RHO_DISRUP(Msun/pc^3) 10.TYPE_DISRUP 11.GC_STATUS', comments = '#')

    f_disrt.close()



##Running the functions
read_property_all()
#cluster_initialization(5000)
#main(5000, 200, 400)
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


