import numpy as np

#import jax.numpy as np

#from jax.scipy.special import gammaincc as gammainc
#from jax.scipy.special import gammaln as gammaln
from scipy.special import gammaincc
from scipy.special import gammaln
from scipy.special import gamma

import astropy.units as u
from astropy.cosmology import Planck15

from statsmodels.stats.weightstats import DescrStatsW

#cosmology helper functions to convert between redshift and lookback time
#zmax is the maximum redshift of star formation 
zmax = 20
zs_i = np.linspace(0, zmax, 1000)
tLs_i = Planck15.lookback_time(zs_i).to(u.Gyr).value
tL_at_z_interp = lambda z: np.interp(z, zs_i, tLs_i)
z_at_tL_interp = lambda t: np.interp(t, tLs_i, zs_i)

#simulated grid parameters
zeta_grid = np.array([0.0002, 0.002, 0.02]) #metallicities (0.01, 0.1, 1 Zsun)
rv_grid = np.array([0.5, 1, 2, 4]) #virial radii in pc
ncl_grid = np.array([2e5, 4e5, 8e5, 1.6e6]) #number of particles. Stellar mass is 0.6 Msun * ncl.

#helper function because jax doesn't have a gamma function defined
def jax_gamma(x):
    #return np.exp(gammaln(x))
    return gamma(x)

def chi_eff(m1, m2, s1, s2, alpha, beta):
    return (m1*s1*np.cos(alpha)+m2*s2*np.cos(beta))/(m1+m2)

def schechter_lower_int(beta, logMstar, logMlo):
    '''
    inputs: power law slope beta, log10 Schechter mass Mstar, log10 minimum integration bound Mlo
    returns the integral M^beta exp(-M/Mstar) dM from Mlo to infinity
    '''
    #change of variables x = M/Mstar 
    #M = x*Mstar, dx = dM/Mstar, dM = dx * Mstar, xlo = Mlo/ Mstar
    # Mstar^(beta + 1) integral [x^beta exp(-x) dx] from xlo to infinity
    lnMstar = logMstar * np.log(10)
    lnMlo = logMlo * np.log(10)
    xlow = np.exp(lnMlo - lnMstar)
    ln_out = (beta + 1) * lnMstar + np.log(gammainc(beta + 1, xlow)) + gammaln(beta + 1) #this last term is because we don't want the normalized version
    return np.exp(ln_out)
    #return gammainc(beta + 1, xlow) #unnormalized

def mean_log10metallicity(z):
    '''Returns the mean log10Z as a function of z
    Assumption is that star-forming gas in GCs has the same metallicity as the rest of the galaxy
    From Madau & Fragos 2017'''
     #this predicts fairly high metallicities
     #need to go beyond z = 7 to get mean metallicity below 0.1 solar, even though GCs in MW have metallicities below 0.1 solar... 
    #probably doesn't matter too much though because metallicity doesn't seem to affect cluster rates as much as other things 
    return 0.153 - 0.074 * z ** 1.34

def metallicity_weights(metals, redshift, sigma_dex = 0.5, Zsun = 0.02):
    '''
    metals: metallicity
    redshift: formation redshift
    sigma_dex: scatter in log10Z
    Returns fraction of star formation in a given metallicity bin at a given redshift
    Assumes the metallicity distribution at each redshift is lognormal, truncated between maximum and minimum simulated metallicity
    assumes metallicity bins are log-spaced
    '''

    log10mean = mean_log10metallicity(redshift) #an array if redshift is an array
   
    x = np.log10(metals/Zsun) 
    x_grid = np.log10(zeta_grid/Zsun)
    
    w = np.exp(-(x - log10mean)**2/(2*sigma_dex**2))
    
    w_grid = np.array([np.exp(-(xg - log10mean)**2/(2*sigma_dex**2)) for xg in x_grid]) #needs to be normalized at every redshift
    
    norm = np.sum(w_grid, axis = 0)

    return w/norm

def mass_weights_powerlaw(cluster_mass, beta = -2, missing_cluster_factor = 4.0):
    '''
    assume cluster mass distribution is a power law with slope beta
    note that Kremer+ 2020 assumes it is lognormal with mean log10M = 5.54 (approximately center of simulated range) and width sigma(log10M) = 0.52
    missing_cluster_factor: contribution from the clusters too big to model directly. Kremer+ 2020 find that this gives a factor of 4 regardless of radius distribution, but they assume a mass distribution much more skewed to heavy systems.
    '''
    w = cluster_mass**(beta + 1) #must take into account that cluster mass is log-spaced, this is dM/dlogM 
    w_grid = (0.6*ncl_grid)**(beta + 1)
    norm = np.sum(w_grid)
    
    return w/norm * missing_cluster_factor

def mass_weights_schechter(cluster_mass, beta = -2, logMstar0 = 6.26):
    '''
    Following section II.B of Antonini & Gieles 2020, this is the *initial* cluster mass function
    not to be confused with present (evolved) MW cluster mass function 
    cluster_mass: initial cluster mass
    beta: power law slope
    logMstar0: in log10(Msun), initial Schechter mass, 2Mc from Antonini & Gieles 2020
    '''
    x = cluster_mass/ 10**logMstar0
    w = x**(beta + 1) * np.exp(-x)
    
    x_grid = 0.6 * ncl_grid/ 10**logMstar0
    w_grid = x_grid**(beta + 1) * np.exp(-x_grid)
    norm = np.sum(w_grid)  
    
    return w/norm

def compute_missing_cluster_factor(beta = -2, logMstar0 = 6.26, logMlo = 2, logMhi = 8, res = 100):
    '''
    beta: power law slope
    logMstar0: in log10(Msun), initial Schechter mass, 2Mc from Antonini & Gieles 2020
    logMlo: in log10(Msun), minimum initial cluster mass
    logMhi: in log10(Msun), maximum initial cluster mass
    returns: factor by which to multiply BBH merger rate to account for cluster masses not simulated. This corresponds to the average number of mergers over the simulated mass range, divided by the average number of mergers over the full mass range from logMlo to logMhi. Note it is usually smaller than 1!
    '''
    
    #asumption is that number of mergers as a function of cluster mass [for a fixed radius] scales as M^1.6 (from Antonini & Gieles 2020)
    
    x_grid_full = np.logspace(logMlo-logMstar0, logMhi-logMstar0, res) #log spaced bins between 100 and 10^8 Msun
    w_grid_full = x_grid_full**(beta + 1) * np.exp(-x_grid_full) #cluster weight according to mass distribution (not normalized)
    norm_full = np.sum(w_grid_full)
    
    average_merge_full = np.sum(x_grid_full**1.6 * w_grid_full/norm_full) #weighted sum of (m/Mstar)**1.6, corresponding to average number of mergers per cluster over the full mass range
    
    x_grid = 0.6 * ncl_grid/ 10**logMstar0
    w_grid = x_grid**(beta + 1) * np.exp(-x_grid)
    norm = np.sum(w_grid)
    average_merge_sim = np.sum(x_grid**1.6 * w_grid/norm) #weighted sum of (m/Mstar)**1.6 in the simulated mass range, corresponding to average number of mergers per simulated cluster

    missing_cluster_factor = average_merge_full/average_merge_sim
    
    return missing_cluster_factor

def compute_disrupted_cluster_factor(beta = -2, logMstar0 = 6.26, logMlo = 2, logMhi = 8, logDelta = 5.33):
    '''
    beta: power law slope of Schecther function describing GC birth mass distribution
    logMstar0: log10 Schechter mass of GC birth mass distribution
    logMlo: log10 minimum cluster mass
    logDelta: log10 of mass lost by clusters between birth and now (excluding stellar mass loss). 
    returns: factor by which to multiply BBH merger density to account for cluster disruption/ evaporation mass loss. 
    '''
    
    #Following Section II of Antonini & Gieles 2020
    #Assumes all clusters lost the same mass Delta (excluding stellar mass loss).
    #Delta is typically inferred by comparing evolved GC mass distribution to birth GC mass distribution. 
    #Integral must be evaluated numerically. 
     
    logMc = logMstar0 - np.log10(2) #log(Mstar0/2)
    
    logm_grid = np.logspace(logMlo, logMhi, 20) #log spaced bins between 100 and 10^8 Msun, preliminary tests suggest 20 is enough
    
    phi_cl0 = 2**(-1-beta) * logm_grid**beta * np.exp(-logm_grid/ 10**logMstar0) #birth mass function 

    phi_cl = (logm_grid + (10**logDelta))**beta * np.exp(-(logm_grid+10**logDelta)/ 10**logMc)

    NBH_initial = np.trapz(phi_cl0 * logm_grid**1.6, logm_grid)
    
    NBH_final = np.trapz(phi_cl * logm_grid**1.6, logm_grid)
    
    K_merge = NBH_initial / NBH_final / 2**1.6 #Divide by M = 2 because factor of 2 just from stellar mass loss so doesn't contribute to BBH rate. check that this is still 2 for arbitrary beta, or does it become 2**(-1-beta). 
    
    return K_merge

def average_mass_schechter(beta = -2, logMstar0 = 6.26, logMlo = 2, logMhi = 8, res = 20):
    
    m_grid = np.logspace(logMlo, logMhi, res)
    
    x = m_grid / 10**logMstar0
    
    pdf_mass = x**(beta + 1) * np.exp(-x)
    
    pdf_mass /= np.trapz(pdf_mass, m_grid)
    
    average_mass = np.trapz(pdf_mass * m_grid, m_grid)
    
    return average_mass

def cluster_number_density_from_mass_density(rho_GC = 7.3e14, beta = -2, logMstar0 = 6.26, logMlo = 2, logMhi = 8, logDelta = 5.33):
    '''
    rho_GC: mass density of GCs *today*, units Msun/ Gpc^3 (e.g. Antonini & Gieles 2020 Sec IIA)
    beta: power law slope of Schecther function describing GC birth mass distribution
    logMstar0: log10 Schechter mass of GC birth mass distribution
    logMlo: log10 minimum cluster mass
    logDelta: log10 of mass lost by clusters between birth and now (excluding stellar mass loss). 
    returns: cluster number density given a mass density, assuming mass distribution follows evolved Schechter function. Units 1/ Gpc^3 (or the units of rho_GC/ Msun)
    '''
    logMc = logMstar0 - np.log10(2) #log(Mstar0/2)
    
    logm_grid = np.logspace(logMlo, logMhi, 20)
    #logm_grid = ncl_grid * 0.6
    
    phi_cl = (logm_grid + (10**logDelta))**beta * np.exp(-(logm_grid+10**logDelta)/ 10**logMc)
    
    #number density = rho/<M> where <M> is average cluster mass \int M p(M) dM 
    average_mass = np.trapz(phi_cl * logm_grid, logm_grid)/ np.trapz(phi_cl, logm_grid)
    
    #average_mass = np.sum(phi_cl * logm_grid)/np.sum(phi_cl)
    
    print('number density from mass density', rho_GC/average_mass)

    return rho_GC/ average_mass 
                             
def radius_weights(cluster_radius, mu_rv = 1, sigma_rv = 1.5):
    '''
    assume cluster size distribution is Gaussian
    cluster_radius: virial radius of given cluster (pc)
    mu_rv: mean radius (pc)
    sigma_rv: standard deviation (pc)
    returns: fractional contribution from the given cluster radius (normalized so that the sum over the radius grid is unity)
    '''
    w = np.exp(-(cluster_radius - mu_rv) ** 2. / (2. * sigma_rv ** 2.)) * cluster_radius #must take into account that cluster radius is log-spaced
    w_grid = np.exp(-(rv_grid - mu_rv) ** 2. / (2. * sigma_rv ** 2.)) * rv_grid
    
    return w/np.sum(w_grid)

def redshift_peak(z, a, b, zp):
    '''
    Madau-like redshift distribution
    a: low redshift is approximately (1 + z)^a
    b: high redshift is approximately (1 + z)^-b
    zp: approximate peak redshift
    '''
    return (1.0+(1.0+zp)**(-a-b))*(1+z)**a/(1.0+((1.0+z)/(1.0+zp))**(a+b))

def sfr_at_z_norm(z, z_gc = 4.5, a = 2.5, b = 2.5):
    '''
    cluster star formation history, normalized to give volumetric number density of 1 Gpc^-3 yr^-1 today
    Assume it is Madau-like with params z_gc, a, b
    z_gc: peak redshift
    a: low redshift power-law slope in (1 + z)
    b: high redshift power-law slope slope in (1 + z)
    '''
    dNdVdt_unnorm = redshift_peak(z, a, b, z_gc) #dN/dVcdt(z) 
    dNdV0_unnorm = np.trapz(redshift_peak(zs_i, a, b, z_gc), tLs_i*1e9) #integrate over lookback time, recall that tLs_i is in Gyr
    dNdVdt = dNdVdt_unnorm/dNdV0_unnorm
    
    return dNdVdt

def sfr_at_z(z, dNdV0 = 2.31e9, z_gc = 4.5, a = 2.5, b = 2.5, disrupted_factor = 1.0): 
    '''
    cluster star formation history (e.g. Fig 5 in Rodriguez & Loeb 2018)
    Assume it is Madau-like with params z_gc, a, b
    z_gc: peak redshift
    a: low redshift power-law slope in (1 + z)
    b: high redshift power-law slope slope in (1 + z)
    dNdV0: number density in comoving Gpc^-3 at z = 0, found by integrating the sfr dN/dVdt over all t. Kremer+ 2020 assumes volumetric number density of 2.31e9 Gpc^-3. In terms of mass density, would be typical cluster mass * 2.31e9 Gpc^-3 yr^-1 or ~5e5 Msun Mpc^-3 yr^-1. If mass density is better known than number density, replace this with dM/dV and then divide by typical cluster mass according to assumed mass distribution.
    disrupted_factor: accounts for contribution from clusters that were disrupted/ evaporated before the present day, which has the same effect as adjusting the cluster number density 
    returns: number density (comoving Gpc^-3 yr^-1) evaluated at z
    '''
    dNdVdt = sfr_at_z_norm(z, z_gc, a, b)
    
    dNdVdt_norm = dNdVdt * dNdV0 * disrupted_factor
    
    return dNdVdt_norm
    

def read_data(sourcepath):
    
    #load in data
    data_gwcap = np.genfromtxt(sourcepath+'GWcap_BBH_maingrid.dat')
    gwc_type = data_gwcap[:,3]
    data_inclu = np.genfromtxt(sourcepath+'Incluster_BBH_maingrid.dat')
    data_esc = np.genfromtxt(sourcepath+'Esc_BBH_maingrid.dat')
    tmer_esc = data_esc[:,2]+data_esc[:,3]

    ###extract spin info since they are in another file###
    allbbh = np.genfromtxt(sourcepath+'All_BBH_with_gen_spin.txt')
    allt_mer = allbbh[:,1]; allm0 = allbbh[:,10]; allmodelno = allbbh[:,0]
    allspin0 = allbbh[:,20]; allspin1 = allbbh[:,21]
    allt_code = allbbh[:,2]
    allid0 = allbbh[:,5]
    
    S0 = []; S1 = []
    for ii in range(len(data_gwcap[:,0])):
        if gwc_type[ii]!=2:
            continue
        S0.append(allspin0[(allt_mer==data_gwcap[:,2][ii]) & (allm0==data_gwcap[:,10][ii]) & (allmodelno==data_gwcap[:,0][ii])][0])
        S1.append(allspin1[(allt_mer==data_gwcap[:,2][ii]) & (allm0==data_gwcap[:,10][ii]) & (allmodelno==data_gwcap[:,0][ii])][0])
        if len(allspin0[(allt_mer==data_gwcap[:,2][ii]) & (allm0==data_gwcap[:,10][ii]) & (allmodelno==data_gwcap[:,0][ii])])>1:
            print('error', data_gwcap[:,2][ii], data_gwcap[:,0][ii])

    for jj in range(len(data_inclu[:,0])):
        S0.append(allspin0[(allt_mer==data_inclu[:,2][jj]) & (allm0==data_inclu[:,7][jj]) & (allmodelno==data_inclu[:,0][jj])][0])
        S1.append(allspin1[(allt_mer==data_inclu[:,2][jj]) & (allm0==data_inclu[:,7][jj]) & (allmodelno==data_inclu[:,0][jj])][0])
        if len(allspin0[(allt_mer==data_inclu[:,2][jj]) & (allm0==data_inclu[:,7][jj]) & (allmodelno==data_inclu[:,0][jj])]) > 1:
            print('error', data_inclu[:,2][jj], data_inclu[:,0][jj])

    for kk in range(len(data_esc[:,0])):
        if tmer_esc[kk]>=14000.: continue
        S0.append(allspin0[(allt_code==data_esc[:,1][kk]) & (allid0==data_esc[:,6][kk]) & (allmodelno==data_esc[:,0][kk])][0])
        S1.append(allspin1[(allt_code==data_esc[:,1][kk]) & (allid0==data_esc[:,6][kk]) & (allmodelno==data_esc[:,0][kk])][0])
        if len(allspin0[(allt_code==data_esc[:,1][kk]) & (allid0==data_esc[:,6][kk]) & (allmodelno==data_esc[:,0][kk])])>1:
            print('error', data_esc[:,1][kk], data_esc[:,0][kk], allspin0[(allt_code==data_esc[:,1][kk]) & (allm0==data_esc[:,4][kk]) & (allmodelno==data_esc[:,0][kk])])

    S0 = np.array(S0); S1 = np.array(S1)    

    numsim = np.array(list(data_gwcap[:,0][gwc_type==2])+list(data_inclu[:,0])+list(data_esc[:,0][tmer_esc<14000.])).astype(int)
    tgw = np.array(list(data_gwcap[:,2][gwc_type==2]/1000.)+list(data_inclu[:,2]/1000.)+list(tmer_esc[tmer_esc<14000.]/1000.))  ##merger times in Gyr
    M0 = np.array(list(data_gwcap[:,10][gwc_type==2])+list(data_inclu[:,7])+list(data_esc[:,4][tmer_esc<14000.]))
    M1 = np.array(list(data_gwcap[:,11][gwc_type==2])+list(data_inclu[:,8])+list(data_esc[:,5][tmer_esc<14000.]))
    
    
    paths = np.genfromtxt(sourcepath+'path_allfinished_newruns_maingrid.dat', dtype='str')
    paths = paths[:,0]
    ncll = []; zb = []; rvv = []
    for xx in range(len(numsim)):

        s=paths[numsim[xx]].split('/')
        n_star=float(s[-2])
        z=float(s[-3][1:])
        rg=int(s[-4][2:])
        rv=float(s[-5][2:])

        ncll.append(n_star); rvv.append(rv); zb.append(z)
    
    ncll = np.array(ncll); rvv = np.array(rvv); zb = np.array(zb)

    #cluster params and merger time for each BBH merger(?) -- total is 4330 mergers
    #numsim = data[:,0] #simulation number
    #rvv = data[:,1] #virial radii
    #zb = data[:,3] #metallicity
    #ncll = data[:,4] #number of particles. multiply by 0.6 Msun to get stellar mass. 
    #tgw = data[:,6] * 1e-3 #merger times in Gyr

    # limit to ncll<2.e6 to have a uniform grid
    #(Idx,) = np.where(ncll < 2.e6) 
    #numsim = numsim[Idx] 
    #rvv = rvv[Idx]
    #zb = zb[Idx]
    #ncll = ncll[Idx]
    #tgw = tgw[Idx]

    #141 different GC simulations -- supposed to be 144 (4 ncl, 4rv, 3 zeta, 3 different galactocentric radii), 
    #but missing the 3 corresponding to (ncll/2e5==8) & (rvv == 0.5) & (zb == 0.0002) (numsim = 3,15,27). 
    #It's probably fine because such low metallicities are very rare at relevant redshifts.
    #To make the grid consistent, add in these missing simulations assuming they are identical to the zb == 0.002 versions. 

    #select sims we want to copy
    copy_sel = (ncll/2e5==8) & (rvv == 0.5) & (zb == 0.002)
    
    ncopy = len(numsim[copy_sel])
    
    #get numsim of missing sims
    #missing_sims = list(set(1+np.arange(143)).difference(set(numsim)))
    
    #make the copies
    rvv_copy = rvv[copy_sel]
    ncll_copy = ncll[copy_sel]
    tgw_copy = tgw[copy_sel]
    M0_copy = M0[copy_sel]
    M1_copy = M1[copy_sel]
    S0_copy = S0[copy_sel]
    S1_copy = S1[copy_sel]
    print(M0_copy, M1_copy, len(M0_copy))
    
    #pretend they correspond to the missing zb
    zb_copy = 0.0002 * np.ones(ncopy)
    
    #label the fake sims with -1* original numsim so we hopefully remember we did something sketchy 
    numsim_copy = -numsim[copy_sel]

    rvv_new = np.concatenate((rvv, rvv_copy))
    ncll_new = np.concatenate((ncll, ncll_copy))
    tgw_new = np.concatenate((tgw, tgw_copy))

    zb_new = np.concatenate((zb, zb_copy))
    
    numsim_new = np.concatenate((numsim, numsim_copy))

    M0_new = np.concatenate((M0, M0_copy))
    M1_new = np.concatenate((M1, M1_copy))

    S0_new = np.concatenate((S0, S0_copy))
    S1_new = np.concatenate((S1, S1_copy))
    
    return numsim_new, rvv_new, zb_new, ncll_new, tgw_new, M0_new, M1_new, S0_new, S1_new

def weighted_quantile(values, quantiles, sample_weight=None, 
                      values_sorted=False, old_style=False):
    """ Very close to numpy.percentile, but supports weights.
    NOTE: quantiles should be in [0, 1]!
    :param values: numpy.array with data
    :param quantiles: array-like with many quantiles needed
    :param sample_weight: array-like of the same length as `array`
    :param values_sorted: bool, if True, then will avoid sorting of
        initial array
    :param old_style: if True, will correct output to be consistent
        with numpy.percentile.
    :return: numpy.array with computed quantiles.
    """
    values = np.array(values)
    quantiles = np.array(quantiles)
    if sample_weight is None:
        sample_weight = np.ones(len(values))
    sample_weight = np.array(sample_weight)
    assert np.all(quantiles >= 0) and np.all(quantiles <= 1), \
        'quantiles should be in [0, 1]'

    if not values_sorted:
        sorter = np.argsort(values)
        values = values[sorter]
        sample_weight = sample_weight[sorter]

    weighted_quantiles = np.cumsum(sample_weight) - 0.5 * sample_weight
    if old_style:
        # To be convenient with numpy.percentile
        weighted_quantiles -= weighted_quantiles[0]
        weighted_quantiles /= weighted_quantiles[-1]
    else:
        weighted_quantiles /= np.sum(sample_weight)
    return np.interp(quantiles, weighted_quantiles, values)

def ave_mass_redshift(data, zmerge, z_gc = 4.5, a = 2.5, b = 2.5, dNdV0 = 2.31e9, logf_disrupted_cluster = 0.5, sigma_dex = 0.5, Zsun = 0.02, mu_rv = 1, sigma_rv = 1.5, beta = -2, logMstar0 = 6.26, logMlo = 2, logMhi = 8):
    #numsim, rvv, zb, ncll, tgw = data[0], data[1], data[2], data[3], data[4]
    m0, m1 = data[5], data[6]
    M0 = np.maximum(m0, m1)
    M1 = np.minimum(m0, m1)

    numsim = data[0]    
    rvv = data[1]
    zb = data[2] 
    ncll = data[3] 
    tgw = data[4]

    #compute mass and radius weights for each simulation based on ncl, rv. 
    mweights = mass_weights_schechter(ncll*0.6, beta, logMstar0)
    rweights = radius_weights(rvv, mu_rv, sigma_rv)
    

    cluster_weight = mweights * rweights * dNdV0 * 10**logf_disrupted_cluster 
    #f_missing_cluster * (remove this and leave the uncertainties in logf_disrupted_cluster)
    #used to be cluster_weight = mweights * rweights * f_missing_cluster * dNdV0 * 10**logf_disrupted_cluster

    tL_merge = tL_at_z_interp(zmerge) #lookback time at merger in Gyr
    tL_form = tL_merge + tgw #lookback time at formation
    z_form = z_at_tL_interp(tL_form) #redshift at formation
    
    metal_weight = metallicity_weights(zb, z_form, sigma_dex, Zsun)

    out = np.sum(M0*cluster_weight * metal_weight * sfr_at_z_norm(z_form, z_gc, a, b))/np.sum(cluster_weight * metal_weight * sfr_at_z_norm(z_form, z_gc, a, b)) #sum over all mergers

    return out

def std_mass_redshift(data, zmerge, z_gc = 4.5, a = 2.5, b = 2.5, dNdV0 = 2.31e9, logf_disrupted_cluster = 0.5, sigma_dex = 0.5, Zsun = 0.02, mu_rv = 1, sigma_rv = 1.5, beta = -2, logMstar0 = 6.26, logMlo = 2, logMhi = 8):
    #numsim, rvv, zb, ncll, tgw = data[0], data[1], data[2], data[3], data[4]
    m0, m1 = data[5], data[6]
    M0 = np.maximum(m0, m1)
    M1 = np.minimum(m0, m1)

    numsim = data[0]
    rvv = data[1]
    zb = data[2]
    ncll = data[3]
    tgw = data[4]

    #compute mass and radius weights for each simulation based on ncl, rv.
    mweights = mass_weights_schechter(ncll*0.6, beta, logMstar0)
    rweights = radius_weights(rvv, mu_rv, sigma_rv)


    cluster_weight = mweights * rweights * dNdV0 * 10**logf_disrupted_cluster
    #f_missing_cluster * (remove this and leave the uncertainties in logf_disrupted_cluster)
    #used to be cluster_weight = mweights * rweights * f_missing_cluster * dNdV0 * 10**logf_disrupted_cluster

    tL_merge = tL_at_z_interp(zmerge) #lookback time at merger in Gyr
    tL_form = tL_merge + tgw #lookback time at formation
    z_form = z_at_tL_interp(tL_form) #redshift at formation

    metal_weight = metallicity_weights(zb, z_form, sigma_dex, Zsun)

    #out = np.sum(M0*cluster_weight * metal_weight * sfr_at_z_norm(z_form, z_gc, a, b))/np.sum(cluster_weight * metal_weight * sfr_at_z_norm(z_form, z_gc, a, b)) #sum over all mergers
    #out_std = math.sqrt(np.sum((M0-out)**2*cluster_weight * metal_weight * sfr_at_z_norm(z_form, z_gc, a, b))/(np.count_nonzero(cluster_weight * metal_weight * sfr_at_z_norm(z_form, z_gc, a, b))-1)*np.count_nonzero(cluster_weight * metal_weight * sfr_at_z_norm(z_form, z_gc, a, b))/np.sum(cluster_weight * metal_weight * sfr_at_z_norm(z_form, z_gc, a, b)))
    
    combined_weights = cluster_weight * metal_weight * sfr_at_z_norm(z_form, z_gc, a, b)
    out_std = DescrStatsW(data=M0, weights=combined_weights).std
    
    return out_std

def percentile_redshift(data, zmerge, perc_input, z_gc = 4.5, a = 2.5, b = 2.5, dNdV0 = 2.31e9, logf_disrupted_cluster = 0.5, sigma_dex = 0.5, Zsun = 0.02, mu_rv = 1, sigma_rv = 1.5, beta = -2, logMstar0 = 6.26, logMlo = 2, logMhi = 8):
    #numsim, rvv, zb, ncll, tgw = data[0], data[1], data[2], data[3], data[4]
    m0, m1 = data[5], data[6]
    M0 = np.maximum(m0, m1)
    M1 = np.minimum(m0, m1)

    numsim = data[0]
    rvv = data[1]
    zb = data[2]
    ncll = data[3]
    tgw = data[4]

    #compute mass and radius weights for each simulation based on ncl, rv.
    mweights = mass_weights_schechter(ncll*0.6, beta, logMstar0)
    rweights = radius_weights(rvv, mu_rv, sigma_rv)


    cluster_weight = mweights * rweights * dNdV0 * 10**logf_disrupted_cluster
    #f_missing_cluster * (remove this and leave the uncertainties in logf_disrupted_cluster)
    #used to be cluster_weight = mweights * rweights * f_missing_cluster * dNdV0 * 10**logf_disrupted_cluster

    tL_merge = tL_at_z_interp(zmerge) #lookback time at merger in Gyr
    tL_form = tL_merge + tgw #lookback time at formation
    z_form = z_at_tL_interp(tL_form) #redshift at formation

    metal_weight = metallicity_weights(zb, z_form, sigma_dex, Zsun)

    combined_weights = cluster_weight * metal_weight * sfr_at_z_norm(z_form, z_gc, a, b)
    #out_quantiles = weighted_quantile(M0, [perc_input], sample_weight=combined_weights, values_sorted=False, old_style=False)  #wikipedia calculation, same as the statsmodels method below.
    
    wq = DescrStatsW(data=M0, weights=combined_weights)
    out_quantiles=wq.quantile(probs=np.array([perc_input]), return_pandas=False)

    return out_quantiles


def spin_redshift(data, zmerge, perc_input, z_gc = 4.5, a = 2.5, b = 2.5, dNdV0 = 2.31e9, logf_disrupted_cluster = 0.5, sigma_dex = 0.5, Zsun = 0.02, mu_rv = 1, sigma_rv = 1.5, beta = -2, logMstar0 = 6.26, logMlo = 2, logMhi = 8):
    m0, m1 = data[5], data[6]
    s0, s1 = data[7], data[8]

    numsim = data[0]
    rvv = data[1]
    zb = data[2]
    ncll = data[3]
    tgw = data[4]

    #compute mass and radius weights for each simulation based on ncl, rv.
    mweights = mass_weights_schechter(ncll*0.6, beta, logMstar0)
    rweights = radius_weights(rvv, mu_rv, sigma_rv)


    cluster_weight = mweights * rweights * dNdV0 * 10**logf_disrupted_cluster
    #f_missing_cluster * (remove this and leave the uncertainties in logf_disrupted_cluster)
    #used to be cluster_weight = mweights * rweights * f_missing_cluster * dNdV0 * 10**logf_disrupted_cluster

    tL_merge = tL_at_z_interp(zmerge) #lookback time at merger in Gyr
    tL_form = tL_merge + tgw #lookback time at formation
    z_form = z_at_tL_interp(tL_form) #redshift at formation

    metal_weight = metallicity_weights(zb, z_form, sigma_dex, Zsun)
    
    combined_weights = cluster_weight * metal_weight * sfr_at_z_norm(z_form, z_gc, a, b)
    
    alpha_list = np.random.uniform(low=0., high=np.pi, size=len(m0))
    beta_list = np.random.uniform(low=0., high=np.pi, size=len(m0))
    
    chieff = chi_eff(m0, m1, s0, s1, alpha_list, beta_list)
   
    signs = np.random.choice([0,1], size=len(chieff[chieff==0]))
    combined_weights_pos = list(combined_weights[chieff>0])+list(combined_weights[chieff==0][signs==1])
    chieff_pos = list(chieff[chieff>0])+list(chieff[chieff==0][signs==1])
    combined_weights_neg = list(combined_weights[chieff<0])+list(combined_weights[chieff==0][signs==0])
    chieff_neg = list(chieff[chieff<0])+list(chieff[chieff==0][signs==0])
    #print(np.mean(combined_weights_pos)/np.std(combined_weights_pos), np.mean(combined_weights_neg)/np.std(combined_weights_neg))
    #print(np.mean(combined_weights)/np.std(combined_weights))

    out = np.sum(chieff*cluster_weight * metal_weight * sfr_at_z_norm(z_form, z_gc, a, b))/np.sum(cluster_weight * metal_weight * sfr_at_z_norm(z_form, z_gc, a, b))
    out_std = DescrStatsW(data=chieff, weights=combined_weights).std

    #wq_pos = DescrStatsW(data=chieff_pos, weights=combined_weights_pos)
    #out_quantiles_pos=wq_pos.quantile(probs=np.array([perc_input]), return_pandas=False)
    #wq_neg = DescrStatsW(data=np.abs(chieff_neg), weights=combined_weights_neg)
    #out_quantiles_neg=-wq_neg.quantile(probs=np.array([perc_input]), return_pandas=False)
    wq_abs = DescrStatsW(data=np.abs(chieff), weights=combined_weights)
    out_quantiles_abs=wq_abs.quantile(probs=np.array([perc_input]), return_pandas=False)

    return out, out_quantiles_abs#out_quantiles_pos, out_quantiles_neg


#each BBH came from a cluster that represents a rate density at some z/time. 

def merger_rate_at_z(zmerge, formation_rate_at_z, tgw, cluster_weight, metal, metal_frac_at_z, sfr_kwargs = {}, metal_kwargs = {}):

    ''' 
    zmerge: desired merger redshift 
    formation_rate_at_z: a function that returns the formation rate (dN/dVcdt) at a given redshift
    tgw: array of delay times between formation and merger (Gyr), each delay time coresponds to one BBH
    cluster_weight: weight assigned to specific cluster, same dimensions as tgw
    metal: metallicity assigned to specific cluster, same dimensions as tgw
    metal_frac_at_z: a function that returns the metallicity fraction at a given formation redshift and metallicity
    sfr_kwargs: other params called by formation_rate_at_z
    metal_kwargs: other params called by metal_pdf_at_z
    returns: merger rate at given zmerge
    '''
   
    tL_merge = tL_at_z_interp(zmerge) #lookback time at merger in Gyr
    tL_form = tL_merge + tgw #lookback time at formation
    z_form = z_at_tL_interp(tL_form) #redshift at formation
    
    metal_weight = metal_frac_at_z(metal, z_form, **metal_kwargs)

    return np.sum(cluster_weight * metal_weight * formation_rate_at_z(z_form, **sfr_kwargs)) #sum over all mergers

def merger_rate_at_z_pop(data, zmerge, mlow, mhigh, z_gc = 4.5, a = 2.5, b = 2.5, dNdV0 = 2.31e9, logf_disrupted_cluster = 0.5, sigma_dex = 0.5, Zsun = 0.02, mu_rv = 1, sigma_rv = 1.5, beta = -2, logMstar0 = 6.26, logMlo = 2, logMhi = 8):
    '''
    data: output of read_data() -- list of numsim, rvv, zb, ncll, tgw
    zmerge: merger redshift
    z_gc: peak formation redshift
    a: formation rate follows (1 + z)^a at low z
    b: formation rate follows (1 + z)^-b at high z
    dNdV0: number density of GCs today in units Gpc^-3
    logf_disrupted_cluster: log10 of the contribution to formation rate at each z from cluster mass lost between formation and today
    sigma_dex: scatter in metallicity-redshift relation
    Zsun: solar metallicity
    mu_rv: mean cluster radius (pc)
    sigma_rv: standard deviation of cluster radius distripution (pc)
    beta: power law slope of birth cluster mass distribution
    logMstar0: log10 Schechter mass of birth cluster mass distribution
    logMlo: log10 minimum GC mass (Msun)
    logMhi: log10 maximum GC mass (Msun)
    '''

    #numsim, rvv, zb, ncll, tgw = data[0], data[1], data[2], data[3], data[4]
    m0, m1 = data[5], data[6]
    M0 = np.maximum(m0, m1)
    M1 = np.minimum(m0, m1)

    numsim = data[0][(M0>=mlow) & (M0<mhigh)] 
    rvv = data[1][(M0>=mlow) & (M0<mhigh)] 
    zb = data[2][(M0>=mlow) & (M0<mhigh)] 
    ncll = data[3][(M0>=mlow) & (M0<mhigh)] 
    tgw = data[4][(M0>=mlow) & (M0<mhigh)]

    #compute mass and radius weights for each simulation based on ncl, rv. 
    mweights = mass_weights_schechter(ncll*0.6, beta, logMstar0)
    rweights = radius_weights(rvv, mu_rv, sigma_rv)
    
    f_missing_cluster = compute_missing_cluster_factor(beta, logMstar0, logMlo, logMhi)
    
    cluster_weight = mweights * rweights * dNdV0 * 10**logf_disrupted_cluster 
    #f_missing_cluster * (remove this and leave the uncertainties in logf_disrupted_cluster)
    #used to be cluster_weight = mweights * rweights * f_missing_cluster * dNdV0 * 10**logf_disrupted_cluster
    
    merger_rate_array = merger_rate_at_z(zmerge, sfr_at_z_norm, tgw, cluster_weight, zb, metallicity_weights, sfr_kwargs = {'z_gc': z_gc, 'a': a, 'b': b}, metal_kwargs = {'sigma_dex': sigma_dex, 'Zsun': Zsun})
    #print('merger_rate_array', merger_rate_array)

    out = np.sum(merger_rate_array)
    
    return out


def merger_rate_at_z_pop_metal(data, zmerge, zmetal, z_gc = 4.5, a = 2.5, b = 2.5, dNdV0 = 2.31e9, logf_disrupted_cluster = 0.5, sigma_dex = 0.5, Zsun = 0.02, mu_rv = 1, sigma_rv = 1.5, beta = -2, logMstar0 = 6.26, logMlo = 2, logMhi = 8):
    '''
    data: output of read_data() -- list of numsim, rvv, zb, ncll, tgw
    zmerge: merger redshift
    z_gc: peak formation redshift
    a: formation rate follows (1 + z)^a at low z
    b: formation rate follows (1 + z)^-b at high z
    dNdV0: number density of GCs today in units Gpc^-3
    logf_disrupted_cluster: log10 of the contribution to formation rate at each z from cluster mass lost between formation and today
    sigma_dex: scatter in metallicity-redshift relation
    Zsun: solar metallicity
    mu_rv: mean cluster radius (pc)
    sigma_rv: standard deviation of cluster radius distripution (pc)
    beta: power law slope of birth cluster mass distribution
    logMstar0: log10 Schechter mass of birth cluster mass distribution
    logMlo: log10 minimum GC mass (Msun)
    logMhi: log10 maximum GC mass (Msun)
    '''

    #numsim, rvv, zb, ncll, tgw = data[0], data[1], data[2], data[3], data[4]
    m0, m1 = data[5][data[2]==zmetal], data[6][data[2]==zmetal]
    M0 = np.maximum(m0, m1)
    M1 = np.minimum(m0, m1)

    numsim = data[0][data[2]==zmetal]
    rvv = data[1][data[2]==zmetal]
    zb = data[2][data[2]==zmetal]
    ncll = data[3][data[2]==zmetal]
    tgw = data[4][data[2]==zmetal]

    
    #compute mass and radius weights for each simulation based on ncl, rv.
    mweights = mass_weights_schechter(ncll*0.6, beta, logMstar0)
    rweights = radius_weights(rvv, mu_rv, sigma_rv)

    f_missing_cluster = compute_missing_cluster_factor(beta, logMstar0, logMlo, logMhi)

    cluster_weight = mweights * rweights * dNdV0 * 10**logf_disrupted_cluster
    #f_missing_cluster * (remove this and leave the uncertainties in logf_disrupted_cluster)
    #used to be cluster_weight = mweights * rweights * f_missing_cluster * dNdV0 * 10**logf_disrupted_cluster

    merger_rate_array = merger_rate_at_z(zmerge, sfr_at_z_norm, tgw, cluster_weight, zb, metallicity_weights, sfr_kwargs = {'z_gc': z_gc, 'a': a, 'b': b}, metal_kwargs = {'sigma_dex': sigma_dex, 'Zsun': Zsun})
    #print('merger_rate_array', merger_rate_array)

    out = np.sum(merger_rate_array)

    return out


def merger_rate_at_z_pop_gen(data, zmerge, ngen, z_gc = 4.5, a = 2.5, b = 2.5, dNdV0 = 2.31e9, logf_disrupted_cluster = 0.5, sigma_dex = 0.5, Zsun = 0.02, mu_rv = 1, sigma_rv = 1.5, beta = -2, logMstar0 = 6.26, logMlo = 2, logMhi = 8):
    '''
    data: output of read_data() -- list of numsim, rvv, zb, ncll, tgw
    zmerge: merger redshift
    z_gc: peak formation redshift
    a: formation rate follows (1 + z)^a at low z
    b: formation rate follows (1 + z)^-b at high z
    dNdV0: number density of GCs today in units Gpc^-3
    logf_disrupted_cluster: log10 of the contribution to formation rate at each z from cluster mass lost between formation and today
    sigma_dex: scatter in metallicity-redshift relation
    Zsun: solar metallicity
    mu_rv: mean cluster radius (pc)
    sigma_rv: standard deviation of cluster radius distripution (pc)
    beta: power law slope of birth cluster mass distribution
    logMstar0: log10 Schechter mass of birth cluster mass distribution
    logMlo: log10 minimum GC mass (Msun)
    logMhi: log10 maximum GC mass (Msun)
    '''

    #numsim, rvv, zb, ncll, tgw = data[0], data[1], data[2], data[3], data[4]
    m0, m1 = data[5], data[6]
    M0 = np.maximum(m0, m1)
    M1 = np.minimum(m0, m1)
 
    s0, s1 = data[7], data[8]

    if ngen=='1G':
        numsim = data[0][(s0==0.) & (s1==0.)]
        rvv = data[1][(s0==0.) & (s1==0.)]
        zb = data[2][(s0==0.) & (s1==0.)]
        ncll = data[3][(s0==0.) & (s1==0.)]
        tgw = data[4][(s0==0.) & (s1==0.)]
    else:
        numsim = data[0][(s0>0.) | (s1>0.)]
        rvv = data[1][(s0>0.) | (s1>0.)]
        zb = data[2][(s0>0.) | (s1>0.)]
        ncll = data[3][(s0>0.) | (s1>0.)]
        tgw = data[4][(s0>0.) | (s1>0.)]


    #compute mass and radius weights for each simulation based on ncl, rv.
    mweights = mass_weights_schechter(ncll*0.6, beta, logMstar0)
    rweights = radius_weights(rvv, mu_rv, sigma_rv)

    f_missing_cluster = compute_missing_cluster_factor(beta, logMstar0, logMlo, logMhi)

    cluster_weight = mweights * rweights * dNdV0 * 10**logf_disrupted_cluster
    #f_missing_cluster * (remove this and leave the uncertainties in logf_disrupted_cluster)
    #used to be cluster_weight = mweights * rweights * f_missing_cluster * dNdV0 * 10**logf_disrupted_cluster

    merger_rate_array = merger_rate_at_z(zmerge, sfr_at_z_norm, tgw, cluster_weight, zb, metallicity_weights, sfr_kwargs = {'z_gc': z_gc, 'a': a, 'b': b}, metal_kwargs = {'sigma_dex': sigma_dex, 'Zsun': Zsun})
    #print('merger_rate_array', merger_rate_array)

    out = np.sum(merger_rate_array)

    return out


def merger_rate_at_z_pop_selfconsistentfactors(data, zmerge, mlow, mhigh, z_gc = 4.5, a = 2.5, b = 2.5, sigma_dex = 0.5, Zsun = 0.02, mu_rv = 1, sigma_rv = 1.5, beta = -2, logMstar0 = 6.26, rho_GC = 7.3e14, logDelta = 5.33, logMlo = 2, logMhi = 8, average_M_evolved = None):
    '''
    data: output of read_data() -- list of numsim, rvv, zb, ncll, tgw
    zmerge: merger redshift
    z_gc: peak formation redshift
    a: formation rate follows (1 + z)^a at low z
    b: formation rate follows (1 + z)^-b at high z
    sigma_dex: scatter in metallicity-redshift relation
    Zsun: solar metallicity
    mu_rv: mean cluster radius (pc)
    sigma_rv: standard deviation of cluster radius distripution (pc)
    beta: power law slope of birth cluster mass distribution
    logMstar0: log10 Schechter mass of birth cluster mass distribution
    rho_GC: mass density of GCs today (Msun/ Gpc^3)
    logDelta: log10 mass (Msun) lost by GCs between formation and today (excluding stellar mass loss)
    logMlo: log10 minimum GC mass (Msun)
    logMhi: log10 maximum GC mass (Msun)
    average_M_evolved: average GC mass of evolved clusters (Msun) used to compute number density dNdV0 from rho_GC. If None, then average evolved mass is computed from other parameters assuming model of Antonini & Gieles 2020. Typical value is 3e5. 
    '''
    
    if average_M_evolved:
        dNdV0 = rho_GC/ average_M_evolved
    else:
        dNdV0 = cluster_number_density_from_mass_density(rho_GC, beta, logMstar0, logMlo, logMhi, logDelta)
    
    f_missing_cluster = compute_missing_cluster_factor(beta, logMstar0, logMlo, logMhi)
    
    f_disrupted_cluster = compute_disrupted_cluster_factor(beta, logMstar0, logMlo, logMhi, logDelta)
    #print('disrupt factor', f_disrupted_cluster)

    M0, M1 = data[5], data[6]
    numsim = data[0][(M0>=mlow) & (M0<mhigh)] 
    rvv = data[1][(M0>=mlow) & (M0<mhigh)] 
    zb = data[2][(M0>=mlow) & (M0<mhigh)] 
    ncll = data[3][(M0>=mlow) & (M0<mhigh)] 
    tgw = data[4][(M0>=mlow) & (M0<mhigh)] 
    
            
    #compute mass and radius weights for each simulation based on ncl, rv. 
    mweights = mass_weights_schechter(ncll*0.6, beta, logMstar0)
    rweights = radius_weights(rvv, mu_rv, sigma_rv) 
    
    cluster_weight = mweights * rweights * dNdV0 * f_missing_cluster * f_disrupted_cluster
    
    merger_rate_array = merger_rate_at_z(zmerge, sfr_at_z_norm, tgw, cluster_weight, zb, metallicity_weights, sfr_kwargs = {'z_gc': z_gc, 'a': a, 'b': b}, metal_kwargs = {'sigma_dex': sigma_dex, 'Zsun': Zsun})
    
    out = np.sum(merger_rate_array)
    
    return out


    
