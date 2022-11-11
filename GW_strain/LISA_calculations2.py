import numpy as np

from scipy import special
from scipy import integrate
from scipy.integrate import ode
from scipy.optimize import brentq
pi = np.pi

import pdb
import matplotlib.pyplot as plt
import scipy.constants as ct

from scipy.integrate import odeint#, solve_ivp
from scipy.special import jv #bessel function of the first kind
from astropy.cosmology import Planck15 as cosmo
from scipy.interpolate import interp1d


Msun = 1.99e30
G = 6.67e-11
c = 3.0e8

def e_factor(e):
    return (1+e)**1.1954 / (1-e**2.)**1.5

def chirp_mass(m1,m2,z=0):
    return (1+z)*(m1*m2)**0.6 / (m1+m2)**0.2

def eta(m1,m2):
    return (m1*m2)/(m1+m2)**2

def deda_peters(a,e):
    num = 12*a*(1+(73./24)*e**2 + (37./96)*e**4)
    denom = 19*e*(1-e**2)*(1+(121./304)*e**2)
    return denom/num

def inspiral_time_peters(a0,e0,m1,m2,af=0):
    """
    Computes the inspiral time, in Gyr, for a binary
    a0 in Au, and masses in solar masses
    
    if different af is given, computes the time from a0,e0
    to that af
    
    for af=0, just returns inspiral time
    for af!=0, returns (t_insp,af,ef)
    """
    coef = 6.086768e-11 #G^3 / c^5 in au, gigayear, solar mass units
    beta = (64./5.) * coef * m1 * m2 * (m1+m2)
    
    if e0 == 0:
        print(e0,a0)
        if not af == 0:
            print("ERROR: doesn't work for circular binaries")
            return 0
        return a0**4 / (4*beta)
    
    c0 = a0 * (1.-e0**2.) * e0**(-12./19.) * (1.+(121./304.)*e0**2.)**(-870./2299.)
    
    if af == 0:
        eFinal = 0.
    else:
        r = ode(deda_peters)
        r.set_integrator('lsoda')
        r.set_initial_value(e0,a0)
        r.integrate(af)
        if not r.successful():
            print("ERROR, Integrator failed!")
        else:
            eFinal = r.y[0]      
    
    time_integrand = lambda e: e**(29./19.)*(1.+(121./304.)*e**2.)**(1181./2299.) / (1.-e**2.)**1.5
    integral,abserr = integrate.quad(time_integrand,eFinal,e0)
    
    if af==0:
        return integral * (12./19.) * c0**4. / beta
    else:
        return (integral * (12./19.) * c0**4. / beta), af, eFinal

def a_at_fLow(m1,m2,fLow = 5):
    """
    Computes the semi-major axis at an orbital frequency of fLow
    Masses in solar masses, fLow in Hz
    """
    G = 3.9652611e-14 # G in au,solar mass, seconds
    quant = G*(m1+m2) / (4*pi**2 * fLow**2)
    return quant**(1./3)

def eccentricity_at_fLow(m1,m2,a_0,e_0,fLow=5):
    """
    Computes the eccentricity at a given fLow
    
    Masses are in solar masses, a_0 in AU
    
    NOTE!!! The frequency here is the ORBITAL frequency.  
    So if you want the eccentricity at a G-wave frequency of 10, set fLow=5
    (divide by 2)
    """
    a_low = a_at_fLow(m1,m2,fLow)
    
    r = ode(deda_peters)
    r.set_integrator('lsoda')
    r.set_initial_value(e_0,a_0)
    
    r.integrate(a_low)
    
    if not r.successful():
        print("ERROR, Integrator failed!")
    else:
        return r.y[0]
    

def comovingDistance(z):
    h = 0.679
    omegaM = 0.306
    omegaK = 0.
    omegaL = 1 - 0.306
    dh = 3000. / h
    e = lambda zp: 1./np.sqrt(omegaM*(1+zp)**3 + omegaL)
    return dh*integrate.quad(e,0,z)[0]

def luminosityDistance(z):
     return (1+z)*comovingDistance(z)

def zAtLuminosityDistance(d):
    zero = lambda z: luminosityDistance(z) - d
    return brentq(zero,0,5)

#zz = linspace(0,3.6,5000)
#cmcm = array([luminosityDistance(z) for z in zz])
#def luminosityDistancePrecomputed(z):
#    idx = int(z/3.6*5000)
#    if idx > 5000:
#        return inf
#    return numpy.interp(z,zz,cmcm)

def lookbackTime(z):
    h = 0.679
    omegaM = 0.306
    omegaK = 0.
    omegaL = 1 - 0.306
    th = 9.78/h
    e = lambda zp: 1./(np.sqrt(omegaM*(1+zp)**3 + omegaL)*(1+zp))
    return th*integrate.quad(e,0,z)[0]

def zAtLookbackTime(t):
    zero = lambda z: lookbackTime(z) - t
    return brentq(zero,0,10)

#def dVcdz(z):
#    h = 0.679
#    omegaM = 0.306
#    omegaK = 0.
#    omegaL = 1 - 0.306
#    dh = 3000. / h
#    e = sqrt(omegaM*(1+z)**3 + omegaL)
#    return 4*pi*dh*comovingDistance(z)**2/e

def g_n(n,m0,m1,Porb,e):
    gn = pow(n,4.)/32.*(pow(special.jn(n-2.,n*e) - 2*e*special.jn(n-1.,n*e) + 2./n*special.jn(n,n*e) + 2*e*special.jn(n+1.,n*e) - special.jn(n+2.,n*e),2.) + (1-e*e)*pow(special.jn(n-2.,n*e) - 2*special.jn(n,n*e) + special.jn(n+2.,n*e),2.) + 4/(3*n*n)*pow(special.jn(n,n*e),2))
    return gn
def fdotn(m0,m1,Porb,e):
    F = (1 + 73./24.*pow(e,2.) + 37./96.*pow(e,4.))*pow(1-e*e,-7./2.)
    nu = 1/Porb*1/c
    Mc = pow(m0*m1,3./5.)/pow(m0+m1,1./5.)*G/(c*c)
    fdotn = 96./(10.*np.pi)*pow(Mc,5./3.)*pow(2*np.pi*nu,11./3.)*F
    return fdotn
def h_o(m0,m1,Porb,d,n):
    forb = 1/Porb
    m0 = m0*Msun
    m1 = m1*Msun
    Mc = pow(m0*m1,3./5.)/pow(m0+m1,1./5.)
    h_o = G/pow(c,2.) * Mc/d * pow(G/pow(c,3.)*np.pi*forb*Mc,2./3.)
    return h_o
def hn(m1,m2,a,e,d,n):
    """
    Takes d in Mpc
    """
    D = d*3.086e16*1.e6
    Porb = np.sqrt(4.*np.pi**2./(G*Msun*(m1+m2))*(a*1.5e11)**3.)
    f_orb = 1./Porb
    #print Porb, f_orb
    F = (1 + 73./24.*pow(e,2.) + 37./96.*pow(e,4.))*pow(1-e*e,-7./2.)
    z = zAtLuminosityDistance(d)  # takes d in Mpc
    #print z
    Mchirp = chirp_mass(m1,m2,z=0)
    gn = g_n(n,m1,m2,Porb,e)
    h_n = 1/(np.pi*D)*pow(2*1./c**3.0 * (G*Mchirp*Msun)**(5./3.) / (3*np.pi**(1./3.)*(1+z)**(1./3.)*(n*f_orb/(1+z))**(1./3.)) * (2./n)**(2./3.) * gn/F , 0.5)
    return h_n

def SNR_calculator(m0,m1,a,e,d):
    Porb = (4.*np.pi**2./(6.67e-11*1.99e30*(m0+m1))*(a*1.5e11)**3.)**0.5
    f_orb = 1/Porb

    d = d*1.e6*3.086e16     # CONVERT D TO METERS
    SNR_tot_sq = 0
    ### IMPORT LISA NOISE CURVE
    data2 = np.loadtxt('/projects/b1095/syr904/MyCodes/PythonCodes/GW_strain/characteristic_noise_strain.dat')
    fn = data2[:,0]
    hn_lisa = np.sqrt(fn)*data2[:,1]
    ###########################
    Tobs = 10.0*3.15e7
    if e < 0.0001:
        h_lisa = np.interp(f_orb,fn[:],hn_lisa[:])
        ho = h_o(m0,m1,Porb,d,1.)
        SNR = np.sqrt(2*pow(ho,2.)*Tobs/pow(h_lisa,2.))
    if e >= 0.0001:
        for n in range(1,100):
            h_lisa = np.interp(f_orb*n,fn[:],hn_lisa[:])
            ho = h_o(m0,m1,Porb,d,n)
            gn = g_n(n,m0,m1,Porb,e)
            SNR_n_sq = pow(ho,2.)*gn/pow(n,2.)*Tobs/pow(h_lisa,2.)
            SNR_tot_sq = SNR_tot_sq + SNR_n_sq

        SNR_tot = np.sqrt(2.*SNR_tot_sq)
        SNR = SNR_tot

    return SNR

def Integrate(m1,m2,a0,e0,t_f):
    #t_inspiral = inspiral_time_peters(a0,e0,m1,m2,af=0)*1000.
    M1 = m1*1.989e30
    M2 = m2*1.989e30
    G = 6.67408e-11
    c = 2.9979e8
    A = a0*1.496e11
    ECC = e0
    seconds = 3.154e7
    g = 6.086768e-11

    def binary(y, t):
        a, e = y
        dadt = -64./5.*G**3*M1*M2*(M1+M2)/(c**5.*a**3.*(1-e**2.)**3.5) * (1 + 73/24.*e**2. + 37/96.*e**4.)
        dedt = -304./15.*e*G**3*M1*M2*(M1+M2)/(c**5*a**4*(1-e**2.)**2.5) * (1 + 121./304*e**2.)

        dydt = [dadt, dedt]
        return dydt
    y0 = [A,ECC]
    t_start = 1.e-2
    delta = np.logspace(np.log10(t_start),np.log10(t_f),1000)
    t = np.zeros(1000)
    for i in range(1,len(delta)):
        #print t_f, delta[i], t_f + t_start - delta[len(delta)-i-1]
        t[i] = (t_f + t_start - delta[len(delta)-i-1])*seconds
    #t = np.linspace(0.0,t_f*seconds, 1000.)
    sol = odeint(binary, y0, t)
    A_F = sol[:,0]
    E_F = sol[:,1]
    #for i in range(len(delta)):
    #   print delta[i]
    #for i in range(len(A_F)):
        #      print t[i]/seconds, A_F[i]/1.5e11, E_F[i]
    #return t/seconds, A_F/1.496e11, E_F
    return A_F[-1]/1.496e11, E_F[-1]

def t_inspiral_2(a0,e0,m1,m2,t_flag=0, array=0, LIGO=0):
    """
    Computes inspiral time in years

    if t_flag > 0 (i.e. if you give it a t_flag value in years), this computes the SMA and ecc at t=t_flag

    if array=1, output the full a and e array from t=0 to t=t_flag

    if LIGO =1, integrate ALL the way to inspiral
    """
    ### DO NOT CHANGE!!! YOU CAN MESS WITH LISA_calculations.py if you want but not this one!!!
    #################

    M1 = m1*1.989e30
    M2 = m2*1.989e30
    G = 6.67408e-11
    c = 2.9979e8
    A = a0*1.496e11
    ECC = e0
    TIME = 0
    seconds = 3.154e7

    a_array = []
    ecc_array = []
    t_array = []

    if LIGO == 1:
        DELTA_T_MIN = 1.e-7
    if LIGO == 2:
        DELTA_T_MIN = 1.e-10
    if LIGO == 0:
        DELTA_T_MIN = 0.1

    t_guess = inspiral_time_peters(a0,e0,m1,m2,af=0)*1.e9  # in years
    if t_guess > 1.e25:
        if t_flag == 0:
            print('Inspiral time too long')
            return t_guess
        else:
            print('Inspiral time too long')
            return a0, e0
    t_max = int(np.round(t_guess*1.2))
    #delta = int(t_max/100000.)
    if t_max == 0:
        t_max = t_guess*1.2
    
    #print delta/seconds
    #if t_flag > 0:
    #   t_start = t_flag*.0001
    #   t_max = t_flag*1.2
    #else:
    #   t_start = t_max*.0001
    t_start = t_max*0.0001
    delta = np.logspace(np.log10(t_start),np.log10(t_max),10000)
    t = np.zeros(10000)
    for i in range(1,len(delta)):
        #print t_f, delta[i], t_f + t_start - delta[len(delta)-i-1]
        t[i] = (t_max + t_start - delta[len(delta)-i-1])*seconds

    a = A
    e = ECC

    for i in range(len(t)):
        DELTA = t[i+1] - t[i]
        #print i, a/1.496e11, e, t[i]/seconds, DELTA/seconds, 'v1'
        a0 = a
        e0 = e
        a_array.append(a/1.496e11)
        ecc_array.append(e)
        t_array.append(t[i]/seconds)
        a = a-64./5.*G**3*M1*M2*(M1+M2)/(c**5.*a**3.*(1-e**2.)**3.5) * (1 + 73/24.*e**2. + 37/96.*e**4.)*DELTA  
        e = e-304./15.*e*G**3*M1*M2*(M1+M2)/(c**5*a**4*(1-e**2.)**2.5) * (1 + 121./304*e**2.)*DELTA
        if t_flag > 0 and t[i]/seconds > t_flag:
            SMA = a0/1.486e11
            ECCENTRICITY = e0
            DELTA_T_FINAL = 1.e-10
            t_FINAL = t[i]/seconds
            break
        if a/1.496e11 < 0 or e < 0:
            t_FINAL = t[i]/seconds
            #DELTA_T_FINAL = (t[i] - t[i-1])/seconds
            DELTA_T_FINAL = inspiral_time_peters(a0/1.486e11,e0,m1,m2,af=0)*1.e9*1.2
            #print 'break1', a0/1.496e11, e0, t_FINAL, DELTA_T_FINAL
            break

    if DELTA_T_FINAL > DELTA_T_MIN:
        flag = 0
        while flag == 0:

            T_START = t_FINAL
            #t_max = t_START + DELTA_T_FINAL
            #t = np.linspace(T_START*seconds,t_max*seconds,100)
            t_max = DELTA_T_FINAL
            #print T_START, t_max, 'look here'
            t_start = DELTA_T_FINAL*0.0001
            delta = np.logspace(np.log10(t_start),np.log10(t_max),1000) 
            t = np.zeros(1000)
            #print T_START, t_max, np.log10(T_START),np.log10(t_max), 'look here'
            for i in range(1,len(delta)):
                #print t_f, delta[i], t_f + t_start - delta[len(delta)-i-1]
                t[i] = (t_max + t_start - delta[len(delta)-i-1])*seconds        
                #print i, t[i]/seconds, delta[i], 'here'
            a = a0
            e = e0
            for i in range(len(t)):
                DELTA = t[i+1] - t[i]
                #print i, a/1.496e11, e, t[i]/seconds, DELTA/seconds, 'v2'
                a0 = a
                e0 = e
                a_array.append(a/1.496e11)
                ecc_array.append(e)
                t_array.append(t[i]/seconds+T_START)
                a = a-64./5.*G**3*M1*M2*(M1+M2)/(c**5.*a**3.*(1-e**2.)**3.5) * (1 + 73/24.*e**2. + 37/96.*e**4.)*DELTA
                e = e-304./15.*e*G**3*M1*M2*(M1+M2)/(c**5*a**4*(1-e**2.)**2.5) * (1 + 121./304*e**2.)*DELTA
                if t_flag > 0 and (t[i]/seconds + T_START) > t_flag:
                    SMA = a0/1.486e11
                    ECCENTRICITY = e0
                    DELTA_T_FINAL = 1.e-10
                    t_FINAL = t_FINAL + T_START
                    break
                if a/1.496e11 < 0 or e < 0:
                    t_FINAL = t[i]/seconds
                    #DELTA_T_FINAL = (t[i] - t[i-1])/seconds
                    t_FINAL = t_FINAL + T_START
                    if t_flag > 0:
                        DELTA_T_FINAL = (t_flag - t_FINAL)*1.2
                    else:
                        DELTA_T_FINAL = inspiral_time_peters(a0/1.486e11,e0,m1,m2,af=0)*1.e9*1.2
                    #print 'break1', a0/1.496e11, e0, t_FINAL, DELTA_T_FINAL
                    break
            if DELTA_T_FINAL < DELTA_T_MIN:
                #print 'done!', t_FINAL, DELTA_T_FINAL
                flag = 1
    if t_flag == 0 and array == 0:
        return t_FINAL
    if t_flag > 0 and array == 0:
        return SMA, ECCENTRICITY
    if array == 1:
        return t_array, a_array, ecc_array

def snr_outside(f, hn, mode_vals, averaging_factor=16/5.):
    snr_squared_per_mode = averaging_factor*np.trapz(1./f*(mode_vals/hn[np.newaxis,:])**2, x=f)
    #check on sqrt(2) in paper also and BOWIE
    return np.sqrt(snr_squared_per_mode.sum())

def snr(m1, m2, a0, e0, d, n_max, fn, hn_lisa):

    z = zAtLuminosityDistance(d) # takes d in Mpc

    #t_insp = LISA.t_inspiral_2(a0,e0,m1,m2)
    t, a, e = t_inspiral_2(a0,e0,m1,m2,t_flag=10.0, array=1, LIGO=0)
    t = np.asarray(t)
    a = np.asarray(a)
    e = np.asarray(e)

    f_orb = np.sqrt((m1+m2)/(4.*np.pi**2*a**3))

    interp_funcs = []
    for n in np.arange(1,n_max):
        hcn_array = []
        #for j in range(len(a)):
            #print m1, m2, a[j], e[j], d, n
        hcn = hn(m1,m2,a,e,d,n)
        #hcn_array.append(hcn)
        #hcn = np.asarray(hcn_array)
        f_orb_new = np.logspace(np.log10(f_orb[0]), np.log10(f_orb[-1]), 100)

        interp_funcs.append(interp1d(n*f_orb_new, np.interp(f_orb_new, f_orb, hcn), bounds_error=False, fill_value=1e-60))

    mode_vals = np.asarray([interp(fn) for interp in interp_funcs])

    #overall_snr_sum_inside = snr_inside(fn, hn, mode_vals)
    overall_snr_sum_outside = snr_outside(fn, hn_lisa, mode_vals)

    """
    plt.plot(t/ct.c/ct.Julian_year, e/e0, label='e')
    plt.plot(t/ct.c/ct.Julian_year, a/a0, label='a')
    plt.legend()
    plt.show()
    """
    return overall_snr_sum_outside

def t_inspiral_Hansen(a0,e0,m1,m2,array=0):
    """
    Computes inspiral time in years

    if t_flag > 0 (i.e. if you give it a t_flag value in years), this computes the SMA and ecc at t=t_flag

    if array=1, output the full a and e array from t=0 to t=t_flag

    if LIGO =1, integrate ALL the way to inspiral
    """
    M1 = m1*1.989e30
    M2 = m2*1.989e30
    G = 6.67408e-11
    c = 2.9979e8
    A0 = a0*1.496e11
    ECC = e0
    TIME = 0
    seconds = 3.154e7

    c0 = A0*(1.-e0**2.)/e0**(12./19.) * (1.+121./304.*e0**2.)**(-870./2299.)
    R_1 = c0/2. * (1.+121./304.)**(870./2299.)
    Delta_E_GW = 85.*np.pi/(12.*np.sqrt(2.)) * G**(7./2.)/c**5.*(M1*M2)**2.*(M1+M2)**(0.5)/R_1**(7./2.)
    a_1 = G*M1*M2/(2.*Delta_E_GW)
    Beta = 64./5.*G**3.*M1*M2*(M1+M2)/c**5.
    T_INSP = 768./425.*a_1**4./(4.*Beta)*(R_1/a_1)**(7./2.)*(2.)**(7./2.)/seconds  # in years

    #print c0, R_1, Delta_E_GW, Beta
    #print a0, a_1/1.5e11
    #print 't_insp =', T_INSP, 'a_1 =', a_1/1.5e11

    a_array = []
    ecc_array = []
    t_array = []

    DELTA_T_MIN = 0.1

    t_max = T_INSP*10.0
    #delta = int(t_max/100000.)

    t_start = t_max*1.e-10
    t_arr = np.logspace(np.log10(t_start),np.log10(t_max),50000)
    t = []
    t.append(0.0)
    for j in range(len(t_arr)):
        t.append(t_arr[j])

    a = A0
    e = ECC

    for i in range(0,len(t)):
        DELTA = (t[i+1] - t[i])*seconds
        #print i, a/1.496e11, e, t[i]/seconds, DELTA/seconds, 'v1'
        a0 = a
        e0 = e
        #print t[i], a/1.5e11, e
        a_array.append(a/1.496e11)
        ecc_array.append(e)
        t_array.append(t[i])
        a = a + 64./5.*G**3*M1*M2*(M1+M2)/(c**5.*a**3.*(1-e**2.)**3.5) * (1 + 73/24.*e**2. + 37/96.*e**4.)*DELTA
        e = e + 304./15.*e*G**3*M1*M2*(M1+M2)/(c**5*a**4*(1-e**2.)**2.5) * (1 + 121./304*e**2.)*DELTA
        if e >= 0.9999:
            t_FINAL = t[i]
            DELTA_T_FINAL = (t[i] - t[i-1])
            break
    #print 't_insp =', T_INSP, 'a_1 =', a_1/1.5e11

    if array == 0:
        return T_INSP
    if array == 1:
        return T_INSP, t_array, a_array, ecc_array

def fdotn_2(m0,m1,Porb,e,n,z):
    F = (1 + 73./24.*pow(e,2.) + 37./96.*pow(e,4.))*pow(1-e*e,-7./2.)
    M1 = m0*Msun
    M2 = m1*Msun
    f_orb = 1./Porb
    Mc = chirp_mass(M1,M2)
    fdotn = n * 96./(10.*np.pi) * (1./c**5.) * (G*Mc)**(5./3.) * (2.*np.pi*f_orb/(1.+z))**(11./3.) * F
    #fdotn = 1./(1+z) * 192./10.* n * f_orb * (G*(M1+M2)/(4.*np.pi**2.)*Porb**2.)**(-4./3.) * G**3./c**5.*M1*M2*(M1+M2)*F
    return fdotn

