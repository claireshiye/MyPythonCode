import numpy as np
import constants
import math
from scipy import integrate

def metallicity(m,s):
    """give the metalicity and tell how to convert
    options are:
        ztofe/h
        fe/htoz"""

    if s=='ztofe/h':
        m1 = log10(m/0.02)
    elif s=='fe/htoz':
        m1 = 0.02*10**m
    return m1

def find_rtidal(mc,vg=220.,rg=7400.):
    """takes 
    the cluster mass (mc) in solar mass, 
    the galactic circular velocity (vg) in km/s, 
    and 
    the galactocentric distance (rg) in pc
    returns: the tidal radius of the cluster in pc"""
    r_tidal = (constants.G * mc*constants.Msun / 2 / (vg*constants.km)**2.)**(1./3.) * (rg*constants.PC)**(2./3.) / constants.PC
    return r_tidal

def pc_to_arcsec(pc, R_sun):
    """takes size in pc and distance from sun (in kpc) and converts to arcsec"""
    arcsec = pc/(R_sun*1000)*206271
    return arcsec

def arcmin_to_pc(arcmin, R_sun):
    """takes size in arcmin and distance from sun (in kpc) and converts to pc"""
    pc = arcmin*(R_sun*1000)/3437.75
    return pc

def arcsec_to_pc(arcsec, R_sun):
    """takes size in arcsec and distance from sun (in kpc) and converts to pc"""
    arcmin = arcsec/60.0
    pc = arcmin*(R_sun*1000)/3437.75
    return pc

def SB_converter(I_v):
    """takes I_v (V-band Luminosity/pc^2) and converts to mu_v in V_mag/arcsec^2"""
    SB_arcsec = (21.572+4.83-2.5*np.log10(I_v))
    #SB_arcsec = (21.572-2.5*np.log10(I_v))
    return SB_arcsec

def SB_converter_tot(I):
        """takes I (Luminosity/pc^2) and converts to mu in mag/arcsec^2"""
        SB_arcsec = (21.572+4.74-2.5*np.log10(I))
        return SB_arcsec

def V_band_Lum_calculator(R,L):
    """calculates v-band luminosity from the radius (in units of RSUN) and total Luminosity (in units of LSUN) by integrating Plank's law"""
    sigma = 5.67*10**(-8) ### Stefan-Boltzmann constant in mks
    h = 6.626*10**(-34)
    c = 3.0*10**8
    k = 1.38*10**(-23)
    Lsun = 3.828*10**(26)
    L_mks = L*Lsun
    Rsun = 696000000
    R_mks = R*Rsun
    T = (L_mks/(4*math.pi*R_mks**2*sigma))**(.25)    ### Temperature of star as function of total Luminosity and radius
    f = lambda x: 8*math.pi*R_mks**2*h*c**2*x**(-5)*(math.exp(h*c/(x*k*T))-1)**(-1) 
    X = integrate.quad(f, 4*10**(-7), 7*10**(-7))
    L_v = X[0]/Lsun
    return L_v   # In units of Solar luminosity

def find_MS_turnoff(t):
    """given the time in Myr it finds the MS turn-off mass in Solar masses.  Very simple now.  Need to make the MS lifetime formula better. """
    t_yr = t*10**6
    lm = (9.921 - log10(t_yr))/3.6648
    m = 10**lm
    return(m)    

def find_t_ms(z, m):
    eta = log10(z/0.02)
    a1 = 1.593890e3+2.053038e3*eta+1.231226e3*eta**2.+2.327785e2*eta**3.
    a2 = 2.706708e3+ 1.483131e3*eta+ 5.772723e2*eta**2.+ 7.411230e1*eta**3.
    a3 = 1.466143e2 - 1.048442e2*eta - 6.795374e1*eta**2. - 1.391127e1*eta**3.
    a4 = 4.141960e-2 + 4.564888e-2*eta + 2.958542e-2*eta**2 + 5.571483e-3*eta**3.
    a5 = 3.426349e-1
    a6 = 1.949814e1 + 1.758178*eta - 6.008212*eta**2. - 4.470533*eta**3.
    a7 = 4.903830
    a8 = 5.212154e-2 + 3.166411e-2*eta - 2.750074e-3*eta**2. - 2.271549e-3*eta**3.
    a9 = 1.312179 - 3.294936e-1*eta + 9.231860e-2*eta**2. + 2.610989e-2*eta**3.
    a10 = 8.073972e-1


    m_hook = 1.0185 + 0.16015*eta + 0.0892*eta**2.
    m_HeF = 1.995 + 0.25*eta + 0.087*eta**2.
    m_FGB = 13.048*(z/0.02)**0.06/(1+0.0012*(0.02/z)**1.27)

    t_BGB = (a1+a2*m**4.+a3*m**5.5+m**7.)/(a4*m**2.+a5*m**7.)
    x = max([0.95,min([0.95-0.03*(eta+0.30103)]),0.99])
    mu = max(0.5, 1.0-0.01* max(a6/(m**a7) , a8+a9/m**a10))
    t_hook = mu*t_BGB

    t_MS = max(t_hook, x*t_BGB)

    return (t_MS)

def find_MS_TO(t, z, mguess):
    tguess = find_t_ms(z, mguess)
    #print mguess, tguess, (t-tguess)/t

    while (t-tguess)/t > 0.00005:
        mguess -= 0.00001
        tguess = find_t_ms(z, mguess)
        #print mguess, tguess, (t-tguess)/t
    mto = mguess
    return mto

def Integrate_SBP(filestring):
    data = np.loadtxt('/projects/b1011/kyle/cmc/branches/cmc-mpi/m22_project/rundir/'+filestring+'/initial.snap0010.2Dproj.dat')
    L = 0
    for i in range(0,587943):
        L = L + data[i,1]
        L_half = 0
    for i in range(1,587943):
        L_half = L_half + data[i,1]
        if L_half > L/2:
            return data[i-1,0], L_half, L
            break

def Lum_to_Mag(I_v):
        """takes I_v (V-band Luminosity/pc^2) and converts to mu_v in V_mag/arcsec^2"""
        SB_arcsec = (21.572-2.5*np.log10(I_v))
        return SB_arcsec
