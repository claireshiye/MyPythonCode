import numpy as np 
import pdb 
import matplotlib.pyplot as plt 
import scipy.constants as ct

from scipy.integrate import odeint, solve_ivp
from scipy.special import jv #bessel function of the first kind
from astropy.cosmology import Planck15 as cosmo
from scipy.interpolate import interp1d

import LISA_calculations as LISA

M_sun = 1.989e30

def c0_func(a0, e0):
	return a0*(1.-e0**2)/e0**(12./19.) * (1.+(121./304.)*e0**2)**(-870./2299.)

def beta_func(m1, m2):
	return 64./5.*m1*m2*(m1+m2)

def integrand(e):
	return e**(29./19.)*(1.+(121./304.)*e**2)**(1181./2299.)/(1.-e**2)**(3./2.)

def a_from_e(e, a0, e0):
	#peters equation 5.11
	c0 = c0_func(a0, e0)
	return c0*e**(12./19.)/(1.-e**2) * (1.+(121./304.)*e**2)**(870./2299.)

def t_merg(a0,e0, m1, m2, t_obs=10.0):
	beta = beta_func(m1, m2)
	c0 = c0_func(a0, e0)

	#G = 6.67e-11
	#M1 = m1/ct.G*ct.c**2.  #put masses in kg
	#M2 = m2/ct.G*ct.c**2.
	#c = 3.0e8
	#delta = 304./15.*e0*G**3*M1*M2*(M1+M2)/(c**5*a0**4*(1-e0**2.)**2.5) * (1 + 121./304*e0**2.)*10.*3.15e7
	#print 'delta=',delta

	#epsilon = delta*1.e-3
	epsilon = 1.e-8

	prefactor = 12./19.*c0**4./beta

	ecc_out, a_out, t_out = [], [], []

	t = 0.0
	num_step = 1
	while t/ct.c/ct.Julian_year<=t_obs:
		if e0-num_step*epsilon<0.0:
			break
		ecc_vals = np.linspace(e0-num_step*epsilon, e0, 1000)
		integrand_vals = integrand(ecc_vals)

		t = prefactor*np.trapz(integrand_vals, x=ecc_vals)

		t_out.append(t)
		a_out.append(a_from_e(e0-num_step*epsilon, a0, e0))
		ecc_out.append(e0-num_step*epsilon)
		#print(num_step, t/ct.c/ct.Julian_year)
		num_step +=1
	return  np.asarray(t_out), np.asarray(ecc_out), np.asarray(a_out)

def chirp_mass(m1, m2):
	return (m1*m2)**(3./5.)/(m1+m2)**(1./5.)

def F_func(e):
	return (1.+(73./24.)*e**2.+(37./96.)*e**4.)/(1.-e**2.)**(7./2.)

def g_func(n,e):
	return n**4./32. * ((jv(n-2., n*e) - 2.*e*jv(n-1., n*e) + 2./n*jv(n, n*e) + 2.*e*jv(n+1., n*e) - jv(n+2., n*e))**2. + (1.-e**2.)*(jv(n-2., n*e) - 2.*jv(n, n*e) + jv(n+2., n*e))**2. + 4./(3.*n**2.)*(jv(n, n*e))**2.)

def dEndfr(m1, m2, z, n, f, e):
	#eq 4 from orazio and samsing
	# takes f in detector frame
	Mc = chirp_mass(m1, m2)/(1.+z)  ## chirp mass in observed frame
	return np.pi**(2./3.)*Mc**(5./3.)/(3.*(f*(1.+z))**(1./3.))*(2./n)**(2./3.)*g_func(n,e)/F_func(e)

def hcn_func(m1, m2, z, n, f, e):
	D = cosmo.luminosity_distance(z).value*1e6*ct.parsec
	return 1./(np.pi*D)*np.sqrt(2.*dEndfr(m1, m2, z, n, f, e))

def snr_outside(f, mode_vals, noise_vals, averaging_factor=1.):
	snr_squared_per_mode = averaging_factor*np.trapz(1./f*(mode_vals/noise_vals)**2., x=f)
	#check on sqrt(2) in paper also and BOWIE
	return np.sqrt(1.)*np.sqrt(snr_squared_per_mode.sum())

def snr_inside(f, hn, mode_vals, averaging_factor=16/5.):
	hc = np.sum(mode_vals, axis=0)
	#check on sqrt(2) in paper also and BOWIE
	return np.sqrt(2.)*np.sqrt(averaging_factor*np.trapz(1./f*(hc/hn)**2, x=f))

def snr(m1, m2, a0, e0, d, n_max, noise_interp, t_final):

	z = LISA.zAtLuminosityDistance(d) # takes d in Mpc

	t, a, e = t_inspiral_2(a0/1.5e11,e0,m1,m2,t_flag=t_final, array=1, LIGO=0) # Calculates a, e evolution until t_final. t_final = LISA observation time in years or t_merger
	t = np.asarray(t)*3.15e7*3.0e8
	a = np.asarray(a)*1.5e11
	e = np.asarray(e)
	#for j in range(len(t)):
	#	print t[j]/3.0e8/3.15e7, a[j]/1.5e11, e[j]

	m1 = m1*M_sun*ct.G/ct.c**2.   # converting to c=G=1 units
        m2 = m2*M_sun*ct.G/ct.c**2.

        #t, e, a = t_merg(a0,e0, m1, m2)
	
	# wrong! f_orb = np.sqrt((m1+m2)/(4.*np.pi**2.*a**3.))/(1.+z)   # divide by (1+z) to convert from source frame to detector frame
	f_orb = np.sqrt((m1+m2)/(4.*np.pi**2.*a**3.))

	mode_vals = []
	freqs = []
	for n in np.arange(1.,n_max):
		hcn = hcn_func(m1, m2, z, n, n*f_orb, e)
		f_orb_new = np.logspace(np.log10(f_orb[0]), np.log10(f_orb[-1]), 100)

		freqs.append(n*f_orb_new*ct.c)
		mode_vals.append(np.interp(f_orb_new, f_orb, hcn))
		#interp_funcs.append(interp1d(n*f_orb_new*ct.c, np.interp(f_orb_new, f_orb, hcn), bounds_error=False, fill_value=1.e-60))

	mode_vals = np.asarray(mode_vals)
	freqs = np.asarray(freqs)
	noise_vals = noise_interp(freqs*(1.+z))
	#mode_vals = np.asarray([interp(fn) for interp in interp_funcs])
	#import pdb
	#pdb.set_trace()

	#overall_snr_sum_inside = snr_inside(fn, hn, mode_vals)
	overall_snr_sum_outside = snr_outside(freqs, mode_vals, noise_vals)

	"""
	plt.plot(t/ct.c/ct.Julian_year, e/e0, label='e')
	plt.plot(t/ct.c/ct.Julian_year, a/a0, label='a')
	plt.legend()
	plt.show()
	"""
	return overall_snr_sum_outside


"""
def dadt(a,e, m1, m2):
	#m1 and m2 are in geometrized units
	return -64./5. * m1*m2*(m1+m2)/(a**3*(1-e**2)**(7/2)) * (1.+(73./24.)*e**2+(37./96.)*e**4)

def dedt(a,e,m1,m2):
	return -304./15. * e*m1*m2*(m1+m2)/(a**4*(1-e**2)**(5/2)) * (1+(121./304.)*e**2)

#def ecc_gw_emission(y, t, m1, m2):
def ecc_gw_emission(t, y, m1, m2):
	#pdb.set_trace()
	a, e = y
	dydt = [dadt(a,e, m1, m2), dedt(a,e, m1, m2)]
	return dydt

m1 = 3e1*M_sun*ct.G/ct.c**2
m2 = 3e1*M_sun*ct.G/ct.c**2

f0_gw = 1e-1 #hz - gw frequency
f0_gw = f0_gw/ct.c #per meters

f0_orb = 0.5*f0_gw

a0 = ((m1+m2)/(4.*np.pi**2*f0_orb**2))**(1/3)

a0_normalizer = a0
e0 = 0.3

m1 = m1/a0
m2 = m2/a0

#t = np.linspace(0.0, .25*ct.Julian_year, 10001)*ct.c
t_span = (0.0, 10.0*ct.c/a0)

a0= a0/a0

y0 = [a0, e0]

#sol = odeint(ecc_gw_emission, y0, t, args=(m1, m2))
sol = solve_ivp(lambda t,y: ecc_gw_emission(t, y, m1, m2), t_span, y0)
pdb.set_trace()
"""

def t_inspiral_2(a0,e0,m1,m2, t_flag=0, array=0, LIGO=0, steps=1000):
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
        A = a0*1.496e11
        ECC = e0
        TIME = 0
        seconds = 3.154e7

        a_array = []
        ecc_array = []
        t_array = []

        if LIGO == 1:
                DELTA_T_MIN = 1.e-7
        else:
                DELTA_T_MIN = 0.1

#########
        t_guess = LISA.inspiral_time_peters(a0,e0,m1,m2,af=0)*1.e9  # in years
        #if t_guess > 1.e12:
                #if t_flag == 0:
                #        return t_guess
                #else:
                #        return a0, e0
	t_max = t_guess*1.2

        if t_flag > 0:
                t_start = t_flag*.0001
                t_max = t_flag*1.2
        else:
                t_start = t_max*.0001
########
	
        delta = np.logspace(np.log10(t_start),np.log10(t_max),steps)
        t = np.zeros(steps)
        for i in range(1,len(delta)):
                #print t_f, delta[i], t_f + t_start - delta[len(delta)-i-1]
                t[i] = (t_max + t_start - delta[len(delta)-i-1])*seconds

        a = A
        e = ECC

        for j in range(len(t)-1):
                DELTA = t[j+1] - t[j]
                a0 = a
                e0 = e
                a_array.append(a/1.496e11)
                ecc_array.append(e)
                t_array.append(t[j]/seconds)
                a = a-64./5.*G**3*M1*M2*(M1+M2)/(c**5.*a**3.*(1-e**2.)**3.5) * (1 + 73/24.*e**2. + 37/96.*e**4.)*DELTA
                e = e-304./15.*e*G**3*M1*M2*(M1+M2)/(c**5*a**4*(1-e**2.)**2.5) * (1 + 121./304*e**2.)*DELTA
                #print a, e
                if t_flag > 0 and t[j]/seconds > t_flag:
                        SMA = a0/1.486e11
                        ECCENTRICITY = e0
                        DELTA_T_FINAL = 1.e-10
                        t_FINAL = t[j]/seconds
                        break
                if a/1.496e11 < 0 or e < 0:
                        t_FINAL = t[j]/seconds
                        #DELTA_T_FINAL = (t[i] - t[i-1])/seconds
                        DELTA_T_FINAL = LISA.inspiral_time_peters(a0/1.486e11,e0,m1,m2,af=0)*1.e9*1.2
                        #print DELTA_T_FINAL
                        break

        #if 'DELTA_T_FINAL' in locals()==False:
                #DELTA_T_FINAL=1.0e-6
                #t_FINAL=t[-1]/seconds

        if DELTA_T_FINAL > DELTA_T_MIN:
                flag = 0
                while flag == 0:

                        T_START = t_FINAL
                        #t_max = t_START + DELTA_T_FINAL
                        #t = np.linspace(T_START*seconds,t_max*seconds,100)
                        t_max = DELTA_T_FINAL
                        t_start = DELTA_T_FINAL*0.0001
                        delta = np.logspace(np.log10(t_start),np.log10(t_max), steps)
                        t = np.zeros(steps)
                        for i in range(1,len(delta)):
                                t[i] = (t_max + t_start - delta[len(delta)-i-1])*seconds
                        a = a0
                        e = e0
                        for k in range(len(t)-1):
                                DELTA = t[k+1] - t[k]
                                a0 = a
                                e0 = e
                                a_array.append(a/1.496e11)
                                ecc_array.append(e)
                                t_array.append(t[k]/seconds+T_START)
                                a = a-64./5.*G**3*M1*M2*(M1+M2)/(c**5.*a**3.*(1-e**2.)**3.5) * (1 + 73/24.*e**2. + 37/96.*e**4.)*DELTA
                                e = e-304./15.*e*G**3*M1*M2*(M1+M2)/(c**5*a**4*(1-e**2.)**2.5) * (1 + 121./304*e**2.)*DELTA
                                if t_flag > 0 and (t[k]/seconds + T_START) > t_flag:
                                        SMA = a0/1.486e11
                                        ECCENTRICITY = e0
                                        DELTA_T_FINAL = 1.e-10
                                        t_FINAL = t_FINAL + T_START
                                        break
                                if a/1.496e11 < 0 or e < 0:
                                        t_FINAL = t[k]/seconds
                                        #DELTA_T_FINAL = (t[i] - t[i-1])/seconds
                                        t_FINAL = t_FINAL + T_START
                                        if t_flag > 0:
                                                DELTA_T_FINAL = (t_flag - t_FINAL)*1.2
                                        else:
                                                DELTA_T_FINAL = LISA.inspiral_time_peters(a0/1.486e11,e0,m1,m2,af=0)*1.e9*1.2
                                        break
                        if DELTA_T_FINAL < DELTA_T_MIN:
                                flag = 1
        if t_flag == 0 and array == 0:
                return t_FINAL
        if t_flag > 0 and array == 0:
                return SMA, ECCENTRICITY
        if array == 1:
                return t_array, a_array, ecc_array


if __name__ == '__main__':

	noise = np.genfromtxt('noise_curves/PL.txt', skip_header=4, names=True)
	fn = noise['f']
	hn = np.sqrt(fn)*noise['ASD']

	e0 = 0.9
	m1 = 30.0
	m2 = 30.0

	m1_trans = m1*M_sun*ct.G/ct.c**2
	m2_trans = m2*M_sun*ct.G/ct.c**2

	z = 0.01
	nmax=1000

	f0_gw = 1e-3 #hz - gw frequency
	f0_gw = f0_gw/ct.c #per meters

	f0_orb = 0.5*f0_gw

	a0 = ((m1_trans+m2_trans)/(4.*np.pi**2*f0_orb**2))**(1/3)

	print(snr(m1, m2, a0, e0, z, nmax, fn, hn))
