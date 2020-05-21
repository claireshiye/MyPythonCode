import numpy as np
import math
import astropy

GMsun=1.327*10**20 ##m3s-2
AU=1.496*10**13  ##cm
twopi=6.283185307179586
daysec=86400.
yearsc=3.1557*10**7
mpctokm=3.086e+19

def find_init_conditions(thepath):
	s=thepath.split('/')
	n_star=float(s[-2])
	z=float(s[-3][1:])
	rg=int(s[-4][2:])
	rv=float(s[-5][2:])

	return n_star, z, rg, rv

def au_to_period(semimajor, m1, m2):
	mtot=m1+m2  ##in solar unit
	period=twopi*math.sqrt((semimajor*AU)**3/(mtot*GMsun*10**6)) ##in seconds
	period=period/(3600.*24.) ##in days

	return period
	

def period_to_au(period, m1, m2):
	mtot=m1+m2
	x1=(period*daysec/twopi)**2
	sma=(x1*mtot*GMsun*10**6)**(1./3.)/AU

	return sma


def chirpmass(m1, m2):
	a=(m1*m2)**(3./5.)
	b=(m1+m2)**(1./5.)
	M=a/b

	return M


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


def psr_deathline(P, B):  ##P in seconds and B in Gauss
	#print B/P**2
	if (B/P**2)>0.17*10**12:
		return "yes"
	else:
		return "no"


def ttoredshift(t):    
##t in Gyr; return the redshift
##Module astropy also has a very handy cosmology calculator for this. Check that the value from my calculations is the same as theirs.
	H0=69.6  ##km/sec/Mpc
	omega_m=0.286
	omega_vac=0.714

	A=(2.*H0**(-1)/3/omega_vac**0.5)*mpctokm/yearsc/10**9   ##in Gyr
	B=(omega_vac/omega_m)**0.5

	Z=pow(np.sinh(t/A)/B, -2./3.)-1

	return Z


def redshifttot(Z):
	##Return the corresponding t in Gyr for a redshift
	from astropy.cosmology import FlatLambdaCDM
	cosmo = FlatLambdaCDM(H0=69.6, Om0=0.286)
	return cosmo.age(Z).value


def ComovingDistance(redshift):
	from astropy.cosmology import FlatLambdaCDM
	cosmo = FlatLambdaCDM(H0=69.6, Om0=0.286)
	return cosmo.comoving_distance(redshift)



