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
	period=twopi*np.sqrt((semimajor*AU)**3/(mtot*GMsun*10**6)) ##in seconds
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

def SD_converter(Epsilon, Rsun):
    """takes the surface brightness (1/pc^2) and converts to surface brightness (1/arcsec^2)"""
    SD_arcsec = Epsilon/pc_to_arcsec(1, Rsun)**2
    return SD_arcsec


def dms2degree(dms_value):
    dms_str = dms_value.split(':')
    if dms_str[0][0]=='-':
        degree = float(dms_str[0])-float(dms_str[1])/60.-float(dms_str[2])/3600.
    else:
    	degree = float(dms_str[0])+float(dms_str[1])/60.+float(dms_str[2])/3600.

    return degree


def pm2vel(distance, pm):  ##pm in mas/yr, vel in km/s, distance in kpc
    return 4.74*distance*pm


def vel2pm(distance, vel):  ##pm in mas/yr, vel in km/s, distance in kpc
    return vel/4.74/distance


def pc2arcsec(distance, rpc):  ##distance in kpc
    return rpc/(distance*1000.)*(180./np.pi)*3600


def arcsec2pc(distance, r_arcsec):  ##distance in kpc
    return (r_arcsec/3600.)*(np.pi/180.)*(distance*1000.)


