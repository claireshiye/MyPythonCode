import numpy as np
import math

GMsun=1.327*10**20 ##m3s-2
AU=1.496*10**13  ##cm
twopi=6.283185307179586
daysec=86400.

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
