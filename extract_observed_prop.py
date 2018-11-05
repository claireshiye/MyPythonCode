import numpy as np
import matplotlib.pyplot as plt
import scripts
import glob
import subprocess
import constants
import scipy.integrate
import json

def read_units(filestring):
	"""reads from the supplied conv file and stores the physical units"""
	s = filestring+'.conv.sh'
	f=open(s,'r')
	count=0
	while count<10:
		a=f.readline()
		b=a.split('=')
		if b[0]=='massunitcgs':
			massunitcgs = float(b[1])
			count+=1
		elif b[0]=='massunitmsun':
			massunitmsun = float(b[1])
			count+=1
		elif b[0]=='mstarunitcgs':
			mstarunitcgs = float(b[1])
			count+=1
		elif b[0]=='mstarunitmsun':
			mstarunitmsun = float(b[1])
			count+=1
		elif b[0]=='lengthunitcgs':
			lengthunitcgs = float(b[1])
			count+=1
		elif b[0]=='lengthunitparsec':
			lengthunitparsec = float(b[1])
			count+=1
		elif b[0]=='timeunitcgs':
			timeunitcgs = float(b[1])
			count+=1
		elif b[0]=='timeunitsmyr':
			timeunitsmyr = float(b[1])
			count+=1
		elif b[0]=='nbtimeunitcgs':
			nbtimeunitcgs = float(b[1])
			count+=1
		elif b[0]=='nbtimeunitsmyr':
			nbtimeunitsmyr = float(b[1])
			count+=1
	f.close()
	units = []
	unittype = [('m_cgs', float), ('m_msun', float), ('mstar_cgs', float), ('mstar_msun', float), ('l_cgs', float), ('l_pc', float), ('t_cgs', float),('t_myr', float), ('nbt_cgs', float), ('nbt_myr', float)]
	units.append((massunitcgs, massunitmsun, mstarunitcgs, mstarunitmsun, lengthunitcgs, lengthunitparsec, timeunitcgs, timeunitsmyr, nbtimeunitcgs, nbtimeunitsmyr))
	units = np.array(units, dtype=unittype)
	#return (massunitcgs, massunitmsun, mstarunitcgs, mstarunitmsun, lengthunitcgs, lengthunitparsec, timeunitcgs, timeunitsmyr, nbtimeunitcgs, nbtimeunitsmyr)
	return units


def convert_to_3d(r, vr, vt, SEEDY=100):
        np.random.seed(SEEDY)
        if np.shape(r)==():  #single r, vr, and vt floats given
                #print 'came here'
                sintheta = np.random.uniform(low=-1., high=1.)
                phi = np.random.uniform(low=0., high=2.*np.pi)
                anglev = np.random.uniform(low=0., high=2.*np.pi)
        else:               #the full list of r, vr, and vt given
                #print 'came here too'
                r = np.array(r)
                vr = np.array(vr)
                vt = np.array(vt)

                sintheta = np.random.uniform(low=-1., high=1., size=len(r))
                phi = np.random.uniform(low=0., high=2.*np.pi, size=len(r))
                anglev = np.random.uniform(low=0., high=2.*np.pi, size=len(r))

        costheta = (1-sintheta**2.)**0.5

        rz = r*sintheta
        rx = r*costheta*np.cos(phi)
        ry = r*costheta*np.sin(phi)

        magv = (vr*vr + vt*vt)**0.5
        thetadot = np.cos(anglev) * vt/r
        phidot = np.sin(anglev)*vt/(r*costheta)

        vx = vr * costheta * np.cos(phi) - r * phidot * costheta * np.sin(phi) - r * thetadot * sintheta * np.cos(phi)
        vy = vr * costheta * np.sin(phi) + r * phidot * costheta * np.cos(phi) - r * thetadot * sintheta * np.sin(phi)
        vz = vr * sintheta + r * thetadot * costheta

        r3d = np.array([rx, ry, rz])
        v3d = np.array([vx, vy, vz])

        return r3d, v3d


def project_and_radially_sort(r3d, PROJ=(0,1)):
	r2d = (r3d[PROJ[0],:]**2. + r3d[PROJ[1],:]**2.)**0.5
	ind = np.argsort(r2d)
	return r2d, ind

def make_2D_projection(filestring, snapno, units, SEEDY=100, PROJ=(0,1)):
	#units = scripts.read_units(filestring)
	lpc = units[0]['l_pc']
	kms = 1e-5 * units[0]['l_cgs']/units[0]['nbt_cgs']
	t_myr = scripts.find_t_myr(filestring, snapno) 

	writefilename=filestring+'.snap'+snapno+'.2Dproj.dat'
	writefile=open(writefilename, 'w')
	writefile.write("#t=%g\n#1.r2D(pc) 2.Ltot(Lsun) 3.binflag 4.startype 5.L(Lsun) 6.startype0 7.startype1 8.L0(Lsun) 9.L1(Lsun) 10.Mtot(Msun) 11.M0(Msun) 12.M1(Msun) 13.id 14.id0 15.id1 16.rx(pc) 17.ry(pc) 18.rz(pc) 19.vx(km/s) 20.vy(km/s) 21.vz(km/s)\n" %(t_myr))

	#read the snapfile
	snapfile = filestring+'.snap'+snapno+'.dat.gz'
	colnos = (2, 7, 14, 15, 17, 18, 19, 20, 1, 8, 9, 3, 4, 0, 10, 11)
	#0-r, 1-binflag, 2-startype, 3-L, 4-startype0, 5-startype1, 6-L0, 7-L1, 8-Mtot, 9-M0, 10-M1, 11-vr, 12-vt, 13-id, 14-id0, 15-id1
	data = np.genfromtxt(snapfile, usecols=colnos)
	#data = np.genfromtxt(snapfile)
	r = data[:,0]*lpc
	vr = data[:,11]*kms
	vt = data[:,12]*kms
	r3d, v3d = convert_to_3d(r, vr, vt, SEEDY=SEEDY)
	r2d, ind = project_and_radially_sort(r3d, PROJ=PROJ)
	valid_line=0
	print 'N:', len(ind)
	for i in range(len(ind)):
		try:
			for j in range(len(data[ind[i]])):
				if str(data[ind[i],j])=='nan' or str(data[ind[i],j])=='inf':
					valid_line = 0
					print data[ind[i],:]
					raise StopIteration()
				else:
					valid_line = 1
		except StopIteration:
			pass
		if valid_line==1:
			if data[ind[i],1]==1.:
				Ltot = data[ind[i],6]+data[ind[i],7]
				Mtot = data[ind[i],9]+data[ind[i],10]
			else:
				Ltot = data[ind[i],3]
				Mtot = data[ind[i],8]
			writefile.write("%g %g %g %g %g %g %g %g %g %g %g %g %d %d %d %g %g %g %g %g %g\n" %(r2d[ind[i]], Ltot, data[ind[i],1], data[ind[i],2], Ltot, data[ind[i],4], data[ind[i],5], data[ind[i],6], data[ind[i],7], Mtot, data[ind[i],9], data[ind[i],10], data[ind[i],13], data[ind[i],14], data[ind[i],15], r3d[0,ind[i]], r3d[1,ind[i]], r3d[2,ind[i]], v3d[0,ind[i]], v3d[1,ind[i]], v3d[2,ind[i]], ))
	writefile.close()
	return r, vr, vt, r3d, v3d


def fit_king_cum_func(x, p0, p1):
	Lcum = np.pi * p0 * p1 * p1 * np.log( np.abs(1 + (x/p1)**2.))
	return Lcum

def fit_king_cum_curvefit(r2d, Lcum, p0guess=[1e5, 1.]):
	import scipy.optimize as opt
	#print 'came here'
	p_opt, p_cov = opt.curve_fit(fit_king_cum_func, r2d, Lcum, p0=p0guess)
	#print 'came here too'
	return p_opt, p_cov

def get_obs_props(filestring, snapno, FAC=1.):
	filename = filestring+'.snap'+snapno+'.2Dproj.dat'
	with open(filename, 'r') as f:
		try:
			for line in f:
				if line.rfind("#t=")>-1:
					t = float(line.split('=')[1].split()[0])
					raise StopIteration()
		except StopIteration:
			pass
	data = np.genfromtxt(filename, usecols=(0,1,9,18,19,20,))
	Mtot = np.sum(data[:,2])
	Mave = np.mean(data[:,2])
	Lcum = np.cumsum(data[:,1])
	Ltot = Lcum[-1]
	halfL = 0.5*Ltot
	try:
		for i in range(len(data)-1):
			L1, L2 = Lcum[i], Lcum[i+1]
			r1, r2 = data[i,0], data[i+1,0]
			if L1<=halfL<=L2:
				rhl = r1 + (r2-r1)*(halfL-L1)/(L2-L1)
				drhl = 0.5*(r2-r1)
				raise StopIteration()
	except StopIteration:
		N = i
		print 'found rhl:', rhl, drhl, i, N

	#guess_sigma = np.sum(data[:N, 1])/ (np.pi*data[N-1,0]**2.)
	#print N
	guess_sigma = Lcum[N]/(np.pi*data[N-1,0]**2.)
	guess_rc = 0.5*rhl
	p0guess = [guess_sigma, guess_rc]
	lim = rhl*FAC
	rselect, Lcumselect, count = [], [], 0
	while count<len(data) and data[count,0]<=lim:
		rselect.append(data[count,0])
		Lcumselect.append(Lcum[count])
		count+=1
		#print count
	rselect = np.array(rselect)
	Lcumselect = np.array(Lcumselect)
	p_opt, p_cov = fit_king_cum_curvefit(rselect, Lcumselect, p0guess=p0guess)
	p_err = np.sqrt(np.diag(p_cov))
	sigmac, sigmacerr, rc, rcerr = np.abs(p_opt[0]), np.abs(p_err[0]), np.abs(p_opt[1]), np.abs(p_err[1])
	print 'sigmac:', sigmac, 'sigmacerr:', sigmacerr, 'rc:', rc, 'rcerr:', rcerr
	vselect1, vselect2, vselect3 = [], [], []
	count = 0
	#print count, data[count,0]
	while data[count,0]<=rc and data[count,0]<=rhl:
		vselect1.append(data[count,3])
		vselect2.append(data[count,4])
		vselect3.append(data[count,5])
		count+=1
		#print count
	vselect1, vselect2, vselect3 = np.array(vselect1), np.array(vselect2), np.array(vselect3)
	vsigma = [np.std(vselect1), np.std(vselect2), np.std(vselect3)]
	vsigma = np.array(vsigma) 
	vsigmac = np.mean(vsigma)
	dvsigmac = np.std(vsigma)

	props = {'t': t,
		'rhl': rhl,
		'drhl': drhl,
		'rc': rc,
		'drc': rcerr,
		'sigmac': sigmac,
		'dsigmac': sigmacerr,
		'Ltot': Ltot,
		'M/L': Mtot/Ltot,
		'Mave': Mave,
		'vsigmac_rv': vsigmac,
		'dvsigmac_rv': dvsigmac
		}

	return Mtot, props
