import numpy as np
from glob import glob
import collections
from collections import Counter
import os,sys
import subprocess
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
#import seaborn as sns
import gzip
import math
import re
import history_maker_full4 as hi4
import history_maker_full5 as hi5
import history_cmc as hic
import dynamics as dyn
import scripts3
import scripts1
import scripts2
from scipy import stats
import ecc_calc as gwcalc



def find_Lgammaray_lastsnapshot(sourcedir, eta_gamma):
	dataobs=np.genfromtxt(sourcedir+'clusterproperty_nondissolved_last.dat')
	Mtot=np.array(dataobs[:,2])
	datamsp=np.genfromtxt(sourcedir+'msp_last.dat')
	#datapsr=np.genfromtxt(sourcedir+'normalpsr_last.dat')
	Bmsp=np.array(datamsp[:,4]); Pmsp=np.array(datamsp[:,5]); modelmsp=datamsp[:,0]
	#Bpsr=np.array(datapsr[:,4]); Ppsr=np.array(datapsr[:,5]); modelpsr=datapsr[:,0]

	Cscale=9.6*10**33  ##in erg/s
	Lgamma_msp=Cscale*(eta_gamma/0.2)*(Bmsp/10**8.5)**2*(3./(Pmsp*1000.))**4
	#Lgamma_psr=Cscale*(eta_gamma/0.2)*(Bpsr/10**8.5)**2*(3./(Ppsr*1000.))**4
	#print Lgamma_msp, Lgamma_psr

	Lgamma_tot=[]
	for i in range(108):
		ltot=0
		for j in range(len(modelmsp)):
			if modelmsp[j]==i:
				ltot+=Lgamma_msp[j]
		#for k in range(len(modelpsr)):
		#	if modelpsr[k]==i:
		#		ltot+=Lgamma_psr[k]

		Lgamma_tot.append(ltot)


	return Lgamma_tot, Mtot


def find_Lgammaray_alltime(pathlist, start, end):
	Cscale=9.6*10**33  ##in erg/s
	sourcedir=np.genfromtxt(pathlist,dtype=str)
	for i in range(start, end):
		fgamma=open('/projects/b1011/syr904/projects/PULSAR/kickgrid_runs/MSPBHinGC/data/Lgamma/Lgamma_'+str(i).zfill(2), 'a+', 0)
		fgamma.write('#1.T(Myr) 2.Ltot 3.Nmsp\n')
		filestr=sourcedir[i]+'/initial'
		snaps=np.sort(glob(filestr+'.snap*.dat.gz'))
		for j in range(len(snaps)):
			Ltot=0; nmsp=0
			t=get_time(snaps[j])
			t_conv=conv('t', filestr+'.conv.sh')
			t=t*t_conv   ##in Myr
			with gzip.open(snaps[j], 'r') as fsnap:
				next(fsnap)
				next(fsnap)
				for line in fsnap:
					datasnap=line.split()
					if int(datasnap[14])==13:
						Pspin=twopi*yearsc/float(datasnap[59])  ##in sec
						B=float(datasnap[60])
						if Pspin<=0.03:
							nmsp+=1
							Ltot+=Cscale*(2./0.2)*(B/10**8.5)**2*(3./(Pspin*1000.))**4

					if int(datasnap[17])==13:
						Pspin=twopi*yearsc/float(datasnap[45])  ##in sec
						B=float(datasnap[47])
						if Pspin<=0.03:
							nmsp+=1
							Ltot+=Cscale*(2./0.2)*(B/10**8.5)**2*(3./(Pspin*1000.))**4

					if int(datasnap[18])==13:
						Pspin=twopi*yearsc/float(datasnap[46])  ##in sec
						B=float(datasnap[48])
						if Pspin<=0.03:
							nmsp+=1
							Ltot+=Cscale*(2./0.2)*(B/10**8.5)**2*(3./(Pspin*1000.))**4


			fgamma.write('%f %e %d\n'%(t, Ltot, nmsp))

			print j

		fgamma.close()

		print 'model=', i