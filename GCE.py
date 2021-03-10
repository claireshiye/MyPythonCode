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
#import history_maker_full4 as hi4
#import history_maker_full5 as hi5
import history_cmc as hic
import dynamics as dyn
import scripts3
import scripts1
import scripts2
from scipy import stats
import ecc_calc as gwcalc



def find_Lgammaray_lastsnapshot(filepath, eta_gamma):
	dataobs=np.genfromtxt(filepath+'clusterproperty_maingrid_last.dat')
	Mtot=np.array(dataobs[:,2]); status = dataobs[:,-1]
	datamsp=np.genfromtxt(filepath+'msp_maingrid_last.dat')
	datapsr=np.genfromtxt(filepath+'normalpsr_maingrid_last.dat')
	Bmsp=np.array(datamsp[:,4]); Pmsp=np.array(datamsp[:,5]); modelmsp=datamsp[:,0]
	Bpsr=np.array(datapsr[:,4]); Ppsr=np.array(datapsr[:,5]); modelpsr=datapsr[:,0]
	print(len(modelpsr))

	Cscale=9.6*10**33  ##in erg/s
	Lgamma_msp=Cscale*(eta_gamma/0.2)*(Bmsp/10**8.5)**2*(3./(Pmsp*1000.))**4
	Lgamma_psr=Cscale*(eta_gamma/0.2)*(Bpsr/10**8.5)**2*(3./(Ppsr*1000.))**4
	#print(len(Lgamma_msp), len(Lgamma_psr))

	L_tot=[]; L_msp = []; L_psr = []
	for ii in range(len(Mtot)):
		ltot = 0; lmsp = 0; lpsr = 0
		for j in range(len(modelmsp)):
			if int(modelmsp[j])==ii:
				lmsp += Lgamma_msp[j]
				ltot += Lgamma_msp[j]

		for k in range(len(modelpsr)):
			if int(modelpsr[k])==ii:
				lpsr += Lgamma_psr[k]
				ltot += Lgamma_psr[k]

		L_tot.append(ltot); L_msp.append(lmsp); L_psr.append(lpsr)


	return L_tot, L_msp, L_psr, Mtot


def find_Lgammaray_alltime(pathlist, start, end, eta_gamma):
    Cscale=9.6*10**33  ##in erg/s
    sourcedir=np.genfromtxt(pathlist,dtype=str)
    paths = sourcedir[:,0]; status = sourcedir[:,1]

    fgamma=open('/projects/b1095/syr904/projects/GCE/catalog/Lgamma_alltime_catalog.dat', 'a+')
    fgamma.write('#1.Model 2.T(Myr) 3.Lmsp\n')

    for ii in range(start, end):
        filestr=paths[ii]+'initial'
        t_conv=dyn.conv('t', filestr+'.conv.sh')

        Lmsp = []; time = []

        with open(filestr+'.morepulsars.dat', 'r') as fpsr:
            next(fpsr)
            for line in fpsr:
            	data=line.split()
            	if int(data[2])!=1:
            		Pspin=float(data[9])  ##in sec
            		B=float(data[7])
            		if Pspin<=0.03:
            		    Lmsp.append(Cscale*(eta_gamma/0.2)*(B/10**8.5)**2*(3./(Pspin*1000.))**4)
            		    time.append(float(data[1])*t_conv)

            	else:
                    if int(data[11])==13:
                    	Pspin=float(data[9])  ##in sec
                    	B=float(data[7])
                    	if Pspin<=0.03:
                    		Lmsp.append(Cscale*(eta_gamma/0.2)*(B/10**8.5)**2*(3./(Pspin*1000.))**4)
                    		time.append(float(data[1])*t_conv)

                    if int(data[12])==13:
                    	Pspin=float(data[10])  ##in sec
                    	B=float(data[8])
                    	if Pspin<=0.03:
                    		Lmsp.append(Cscale*(eta_gamma/0.2)*(B/10**8.5)**2*(3./(Pspin*1000.))**4)
                    		time.append(float(data[1])*t_conv)
        
        #print(time)

        #sums = {}
        #for key, value in zip(time,Lmsp):
        #    try:
        #        sums[key] += value
        #    except KeyError:
        #        sums[key] = value

        #print(sums, sums.keys())
        #allkey = sums.keys()
        #for key in sums:
        #	thetime = float(key)
        #	#print(thetime)
        #	theL = sums[key]
        #	fgamma.write('%d %f %e\n'%(ii, thetime, theL))
        
        allkey = list(Counter(time).keys())
        print(allkey)
        for x in range(len(allkey)):
        	theL = 0
        	for y in range(len(time)):
        		if time[y]==allkey[x]:
        			theL+=Lmsp[y]

        	fgamma.write('%d %f %e\n'%(ii, allkey[x], theL))

        print(ii)

    fgamma.close()
