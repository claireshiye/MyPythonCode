import numpy as np
from glob import glob
import collections
from collections import Counter
import os,sys
import subprocess
import matplotlib
import matplotlib.lines as mlines
#matplotlib.use('PDF')
import matplotlib.pyplot as plt
#import seaborn as sns
import gzip
import math
import re
import history_maker_full5
import history_cmc as hic

##Find observed rc, rh in models
def find_obsrcrh(snapshotobs):
        #snapobs=np.sort(glob(filepath+'/'+'initial.snap*.obs_params.dat'))
        dataobs=np.genfromtxt(snapshotobs)
        rc=dataobs[0, 7]; rhl=dataobs[0, 8]; t_Gyr=dataobs[0, 10]/1000.0; mass=dataobs[0, 12]

        return rc, rhl, t_Gyr, mass


def find_obsrcrh_lastsnap_allmodels(pathlist, start, end):
	pref='initial'
	sourcedir=np.genfromtxt(pathlist, dtype='|S')
	RC=[]; RHL=[]; T=[]; NBH=[]; NTOT=[]; M=[]
	for i in range(start, end):
		filepath=sourcedir[i]
		pref='initial'
		filestr=filepath+'/'+pref
		snapobs=np.sort(glob(filestr+'.snap*.obs_params.dat'))
		lastsnapobs=snapobs[-1]
		Rc, Rhl, T_Gyr, m=find_obsrcrh(lastsnapobs)
		RC.append(Rc); RHL.append(Rhl); T.append(T_Gyr); M.append(m)
        	Nbh, Ntot=find_NBH_NTOT_last(filestr)
        	NBH.append(float(Nbh)); NTOT.append(float(Ntot))
		
		print i

   	np.savetxt('/projects/b1011/syr904/projects/PULSAR/kickgrid_runs/paper/data/kickgrid_obsproperty_newmodel.dat', np.c_[T, NTOT, NBH, RC, RHL, M], fmt='%f %d %d %f %f %f', delimiter=' ', header='t_Gyr Ntot Nbh rc rhl M', comments='#')


##Find BH in the last timestep
def find_NBH_NTOT_last(filestring):
	filebh=filestring+'.bh.dat'
    	filedyn=filestring+'.dyn.dat'
    	with open(filebh, 'r') as fbh:
        	for line in fbh:pass
        	lastbh=line
    	databh=lastbh.split()
   	nbhlast=float(databh[2])

    	with open(filedyn, 'r') as fdyn:
        	for line in fdyn:pass
        	lastdyn=line
    	datadyn=lastdyn.split()
    	ntotlast=float(datadyn[3])

    	return nbhlast, ntotlast


##Find BH in random timestep
def find_NBH_NTOT(filestring, time):
    nbh=0; ntot=0; mass=0

    filebh=filestring+'.bh.dat'
    databh=np.genfromtxt(filebh)
    for i in range(len(databh[:,1])):
        if databh[:,1][i]==time:
            nbh=databh[:,2][i]

    filedyn=filestring+'.dyn.dat'
    datadyn=np.genfromtxt(filedyn)
    for j in range(len(datadyn[:,0])):
        if datadyn[:,0][j]==time:
            ntot=datadyn[:,3][j]; mass=datadyn[:,4][j] ##mass in code unit

    return nbh, ntot, mass



