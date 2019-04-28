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


def conv_dict(): return {'l':15, 't':19, 'm':7}    #?


def conv(unit,filepath):   # Returns the unit conversion multiplier given a simulation's *.conv.sh file and a unit (either 'l' or 't' or 'm')
    dict = conv_dict()
    from re import findall
    with open(filepath,'r') as f:
        head = [next(f) for x in xrange(24)]
    return float(findall('\d+[\.]?\d*',head[dict[unit]])[0])


def get_time(filepath):      # Returns the cluster's age for a given snapshot
    import gzip
    from re import findall
    with gzip.open(filepath,'r') as f: contents = f.readline()
    if not findall('\d+[\.]?\d*',contents):        # Returns time = 0 for snapshot files without a time header
        print 'snapshot empty'; return float(0)
    else: return float(findall('\d+[\.]?\d*',contents)[0])



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
    nbhlast=int(databh[2])

    with open(filedyn, 'r') as fdyn:
        for line in fdyn:pass
        lastdyn=line
    datadyn=lastdyn.split()
    ntotlast=int(datadyn[3])

    return nbhlast, ntotlast


##Find BH in random timestep
def find_NBH_NTOT(filestring, time, tconv):
    nbh=0; ntot=0; mass=0

    filebh=filestring+'.bh.dat'
    databh=np.genfromtxt(filebh)
    for i in range(len(databh[:,1])):
        if databh[:,1][i]*tconv>=time:
            nbh=databh[:,2][i]
            break

    filedyn=filestring+'.dyn.dat'
    datadyn=np.genfromtxt(filedyn)
    for j in range(len(datadyn[:,0])):
        if datadyn[:,0][j]*tconv>=time:
            ntot=datadyn[:,3][j]; mass=datadyn[:,4][j] ##mass in code unit
            break

    return nbh, ntot, mass


##Find Mtot, rc, rh, rho0 in the last timestep
def find_rcrh_mtotrho0_last(filestring):
    filedyn=filestring+'.dyn.dat'
    with open(filedyn, 'r') as fdyn:
        for line in fdyn: pass
        lastdyn=line
    datadyn=lastdyn.split()
    mtot=float(datadyn[4]); rc=float(datadyn[7]); rh=float(datadyn[20]); rhoc=float(datadyn[21])

    return mtot, rc, rh, rhoc    ##in code unit


##Find Mtot, rc, rh, rho0 at a random timestep
def find_rcrh_mtotrho0(filestring, time, tconv):
    filedyn=filestring+'.dyn.dat'
    with open(filedyn, 'r') as fdyn:
        next(fdyn)
        next(fdyn)
        for line in fdyn:
            datadyn=line.split()
            if float(datadyn[0])*tconv>=time: 
                mtot=float(datadyn[4]); rc=float(datadyn[7]); rh=float(datadyn[20]); rhoc=float(datadyn[21])
                break

    return mtot, rc, rh, rhoc    ##in code unit


##Find core mass in the last snapshot
def find_mc(filestring, r_c):   ##r_c in code units
    mc=0; mc_wobh=0
    snaps=np.sort(glob(filestring+'.snap*.dat.gz'))
    with gzip.open(snaps[-1], 'r') as flastsnap:
        for _ in xrange(2): next(flastsnap)
        for line in flastsnap:
            datalast=line.split()
            r=float(datalast[2]); m=float(datalast[1])
            binflag=int(datalast[7]); k=int(datalast[14]); k0=int(datalast[17]); k1=int(datalast[18]); m0=float(datalast[8]); m1=float(datalast[9])

            if r<=r_c: 
                mc+=m
                if binflag!=1 and k!=14:
                    mc_wobh+=m
                if binflag==1 and (k0!=14 and k1!=14):
                    mc_wobh+=m
                if binflag==1 and (k0!=14 and k1==14):
                    mc_wobh+=m0
                if binflag==1 and (k0==14 and k1!=14):
                    mc_wobh+=m1
            if r>r_c: break

    return mc, mc_wobh


##Find the number of NSs and NS binaries in the initial.ns.dat file at a random timestep
def find_Nns(filestring, timemin, timemax, tconv):
    datans=np.genfromtxt(filestring+'.ns.dat')
    for i in range(len(datans[:,0])):
        if datans[:,0][i]*tconv>=timemin and datans[:,0][i]*tconv<=timemax:
            nns=datans[:,1][i]; ndns=datans[:,7][i]; nnsbh=datans[:,8][i]
            break

    return nns, ndns, nnsbh


def find_Nns_last(filestring):
    datans=np.genfromtxt(filestring+'.ns.dat')
    nns=datans[-1,1]; ndns=datans[-1,7]; nnsbh=datans[-1,8]

    return nns, ndns, nnsbh


##Find Nbh, Mtot, Nns, Ndns
def find_clusterparameter_allmodel(pathlist, start, end):
    model=[]; NBH=[]; MTOT=[]; RC=[]; RH=[]; NNS=[]; NDNS=[]; NNSBH=[]; status=[]
    sourcedir=np.genfromtxt(pathlist, dtype=str)
    for i in range(start, end):
        filestr=sourcedir[i]+'/initial'
        t_conv=dyn.conv('t', filestr+'.conv.sh')
        l_conv=dyn.conv('l', filestr+'.conv.sh')
        m_conv=dyn.conv('m', filestr+'.conv.sh')

        ##Numbers at a certain time
        #Nbh, Ntot, Mtot=dyn.find_NBH_NTOT(filestr, 12000., t_conv)
        #Mtot, Rc, Rh, Rhoc=dyn.find_rcrh_mtotrho0(filestr, 12000., t_conv)
        #Nns, Ndns, Nnsbh=dyn.find_Nns(filestr, 12000., t_conv)

        ##Numbers at the last snapshot
        Nbh, Ntot=dyn.find_NBH_NTOT_last(filestr)
        Mtot, Rc, Rh, Rhoc=dyn.find_rcrh_mtotrho0_last(filestr)

        filedyn=filestr+'.dyn.dat'
        with open(filedyn, 'r') as fdyn:
            for line in fdyn: pass
            lastdyn=line
        datadyn=lastdyn.split()
        t_last=float(datadyn[0])*t_conv

        ##If the cluster dissolved, get the numbers for the last snapshot; if not, get the numbers within the last 3Gyr
        if t_last<12000.:
            Nns, Ndns, Nnsbh=dyn.find_Nns_last(filestr)
            status.append(1)
        else:
            Nns, Ndns, Nnsbh=dyn.find_Nns(filestr, 9000., 12000., t_conv)
            status.append(0)


        model.append(i); NBH.append(Nbh); MTOT.append(Mtot*m_conv); RC.append(Rc*l_conv); RH.append(Rh*l_conv)
        NNS.append(Nns); NDNS.append(Ndns); NNSBH.append(Nnsbh)

        print i


    np.savetxt(savepath+'/newruns/clusterproperty.dat', np.c_[model, NBH, MTOT, RC, RH, NNS, NDNS, NNSBH, status], fmt='%d %d %f %f %f %d %d %d %d', header='1.Model 2.Nbh 3.Mtot(Msun) 4.rc(pc) 5.rh(pc) 6.Nns 7.Ndns 8.Nnsbh 9.Dissolved?', delimiter='', comments='#')


