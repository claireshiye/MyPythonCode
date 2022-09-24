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
import scipy
from scipy.interpolate import interp1d
from scipy import stats
import scipy.integrate as integrate
import ecc_calc as gwcalc
import ns
import pandas as pd

yearsc=31557600.
twopi=6.283185307179586
Gconst=6.674*10**-8 ##cm3*g-1*s-2
clight=3*10**10 ##cm/s
Msun=2*10**33 ##gram
AU=1.496*10**13  ##cm
PC=3.086*10**18  ##cm
Kconst=9.87*10**-48 ##yr/G^2
Lsun=4.02*10**16 ##mJy*kpc^2

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


def find_Lgammaray_alltime(pathlist, start, end, eta_gamma=0.1):
    Cscale=9.6*10**33  ##in erg/s
    #sourcedir=np.genfromtxt(pathlist,dtype=str)
    #paths = sourcedir[:,0]; status = sourcedir[:,1]
    paths = ['/projects/b1091/CMC_Grid_March2019/rundir/rv1/rg2/z0.02/4e5/']

    #fgamma=open('/projects/b1095/syr904/N1e7_rg8_f5_z0.002/Lgamma_alltime_N1e7.dat', 'a+')
    #fgamma.write('#1.Model 2.T(Myr) 3.Lmsp\n')

    for ii in range(start, end):
        print(paths[ii])
        filestr=paths[ii]+'initial'
        t_conv=dyn.conv('t', filestr+'.conv.sh')

        s=paths[ii].split('/')
        n_star=s[-2]
        z=s[-3][1:]
        rg=s[-4][2:]
        rv=s[-5][2:]

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
                        print(Pspin,B)

                else:
                    if int(data[11])==13:
                        Pspin=float(data[9])  ##in sec
                        B=float(data[7])
                        if Pspin<=0.03:
                            Lmsp.append(Cscale*(eta_gamma/0.2)*(B/10**8.5)**2*(3./(Pspin*1000.))**4)
                            time.append(float(data[1])*t_conv)
                            print(Pspin,B)

                    if int(data[12])==13:
                        Pspin=float(data[10])  ##in sec
                        B=float(data[8])
                        if Pspin<=0.03:
                            Lmsp.append(Cscale*(eta_gamma/0.2)*(B/10**8.5)**2*(3./(Pspin*1000.))**4)
                            time.append(float(data[1])*t_conv)
                            print(Pspin,B)
        
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

        #fgamma=open('/projects/b1095/syr904/projects/GCE/catalog/data_lgamma/model_rv'+rv+'_rg'+rg+'_z'+z+'_'+n_star+'.txt', 'a+')
        #fgamma.write('#1.Model 2.T(Myr) 3.Lmsp\n')
        
        allkey = list(Counter(time).keys())
        #print(allkey)
        for x in range(len(allkey)):
        	theL = 0
        	for y in range(len(time)):
        		if time[y]==allkey[x]:
        			theL+=Lmsp[y]

        	#fgamma.write('%d %f %e\n'%(ii, allkey[x], theL))

        print(ii)

    #fgamma.close()


##Need to work on this
def find_Lgammaray_1013Gyr(filepath, eta_gamma):
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


def find_Lave_alltime_allmsp(pathlist, start, end, eta_gamma=0.1):
    Cscale=9.6*10**33  ##in erg/s
    sourcedir=np.genfromtxt(pathlist,dtype=str)
    paths = sourcedir[:,0]; status = sourcedir[:,1]

    #fgamma=open('/projects/b1095/syr904/N1e7_rg8_f5_z0.002/Lgamma_alltime_N1e7.dat', 'a+')
    #fgamma.write('#1.Model 2.T(Myr) 3.Lmsp\n')

    Lmsp = []; time = []
    for ii in range(start, end):
        print(paths[ii])
        filestr=paths[ii]+'initial'
        t_conv=dyn.conv('t', filestr+'.conv.sh')

        s=paths[ii].split('/')
        n_star=s[-2]
        z=s[-3][1:]
        rg=s[-4][2:]
        rv=s[-5][2:]
        
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
                        #print(Pspin,B)

                else:
                    if int(data[11])==13:
                        Pspin=float(data[9])  ##in sec
                        B=float(data[7])
                        if Pspin<=0.03:
                            Lmsp.append(Cscale*(eta_gamma/0.2)*(B/10**8.5)**2*(3./(Pspin*1000.))**4)
                            time.append(float(data[1])*t_conv)
                            #print(Pspin,B)

                    if int(data[12])==13:
                        Pspin=float(data[10])  ##in sec
                        B=float(data[8])
                        if Pspin<=0.03:
                            Lmsp.append(Cscale*(eta_gamma/0.2)*(B/10**8.5)**2*(3./(Pspin*1000.))**4)
                            time.append(float(data[1])*t_conv)
                            #print(Pspin,B)
        
        
        print(ii)
    print(np.average(Lmsp))


##Find the gamma-ray luminosity for individual MSPs for models with rg=2 in the catalog model for all time steps
def find_Lgamma_allmsp_alltime(pathlist, start, end, eta_gamma=0.1):
    Cscale=9.6*10**33  ##in erg/s
    sourcedir=np.genfromtxt(pathlist,dtype=str)
    paths = sourcedir[:,0]; status = sourcedir[:,1]

    Lmsp = []; time = []
    for ii in range(start, end):
        filestr=paths[ii]+'initial'
        t_conv=dyn.conv('t', filestr+'.conv.sh')

        s=paths[ii].split('/')
        n_star=float(s[-2])
        z=float(s[-3][1:])
        rg=float(s[-4][2:])
        rv=float(s[-5][2:])

        if rg>2.: continue

        print(paths[ii])

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

    
    f=open('/projects/b1095/syr904/projects/GCE/catalog/lgamma_allmsp_alltime_rg2.txt', 'a+')
    t_array = np.linspace(0, 11600, 117)
    for xx in range(len(t_array)-1):
        f.write("%f " % t_array[xx])
        for yy in range(len(time)):
            if t_array[xx]<=time[yy]<t_array[xx+1]:
                #Lmsp_timestep.append(Lmsp[yy])
                f.write("%e " % Lmsp[yy])
 
        f.write("\n")

    f.close()
    #print([t_array[xx]]+Lmsp_timestep)


##Find the ejected MSPs for models with rg=2 in the catalog model for all time steps
def find_allmsp_alltime_ejected(pathlist, start, end):
    sourcedir=np.genfromtxt(pathlist,dtype=str)
    paths = sourcedir[:,0]; status = sourcedir[:,1]

    t_array = np.linspace(0, 11500, 116)

    N_label = ['2e5','4e5','8e5','16e5']
    nmsp_tot = []
    for kk in range(len(N_label)):
        n_model = 0
        nmsp_tot.append([])
        nmsp = []; alltime = []
        for ii in range(start, end):
            filestr=paths[ii]+'initial'
            t_conv=dyn.conv('t', filestr+'.conv.sh')

            s=paths[ii].split('/')
            n_star=float(s[-2])
            z=float(s[-3][1:])
            rg=float(s[-4][2:])
            rv=float(s[-5][2:])

            if rg>2.: 
                continue
            if n_star!=float(N_label[kk]):
                continue

            print(paths[ii])
            n_model+=1
            print(n_model)
        
            datans = np.genfromtxt(filestr+'.esc_ns.dat')
            if len(datans)==0:
                #print('no ejection')
                continue
            elif len(datans)==1:
                #print('one ejection')
                continue
            else:
                #print('multiple ejections')
                times = datans[:,0]
                n_psr = datans[:,5]; n_msp = datans[:,6]
                n_nswd = datans[:,9]; n_nsms = datans[:,10]
                n_nsg = datans[:,11]; n_mtb = datans[:,4]
                n_ej = n_msp#+n_mtb#+n_nswd+n_nsms+n_nsg
                #print(times)

            nmsp=nmsp+list(n_ej);  alltime=alltime+list(times)

        print(np.max(nmsp), np.sum(nmsp))
        n_count = 0
        for aa in range(len(nmsp)):
            if nmsp[aa]>0: 
                print(alltime[aa])
                n_count+=1
        print(n_count)

        for xx in range(len(t_array)-1):
            nmsp_n = 0
            for yy in range(len(alltime)):
                if t_array[xx]<=alltime[yy]<t_array[xx+1]:
                    nmsp_n+=nmsp[yy]

            if nmsp_n>0: print(nmsp_n)
            nmsp_n_temp = float(nmsp_n)/float(n_model)
            nmsp_tot[kk].append(nmsp_n_temp)
   
    np.savetxt('/projects/b1095/syr904/projects/GCE/catalog/allmsp_eject_rg2.txt',np.c_[t_array[:-1],nmsp_tot[0],nmsp_tot[1],nmsp_tot[2],nmsp_tot[3]], fmt = '%f %f %f %f %f', header = '1.T(Myr) 2.N_2e5 3.N_4e5 4.N_8e5 5.N_16e5', comments = '#', delimiter = '')



def find_Lgammaray_lastsnapshot_onerun(filepath, eta_gamma = 0.1):
    datamsp=np.genfromtxt(filepath+'msp_last.dat')
    datapsr=np.genfromtxt(filepath+'normalpsr_last.dat')
    Bmsp=np.array(datamsp[:,4]); Pmsp=np.array(datamsp[:,5])
    Bpsr=np.array(datapsr[:,4]); Ppsr=np.array(datapsr[:,5])

    Cscale=9.6*10**33  ##in erg/s
    Lgamma_msp=Cscale*(eta_gamma/0.2)*(Bmsp/10**8.5)**2*(3./(Pmsp*1000.))**4
    Lgamma_psr=Cscale*(eta_gamma/0.2)*(Bpsr/10**8.5)**2*(3./(Ppsr*1000.))**4

    lmsp = sum(Lgamma_msp); lpsr = sum(Lgamma_psr); ltot = lmsp+lpsr

    return ltot, lmsp, lpsr


def find_Lgammaray_alltime_onerun(modelpath, eta_gamma = 0.1): 
    Cscale=9.6*10**33  ##in erg/s

    fgamma=open('/jet/home/yec/behemoth/Lgamma_alltime_behemoth.dat', 'a+')
    fgamma.write('#1.T(Myr) 2.Ltot(erg/s) 3.Lmsp(erg/s)\n')

    snaps_raw=glob(modelpath+'*.snap*.dat.gz')
    snapno = []
    for xx in range(len(snaps_raw)):
        strings = snaps_raw[xx].split('.')
        no = re.findall(r'\d+', strings[1])[0]
        if no != '0000':
            no = no.lstrip('0')
        else:
            no = '0'
        snapno.append(int(no))

    snapnos, snaps = (list(t) for t in zip(*sorted(zip(snapno, snaps_raw))))
    
    t_conv = ns.conv('t', '/ocean/projects/ast190019p/shared/FIRE_balls/1e6-1e7/26859/26859-bigger.conv.sh')
    for ii in range(len(snaps)):
        Ltot = 0; Lmsp = 0
        time = ns.get_time(snaps[ii])
        t_myr = time*t_conv
        with gzip.open(snaps[ii], 'r') as fsnap:
            next(fsnap)
            next(fsnap)

            for line in fsnap:
                datasnap = line.split()
                if int(datasnap[7])!=1:
                    if int(datasnap[14])==13:
                        Pspin=twopi*yearsc/float(datasnap[59])
                        B = float(datasnap[60])
                        deathcut=(Pspin**2)*(0.17*10**12)
                        if deathcut<B:
                            Ltot+=Cscale*(eta_gamma/0.2)*(B/10**8.5)**2*(3./(Pspin*1000.))**4
                            if Pspin<=0.03:
                                Lmsp+=Cscale*(eta_gamma/0.2)*(B/10**8.5)**2*(3./(Pspin*1000.))**4

                else:
                    if int(datasnap[17])==13:
                        Pspin0=twopi*yearsc/float(datasnap[45])
                        B0 = float(datasnap[47])
                        deathcut0=(Pspin0**2)*(0.17*10**12)
                        if deathcut0<B0:
                            Ltot+=Cscale*(eta_gamma/0.2)*(B0/10**8.5)**2*(3./(Pspin0*1000.))**4
                            if Pspin0<=0.03:
                                Lmsp+=Cscale*(eta_gamma/0.2)*(B0/10**8.5)**2*(3./(Pspin0*1000.))**4
                    if int(datasnap[18])==13:
                        Pspin1=twopi*yearsc/float(datasnap[46])
                        B1 = float(datasnap[48])
                        deathcut1=(Pspin1**2)*(0.17*10**12)
                        if deathcut1<B1:
                            Ltot+=Cscale*(eta_gamma/0.2)*(B1/10**8.5)**2*(3./(Pspin1*1000.))**4
                            if Pspin1<=0.03:
                                Lmsp+=Cscale*(eta_gamma/0.2)*(B1/10**8.5)**2*(3./(Pspin1*1000.))**4
        
        fgamma.write('%f %e %e\n'%(t_myr, Ltot, Lmsp))
        print(snaps[ii])

    fgamma.close()



##Find MSPs at the time of their formation in models without morepulsar.dat output
def find_msp_atbirth_behemoth(modelpath):
    t_conv = ns.conv('t', modelpath+'26859-bigger.conv.sh')
    
    snaps_raw=glob(modelpath+'*.snap*.dat.gz')
    snapno = []
    for xx in range(len(snaps_raw)):
        strings = snaps_raw[xx].split('.')
        no = re.findall(r'\d+', strings[1])[0]
        if no != '0000':
            no = no.lstrip('0')
        else:
            no = '0'
        snapno.append(int(no))
    snapnos, snaps = (list(t) for t in zip(*sorted(zip(snapno, snaps_raw))))

    id0 = []; id1 = []; ps = []; mfield = []; tmyr = []; formation = [] 
    m0 = []; m1 = []; sma = []; ecc = []; k0 = []; k1 = []

    for ii in range(len(snaps)):
        time = ns.get_time(snaps[ii])*t_conv
        with gzip.open(snaps[ii], 'r') as fsnap:
            next(fsnap)
            next(fsnap)
            for line in fsnap:
                datasnap = line.split()
                if int(datasnap[7])!=1:
                    spin=twopi*yearsc/float(datasnap[59])
                    bfield = float(datasnap[60])
                    deathcut=(spin**2)*(0.17*10**12)
                    if int(datasnap[14])==13 and deathcut<bfield and spin<=0.03 and int(datasnap[0]) not in id0:
                        id0.append(int(datasnap[0])); id1.append(-100)
                        m0.append(float(datasnap[1])); m1.append(-100)
                        k0.append(int(datasnap[14])); k1.append(-100)
                        ps.append(spin); mfield.append(bfield)
                        sma.append(-100); ecc.append(-100)
                        tmyr.append(time); formation.append(int(datasnap[61]))
                else:
                    if int(datasnap[17])==13:
                        spin=twopi*yearsc/float(datasnap[45])
                        bfield = float(datasnap[47])
                        deathcut=(spin**2)*(0.17*10**12)
                        if deathcut<bfield and spin<=0.03 and int(datasnap[10]) not in id0:
                            id0.append(int(datasnap[10])); id1.append(int(datasnap[11]))
                            m0.append(float(datasnap[8])); m1.append(float(datasnap[9]))
                            k0.append(int(datasnap[17])); k1.append(int(datasnap[18]))
                            ps.append(spin); mfield.append(bfield)
                            sma.append(float(datasnap[12])); ecc.append(float(datasnap[13]))
                            tmyr.append(time); formation.append(int(datasnap[49]))

                    if int(datasnap[18])==13:
                        spin=twopi*yearsc/float(datasnap[46])
                        bfield = float(datasnap[48])
                        deathcut=(spin**2)*(0.17*10**12)
                        if deathcut<bfield and spin<=0.03 and int(datasnap[11]) not in id0:
                            id0.append(int(datasnap[11])); id1.append(int(datasnap[10]))
                            m0.append(float(datasnap[9])); m1.append(float(datasnap[8]))
                            k0.append(int(datasnap[18])); k1.append(int(datasnap[17]))
                            ps.append(spin); mfield.append(bfield)
                            sma.append(float(datasnap[12])); ecc.append(float(datasnap[13]))            
                            tmyr.append(time); formation.append(int(datasnap[50]))

        print(snaps[ii])

    np.savetxt('/jet/home/yec/behemoth/msp_behemoth_atbirth.dat', np.c_[tmyr, id0, id1, m0, m1, k0, k1, mfield, ps, sma, ecc, formation], fmt = '%f %d %d %f %f %d %d %e %f %f %f %d', header = '#1.Time(Myr) 2.ID0 3.ID1 4.M0 5.M1 6.K0 7.K1 8.B(G) 9.Ps(sec) 10.SMA(AU) 11.ECC 12.Formation')



##Find the mass density at half-mass radius
def find_rh_density(modelpath, thesnap):
    filestr = modelpath+'initial'
    snap = modelpath+thesnap
    dynfile = np.genfromtxt(filestr+'.dyn.dat')
    tdyn = dynfile[:,0]; r_h = dynfile[:,20]

    m = []; radius = []

    t_conv = dyn.conv('t', filestr+'.conv.sh')
    l_conv = dyn.conv('l', filestr+'.conv.sh')
    
    print(l_conv)

    tsnap = dyn.get_time(snap)
    for xx in range(len(tdyn)):
        if round(tdyn[xx], 6) == tsnap:
            rh_snap = r_h[xx]

    print(rh_snap)

    with gzip.open(snap, 'r') as fsnap:
        next(fsnap); next(fsnap)
        for line in fsnap:
            data = line.split()
            m.append(float(data[1])); radius.append(float(data[2])*l_conv)
            if float(data[2])>rh_snap:
                break

    #vol = (2*twopi/3)*(rh_snap*l_conv)**3
    vol = (2*twopi/3)*(radius[-1])**3
    averho_rh = sum(m)/vol


    return averho_rh, m, radius



##Estimate the final num of MSP and L_gamma within 3 kpc of the Galactic Center
def estimate_Nmsp_Lgamma_atpresent():
    sample_init = np.genfromtxt('/projects/b1095/syr904/projects/GCE/catalog/cluster_sample_initial.dat', dtype = str)
    N_star = sample_init[:,0]; rv = sample_init[:,1]; rg = sample_init[:,2]; z = sample_init[:,3]
    rg_init = sample_init[:,5].astype(float)

    num_msp = []; L_gamma = []; Mass = []
    for ii in range(len(N_star)):
        if rg_init[ii]<=3.:
            modelpath = '/projects/b1091/CMC_Grid_March2019/rundir/rv'+rv[ii]+'/rg'+rg[ii]+'/z'+z[ii]+'/'+N_star[ii]+'/'

            m_conv = dyn.conv('m', modelpath+'initial.conv.sh')

            nsfile = modelpath+'initial.ns.dat'
            with open(nsfile, 'r') as fns:
                for line in fns:
                    pass
                ns_last_line = line
                datans = ns_last_line.split()

            num_msp.append(int(datans[6]))

            dynfile = modelpath+'initial.dyn.dat'
            with open(dynfile, 'r') as fdyn:
                for line in fdyn:
                    pass
                dyn_last_line = line
                datadyn = dyn_last_line.split()

            Mass.append(float(datadyn[4])*m_conv)

            lgammafile = '/projects/b1095/syr904/projects/GCE/catalog/Lgamma_alltime_maingrid.dat'
            pathfile = '/projects/b1095/syr904/projects/GCE/catalog/path_allfinished_newruns_maingrid.dat'
            datalg = np.genfromtxt(lgammafile)
            modellg = datalg[:,0]; tlg = datalg[:,1]; lmsp = datalg[:,2]

            sourcedir = np.genfromtxt(pathfile, dtype=str)
            paths = list(sourcedir[:,0])
            model_no = paths.index(modelpath)
            print(model_no, modelpath)
            
            lg_last = lmsp[modellg==model_no]
            if not len(lg_last): continue

            lg_last = lg_last[-1]
            print(lg_last, model_no)
            L_gamma.append(lg_last)


    print(sum(num_msp), sum(L_gamma), sum(Mass))


##Printout ejected MSPs
def print_esc_Nns(pathlist, start, end):
    sourcedir=np.genfromtxt(pathlist, dtype='str')
    status=sourcedir[:,1]
    sourcedir=sourcedir[:,0]

    #sourcedir=['/projects/b1091/CMC_Grid_March2019/rundir/rv0.5/rg8/z0.002/1.6e6/']
    
    for ii in range(start, end):
        pref='initial'
        filestr=sourcedir[ii]+pref

        escfile = filestr+'.esc.dat'
        t_conv = dyn.conv('t', filestr+'.conv.sh')
        os.system('rm '+filestr+'.esc_ns.dat')
    
        try:
           fh = open(filestr+'.esc_ns.dat', 'r')
           print('is done')
        except:
        #if True:
            print(sourcedir[ii])
            #fhandle=open('/projects/b1095/syr904/projects/WDs/initial.ns'+str(i)+'.dat', 'a+')
            fhandle=open(filestr+'.esc_ns.dat', 'a+')
            fhandle.write('#1.Totaltime, 2.Nns,tot, 3.Nns,single, 4.Nns,binary, 5.Nns,mtb, 6.Npulsar, 7.Nmsp, 8.Nns-ns, 9.Nns-bh, 10.Nns-wd, 11.Nns-ms, 12.Nns-postms\n')
            
            N_NS=0; N_NS_SIN=0; N_NS_BIN=0; N_NS_MTB=0; N_PULS=0; N_MSP=0; N_NSNS=0; N_NSBH=0; N_NSWD=0; N_NSMS=0; N_NSPOSTMS=0
            
            thold = 10000
            with open(escfile, 'r') as fesc:
                next(fesc)
                for line in fesc:
                    dataesc=line.split()
                    t = dataesc[1]
                    if float(t) > float(thold):
                        fhandle.write('%f %d %d %d %d %d %d %d %d %d %d %d\n'%(float(thold)*t_conv, N_NS, N_NS_SIN, N_NS_BIN, N_NS_MTB, N_PULS, N_MSP, N_NSNS, N_NSBH, N_NSWD, N_NSMS, N_NSPOSTMS))
                        N_NS=0; N_NS_SIN=0; N_NS_BIN=0; N_NS_MTB=0; N_PULS=0; N_MSP=0; N_NSNS=0; N_NSBH=0; N_NSWD=0; N_NSMS=0; N_NSPOSTMS=0

                    if int(dataesc[14])!=1:
                        if int(dataesc[21])==13: 
                            N_NS+=1; N_NS_SIN+=1
                            spin=twopi*yearsc/float(dataesc[60])
                            deathcut=(spin**2)*(0.17*10**12)
                            if deathcut<float(dataesc[61]): 
                                N_PULS+=1
                                if spin<=0.03: 
                                    N_MSP+=1

                    if int(dataesc[14])==1:
                        if int(dataesc[22])==13:
                            N_NS+=1; N_NS_BIN+=1
                            spin0=twopi*yearsc/float(dataesc[43])
                            deathcut0=(spin0**2)*(0.17*10**12)
                            if float(dataesc[42])>=1: N_NS_MTB+=1
                            if deathcut0<float(dataesc[45]): 
                                N_PULS+=1
                                if spin0<=0.03: 
                                    N_MSP+=1

                            if deathcut0<float(dataesc[45]) and float(dataesc[42])>=1:
                                N_NS_MTB-=1
                                    
                            if int(dataesc[23])<2: N_NSMS+=1
                            elif int(dataesc[23])>=10 and int(dataesc[23])<=12: N_NSWD+=1
                            elif int(dataesc[23])==13: N_NSNS+=1
                            elif int(dataesc[23])==14: N_NSBH+=1
                            else: N_NSPOSTMS+=1

                        if int(dataesc[23])==13:
                            N_NS+=1; N_NS_BIN+=1
                            spin1=twopi*yearsc/float(dataesc[44])
                            deathcut1=(spin1**2)*(0.17*10**12)
                            if float(dataesc[41])>=1: N_NS_MTB+=1
                            if deathcut1<float(dataesc[46]): 
                                N_PULS+=1
                                if spin1<=0.03: 
                                    N_MSP+=1

                            if deathcut1<float(dataesc[46]) and float(dataesc[41])>=1:
                                N_NS_MTB-=1

                            if int(dataesc[22])<2: N_NSMS+=1
                            elif int(dataesc[22])>=10 and int(dataesc[22])<=12: N_NSWD+=1
                            elif int(dataesc[22])==13: print('already counted')
                            elif int(dataesc[22])==14: N_NSBH+=1
                            else: N_NSPOSTMS+=1

                    thold = t

                fhandle.write('%f %d %d %d %d %d %d %d %d %d %d %d\n'%(float(thold)*t_conv, N_NS, N_NS_SIN, N_NS_BIN, N_NS_MTB, N_PULS, N_MSP, N_NSNS, N_NSBH, N_NSWD, N_NSMS, N_NSPOSTMS))

            print(ii)

            fhandle.close()


def mean_confidence_interval(data, confidence=0.95):
    a = 1.0 * np.array(data)
    n = len(a)
    m, se = np.mean(a), scipy.stats.sem(a)
    h = se * scipy.stats.t.ppf((1 + confidence) / 2., n-1)
    return m, m-h, m+h


def rm_max_min_list(thelist):
    thelist  = np.array(thelist)
    max_value = max(thelist)
    min_value = min(thelist)

    #print(max_value, min_value)
    if max_value == min_value:
        return thelist

    else:
        newlist = thelist[thelist<max_value]
        #newnewlist = newlist[newlist>min_value]

        #return newnewlist
        return newlist

##Extract and grouping the number of MSPs from the catalog models
def extract_nmsp_nondissolved():
    pathlist = np.genfromtxt('/projects/b1095/syr904/projects/GCE/catalog/path_allfinished_newruns_maingrid.dat', dtype=str)
    paths = pathlist[:,0]; status = pathlist[:,1]

    bin_size = 400

    Nstar = [200000.,400000.,800000.,1600000.]
    RV = [4.,2.,1.,0.5]
    Zmetal = [0.0002,0.002,0.02]
    RG = [2.,8.,20.]

    n_model_mass = [0,0,0,0]; n_model_rv = [0,0,0,0]; n_model_z = [0,0,0]; n_model_rg = [0,0,0]
    n_model_mass_rv = [[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]
    for ii in range(len(paths)):

        ##Initial Conditions
        s=paths[ii].split('/')
        n_star=float(s[-2])
        z=float(s[-3][1:])
        rg=int(s[-4][2:])
        rv=float(s[-5][2:])
        
        for xx in range(4):
            if n_star==Nstar[xx] and status[ii]=='1':
                n_model_mass[xx]+=1
                for yy in range(4):
                    if rv==RV[yy]:
                        n_model_mass_rv[xx][yy]+=1

            if rv==RV[xx] and status[ii]=='1':
                n_model_rv[xx]+=1
                
            
        if z==0.0002 and status[ii]=='1': 
            n_model_z[0]+=1
        if z==0.002 and status[ii]=='1': 
            n_model_z[1]+=1
        if z==0.02 and status[ii]=='1': 
            n_model_z[2]+=1
            
            
        if rg==2 and status[ii]=='1': 
            n_model_rg[0]+=1
        if rg==8 and status[ii]=='1': 
            n_model_rg[1]+=1
        if rg==20 and status[ii]=='1': 
            n_model_rg[2]+=1
        

    ##Grouping models        
    n_msp_rv = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_msp_rv_average = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_msp_rv_average_std = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_msp_rv_median = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_msp_rv_range = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_msp_rv_mave = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_msp_rv_mave_std = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_msp_rv_mmed = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]


    n_msp_mass = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_msp_mass_average = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_msp_mass_average_std = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_msp_mass_median = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_msp_mass_range = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_msp_mass_mave = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_msp_mass_mave_std = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_msp_mass_mmed = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]

    n_msp_mass_rv = [[np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)],[np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)],[np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)],[np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]]
    n_msp_mass_rv_average = [[np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)],[np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)],[np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)],[np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]]
    n_msp_mass_rv_average_std = [[np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)],[np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)],[np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)],[np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]]
    n_msp_mass_rv_median = [[np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)],[np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)],[np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)],[np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]]

    n_msp_z = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_msp_z_average = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_msp_z_average_std = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_msp_z_median = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_msp_z_range = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_msp_z_mave = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_msp_z_mave_std = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_msp_z_mmed = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]

    n_msp_rg = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_msp_rg_average = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_msp_rg_average_std = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_msp_rg_median = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_msp_rg_range = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_msp_rg_mave = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_msp_rg_mave_std = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_msp_rg_mmed = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    
    nmsp_scatter_n = [[[] for _ in range(bin_size)],[[] for _ in range(bin_size)],[[] for _ in range(bin_size)],[[] for _ in range(bin_size)]]

    nmsp_scatter_rv = [[[] for _ in range(bin_size)],[[] for _ in range(bin_size)],[[] for _ in range(bin_size)],[[] for _ in range(bin_size)]]

    nmsp_scatter_z = [[[] for _ in range(bin_size)],[[] for _ in range(bin_size)],[[] for _ in range(bin_size)]]

    nmsp_scatter_rg = [[[] for _ in range(bin_size)],[[] for _ in range(bin_size)],[[] for _ in range(bin_size)]]

    nmsp_scatter_n_mave = [[[] for _ in range(bin_size)],[[] for _ in range(bin_size)],[[] for _ in range(bin_size)],[[] for _ in range(bin_size)]]

    nmsp_scatter_rv_mave = [[[] for _ in range(bin_size)],[[] for _ in range(bin_size)],[[] for _ in range(bin_size)],[[] for _ in range(bin_size)]]

    nmsp_scatter_z_mave = [[[] for _ in range(bin_size)],[[] for _ in range(bin_size)],[[] for _ in range(bin_size)]]

    nmsp_scatter_rg_mave = [[[] for _ in range(bin_size)],[[] for _ in range(bin_size)],[[] for _ in range(bin_size)]]

    nmsp_scatter_mass_rv = [[[[] for _ in range(bin_size)],[[] for _ in range(bin_size)],[[] for _ in range(bin_size)],[[] for _ in range(bin_size)]], [[[] for _ in range(bin_size)],[[] for _ in range(bin_size)],[[] for _ in range(bin_size)],[[] for _ in range(bin_size)]],[[[] for _ in range(bin_size)],[[] for _ in range(bin_size)],[[] for _ in range(bin_size)],[[] for _ in range(bin_size)]],[[[] for _ in range(bin_size)],[[] for _ in range(bin_size)],[[] for _ in range(bin_size)],[[] for _ in range(bin_size)]]]

    t_all = np.linspace(0, 13000., bin_size+1)
    for kk in range(0, 144):
        print(paths[kk])

        ##Initial Conditions
        s=paths[kk].split('/')
        n_star=float(s[-2])
        z=float(s[-3][1:])
        rg=int(s[-4][2:])
        rv=float(s[-5][2:])
        
        t_conv = dyn.conv('t', paths[kk]+'initial.conv.sh')
        m_conv = dyn.conv('m', paths[kk]+'initial.conv.sh')

        if status[kk]=='1':
            datans = np.genfromtxt(paths[kk]+'initial.ns.dat')
            times = datans[:,0]*t_conv
            n_psr = datans[:,5]; n_msp = datans[:,6]
         
            ##Interpolate the number of NS data
            f = interp1d(times, n_msp, kind='nearest')
            t_interpld = np.linspace(0, np.max(times), 3*bin_size)
            n_msp_new = f(t_interpld)
            #print(n_msp_new)

            t_dyn = []; m_dyn = []
            with open(paths[kk]+'initial.dyn.dat', 'r') as fdyn:
                next(fdyn); next(fdyn)
                for line in fdyn:
                    datadyn = line.split()
                    t_dyn.append(t_conv*float(datadyn[0]))
                    m_dyn.append(m_conv*float(datadyn[4]))

            if len(np.sort(glob(paths[kk]+'*.dyn.dat')))>1:
                with open(paths[kk]+'initial2.dyn.dat', 'r') as fdyn:
                    next(fdyn); next(fdyn)
                    for line in fdyn:
                        datadyn = line.split()
                        t_dyn.append(t_conv*float(datadyn[0]))
                        m_dyn.append(m_conv*float(datadyn[4]))

            #f_dyn = interp1d(t_dyn, m_dyn, kind='nearest')
            #tdyn_interpld = np.linspace(0, np.max(t_dyn), 2*len(t_dyn))
            #m_dyn_interpld = f_dyn(tdyn_interpld)
            tdyn_index = np.digitize(t_dyn, t_all)
            df = pd.DataFrame(m_dyn)
            dyn_mass = (df.groupby(tdyn_index).mean()).values
            dyn_mass = np.delete(dyn_mass,-1)
            print(len(dyn_mass))


    
        n_mass = [[],[],[],[]]; n_rv = [[],[],[],[]]; n_z = [[],[],[]]; n_rg = [[],[],[]]
        n_mass_rv = [[[],[],[],[]],[[],[],[],[]],[[],[],[],[]],[[],[],[],[]]]

        for jj in range(len(t_all)-1):
            #dyn_mass_temp = 0
            #count_dyn = 0
            #if status[kk]=='1':
            #    for i in range(len(t_dyn)):
            #        if t_all[jj] <= t_dyn[i] < t_all[jj+1]:
            #            dyn_mass_temp+=m_dyn[i]
            #            count_dyn+=1

            #    if count_dyn!=0:
            #        dyn_mass.append(dyn_mass_temp/count_dyn)
            #    else:
            #        print('no mass')


            #print(jj)
            n_mass_temp = [0,0,0,0]; n_rv_temp = [0,0,0,0]; n_z_temp = [0,0,0]; n_rg_temp = [0,0,0]
            count_mass = [0,0,0,0]; count_rv = [0,0,0,0]; count_z = [0,0,0]; count_rg = [0,0,0]
            n_mass_rv_temp = [[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]
            count_mass_rv = [[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]
            
            ##Group by initial mass and initial rv
            for xx in range(4):
                if n_star==Nstar[xx] and status[kk]=='1':
                    for i in range(len(t_interpld)):
                        if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                            n_mass_temp[xx]+=n_msp_new[i]
                            count_mass[xx]+=1  ##multiple time steps may belong to the same bin

                    for yy in range(4):
                        if rv==RV[yy]:
                            for i in range(len(t_interpld)):
                                if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                                    n_mass_rv_temp[xx][yy]+=n_msp_new[i]
                                    count_mass_rv[xx][yy]+=1

                if rv==RV[xx] and status[kk]=='1':
                    for i in range(len(t_interpld)):
                        if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                            n_rv_temp[xx]+=n_msp_new[i]
                            count_rv[xx]+=1

        
            ##Group by metallicity
            if z==0.0002 and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_z_temp[0]+=n_msp_new[i]
                        count_z[0]+=1
        
            if z==0.002 and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_z_temp[1]+=n_msp_new[i]
                        count_z[1]+=1
        
            if z==0.02 and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_z_temp[2]+=n_msp_new[i]
                        count_z[2]+=1
                    
                    
            ##Group by galactocentric distance
            if rg==2 and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_rg_temp[0]+=n_msp_new[i]
                        count_rg[0]+=1
        
            if rg==8 and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_rg_temp[1]+=n_msp_new[i]
                        count_rg[1]+=1
        
            if rg==20 and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_rg_temp[2]+=n_msp_new[i]
                        count_rg[2]+=1
    
            #print(count_rv[0])
        
            for x in range(4):
                if count_rv[x]!=0:
                    n_rv_temp[x] = n_rv_temp[x]/count_rv[x]
                if count_mass[x]!=0:
                    n_mass_temp[x] = n_mass_temp[x]/count_mass[x]
                    
                n_rv[x].append(n_rv_temp[x])
                n_mass[x].append(n_mass_temp[x])

                for y in range(4):
                    if count_mass_rv[x][y]!=0:
                        n_mass_rv_temp[x][y] = n_mass_rv_temp[x][y]/count_mass_rv[x][y]

                    n_mass_rv[x][y].append(n_mass_rv_temp[x][y])
        
            for x in range(3):
                if count_z[x]!=0:
                    n_z_temp[x] = n_z_temp[x]/count_z[x]
                
                n_z[x].append(n_z_temp[x])
                
                
                if count_rg[x]!=0:
                    n_rg_temp[x] = n_rg_temp[x]/count_rg[x]
                
                n_rg[x].append(n_rg_temp[x])
                            
            
        for y in range(4):
            n_msp_rv[y] = n_msp_rv[y]+np.array(n_rv[y])
            n_msp_mass[y] = n_msp_mass[y]+np.array(n_mass[y])
            n_msp_rv_average[y] = n_msp_rv_average[y] + np.array(n_rv[y])/n_model_rv[y]
            n_msp_mass_average[y] = n_msp_mass_average[y] + np.array(n_mass[y])/n_model_mass[y]

            
        for y in range(3):
            n_msp_z[y] = n_msp_z[y]+np.array(n_z[y])
            n_msp_z_average[y] = n_msp_z_average[y] + np.array(n_z[y])/n_model_z[y]
            n_msp_rg[y] = n_msp_rg[y]+np.array(n_rg[y])
            n_msp_rg_average[y] = n_msp_rg_average[y] + np.array(n_rg[y])/n_model_rg[y]



        for xx in range(4):
            for yy in range(4):
                if n_star==Nstar[xx] and status[kk]=='1':
                    if rv==RV[yy]:
                        nmsp_scatter_mass_rv[xx][yy] = np.hstack((nmsp_scatter_mass_rv[xx][yy], np.split(np.array(n_mass_rv[xx][yy]),len(n_mass_rv[xx][yy]))))

            ##Group by initial mass
            if n_star==Nstar[xx] and status[kk]=='1':
                #print(len(n_mass[0]))
                nmsp_scatter_n[xx] = np.hstack((nmsp_scatter_n[xx], np.split(np.array(n_mass[xx]),len(n_mass[xx]))))
                #print(nmsp_scatter_n2e5)
                nmsp_scatter_n_mave[xx] = np.hstack((nmsp_scatter_n_mave[xx], np.split(np.array(n_mass[xx])/np.array(dyn_mass),len(n_mass[xx]))))
        

            ##Group by initial rv   
            if rv==RV[xx] and status[kk]=='1':
                nmsp_scatter_rv[xx] = np.hstack((nmsp_scatter_rv[xx], np.split(np.array(n_rv[xx]), len(n_rv[xx]))))
                nmsp_scatter_rv_mave[xx] = np.hstack((nmsp_scatter_rv_mave[xx], np.split(np.array(n_rv[xx])/np.array(dyn_mass), len(n_rv[xx]))))
        
                    
        for xx in range(3):
            ##Group by metallicity
            if z==Zmetal[xx] and status[kk]=='1':
                nmsp_scatter_z[xx] = np.hstack((nmsp_scatter_z[xx], np.split(np.array(n_z[xx]), len(n_z[xx]))))
                nmsp_scatter_z_mave[xx] = np.hstack((nmsp_scatter_z_mave[xx], np.split(np.array(n_z[xx])/np.array(dyn_mass), len(n_z[xx]))))
                    
            ##Group by galactocentric distance
            if rg==RG[xx] and status[kk]=='1':
                nmsp_scatter_rg[xx] = np.hstack((nmsp_scatter_rg[xx], np.split(np.array(n_rg[xx]), len(n_rg[xx]))))
                nmsp_scatter_rg_mave[xx] = np.hstack((nmsp_scatter_rg_mave[xx], np.split(np.array(n_rg[xx])/np.array(dyn_mass), len(n_rg[xx]))))


    for ii in range(4):
        for xx in range(bin_size):
            #print(nmsp_scatter_n[ii][xx])
            #print(rm_max_min_list(nmsp_scatter_n[ii][xx]))
            nmsp_scatter_n_temp = rm_max_min_list(nmsp_scatter_n[ii][xx])
            nmsp_scatter_rv_temp = rm_max_min_list(nmsp_scatter_rv[ii][xx])
            nmsp_scatter_rv_mave_temp = rm_max_min_list(nmsp_scatter_rv_mave[ii][xx])
            nmsp_scatter_n_mave_temp = rm_max_min_list(nmsp_scatter_n_mave[ii][xx])

            n_msp_mass_average_std[ii][xx]+=np.std(nmsp_scatter_n_temp)
            n_msp_rv_average_std[ii][xx]+=np.std(nmsp_scatter_rv_temp)

            n_msp_mass_median[ii][xx]+=np.median(nmsp_scatter_n_temp)
            n_msp_rv_median[ii][xx]+=np.median(nmsp_scatter_rv_temp)
            n_msp_mass_range[ii][xx]+=np.ptp(nmsp_scatter_n_temp)
            n_msp_rv_range[ii][xx]+=np.ptp(nmsp_scatter_rv_temp)

            n_msp_rv_mave[ii][xx]+= np.mean(nmsp_scatter_rv_mave_temp)
            n_msp_rv_mave_std[ii][xx]+= np.std(nmsp_scatter_rv_mave_temp)
            n_msp_rv_mmed[ii][xx]+= np.median(nmsp_scatter_rv_mave_temp)

            n_msp_mass_mave[ii][xx]+= np.mean(nmsp_scatter_n_mave_temp)
            n_msp_mass_mave_std[ii][xx]+= np.std(nmsp_scatter_n_mave_temp)
            n_msp_mass_mmed[ii][xx]+= np.median(nmsp_scatter_n_mave_temp)

        for jj in range(4):
            for yy in range(bin_size):
                n_msp_mass_rv[ii][jj][yy]+=np.sum(nmsp_scatter_mass_rv[ii][jj][yy])
                n_msp_mass_rv_average[ii][jj][yy]+=np.mean(nmsp_scatter_mass_rv[ii][jj][yy])
                n_msp_mass_rv_average_std[ii][jj][yy]+=np.std(nmsp_scatter_mass_rv[ii][jj][yy])
                n_msp_mass_rv_median[ii][jj][yy]+=np.median(nmsp_scatter_mass_rv[ii][jj][yy])

    for ii in range(3):
        for xx in range(bin_size):
            nmsp_scatter_z_temp = rm_max_min_list(nmsp_scatter_z[ii][xx])
            nmsp_scatter_rg_temp = rm_max_min_list(nmsp_scatter_rg[ii][xx])
            nmsp_scatter_z_mave_temp = rm_max_min_list(nmsp_scatter_z_mave[ii][xx])
            nmsp_scatter_rg_mave_temp = rm_max_min_list(nmsp_scatter_rg_mave[ii][xx])

            n_msp_z_average_std[ii][xx]+=np.std(nmsp_scatter_z_temp)
            n_msp_rg_average_std[ii][xx]+=np.std(nmsp_scatter_rg_temp)

            n_msp_z_median[ii][xx]+=np.median(nmsp_scatter_z_temp)
            n_msp_rg_median[ii][xx]+=np.median(nmsp_scatter_rg_temp)
            n_msp_z_range[ii][xx]+=np.ptp(nmsp_scatter_z_temp)
            n_msp_rg_range[ii][xx]+=np.ptp(nmsp_scatter_rg_temp)


            n_msp_z_mave[ii][xx]+=np.mean(nmsp_scatter_z_mave_temp)
            n_msp_z_mave_std[ii][xx]+=np.std(nmsp_scatter_z_mave_temp)
            n_msp_z_mmed[ii][xx]+=np.median(nmsp_scatter_z_mave_temp)

            n_msp_rg_mave[ii][xx]+=np.mean(nmsp_scatter_rg_mave_temp)
            n_msp_rg_mave_std[ii][xx]+=np.std(nmsp_scatter_rg_mave_temp)
            n_msp_rg_mmed[ii][xx]+=np.median(nmsp_scatter_rg_mave_temp)

    
    print(n_msp_mass_average_std[0])

    t_all = np.delete(t_all , 0)

    #for z in range(4):
    #    n_msp_rv[z] = np.insert(n_msp_rv[z], 0, 0.); n_msp_rv_average[z] = np.insert(n_msp_rv_average[z], 0, 0.)
    #    n_msp_rv_average_std[z] = np.insert(n_msp_rv_average_std[z], 0, 0.)
    #    n_msp_rv_median[z] = np.insert(n_msp_rv_median[z], 0, 0.)
#
    #    n_msp_mass[z] = np.insert(n_msp_mass[z], 0, 0.); n_msp_mass_average[z] = np.insert(n_msp_mass_average[z], 0, 0.)
    #    n_msp_mass_average_std[z] = np.insert(n_msp_mass_average_std[z], 0, 0.)
    #    n_msp_mass_median[z] = np.insert(n_msp_mass_median[z], 0, 0.)
#
    #    for z1 in range(4):
    #        n_msp_mass_rv[z][z1] = np.insert(n_msp_mass_rv[z][z1], 0, 0.)
    #        n_msp_mass_rv_average[z][z1] = np.insert(n_msp_mass_rv_average[z][z1], 0, 0.)
    #        n_msp_mass_rv_average_std[z][z1] = np.insert(n_msp_mass_rv_average_std[z][z1], 0, 0.)
    #        n_msp_mass_rv_median[z][z1] = np.insert(n_msp_mass_rv_median[z][z1], 0, 0.)


    #for z in range(3):
    #    n_msp_z[z] = np.insert(n_msp_z[z], 0, 0.); n_msp_z_average[z] = np.insert(n_msp_z_average[z], 0, 0.)
    #    n_msp_z_average_std[z] = np.insert(n_msp_z_average_std[z], 0, 0.)
    #    n_msp_z_median[z] = np.insert(n_msp_z_median[z], 0, 0.)
#
#
    #    n_msp_rg[z] = np.insert(n_msp_rg[z], 0, 0.); n_msp_rg_average[z] = np.insert(n_msp_rg_average[z], 0, 0.)
    #    n_msp_rg_average_std[z] = np.insert(n_msp_rg_average_std[z], 0, 0.)
    #    n_msp_rg_median[z] = np.insert(n_msp_rg_median[z], 0, 0.)

    
    filenames = ['nmsp_mass_age_nondissolved_rmmax.dat', 'nmsp_rv_age_nondissolved_rmmax.dat', 'nmsp_z_age_nondissolved_rmmax.dat', 'nmsp_rg_age_nondissolved_rmmax.dat']
    np.savetxt('/projects/b1095/syr904/projects/GCE/catalog/'+filenames[0], np.c_[t_all, n_msp_mass[0], n_msp_mass[1], n_msp_mass[2], n_msp_mass[3], n_msp_mass_average[0], n_msp_mass_average[1], n_msp_mass_average[2], n_msp_mass_average[3], n_msp_mass_average_std[0], n_msp_mass_average_std[1], n_msp_mass_average_std[2], n_msp_mass_average_std[3], n_msp_mass_median[0], n_msp_mass_median[1], n_msp_mass_median[2], n_msp_mass_median[3], n_msp_mass_range[0], n_msp_mass_range[1], n_msp_mass_range[2], n_msp_mass_range[3],n_msp_mass_mave[0], n_msp_mass_mave[1], n_msp_mass_mave[2], n_msp_mass_mave[3], n_msp_mass_mave_std[0], n_msp_mass_mave_std[1], n_msp_mass_mave_std[2], n_msp_mass_mave_std[3], n_msp_mass_mmed[0], n_msp_mass_mmed[1], n_msp_mass_mmed[2], n_msp_mass_mmed[3]], fmt = '%f %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g', header = '1.time(Myr) 2.N_2e5 3.N_4e5 4.N_8e5 5.N_16e5 6.N_2e5_ave 7.N_4e5_ave 8.N_8e5_ave 9.N_16e5_ave 10.N_2e5_ave_std 11.N_4e5_ave_std 12.N_8e5_ave_std 13.N_16e5_ave_std 14.N_2e5_med 15.N_4e5_med 16.N_8e5_med 17.N_16e5_med 18.N_2e5_range 19.N_4e5_range 20.N_8e5_range 21.N_16e5_range 22.N_2e5_mave 23.N_4e5_mave 24.N_8e5_mave 25.N_16e5_mave 26.N_2e5_mave_std 27.N_4e5_mave_std 28.N_8e5_mave_std 29.N_16e5_mave_std 30.N_2e5_mave_med 31.N_4e5_mave_med 32.N_8e5_mave_med 33.N_16e5_mave_med', comments = '#', delimiter = ' ')

    np.savetxt('/projects/b1095/syr904/projects/GCE/catalog/'+filenames[1], np.c_[t_all, n_msp_rv[0], n_msp_rv[1], n_msp_rv[2], n_msp_rv[3], n_msp_rv_average[0], n_msp_rv_average[1], n_msp_rv_average[2], n_msp_rv_average[3], n_msp_rv_average_std[0], n_msp_rv_average_std[1], n_msp_rv_average_std[2], n_msp_rv_average_std[3], n_msp_rv_median[0], n_msp_rv_median[1], n_msp_rv_median[2], n_msp_rv_median[3], n_msp_rv_range[0], n_msp_rv_range[1], n_msp_rv_range[2], n_msp_rv_range[3], n_msp_rv_mave[0], n_msp_rv_mave[1], n_msp_rv_mave[2], n_msp_rv_mave[3], n_msp_rv_mave_std[0], n_msp_rv_mave_std[1], n_msp_rv_mave_std[2], n_msp_rv_mave_std[3], n_msp_rv_mmed[0], n_msp_rv_mmed[1], n_msp_rv_mmed[2], n_msp_rv_mmed[3]], fmt = '%f %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g', header = '1.time(Myr) 2.rv_4 3.rv_2 4.rv_1 5.rv_0.5 6.rv_4_ave 7.rv_2_ave 8.rv_1_ave 9.rv_0.5_ave 10.rv_4_ave_std 11.rv_2_ave_std 12.rv_1_ave_std 13.rv_0.5_ave_std 14.rv_4_med 15.rv_2_med 16.rv_1_med 17.rv_0.5_med 18.rv_4_range 19.rv_2_range 20.rv_1_range 21.rv_0.5_range 22.rv_4_mave 23.rv_2_mave 24.rv_1_mave 25.rv_0.5_mave 26.rv_4_mave_std 27.rv_2_mave_std 28.rv_1_mave_std 29.rv_0.5_mave_std 30.rv_4_mave_med 31.rv_2_mave_med 32.rv_1_mave_med 33.rv_0.5_mave_med', comments = '#', delimiter = ' ')

    np.savetxt('/projects/b1095/syr904/projects/GCE/catalog/'+filenames[2], np.c_[t_all, n_msp_z[0], n_msp_z[1], n_msp_z[2], n_msp_z_average[0], n_msp_z_average[1], n_msp_z_average[2], n_msp_z_average_std[0], n_msp_z_average_std[1], n_msp_z_average_std[2], n_msp_z_median[0], n_msp_z_median[1], n_msp_z_median[2], n_msp_z_range[0], n_msp_z_range[1], n_msp_z_range[2], n_msp_z_mave[0], n_msp_z_mave[1], n_msp_z_mave[2], n_msp_z_mave_std[0], n_msp_z_mave_std[1], n_msp_z_mave_std[2], n_msp_z_mmed[0], n_msp_z_mmed[1], n_msp_z_mmed[2]], fmt = '%f %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g', header = '1.time(Myr) 2.z_0.0002 3.z_0.002 4.z_0.02 5.z_0.0002_ave 6.z_0.002_ave 7.z_0.02_ave 8.z_0.0002_ave_std 9.z_0.002_ave_std 10.z_0.02_ave_std 11.z_0.0002_med 12.z_0.002_med 13.z_0.02_med 14.z_0.0002_range 15.z_0.002_range 16.z_0.02_range 17.z_0.0002_mave 18.z_0.002_mave 19.z_0.02_mave 20.z_0.0002_mave_std 21.z_0.002_mave_std 22.z_0.02_mave_std 23.z_0.0002_mave_med 24.z_0.002_mave_med 25.z_0.02_mave_med', comments = '#', delimiter = ' ')

    np.savetxt('/projects/b1095/syr904/projects/GCE/catalog/'+filenames[3], np.c_[t_all, n_msp_rg[0], n_msp_rg[1], n_msp_rg[2], n_msp_rg_average[0], n_msp_rg_average[1], n_msp_rg_average[2], n_msp_rg_average_std[0], n_msp_rg_average_std[1], n_msp_rg_average_std[2], n_msp_rg_median[0], n_msp_rg_median[1], n_msp_rg_median[2], n_msp_rg_range[0], n_msp_rg_range[1], n_msp_rg_range[2], n_msp_rg_mave[0], n_msp_rg_mave[1], n_msp_rg_mave[2], n_msp_rg_mave_std[0], n_msp_rg_mave_std[1], n_msp_rg_mave_std[2],n_msp_rg_mmed[0], n_msp_rg_mmed[1], n_msp_rg_mmed[2]], fmt = '%f %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g', header = '1.time(Myr) 2.rg_2 3.rg_8 4.rg_20 5.rg_2_ave 6.rg_8_ave 7.rg_20_ave 8.rg_2_ave_std 9.rg_8_ave_std 10.rg_20_ave_std 11.rg_2_med 12.rg_8_med 13.rg_20_med 14.rg_2_range 15.rg_8_range 16.rg_20_range 17.rg_2_mave 18.rg_8_mave 19.rg_20_mave 20.rg_2_mave_std 21.rg_8_mave_std 22.rg_20_mave_std 23.rg_2_mave_med 24.rg_8_mave_med 25.rg_20_mave_med', comments = '#', delimiter = ' ')
    #np.savetxt('/projects/b1095/syr904/projects/GCE/catalog/nmsp_mass_rv_age_nondissolved.dat', np.c_[t_all, n_msp_mass_rv[0][0], n_msp_mass_rv[0][1], n_msp_mass_rv[0][2], n_msp_mass_rv[0][3], n_msp_mass_rv[1][0], n_msp_mass_rv[1][1], n_msp_mass_rv[1][2], n_msp_mass_rv[1][3], n_msp_mass_rv[2][0], n_msp_mass_rv[2][1], n_msp_mass_rv[2][2], n_msp_mass_rv[2][3], n_msp_mass_rv[3][0], n_msp_mass_rv[3][1], n_msp_mass_rv[3][2], n_msp_mass_rv[3][3],n_msp_mass_rv_average[0][0], n_msp_mass_rv_average[0][1], n_msp_mass_rv_average[0][2], n_msp_mass_rv_average[0][3], n_msp_mass_rv_average[1][0], n_msp_mass_rv_average[1][1], n_msp_mass_rv_average[1][2], n_msp_mass_rv_average[1][3], n_msp_mass_rv_average[2][0], n_msp_mass_rv_average[2][1], n_msp_mass_rv_average[2][2], n_msp_mass_rv_average[2][3], n_msp_mass_rv_average[3][0], n_msp_mass_rv_average[3][1], n_msp_mass_rv_average[3][2], n_msp_mass_rv_average[3][3], n_msp_mass_rv_average_std[0][0], n_msp_mass_rv_average_std[0][1], n_msp_mass_rv_average_std[0][2], n_msp_mass_rv_average_std[0][3], n_msp_mass_rv_average_std[1][0], n_msp_mass_rv_average_std[1][1], n_msp_mass_rv_average_std[1][2], n_msp_mass_rv_average_std[1][3], n_msp_mass_rv_average_std[2][0], n_msp_mass_rv_average_std[2][1], n_msp_mass_rv_average_std[2][2], n_msp_mass_rv_average_std[2][3], n_msp_mass_rv_average_std[3][0], n_msp_mass_rv_average_std[3][1], n_msp_mass_rv_average_std[3][2], n_msp_mass_rv_average_std[3][3], n_msp_mass_rv_median[0][0], n_msp_mass_rv_median[0][1], n_msp_mass_rv_median[0][2], n_msp_mass_rv_median[0][3], n_msp_mass_rv_median[1][0], n_msp_mass_rv_median[1][1], n_msp_mass_rv_median[1][2], n_msp_mass_rv_median[1][3], n_msp_mass_rv_median[2][0], n_msp_mass_rv_median[2][1], n_msp_mass_rv_median[2][2], n_msp_mass_rv_median[2][3], n_msp_mass_rv_median[3][0], n_msp_mass_rv_median[3][1], n_msp_mass_rv_median[3][2], n_msp_mass_rv_median[3][3]], fmt = '%f %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g', header = '1.time(Myr) 2.N2e5_rv4 3.N2e5_rv2 4.N2e5_rv1 5.N2e5_rv0.5 6.N4e5_rv4 7.N4e5_rv2 8.N4e5_rv1 9.N4e5_rv0.5 10.N8e5_rv4 11.N8e5_rv2 12.N8e5_rv1 13.N8e5_rv0.5 14.N16e5_rv4 15.N16e5_rv2 16.N16e5_rv1 17.N16e5_rv0.5 18.N2e5_rv4_ave 19.N2e5_rv2_ave 20.N2e5_rv1_ave 21.N2e5_rv0.5_ave 22.N4e5_rv4_ave 23.N4e5_rv2_ave 24.N4e5_rv1_ave 25.N4e5_rv0.5_ave 26.N8e5_rv4_ave 27.N8e5_rv2_ave 28.N8e5_rv1_ave 29.N8e5_rv0.5_ave 30.N16e5_rv4_ave 31.N16e5_rv2_ave 32.N16e5_rv1_ave 33.N16e5_rv0.5_ave 34.N2e5_rv4_ave_std 35.N2e5_rv2_ave_std 36.N2e5_rv1_ave_std 37.N2e5_rv0.5_ave_std 38.N4e5_rv4_ave_std 39.N4e5_rv2_ave_std 40.N4e5_rv1_ave_std 41.N4e5_rv0.5_ave_std 42.N8e5_rv4_ave_std 43.N8e5_rv2_ave_std 44.N8e5_rv1_ave_std 45.N8e5_rv0.5_ave_std 46.N16e5_rv4_ave_std 47.N16e5_rv2_ave_std 48.N16e5_rv1_ave_std 49.N16e5_rv0.5_ave_std 50.N2e5_rv4_med 51.N2e5_rv2_med 52.N2e5_rv1_med 53.N2e5_rv0.5_med 54.N4e5_rv4_med 55.N4e5_rv2_med 56.N4e5_rv1_med 57.N4e5_rv0.5_med 58.N8e5_rv4_med 59.N8e5_rv2_med 60.N8e5_rv1_med 61.N8e5_rv0.5_med 62.N16e5_rv4_med 63.N16e5_rv2_med 64.N16e5_rv1_med 65.N16e5_rv0.5_med', comments = '#', delimiter = ' ')


##Extract and grouping the number of MSPs from the catalog models
def extract_nmsp():
    pathlist = np.genfromtxt('/projects/b1095/syr904/projects/GCE/catalog/path_allfinished_newruns_maingrid.dat', dtype=str)
    paths = pathlist[:,0]; status = pathlist[:,1]

    bin_size = 400

    n_model_mass = [0,0,0,0]; n_model_rv = [0,0,0,0]; n_model_z = [0,0,0]; n_model_rg = [0,0,0]
    for ii in range(len(paths)):

        ##Initial Conditions
        s=paths[ii].split('/')
        n_star=float(s[-2])
        z=float(s[-3][1:])
        rg=int(s[-4][2:])
        rv=float(s[-5][2:])

        if rg>2:
            continue
        
        if n_star==200000.:# and status[ii]=='1': 
            n_model_mass[0]+=1
        if n_star==400000.:# and status[ii]=='1': 
            n_model_mass[1]+=1
        if n_star==800000.:# and status[ii]=='1': 
            n_model_mass[2]+=1
        if n_star==1600000.:# and status[ii]=='1': 
            n_model_mass[3]+=1
            
        if rv==4.:# and status[ii]=='1': 
            n_model_rv[0]+=1
        if rv==2.:# and status[ii]=='1': 
            n_model_rv[1]+=1
        if rv==1.:# and status[ii]=='1': 
            n_model_rv[2]+=1
        if rv==0.5:# and status[ii]=='1': 
            n_model_rv[3]+=1
            
            
        if z==0.0002:# and status[ii]=='1': 
            n_model_z[0]+=1
        if z==0.002:# and status[ii]=='1': 
            n_model_z[1]+=1
        if z==0.02:# and status[ii]=='1': 
            n_model_z[2]+=1
            
            
        if rg==2:# and status[ii]=='1': 
            n_model_rg[0]+=1
        if rg==8:# and status[ii]=='1': 
            n_model_rg[1]+=1
        if rg==20:# and status[ii]=='1': 
            n_model_rg[2]+=1
        

    ##Grouping models        
    n_msp_rv = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_msp_rv_average = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_msp_rv_average_std = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_msp_rv_median = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]

    n_msp_mass = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_msp_mass_average = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_msp_mass_average_std = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_msp_mass_median = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]

    n_msp_z = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_msp_z_average = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_msp_z_average_std = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_msp_z_median = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]


    n_msp_rg = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_msp_rg_average = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_msp_rg_average_std = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_msp_rg_median = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    
    nmsp_scatter_n2e5 = [[] for _ in range(bin_size)]
    nmsp_scatter_n4e5 = [[] for _ in range(bin_size)]
    nmsp_scatter_n8e5 = [[] for _ in range(bin_size)]
    nmsp_scatter_n16e5 = [[] for _ in range(bin_size)]

    nmsp_scatter_rv4 = [[] for _ in range(bin_size)]
    nmsp_scatter_rv2 = [[] for _ in range(bin_size)]
    nmsp_scatter_rv1 = [[] for _ in range(bin_size)]
    nmsp_scatter_rv05 = [[] for _ in range(bin_size)]

    nmsp_scatter_z00002 = [[] for _ in range(bin_size)]
    nmsp_scatter_z0002 = [[] for _ in range(bin_size)]
    nmsp_scatter_z002= [[] for _ in range(bin_size)]

    nmsp_scatter_rg2 = [[] for _ in range(bin_size)]
    nmsp_scatter_rg8 = [[] for _ in range(bin_size)]
    nmsp_scatter_rg20= [[] for _ in range(bin_size)]


    t_all = np.linspace(0, 13000., bin_size+1)
    for kk in range(0, 144):
        print(paths[kk])

        ##Initial Conditions
        s=paths[kk].split('/')
        n_star=float(s[-2])
        z=float(s[-3][1:])
        rg=int(s[-4][2:])
        rv=float(s[-5][2:])

        if rg>2:
            continue
        
        t_conv = dyn.conv('t', paths[kk]+'initial.conv.sh')
        datans = np.genfromtxt(paths[kk]+'initial.ns.dat')
        times = np.array(datans[:,0])*t_conv
        n_psr = datans[:,5]; n_msp = datans[:,6]
        
        ##Interpolate the number of NS data
        f = interp1d(times, n_msp, kind='nearest')
        t_interpld = np.linspace(0, np.max(times), 3*bin_size)
        n_msp_new = f(t_interpld)
        #print(n_msp_new)
    
        n_mass = [[],[],[],[]]; n_rv = [[],[],[],[]]; n_z = [[],[],[]]; n_rg = [[],[],[]]
        for jj in range(len(t_all)-1):
            #print(jj)
            n_mass_temp = [0,0,0,0]; n_rv_temp = [0,0,0,0]; n_z_temp = [0,0,0]; n_rg_temp = [0,0,0]
            count_mass = [0,0,0,0]; count_rv = [0,0,0,0]; count_z = [0,0,0]; count_rg = [0,0,0]
            
            ##Group by initial mass
            if n_star==200000.:# and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_mass_temp[0]+=n_msp_new[i]
                        count_mass[0]+=1  ##multiple time steps may belong to the same bin
        
            if n_star==400000.:# and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_mass_temp[1]+=n_msp_new[i]
                        count_mass[1]+=1
        
            if n_star==800000.:# and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_mass_temp[2]+=n_msp_new[i]
                        count_mass[2]+=1
        
            if n_star==1600000.:# and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_mass_temp[3]+=n_msp_new[i]
                        count_mass[3]+=1
            
            ##Group by initial rv   
            if rv==4.:# and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_rv_temp[0]+=n_msp_new[i]
                        count_rv[0]+=1
        
            if rv==2.:# and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_rv_temp[1]+=n_msp_new[i]
                        count_rv[1]+=1
        
            if rv==1.:# and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_rv_temp[2]+=n_msp_new[i]
                        count_rv[2]+=1
        
            if rv==0.5:# and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_rv_temp[3]+=n_msp_new[i]
                        count_rv[3]+=1
                    
        
            ##Group by metallicity
            if z==0.0002:# and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_z_temp[0]+=n_msp_new[i]
                        count_z[0]+=1
        
            if z==0.002:# and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_z_temp[1]+=n_msp_new[i]
                        count_z[1]+=1
        
            if z==0.02:# and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_z_temp[2]+=n_msp_new[i]
                        count_z[2]+=1
                    
                    
            ##Group by galactocentric distance
            if rg==2:# and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_rg_temp[0]+=n_msp_new[i]
                        count_rg[0]+=1
        
            if rg==8:# and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_rg_temp[1]+=n_msp_new[i]
                        count_rg[1]+=1
        
            if rg==20:# and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_rg_temp[2]+=n_msp_new[i]
                        count_rg[2]+=1
    
            #print(count_rv[0])
        
            for x in range(4):
                if count_rv[x]!=0:
                    n_rv_temp[x] = n_rv_temp[x]/count_rv[x]
                if count_mass[x]!=0:
                    n_mass_temp[x] = n_mass_temp[x]/count_mass[x]
                    
                n_rv[x].append(n_rv_temp[x])
                n_mass[x].append(n_mass_temp[x])
        
            for x in range(3):
                if count_z[x]!=0:
                    n_z_temp[x] = n_z_temp[x]/count_z[x]
                
                n_z[x].append(n_z_temp[x])     
                
                if count_rg[x]!=0:
                    n_rg_temp[x] = n_rg_temp[x]/count_rg[x]
                
                n_rg[x].append(n_rg_temp[x])
                            
            
        for y in range(4):
            n_msp_rv[y] = n_msp_rv[y]+np.array(n_rv[y])
            n_msp_mass[y] = n_msp_mass[y]+np.array(n_mass[y])
            n_msp_rv_average[y] = n_msp_rv_average[y] + np.array(n_rv[y])/n_model_rv[y]
            n_msp_mass_average[y] = n_msp_mass_average[y] + np.array(n_mass[y])/n_model_mass[y]
            
        for y in range(3):
            n_msp_z[y] = n_msp_z[y]+np.array(n_z[y])
            n_msp_z_average[y] = n_msp_z_average[y] + np.array(n_z[y])/n_model_z[y]
            n_msp_rg[y] = n_msp_rg[y]+np.array(n_rg[y])
            n_msp_rg_average[y] = n_msp_rg_average[y] + np.array(n_rg[y])/n_model_rg[y]


        ##Group by initial mass
        if n_star==200000.:# and status[kk]=='1':
            #print(len(n_mass[0]))
            nmsp_scatter_n2e5 = np.hstack((nmsp_scatter_n2e5, np.split(np.array(n_mass[0]),len(n_mass[0]))))
            #print(nmsp_scatter_n2e5)
        
        if n_star==400000.:# and status[kk]=='1':
            nmsp_scatter_n4e5 = np.hstack((nmsp_scatter_n4e5, np.split(np.array(n_mass[1]), len(n_mass[1]))))
        
        if n_star==800000.:# and status[kk]=='1':
            nmsp_scatter_n8e5 = np.hstack((nmsp_scatter_n8e5, np.split(np.array(n_mass[2]), len(n_mass[2]))))

        if n_star==1600000.:# and status[kk]=='1':
            nmsp_scatter_n16e5 = np.hstack((nmsp_scatter_n16e5, np.split(np.array(n_mass[3]), len(n_mass[3]))))
            

        ##Group by initial rv   
        if rv==4.:# and status[kk]=='1':
            nmsp_scatter_rv4 = np.hstack((nmsp_scatter_rv4, np.split(np.array(n_rv[0]), len(n_rv[0]))))
        
        if rv==2.:# and status[kk]=='1':
            nmsp_scatter_rv2 = np.hstack((nmsp_scatter_rv2, np.split(np.array(n_rv[1]), len(n_rv[1]))))
        
        if rv==1.:# and status[kk]=='1':
            nmsp_scatter_rv1 = np.hstack((nmsp_scatter_rv1, np.split(np.array(n_rv[2]), len(n_rv[2]))))
        
        if rv==0.5:# and status[kk]=='1':
            nmsp_scatter_rv05 = np.hstack((nmsp_scatter_rv05, np.split(np.array(n_rv[3]), len(n_rv[3]))))
                    
    
        ##Group by metallicity
        if z==0.0002:# and status[kk]=='1':
            nmsp_scatter_z00002 = np.hstack((nmsp_scatter_z00002, np.split(np.array(n_z[0]), len(n_z[0]))))
        
        if z==0.002:# and status[kk]=='1':
            nmsp_scatter_z0002 = np.hstack((nmsp_scatter_z0002, np.split(np.array(n_z[1]), len(n_z[1]))))
        
        if z==0.02:# and status[kk]=='1':
            nmsp_scatter_z002 = np.hstack((nmsp_scatter_z002, np.split(np.array(n_z[2]), len(n_z[2]))))
                    
                    
        ##Group by galactocentric distance
        if rg==2:# and status[kk]=='1':
            nmsp_scatter_rg2 = np.hstack((nmsp_scatter_rg2, np.split(np.array(n_rg[0]), len(n_rg[0]))))
        
        if rg==8:# and status[kk]=='1':
            nmsp_scatter_rg8 = np.hstack((nmsp_scatter_rg8, np.split(np.array(n_rg[1]), len(n_rg[1]))))
        
        if rg==20:# and status[kk]=='1':
            nmsp_scatter_rg20 = np.hstack((nmsp_scatter_rg20, np.split(np.array(n_rg[2]), len(n_rg[2]))))

    for ii in range(4):
        for xx in range(bin_size):
            if ii == 0:
                n_msp_mass_average_std[ii][xx]+=np.std(nmsp_scatter_n2e5[xx])
                n_msp_rv_average_std[ii][xx]+=np.std(nmsp_scatter_rv4[xx])

                n_msp_mass_median[ii][xx]+=np.median(nmsp_scatter_n2e5[xx])
                n_msp_rv_median[ii][xx]+=np.median(nmsp_scatter_rv4[xx])
            if ii == 1:
                n_msp_mass_average_std[ii][xx]+=np.std(nmsp_scatter_n4e5[xx])
                n_msp_rv_average_std[ii][xx]+=np.std(nmsp_scatter_rv2[xx])

                n_msp_mass_median[ii][xx]+=np.median(nmsp_scatter_n4e5[xx])
                n_msp_rv_median[ii][xx]+=np.median(nmsp_scatter_rv2[xx])
            if ii == 2:
                n_msp_mass_average_std[ii][xx]+=np.std(nmsp_scatter_n8e5[xx])
                n_msp_rv_average_std[ii][xx]+=np.std(nmsp_scatter_rv1[xx])

                n_msp_mass_median[ii][xx]+=np.median(nmsp_scatter_n8e5[xx])
                n_msp_rv_median[ii][xx]+=np.median(nmsp_scatter_rv1[xx])
            if ii == 3:
                n_msp_mass_average_std[ii][xx]+=np.std(nmsp_scatter_n16e5[xx])
                n_msp_rv_average_std[ii][xx]+=np.std(nmsp_scatter_rv05[xx])

                n_msp_mass_median[ii][xx]+=np.median(nmsp_scatter_n16e5[xx])
                n_msp_rv_median[ii][xx]+=np.median(nmsp_scatter_rv05[xx])

    for ii in range(3):
        for xx in range(bin_size):
            if ii == 0:
                n_msp_z_average_std[ii][xx]+=np.std(nmsp_scatter_z00002[xx])
                n_msp_rg_average_std[ii][xx]+=np.std(nmsp_scatter_rg2[xx])

                n_msp_z_median[ii][xx]+=np.median(nmsp_scatter_z00002[xx])
                n_msp_rg_median[ii][xx]+=np.median(nmsp_scatter_rg2[xx])
            if ii == 1:
                n_msp_z_average_std[ii][xx]+=np.std(nmsp_scatter_z0002[xx])
                n_msp_rg_average_std[ii][xx]+=np.std(nmsp_scatter_rg8[xx])

                n_msp_z_median[ii][xx]+=np.median(nmsp_scatter_z0002[xx])
                n_msp_rg_median[ii][xx]+=np.median(nmsp_scatter_rg8[xx])
            if ii == 2:
                n_msp_z_average_std[ii][xx]+=np.std(nmsp_scatter_z002[xx])
                n_msp_rg_average_std[ii][xx]+=np.std(nmsp_scatter_rg20[xx])

                n_msp_z_median[ii][xx]+=np.median(nmsp_scatter_z002[xx])
                n_msp_rg_median[ii][xx]+=np.median(nmsp_scatter_rg20[xx])
    
    print(n_msp_mass_average_std[0])

    for z in range(4):
        n_msp_rv[z] = np.insert(n_msp_rv[z], 0, 0.); n_msp_rv_average[z] = np.insert(n_msp_rv_average[z], 0, 0.)
        n_msp_rv_average_std[z] = np.insert(n_msp_rv_average_std[z], 0, 0.)
        n_msp_rv_median[z] = np.insert(n_msp_rv_median[z], 0, 0.)

        n_msp_mass[z] = np.insert(n_msp_mass[z], 0, 0.); n_msp_mass_average[z] = np.insert(n_msp_mass_average[z], 0, 0.)
        n_msp_mass_average_std[z] = np.insert(n_msp_mass_average_std[z], 0, 0.)
        n_msp_mass_median[z] = np.insert(n_msp_mass_median[z], 0, 0.)

    for z in range(3):
        n_msp_z[z] = np.insert(n_msp_z[z], 0, 0.); n_msp_z_average[z] = np.insert(n_msp_z_average[z], 0, 0.)
        n_msp_z_average_std[z] = np.insert(n_msp_z_average_std[z], 0, 0.)
        n_msp_z_median[z] = np.insert(n_msp_z_median[z], 0, 0.)


        n_msp_rg[z] = np.insert(n_msp_rg[z], 0, 0.); n_msp_rg_average[z] = np.insert(n_msp_rg_average[z], 0, 0.)
        n_msp_rg_average_std[z] = np.insert(n_msp_rg_average_std[z], 0, 0.)
        n_msp_rg_median[z] = np.insert(n_msp_rg_median[z], 0, 0.)

    
    filenames = ['nmsp_mass_age_rg2.dat', 'nmsp_rv_age.dat', 'nmsp_z_age.dat', 'nmsp_rg_age.dat']
    np.savetxt('/projects/b1095/syr904/projects/GCE/catalog/'+filenames[0], np.c_[t_all, n_msp_mass[0], n_msp_mass[1], n_msp_mass[2], n_msp_mass[3], n_msp_mass_average[0], n_msp_mass_average[1], n_msp_mass_average[2], n_msp_mass_average[3], n_msp_mass_average_std[0], n_msp_mass_average_std[1], n_msp_mass_average_std[2], n_msp_mass_average_std[3], n_msp_mass_median[0], n_msp_mass_median[1], n_msp_mass_median[2], n_msp_mass_median[3]], fmt = '%f %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g', header = '1.time(Myr) 2.N_2e5 3.N_4e5 4.N_8e5 5.N_16e5 6.N_2e5_ave 7.N_4e5_ave 8.N_8e5_ave 9.N_16e5_ave 10.N_2e5_ave_std 11.N_4e5_ave_std 12.N_8e5_ave_std 13.N_16e5_ave_std 14.N_2e5_med 15.N_4e5_med 16.N_8e5_med 17.N_16e5_med', comments = '#', delimiter = ' ')
    #np.savetxt('/projects/b1095/syr904/projects/GCE/catalog/'+filenames[1], np.c_[t_all, n_msp_rv[0], n_msp_rv[1], n_msp_rv[2], n_msp_rv[3], n_msp_rv_average[0], n_msp_rv_average[1], n_msp_rv_average[2], n_msp_rv_average[3], n_msp_rv_average_std[0], n_msp_rv_average_std[1], n_msp_rv_average_std[2], n_msp_rv_average_std[3], n_msp_rv_median[0], n_msp_rv_median[1], n_msp_rv_median[2], n_msp_rv_median[3]], fmt = '%f %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g', header = '1.time(Myr) 2.rv_4 3.rv_2 4.rv_1 5.rv_0.5 6.rv_4_ave 7.rv_2_ave 8.rv_1_ave 9.rv_0.5_ave 10.rv_4_ave_std 11.rv_2_ave_std 12.rv_1_ave_std 13.rv_0.5_ave_std 14.rv_4_med 15.rv_2_med 16.rv_1_med 17.rv_0.5_med', comments = '#', delimiter = ' ')
    #np.savetxt('/projects/b1095/syr904/projects/GCE/catalog/'+filenames[2], np.c_[t_all, n_msp_z[0], n_msp_z[1], n_msp_z[2], n_msp_z_average[0], n_msp_z_average[1], n_msp_z_average[2], n_msp_z_average_std[0], n_msp_z_average_std[1], n_msp_z_average_std[2], n_msp_z_median[0], n_msp_z_median[1], n_msp_z_median[2]], fmt = '%f %g %g %g %g %g %g %g %g %g %g %g %g', header = '1.time(Myr) 2.z_0.0002 3.z_0.002 4.z_0.02 5.z_0.0002_ave 6.z_0.002_ave 7.z_0.02_ave 8.z_0.0002_ave_std 9.z_0.002_ave_std 10.z_0.02_ave_std 11.z_0.0002_med 12.z_0.002_med 13.z_0.02_med', comments = '#', delimiter = ' ')
    #np.savetxt('/projects/b1095/syr904/projects/GCE/catalog/'+filenames[3], np.c_[t_all, n_msp_rg[0], n_msp_rg[1], n_msp_rg[2], n_msp_rg_average[0], n_msp_rg_average[1], n_msp_rg_average[2], n_msp_rg_average_std[0], n_msp_rg_average_std[1], n_msp_rg_average_std[2],n_msp_rg_median[0], n_msp_rg_median[1], n_msp_rg_median[2] ], fmt = '%f %g %g %g %g %g %g %g %g %g %g %g %g', header = '1.time(Myr) 2.rg_2 3.rg_8 4.rg_20 5.rg_2_ave 6.rg_8_ave 7.rg_20_ave 8.rg_2_ave_std 9.rg_8_ave_std 10.rg_20_ave_std 11.rg_2_med 12.rg_8_med 13.rg_20_med', comments = '#', delimiter = ' ')



##Extract the number of MSPs from the catalog models
def extract_nmsp_scatter_points():
    pathlist = np.genfromtxt('/projects/b1095/syr904/projects/GCE/catalog/path_allfinished_newruns_maingrid.dat', dtype=str)
    paths = pathlist[:,0]; status = pathlist[:,1]

    ##Grouping models        
    scatter_nmsp_mass = [[],[],[],[],[],[],[],[]]
    scatter_nmsp_rv = [[],[],[],[],[],[],[],[]]
    scatter_nmsp_z = [[],[],[],[],[],[]]
    scatter_nmsp_rg = [[],[],[],[],[],[]]

    for kk in range(0, 144):
        print(paths[kk])

        ##Initial Conditions
        s=paths[kk].split('/')
        n_star=float(s[-2])
        z=float(s[-3][1:])
        rg=int(s[-4][2:])
        rv=float(s[-5][2:])
        
        t_conv = dyn.conv('t', paths[kk]+'initial.conv.sh')
        datans = np.genfromtxt(paths[kk]+'initial.ns.dat')
        times = datans[:,0]*t_conv
        n_psr = datans[:,5]; n_msp = datans[:,6]
        
            
        ##Group by initial mass
        if n_star==200000.:# and status[kk]=='1':
            scatter_nmsp_mass[0]+=list(times)
            scatter_nmsp_mass[1]+=list(n_msp)
        
        if n_star==400000.:# and status[kk]=='1':
            scatter_nmsp_mass[2]+=list(times)
            scatter_nmsp_mass[3]+=list(n_msp)
        
        if n_star==800000.:# and status[kk]=='1':
            scatter_nmsp_mass[4]+=list(times)
            scatter_nmsp_mass[5]+=list(n_msp)

        if n_star==1600000.:# and status[kk]=='1':
            scatter_nmsp_mass[6]+=list(times)
            scatter_nmsp_mass[7]+=list(n_msp)

        ##Group by initial rv   
        if rv==4.:# and status[kk]=='1':
            scatter_nmsp_rv[0]+=list(times)
            scatter_nmsp_rv[1]+=list(n_msp)
        
        if rv==2.:# and status[kk]=='1':
            scatter_nmsp_rv[2]+=list(times)
            scatter_nmsp_rv[3]+=list(n_msp)
        
        if rv==1.:# and status[kk]=='1':
            scatter_nmsp_rv[4]+=list(times)
            scatter_nmsp_rv[5]+=list(n_msp)
        
        if rv==0.5:# and status[kk]=='1':
            scatter_nmsp_rv[6]+=list(times)
            scatter_nmsp_rv[7]+=list(n_msp)
                    
        
        ##Group by metallicity
        if z==0.0002:# and status[kk]=='1':
            scatter_nmsp_z[0]+=list(times)
            scatter_nmsp_z[1]+=list(n_msp)
        
        if z==0.002:# and status[kk]=='1':
            scatter_nmsp_z[2]+=list(times)
            scatter_nmsp_z[3]+=list(n_msp)
        
        if z==0.02:# and status[kk]=='1':
            scatter_nmsp_z[4]+=list(times)
            scatter_nmsp_z[5]+=list(n_msp)
                    
                    
        ##Group by galactocentric distance
        if rg==2:# and status[kk]=='1':
            scatter_nmsp_rg[0]+=list(times)
            scatter_nmsp_rg[1]+=list(n_msp)
        
        if rg==8:# and status[kk]=='1':
            scatter_nmsp_rg[2]+=list(times)
            scatter_nmsp_rg[3]+=list(n_msp)
        
        if rg==20:# and status[kk]=='1':
            scatter_nmsp_rg[4]+=list(times)
            scatter_nmsp_rg[5]+=list(n_msp)


        nlist = ['2e5', '4e5', '8e5', '16e5']
        rvlist = ['4', '2', '1', '0.5']
        rglist = ['2', '8', '20']
        zlist = ['0.0002', '0.002', '0.02']
    for ii in range(4):
        t_mass_sort, msp_mass_sort = (np.array(x) for x in zip(*sorted(zip(scatter_nmsp_mass[ii*2], scatter_nmsp_mass[ii*2+1]))))
        t_rv_sort, msp_rv_sort = (np.array(x) for x in zip(*sorted(zip(scatter_nmsp_rv[ii*2], scatter_nmsp_rv[ii*2+1]))))
        #np.savetxt('/projects/b1095/syr904/projects/GCE/catalog/nmsp_time_scatter_mass_'+nlist[ii]+'.dat', np.c_[t_mass_sort, msp_mass_sort], fmt = '%f %d', comments = '#', header = '1.Time_n'+nlist[ii]+'[Myr] 2.Nmsp_n'+nlist[ii], delimiter = ' ')
        #np.savetxt('/projects/b1095/syr904/projects/GCE/catalog/nmsp_time_scatter_rv_'+rvlist[ii]+'.dat', np.c_[t_rv_sort, msp_rv_sort], fmt = '%f %d', comments = '#', header = '1.Time_rv'+rvlist[ii]+'[Myr] 2.Nmsp_rv'+rvlist[ii], delimiter = ' ')


    for ii in range(3):
        t_z_sort, msp_z_sort = (np.array(x) for x in zip(*sorted(zip(scatter_nmsp_z[ii*2], scatter_nmsp_z[ii*2+1]))))
        t_rg_sort, msp_rg_sort = (np.array(x) for x in zip(*sorted(zip(scatter_nmsp_rg[ii*2], scatter_nmsp_rg[ii*2+1]))))
        #np.savetxt('/projects/b1095/syr904/projects/GCE/catalog/nmsp_time_scatter_z_'+zlist[ii]+'.dat', np.c_[t_z_sort, msp_z_sort], fmt = '%f %d', comments = '#', header = '1.Time_z'+zlist[ii]+'[Myr] 2.Nmsp_z'+zlist[ii], delimiter = ' ')
        np.savetxt('/projects/b1095/syr904/projects/GCE/catalog/nmsp_time_scatter_rg_'+rglist[ii]+'.dat', np.c_[t_rg_sort, msp_rg_sort], fmt = '%f %d', comments = '#', header = '1.Time_rg'+rglist[ii]+'[Myr] 2.Nmsp_rg'+rglist[ii], delimiter = ' ')


##Extract and grouping the number of MSPs from the catalog models
def extract_escaped_nmsp():
    pathlist = np.genfromtxt('/projects/b1095/syr904/projects/GCE/catalog/path_allfinished_newruns_maingrid.dat', dtype=str)
    paths = pathlist[:,0]; status = pathlist[:,1]

    bin_size = 400

    n_model_mass = [0,0,0,0]; n_model_rv = [0,0,0,0]; n_model_z = [0,0,0]; n_model_rg = [0,0,0]
    for ii in range(len(paths)):

        ##Initial Conditions
        s=paths[ii].split('/')
        n_star=float(s[-2])
        z=float(s[-3][1:])
        rg=int(s[-4][2:])
        rv=float(s[-5][2:])

        if n_star==200000.:# and status[ii]=='1': 
            n_model_mass[0]+=1.
        if n_star==400000.:# and status[ii]=='1': 
            n_model_mass[1]+=1.
        if n_star==800000.:# and status[ii]=='1': 
            n_model_mass[2]+=1.
        if n_star==1600000.:# and status[ii]=='1': 
            n_model_mass[3]+=1.
            
        if rv==4.:# and status[ii]=='1': 
            n_model_rv[0]+=1.
        if rv==2.:# and status[ii]=='1': 
            n_model_rv[1]+=1.
        if rv==1.:# and status[ii]=='1': 
            n_model_rv[2]+=1.
        if rv==0.5:# and status[ii]=='1': 
            n_model_rv[3]+=1.
            
            
        if z==0.0002:# and status[ii]=='1': 
            n_model_z[0]+=1.
        if z==0.002:# and status[ii]=='1': 
            n_model_z[1]+=1.
        if z==0.02:# and status[ii]=='1': 
            n_model_z[2]+=1.
            
            
        if rg==2:# and status[ii]=='1': 
            n_model_rg[0]+=1.
        if rg==8:# and status[ii]=='1': 
            n_model_rg[1]+=1.
        if rg==20:# and status[ii]=='1': 
            n_model_rg[2]+=1.
        

    ##Grouping models        
    n_msp_rv = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_msp_rv_average = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_msp_rv_average_std = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_msp_rv_median = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]

    n_msp_mass = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_msp_mass_average = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_msp_mass_average_std = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_msp_mass_median = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]

    n_msp_z = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_msp_z_average = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_msp_z_average_std = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_msp_z_median = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]


    n_msp_rg = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_msp_rg_average = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_msp_rg_average_std = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_msp_rg_median = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]

    nmsp_scatter_n2e5 = [[] for _ in range(bin_size)]
    nmsp_scatter_n4e5 = [[] for _ in range(bin_size)]
    nmsp_scatter_n8e5 = [[] for _ in range(bin_size)]
    nmsp_scatter_n16e5 = [[] for _ in range(bin_size)]

    nmsp_scatter_rv4 = [[] for _ in range(bin_size)]
    nmsp_scatter_rv2 = [[] for _ in range(bin_size)]
    nmsp_scatter_rv1 = [[] for _ in range(bin_size)]
    nmsp_scatter_rv05 = [[] for _ in range(bin_size)]

    nmsp_scatter_z00002 = [[] for _ in range(bin_size)]
    nmsp_scatter_z0002 = [[] for _ in range(bin_size)]
    nmsp_scatter_z002= [[] for _ in range(bin_size)]

    nmsp_scatter_rg2 = [[] for _ in range(bin_size)]
    nmsp_scatter_rg8 = [[] for _ in range(bin_size)]
    nmsp_scatter_rg20= [[] for _ in range(bin_size)]



    t_all = np.linspace(0, 13000., bin_size+1)
    n_sum_ej = []
    for kk in range(0, 144):
        print(paths[kk])

        ##Initial Conditions
        s=paths[kk].split('/')
        n_star=float(s[-2])
        z=float(s[-3][1:])
        rg=int(s[-4][2:])
        rv=float(s[-5][2:])
        
        t_conv = dyn.conv('t', paths[kk]+'initial.conv.sh')
        datans = np.genfromtxt(paths[kk]+'initial.esc_ns.dat')
        times = datans[:,0]#*t_conv
        n_msp = datans[:,6]; nns_mtb = datans[:,4]
        n_ej = n_msp+nns_mtb
        print(len(n_msp), len(n_ej))
        n_sum_ej.append(np.sum(n_ej))
        
        ##Interpolate the number of NS data
        #f = interp1d(times, n_ej, kind='nearest')
        #t_interpld = np.linspace(np.min(times), np.max(times), 3*bin_size)
        #n_ej_new = f(t_interpld)
        t_interpld = times
        n_ej_new = n_ej
        
    
        n_mass = [[],[],[],[]]; n_rv = [[],[],[],[]]; n_z = [[],[],[]]; n_rg = [[],[],[]]
        for jj in range(len(t_all)-1):
            #print(jj)
            n_mass_temp = [0,0,0,0]; n_rv_temp = [0,0,0,0]; n_z_temp = [0,0,0]; n_rg_temp = [0,0,0]
            count_mass = [0,0,0,0]; count_rv = [0,0,0,0]; count_z = [0,0,0]; count_rg = [0,0,0]
            
            ##Group by initial mass
            if n_star==200000.:# and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_mass_temp[0]+=n_ej_new[i]
                        count_mass[0]+=1  ##multiple time steps may belong to the same bin
        
            if n_star==400000.:# and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_mass_temp[1]+=n_ej_new[i]
                        count_mass[1]+=1
        
            if n_star==800000.:# and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_mass_temp[2]+=n_ej_new[i]
                        count_mass[2]+=1
        
            if n_star==1600000.:# and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_mass_temp[3]+=n_ej_new[i]
                        count_mass[3]+=1
            
            ##Group by initial rv   
            if rv==4.:# and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_rv_temp[0]+=n_ej_new[i]
                        count_rv[0]+=1
        
            if rv==2.:# and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_rv_temp[1]+=n_ej_new[i]
                        count_rv[1]+=1
        
            if rv==1.:# and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_rv_temp[2]+=n_ej_new[i]
                        count_rv[2]+=1
        
            if rv==0.5:# and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_rv_temp[3]+=n_ej_new[i]
                        count_rv[3]+=1
                    
        
            ##Group by metallicity
            if z==0.0002:# and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_z_temp[0]+=n_ej_new[i]
                        count_z[0]+=1
        
            if z==0.002:# and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_z_temp[1]+=n_ej_new[i]
                        count_z[1]+=1
        
            if z==0.02:# and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_z_temp[2]+=n_ej_new[i]
                        count_z[2]+=1
                    
                    
            ##Group by galactocentric distance
            if rg==2:# and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_rg_temp[0]+=n_ej_new[i]
                        count_rg[0]+=1
        
            if rg==8:# and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_rg_temp[1]+=n_ej_new[i]
                        count_rg[1]+=1
        
            if rg==20:# and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_rg_temp[2]+=n_ej_new[i]
                        count_rg[2]+=1
    
            #print(count_rv[0])
        
            for x in range(4):
                if count_rv[x]!=0:
                    n_rv_temp[x] = n_rv_temp[x]/count_rv[x]
                if count_mass[x]!=0:
                    n_mass_temp[x] = n_mass_temp[x]/count_mass[x]
                    
                n_rv[x].append(n_rv_temp[x])
                n_mass[x].append(n_mass_temp[x])
        
            for x in range(3):
                if count_z[x]!=0:
                    n_z_temp[x] = n_z_temp[x]/count_z[x]
                
                n_z[x].append(n_z_temp[x])     
                
                if count_rg[x]!=0:
                    n_rg_temp[x] = n_rg_temp[x]/count_rg[x]
                
                n_rg[x].append(n_rg_temp[x])
                            
            
        for y in range(4):
            n_msp_rv[y] = n_msp_rv[y]+np.array(n_rv[y])
            n_msp_mass[y] = n_msp_mass[y]+np.array(n_mass[y])
            n_msp_rv_average[y] = n_msp_rv_average[y] + np.array(n_rv[y])/n_model_rv[y]
            n_msp_mass_average[y] = n_msp_mass_average[y] + np.array(n_mass[y])/n_model_mass[y]
            
        for y in range(3):
            n_msp_z[y] = n_msp_z[y]+np.array(n_z[y])
            n_msp_z_average[y] = n_msp_z_average[y] + np.array(n_z[y])/n_model_z[y]
            n_msp_rg[y] = n_msp_rg[y]+np.array(n_rg[y])
            n_msp_rg_average[y] = n_msp_rg_average[y] + np.array(n_rg[y])/n_model_rg[y]


        ##Group by initial mass
        if n_star==200000.:# and status[kk]=='1':
            #print(len(n_mass[0]))
            nmsp_scatter_n2e5 = np.hstack((nmsp_scatter_n2e5, np.split(np.array(n_mass[0]),len(n_mass[0]))))
            #print(nmsp_scatter_n2e5)
        
        if n_star==400000.:# and status[kk]=='1':
            nmsp_scatter_n4e5 = np.hstack((nmsp_scatter_n4e5, np.split(np.array(n_mass[1]), len(n_mass[1]))))
        
        if n_star==800000.:# and status[kk]=='1':
            nmsp_scatter_n8e5 = np.hstack((nmsp_scatter_n8e5, np.split(np.array(n_mass[2]), len(n_mass[2]))))

        if n_star==1600000.:# and status[kk]=='1':
            nmsp_scatter_n16e5 = np.hstack((nmsp_scatter_n16e5, np.split(np.array(n_mass[3]), len(n_mass[3]))))
            

        ##Group by initial rv   
        if rv==4.:# and status[kk]=='1':
            nmsp_scatter_rv4 = np.hstack((nmsp_scatter_rv4, np.split(np.array(n_rv[0]), len(n_rv[0]))))
        
        if rv==2.:# and status[kk]=='1':
            nmsp_scatter_rv2 = np.hstack((nmsp_scatter_rv2, np.split(np.array(n_rv[1]), len(n_rv[1]))))
        
        if rv==1.:# and status[kk]=='1':
            nmsp_scatter_rv1 = np.hstack((nmsp_scatter_rv1, np.split(np.array(n_rv[2]), len(n_rv[2]))))
        
        if rv==0.5:# and status[kk]=='1':
            nmsp_scatter_rv05 = np.hstack((nmsp_scatter_rv05, np.split(np.array(n_rv[3]), len(n_rv[3]))))
                    
    
        ##Group by metallicity
        if z==0.0002:# and status[kk]=='1':
            nmsp_scatter_z00002 = np.hstack((nmsp_scatter_z00002, np.split(np.array(n_z[0]), len(n_z[0]))))
        
        if z==0.002:# and status[kk]=='1':
            nmsp_scatter_z0002 = np.hstack((nmsp_scatter_z0002, np.split(np.array(n_z[1]), len(n_z[1]))))
        
        if z==0.02:# and status[kk]=='1':
            nmsp_scatter_z002 = np.hstack((nmsp_scatter_z002, np.split(np.array(n_z[2]), len(n_z[2]))))
                    
                    
        ##Group by galactocentric distance
        if rg==2:# and status[kk]=='1':
            nmsp_scatter_rg2 = np.hstack((nmsp_scatter_rg2, np.split(np.array(n_rg[0]), len(n_rg[0]))))
        
        if rg==8:# and status[kk]=='1':
            nmsp_scatter_rg8 = np.hstack((nmsp_scatter_rg8, np.split(np.array(n_rg[1]), len(n_rg[1]))))
        
        if rg==20:# and status[kk]=='1':
            nmsp_scatter_rg20 = np.hstack((nmsp_scatter_rg20, np.split(np.array(n_rg[2]), len(n_rg[2]))))

    for ii in range(4):
        for xx in range(bin_size):
            if ii == 0:
                n_msp_mass_average_std[ii][xx]+=np.std(nmsp_scatter_n2e5[xx])
                n_msp_rv_average_std[ii][xx]+=np.std(nmsp_scatter_rv4[xx])

                n_msp_mass_median[ii][xx]+=np.median(nmsp_scatter_n2e5[xx])
                n_msp_rv_median[ii][xx]+=np.median(nmsp_scatter_rv4[xx])
            if ii == 1:
                n_msp_mass_average_std[ii][xx]+=np.std(nmsp_scatter_n4e5[xx])
                n_msp_rv_average_std[ii][xx]+=np.std(nmsp_scatter_rv2[xx])

                n_msp_mass_median[ii][xx]+=np.median(nmsp_scatter_n4e5[xx])
                n_msp_rv_median[ii][xx]+=np.median(nmsp_scatter_rv2[xx])
            if ii == 2:
                n_msp_mass_average_std[ii][xx]+=np.std(nmsp_scatter_n8e5[xx])
                n_msp_rv_average_std[ii][xx]+=np.std(nmsp_scatter_rv1[xx])

                n_msp_mass_median[ii][xx]+=np.median(nmsp_scatter_n8e5[xx])
                n_msp_rv_median[ii][xx]+=np.median(nmsp_scatter_rv1[xx])
            if ii == 3:
                n_msp_mass_average_std[ii][xx]+=np.std(nmsp_scatter_n16e5[xx])
                n_msp_rv_average_std[ii][xx]+=np.std(nmsp_scatter_rv05[xx])

                n_msp_mass_median[ii][xx]+=np.median(nmsp_scatter_n16e5[xx])
                n_msp_rv_median[ii][xx]+=np.median(nmsp_scatter_rv05[xx])

    for ii in range(3):
        for xx in range(bin_size):
            if ii == 0:
                n_msp_z_average_std[ii][xx]+=np.std(nmsp_scatter_z00002[xx])
                n_msp_rg_average_std[ii][xx]+=np.std(nmsp_scatter_rg2[xx])

                n_msp_z_median[ii][xx]+=np.median(nmsp_scatter_z00002[xx])
                n_msp_rg_median[ii][xx]+=np.median(nmsp_scatter_rg2[xx])
            if ii == 1:
                n_msp_z_average_std[ii][xx]+=np.std(nmsp_scatter_z0002[xx])
                n_msp_rg_average_std[ii][xx]+=np.std(nmsp_scatter_rg8[xx])

                n_msp_z_median[ii][xx]+=np.median(nmsp_scatter_z0002[xx])
                n_msp_rg_median[ii][xx]+=np.median(nmsp_scatter_rg8[xx])
            if ii == 2:
                n_msp_z_average_std[ii][xx]+=np.std(nmsp_scatter_z002[xx])
                n_msp_rg_average_std[ii][xx]+=np.std(nmsp_scatter_rg20[xx])

                n_msp_z_median[ii][xx]+=np.median(nmsp_scatter_z002[xx])
                n_msp_rg_median[ii][xx]+=np.median(nmsp_scatter_rg20[xx])
    
    print(n_msp_mass_average_std[0])
    print(np.sum(n_sum_ej))

    for z in range(4):
        n_msp_rv[z] = np.insert(n_msp_rv[z], 0, 0.); n_msp_rv_average[z] = np.insert(n_msp_rv_average[z], 0, 0.)
        n_msp_rv_average_std[z] = np.insert(n_msp_rv_average_std[z], 0, 0.)
        n_msp_rv_median[z] = np.insert(n_msp_rv_median[z], 0, 0.)

        n_msp_mass[z] = np.insert(n_msp_mass[z], 0, 0.); n_msp_mass_average[z] = np.insert(n_msp_mass_average[z], 0, 0.)
        n_msp_mass_average_std[z] = np.insert(n_msp_mass_average_std[z], 0, 0.)
        n_msp_mass_median[z] = np.insert(n_msp_mass_median[z], 0, 0.)

    for z in range(3):
        n_msp_z[z] = np.insert(n_msp_z[z], 0, 0.); n_msp_z_average[z] = np.insert(n_msp_z_average[z], 0, 0.)
        n_msp_z_average_std[z] = np.insert(n_msp_z_average_std[z], 0, 0.)
        n_msp_z_median[z] = np.insert(n_msp_z_median[z], 0, 0.)


        n_msp_rg[z] = np.insert(n_msp_rg[z], 0, 0.); n_msp_rg_average[z] = np.insert(n_msp_rg_average[z], 0, 0.)
        n_msp_rg_average_std[z] = np.insert(n_msp_rg_average_std[z], 0, 0.)
        n_msp_rg_median[z] = np.insert(n_msp_rg_median[z], 0, 0.)

    
    filenames = ['nmsp_mass_age_ejected.dat', 'nmsp_rv_age_ejected.dat', 'nmsp_z_age_ejected.dat', 'nmsp_rg_age_ejected.dat']
    np.savetxt('/projects/b1095/syr904/projects/GCE/catalog/'+filenames[0], np.c_[t_all, n_msp_mass[0], n_msp_mass[1], n_msp_mass[2], n_msp_mass[3], n_msp_mass_average[0], n_msp_mass_average[1], n_msp_mass_average[2], n_msp_mass_average[3], n_msp_mass_average_std[0], n_msp_mass_average_std[1], n_msp_mass_average_std[2], n_msp_mass_average_std[3], n_msp_mass_median[0], n_msp_mass_median[1], n_msp_mass_median[2], n_msp_mass_median[3]], fmt = '%f %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g', header = '1.time(Myr) 2.N_2e5 3.N_4e5 4.N_8e5 5.N_16e5 6.N_2e5_ave 7.N_4e5_ave 8.N_8e5_ave 9.N_16e5_ave 10.N_2e5_ave_std 11.N_4e5_ave_std 12.N_8e5_ave_std 13.N_16e5_ave_std 14.N_2e5_med 15.N_4e5_med 16.N_8e5_med 17.N_16e5_med', comments = '#', delimiter = ' ')
    np.savetxt('/projects/b1095/syr904/projects/GCE/catalog/'+filenames[1], np.c_[t_all, n_msp_rv[0], n_msp_rv[1], n_msp_rv[2], n_msp_rv[3], n_msp_rv_average[0], n_msp_rv_average[1], n_msp_rv_average[2], n_msp_rv_average[3], n_msp_rv_average_std[0], n_msp_rv_average_std[1], n_msp_rv_average_std[2], n_msp_rv_average_std[3], n_msp_rv_median[0], n_msp_rv_median[1], n_msp_rv_median[2], n_msp_rv_median[3]], fmt = '%f %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g', header = '1.time(Myr) 2.rv_4 3.rv_2 4.rv_1 5.rv_0.5 6.rv_4_ave 7.rv_2_ave 8.rv_1_ave 9.rv_0.5_ave 10.rv_4_ave_std 11.rv_2_ave_std 12.rv_1_ave_std 13.rv_0.5_ave_std 14.rv_4_med 15.rv_2_med 16.rv_1_med 17.rv_0.5_med', comments = '#', delimiter = ' ')
    np.savetxt('/projects/b1095/syr904/projects/GCE/catalog/'+filenames[2], np.c_[t_all, n_msp_z[0], n_msp_z[1], n_msp_z[2], n_msp_z_average[0], n_msp_z_average[1], n_msp_z_average[2], n_msp_z_average_std[0], n_msp_z_average_std[1], n_msp_z_average_std[2], n_msp_z_median[0], n_msp_z_median[1], n_msp_z_median[2]], fmt = '%f %g %g %g %g %g %g %g %g %g %g %g %g', header = '1.time(Myr) 2.z_0.0002 3.z_0.002 4.z_0.02 5.z_0.0002_ave 6.z_0.002_ave 7.z_0.02_ave 8.z_0.0002_ave_std 9.z_0.002_ave_std 10.z_0.02_ave_std 11.z_0.0002_med 12.z_0.002_med 13.z_0.02_med', comments = '#', delimiter = ' ')
    np.savetxt('/projects/b1095/syr904/projects/GCE/catalog/'+filenames[3], np.c_[t_all, n_msp_rg[0], n_msp_rg[1], n_msp_rg[2], n_msp_rg_average[0], n_msp_rg_average[1], n_msp_rg_average[2], n_msp_rg_average_std[0], n_msp_rg_average_std[1], n_msp_rg_average_std[2],n_msp_rg_median[0], n_msp_rg_median[1], n_msp_rg_median[2] ], fmt = '%f %g %g %g %g %g %g %g %g %g %g %g %g', header = '1.time(Myr) 2.rg_2 3.rg_8 4.rg_20 5.rg_2_ave 6.rg_8_ave 7.rg_20_ave 8.rg_2_ave_std 9.rg_8_ave_std 10.rg_20_ave_std 11.rg_2_med 12.rg_8_med 13.rg_20_med', comments = '#', delimiter = ' ')


##Extract and grouping the gamma luminosity from the catalog models
def extract_lgamma_msp():
    pathlist = np.genfromtxt('/projects/b1095/syr904/projects/GCE/catalog/path_allfinished_newruns_maingrid.dat', dtype=str)
    paths = pathlist[:,0]; status = pathlist[:,1]

    bin_size = 400

    n_model_mass = [0,0,0,0]; n_model_rv = [0,0,0,0]; n_model_z = [0,0,0]; n_model_rg = [0,0,0]
    for ii in range(len(paths)):

        ##Initial Conditions
        s=paths[ii].split('/')
        n_star=float(s[-2])
        z=float(s[-3][1:])
        rg=int(s[-4][2:])
        rv=float(s[-5][2:])

        if rg>2:
            continue
        
        if n_star==200000.:# and status[ii]=='1': 
            n_model_mass[0]+=1
        if n_star==400000.:# and status[ii]=='1': 
            n_model_mass[1]+=1
        if n_star==800000.:# and status[ii]=='1': 
            n_model_mass[2]+=1
        if n_star==1600000.:# and status[ii]=='1': 
            n_model_mass[3]+=1
            
        if rv==4.:# and status[ii]=='1': 
            n_model_rv[0]+=1
        if rv==2.:# and status[ii]=='1': 
            n_model_rv[1]+=1
        if rv==1.:# and status[ii]=='1': 
            n_model_rv[2]+=1
        if rv==0.5:# and status[ii]=='1': 
            n_model_rv[3]+=1
            
            
        if z==0.0002:# and status[ii]=='1': 
            n_model_z[0]+=1
        if z==0.002:# and status[ii]=='1': 
            n_model_z[1]+=1
        if z==0.02:# and status[ii]=='1': 
            n_model_z[2]+=1
            
            
        if rg==2:# and status[ii]=='1': 
            n_model_rg[0]+=1
        if rg==8:# and status[ii]=='1': 
            n_model_rg[1]+=1
        if rg==20:# and status[ii]=='1': 
            n_model_rg[2]+=1
        

    ##Grouping models        
    n_msp_rv = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_msp_rv_average = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_msp_rv_average_std = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_msp_rv_median = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]

    n_msp_mass = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_msp_mass_average = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_msp_mass_average_std = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_msp_mass_median = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]

    n_msp_z = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_msp_z_average = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_msp_z_average_std = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_msp_z_median = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]


    n_msp_rg = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_msp_rg_average = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_msp_rg_average_std = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_msp_rg_median = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    
    nmsp_scatter_n2e5 = [[] for _ in range(bin_size)]
    nmsp_scatter_n4e5 = [[] for _ in range(bin_size)]
    nmsp_scatter_n8e5 = [[] for _ in range(bin_size)]
    nmsp_scatter_n16e5 = [[] for _ in range(bin_size)]

    nmsp_scatter_rv4 = [[] for _ in range(bin_size)]
    nmsp_scatter_rv2 = [[] for _ in range(bin_size)]
    nmsp_scatter_rv1 = [[] for _ in range(bin_size)]
    nmsp_scatter_rv05 = [[] for _ in range(bin_size)]

    nmsp_scatter_z00002 = [[] for _ in range(bin_size)]
    nmsp_scatter_z0002 = [[] for _ in range(bin_size)]
    nmsp_scatter_z002= [[] for _ in range(bin_size)]

    nmsp_scatter_rg2 = [[] for _ in range(bin_size)]
    nmsp_scatter_rg8 = [[] for _ in range(bin_size)]
    nmsp_scatter_rg20= [[] for _ in range(bin_size)]


    t_all = np.linspace(0, 13000., bin_size+1)
    for kk in range(0, 144):
        print(paths[kk])

        ##Initial Conditions
        s=paths[kk].split('/')
        n_star=float(s[-2])
        z=float(s[-3][1:])
        rg=int(s[-4][2:])
        rv=float(s[-5][2:])

        N_star=s[-2]
        Z=s[-3][1:]
        Rg=s[-4][2:]
        Rv=s[-5][2:]


        if rg>2:
            continue
        
        t_conv = dyn.conv('t', paths[kk]+'initial.conv.sh')
        filename = '/projects/b1095/syr904/projects/GCE/catalog/data_lgamma/model_rv'+Rv+'_rg'+Rg+'_z'+Z+'_'+N_star+'.txt'
        datal = np.genfromtxt(filename)
        if len(datal) < 1:
            continue
        elif len(datal.shape)==1:
            #times = datal[1]
            #l_msp = datal[2]
            #t_interpld = [times,times]
            #l_msp_new = [l_msp,l_msp]
            continue
        else:
            times = datal[:,1]
            l_msp = datal[:,2]
        
            if len(times) < 3*bin_size:
                ##Interpolate the Lgamma data
                f = interp1d(times, l_msp, kind='nearest')
                t_interpld = np.linspace(np.min(times), np.max(times), 3*bin_size)
                l_msp_new = f(t_interpld)
                #print(n_msp_new)
            else:
                t_interpld = times
                l_msp_new = l_msp
    
        n_mass = [[],[],[],[]]; n_rv = [[],[],[],[]]; n_z = [[],[],[]]; n_rg = [[],[],[]]
        for jj in range(len(t_all)-1):
            #print(jj)
            n_mass_temp = [0,0,0,0]; n_rv_temp = [0,0,0,0]; n_z_temp = [0,0,0]; n_rg_temp = [0,0,0]
            count_mass = [0,0,0,0]; count_rv = [0,0,0,0]; count_z = [0,0,0]; count_rg = [0,0,0]
            
            ##Group by initial mass
            if n_star==200000.:# and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_mass_temp[0]+=l_msp_new[i]
                        count_mass[0]+=1  ##multiple time steps may belong to the same bin
        
            if n_star==400000.:# and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_mass_temp[1]+=l_msp_new[i]
                        count_mass[1]+=1
        
            if n_star==800000.:# and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_mass_temp[2]+=l_msp_new[i]
                        count_mass[2]+=1
        
            if n_star==1600000.:# and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_mass_temp[3]+=l_msp_new[i]
                        count_mass[3]+=1
            
            ##Group by initial rv   
            if rv==4.:# and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_rv_temp[0]+=l_msp_new[i]
                        count_rv[0]+=1
        
            if rv==2.:# and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_rv_temp[1]+=l_msp_new[i]
                        count_rv[1]+=1
        
            if rv==1.:# and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_rv_temp[2]+=l_msp_new[i]
                        count_rv[2]+=1
        
            if rv==0.5:# and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_rv_temp[3]+=l_msp_new[i]
                        count_rv[3]+=1
                    
        
            ##Group by metallicity
            if z==0.0002:# and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_z_temp[0]+=l_msp_new[i]
                        count_z[0]+=1
        
            if z==0.002:# and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_z_temp[1]+=l_msp_new[i]
                        count_z[1]+=1
        
            if z==0.02:# and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_z_temp[2]+=l_msp_new[i]
                        count_z[2]+=1
                    
                    
            ##Group by galactocentric distance
            if rg==2:# and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_rg_temp[0]+=l_msp_new[i]
                        count_rg[0]+=1
        
            if rg==8:# and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_rg_temp[1]+=l_msp_new[i]
                        count_rg[1]+=1
        
            if rg==20:# and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_rg_temp[2]+=l_msp_new[i]
                        count_rg[2]+=1
    
            #print(count_rv[0])
        
            for x in range(4):
                if count_rv[x]!=0:
                    n_rv_temp[x] = n_rv_temp[x]/count_rv[x]
                if count_mass[x]!=0:
                    n_mass_temp[x] = n_mass_temp[x]/count_mass[x]
                    
                n_rv[x].append(n_rv_temp[x])
                n_mass[x].append(n_mass_temp[x])
        
            for x in range(3):
                if count_z[x]!=0:
                    n_z_temp[x] = n_z_temp[x]/count_z[x]
                
                n_z[x].append(n_z_temp[x])     
                
                if count_rg[x]!=0:
                    n_rg_temp[x] = n_rg_temp[x]/count_rg[x]
                
                n_rg[x].append(n_rg_temp[x])
                            
            
        for y in range(4):
            n_msp_rv[y] = n_msp_rv[y]+np.array(n_rv[y])
            n_msp_mass[y] = n_msp_mass[y]+np.array(n_mass[y])
            n_msp_rv_average[y] = n_msp_rv_average[y] + np.array(n_rv[y])/n_model_rv[y]
            n_msp_mass_average[y] = n_msp_mass_average[y] + np.array(n_mass[y])/n_model_mass[y]
            
        for y in range(3):
            n_msp_z[y] = n_msp_z[y]+np.array(n_z[y])
            n_msp_z_average[y] = n_msp_z_average[y] + np.array(n_z[y])/n_model_z[y]
            n_msp_rg[y] = n_msp_rg[y]+np.array(n_rg[y])
            n_msp_rg_average[y] = n_msp_rg_average[y] + np.array(n_rg[y])/n_model_rg[y]


        ##Group by initial mass
        if n_star==200000.:# and status[kk]=='1':
            #print(len(n_mass[0]))
            nmsp_scatter_n2e5 = np.hstack((nmsp_scatter_n2e5, np.split(np.array(n_mass[0]),len(n_mass[0]))))
            #print(nmsp_scatter_n2e5)
        
        if n_star==400000.:# and status[kk]=='1':
            nmsp_scatter_n4e5 = np.hstack((nmsp_scatter_n4e5, np.split(np.array(n_mass[1]), len(n_mass[1]))))
        
        if n_star==800000.:# and status[kk]=='1':
            nmsp_scatter_n8e5 = np.hstack((nmsp_scatter_n8e5, np.split(np.array(n_mass[2]), len(n_mass[2]))))

        if n_star==1600000.:# and status[kk]=='1':
            nmsp_scatter_n16e5 = np.hstack((nmsp_scatter_n16e5, np.split(np.array(n_mass[3]), len(n_mass[3]))))
            

        ##Group by initial rv   
        if rv==4.:# and status[kk]=='1':
            nmsp_scatter_rv4 = np.hstack((nmsp_scatter_rv4, np.split(np.array(n_rv[0]), len(n_rv[0]))))
        
        if rv==2.:# and status[kk]=='1':
            nmsp_scatter_rv2 = np.hstack((nmsp_scatter_rv2, np.split(np.array(n_rv[1]), len(n_rv[1]))))
        
        if rv==1.:# and status[kk]=='1':
            nmsp_scatter_rv1 = np.hstack((nmsp_scatter_rv1, np.split(np.array(n_rv[2]), len(n_rv[2]))))
        
        if rv==0.5:# and status[kk]=='1':
            nmsp_scatter_rv05 = np.hstack((nmsp_scatter_rv05, np.split(np.array(n_rv[3]), len(n_rv[3]))))
                    
    
        ##Group by metallicity
        if z==0.0002:# and status[kk]=='1':
            nmsp_scatter_z00002 = np.hstack((nmsp_scatter_z00002, np.split(np.array(n_z[0]), len(n_z[0]))))
        
        if z==0.002:# and status[kk]=='1':
            nmsp_scatter_z0002 = np.hstack((nmsp_scatter_z0002, np.split(np.array(n_z[1]), len(n_z[1]))))
        
        if z==0.02:# and status[kk]=='1':
            nmsp_scatter_z002 = np.hstack((nmsp_scatter_z002, np.split(np.array(n_z[2]), len(n_z[2]))))
                    
                    
        ##Group by galactocentric distance
        if rg==2:# and status[kk]=='1':
            nmsp_scatter_rg2 = np.hstack((nmsp_scatter_rg2, np.split(np.array(n_rg[0]), len(n_rg[0]))))
        
        if rg==8:# and status[kk]=='1':
            nmsp_scatter_rg8 = np.hstack((nmsp_scatter_rg8, np.split(np.array(n_rg[1]), len(n_rg[1]))))
        
        if rg==20:# and status[kk]=='1':
            nmsp_scatter_rg20 = np.hstack((nmsp_scatter_rg20, np.split(np.array(n_rg[2]), len(n_rg[2]))))

    for ii in range(4):
        for xx in range(bin_size):
            if ii == 0:
                n_msp_mass_average_std[ii][xx]+=np.std(nmsp_scatter_n2e5[xx])
                n_msp_rv_average_std[ii][xx]+=np.std(nmsp_scatter_rv4[xx])

                n_msp_mass_median[ii][xx]+=np.median(nmsp_scatter_n2e5[xx])
                n_msp_rv_median[ii][xx]+=np.median(nmsp_scatter_rv4[xx])
            if ii == 1:
                n_msp_mass_average_std[ii][xx]+=np.std(nmsp_scatter_n4e5[xx])
                n_msp_rv_average_std[ii][xx]+=np.std(nmsp_scatter_rv2[xx])

                n_msp_mass_median[ii][xx]+=np.median(nmsp_scatter_n4e5[xx])
                n_msp_rv_median[ii][xx]+=np.median(nmsp_scatter_rv2[xx])
            if ii == 2:
                n_msp_mass_average_std[ii][xx]+=np.std(nmsp_scatter_n8e5[xx])
                n_msp_rv_average_std[ii][xx]+=np.std(nmsp_scatter_rv1[xx])

                n_msp_mass_median[ii][xx]+=np.median(nmsp_scatter_n8e5[xx])
                n_msp_rv_median[ii][xx]+=np.median(nmsp_scatter_rv1[xx])
            if ii == 3:
                n_msp_mass_average_std[ii][xx]+=np.std(nmsp_scatter_n16e5[xx])
                n_msp_rv_average_std[ii][xx]+=np.std(nmsp_scatter_rv05[xx])

                n_msp_mass_median[ii][xx]+=np.median(nmsp_scatter_n16e5[xx])
                n_msp_rv_median[ii][xx]+=np.median(nmsp_scatter_rv05[xx])

    for ii in range(3):
        for xx in range(bin_size):
            if ii == 0:
                n_msp_z_average_std[ii][xx]+=np.std(nmsp_scatter_z00002[xx])
                n_msp_rg_average_std[ii][xx]+=np.std(nmsp_scatter_rg2[xx])

                n_msp_z_median[ii][xx]+=np.median(nmsp_scatter_z00002[xx])
                n_msp_rg_median[ii][xx]+=np.median(nmsp_scatter_rg2[xx])
            if ii == 1:
                n_msp_z_average_std[ii][xx]+=np.std(nmsp_scatter_z0002[xx])
                n_msp_rg_average_std[ii][xx]+=np.std(nmsp_scatter_rg8[xx])

                n_msp_z_median[ii][xx]+=np.median(nmsp_scatter_z0002[xx])
                n_msp_rg_median[ii][xx]+=np.median(nmsp_scatter_rg8[xx])
            if ii == 2:
                n_msp_z_average_std[ii][xx]+=np.std(nmsp_scatter_z002[xx])
                n_msp_rg_average_std[ii][xx]+=np.std(nmsp_scatter_rg20[xx])

                n_msp_z_median[ii][xx]+=np.median(nmsp_scatter_z002[xx])
                n_msp_rg_median[ii][xx]+=np.median(nmsp_scatter_rg20[xx])
    
    print(n_msp_mass_average_std[0])

    for z in range(4):
        n_msp_rv[z] = np.insert(n_msp_rv[z], 0, 0.); n_msp_rv_average[z] = np.insert(n_msp_rv_average[z], 0, 0.)
        n_msp_rv_average_std[z] = np.insert(n_msp_rv_average_std[z], 0, 0.)
        n_msp_rv_median[z] = np.insert(n_msp_rv_median[z], 0, 0.)

        n_msp_mass[z] = np.insert(n_msp_mass[z], 0, 0.); n_msp_mass_average[z] = np.insert(n_msp_mass_average[z], 0, 0.)
        n_msp_mass_average_std[z] = np.insert(n_msp_mass_average_std[z], 0, 0.)
        n_msp_mass_median[z] = np.insert(n_msp_mass_median[z], 0, 0.)

    for z in range(3):
        n_msp_z[z] = np.insert(n_msp_z[z], 0, 0.); n_msp_z_average[z] = np.insert(n_msp_z_average[z], 0, 0.)
        n_msp_z_average_std[z] = np.insert(n_msp_z_average_std[z], 0, 0.)
        n_msp_z_median[z] = np.insert(n_msp_z_median[z], 0, 0.)


        n_msp_rg[z] = np.insert(n_msp_rg[z], 0, 0.); n_msp_rg_average[z] = np.insert(n_msp_rg_average[z], 0, 0.)
        n_msp_rg_average_std[z] = np.insert(n_msp_rg_average_std[z], 0, 0.)
        n_msp_rg_median[z] = np.insert(n_msp_rg_median[z], 0, 0.)

    
    filenames = ['lgamma_msp_mass_age_rg2.dat', 'nmsp_rv_age.dat', 'nmsp_z_age.dat', 'nmsp_rg_age.dat']
    np.savetxt('/projects/b1095/syr904/projects/GCE/catalog/'+filenames[0], np.c_[t_all, n_msp_mass[0], n_msp_mass[1], n_msp_mass[2], n_msp_mass[3], n_msp_mass_average[0], n_msp_mass_average[1], n_msp_mass_average[2], n_msp_mass_average[3], n_msp_mass_average_std[0], n_msp_mass_average_std[1], n_msp_mass_average_std[2], n_msp_mass_average_std[3], n_msp_mass_median[0], n_msp_mass_median[1], n_msp_mass_median[2], n_msp_mass_median[3]], fmt = '%f %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g', header = '1.time(Myr) 2.L_2e5 3.L_4e5 4.L_8e5 5.L_16e5 6.L_2e5_ave 7.L_4e5_ave 8.L_8e5_ave 9.L_16e5_ave 10.L_2e5_ave_std 11.L_4e5_ave_std 12.L_8e5_ave_std 13.L_16e5_ave_std 14.L_2e5_med 15.L_4e5_med 16.L_8e5_med 17.L_16e5_med', comments = '#', delimiter = ' ')
    #np.savetxt('/projects/b1095/syr904/projects/GCE/catalog/'+filenames[1], np.c_[t_all, n_msp_rv[0], n_msp_rv[1], n_msp_rv[2], n_msp_rv[3], n_msp_rv_average[0], n_msp_rv_average[1], n_msp_rv_average[2], n_msp_rv_average[3], n_msp_rv_average_std[0], n_msp_rv_average_std[1], n_msp_rv_average_std[2], n_msp_rv_average_std[3], n_msp_rv_median[0], n_msp_rv_median[1], n_msp_rv_median[2], n_msp_rv_median[3]], fmt = '%f %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g', header = '1.time(Myr) 2.rv_4 3.rv_2 4.rv_1 5.rv_0.5 6.rv_4_ave 7.rv_2_ave 8.rv_1_ave 9.rv_0.5_ave 10.rv_4_ave_std 11.rv_2_ave_std 12.rv_1_ave_std 13.rv_0.5_ave_std 14.rv_4_med 15.rv_2_med 16.rv_1_med 17.rv_0.5_med', comments = '#', delimiter = ' ')
    #np.savetxt('/projects/b1095/syr904/projects/GCE/catalog/'+filenames[2], np.c_[t_all, n_msp_z[0], n_msp_z[1], n_msp_z[2], n_msp_z_average[0], n_msp_z_average[1], n_msp_z_average[2], n_msp_z_average_std[0], n_msp_z_average_std[1], n_msp_z_average_std[2], n_msp_z_median[0], n_msp_z_median[1], n_msp_z_median[2]], fmt = '%f %g %g %g %g %g %g %g %g %g %g %g %g', header = '1.time(Myr) 2.z_0.0002 3.z_0.002 4.z_0.02 5.z_0.0002_ave 6.z_0.002_ave 7.z_0.02_ave 8.z_0.0002_ave_std 9.z_0.002_ave_std 10.z_0.02_ave_std 11.z_0.0002_med 12.z_0.002_med 13.z_0.02_med', comments = '#', delimiter = ' ')
    #np.savetxt('/projects/b1095/syr904/projects/GCE/catalog/'+filenames[3], np.c_[t_all, n_msp_rg[0], n_msp_rg[1], n_msp_rg[2], n_msp_rg_average[0], n_msp_rg_average[1], n_msp_rg_average[2], n_msp_rg_average_std[0], n_msp_rg_average_std[1], n_msp_rg_average_std[2],n_msp_rg_median[0], n_msp_rg_median[1], n_msp_rg_median[2] ], fmt = '%f %g %g %g %g %g %g %g %g %g %g %g %g', header = '1.time(Myr) 2.rg_2 3.rg_8 4.rg_20 5.rg_2_ave 6.rg_8_ave 7.rg_20_ave 8.rg_2_ave_std 9.rg_8_ave_std 10.rg_20_ave_std 11.rg_2_med 12.rg_8_med 13.rg_20_med', comments = '#', delimiter = ' ')

    
