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
import ns

twopi = 2*np.pi


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
    sourcedir=np.genfromtxt(pathlist,dtype=str)
    paths = sourcedir[:,0]; status = sourcedir[:,1]

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

        fgamma=open('/projects/b1095/syr904/projects/GCE/catalog/data_lgamma/model_rv'+rv+'_rg'+rg+'_z'+z+'_'+n_star+'.txt', 'a+')
        fgamma.write('#1.Model 2.T(Myr) 3.Lmsp\n')
        
        allkey = list(Counter(time).keys())
        #print(allkey)
        for x in range(len(allkey)):
        	theL = 0
        	for y in range(len(time)):
        		if time[y]==allkey[x]:
        			theL+=Lmsp[y]

        	fgamma.write('%d %f %e\n'%(ii, allkey[x], theL))

        print(ii)

    fgamma.close()
        

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











