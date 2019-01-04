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
import history_maker_full4 as hi4
import history_maker_full5 as hi5
import history_cmc as hic
import dynamics as dyn
import scripts3


#path = '/projects/b1011/kyle/cmc/branches/cmc-mpi/m22_project/rundir/m10'
#fo = 'm10_400000e5_4.5_1.0_0.05_FULL'
#fname = 'initial.pulsars.dat'
yearsc=31557600.
twopi=6.283185307179586
Gconst=6.674*10**-8 ##cm3*g-1*s-2
clight=3*10**10 ##cm/s
Msun=2*10**33 ##gram
AU=1.496*10**13  ##cm
PC=3.086*10**18  ##cm
Kconst=9.87*10**-48 ##yr/G^2
Lsun=4.02*10**16 ##mJy*kpc^2



def readdata_freire():
	#from astropy.extern.six.moves.urllib import request
	#url = 'http://www.naic.edu/~pfreire/GCpsr.txt'
	#open('GCpsr.txt', 'wb').write(request.urlopen(url).read())
	
	K=5.7*10**19
	Ps=[]; Pdots=[]; Bs=[]; Pb=[]; Pdotb=[]; Bb=[]  ##P unit ms, dpdt has a factor of 10^-20
	with open('/projects/b1011/syr904/projects/PULSAR/GCpsr.txt', 'rb') as f:
	    for _ in xrange(4):
	        next(f)
	    for line in f:
	        data=line.split()
	        if not data: continue
	        if str(data[0][0])=='J' or data[0][0]=='B': 
	            if str(data[3])!='*':
	                if str(data[3][0])=='<':continue    
	                if str(data[-1])=='i': 
	                    Ps.append(float(data[2]))
	                    dpdts=data[3].split('(')
	                    if len(dpdts)>1:
	                        errs=dpdts[1].split(')')
	                        if errs[1]!='':
	                            Pdots.append(float(dpdts[0])*10**-15)
	                        else: 
	                            Pdots.append(float(dpdts[0])*10**-20)
	                    else: 
	                        Pdots.append(float(dpdts[0])*10**-20)
	                else: 
	                    Pb.append(float(data[2]))
	                    dpdtb=data[3].split('(')
	                    if len(dpdtb)>1:
	                        errb=dpdtb[1].split(')')
	                        if errb[1]!='':
	                            Pdotb.append(float(dpdtb[0])*10**-15)
	                        else:
	                            Pdotb.append(float(dpdtb[0])*10**-20)
	                    else:
	                        Pdotb.append(float(dpdtb[0])*10**-20)
    
	#print Pdots, Pdotb, Ps, Pb
	#Bs=K*np.sqrt(np.abs(Pdots)*np.array(Ps)*0.001)
	#Bb=K*np.sqrt(np.abs(Pdotb)*np.array(Pb)*0.001)
	Ps=np.array(Ps); Pb=np.array(Pb); Pdots=np.array(Pdots); Pdotb=np.array(Pdotb)
	#Bs=np.array(Bs); Bb=np.array(Bb)
	#print Bs, Bb
	return Ps, Pb, Pdots, Pdotb


def conv_dict(): return {'l':15, 't':19, 'm':7}    #?


def conv(unit,filepath):   # Returns the unit conversion multiplier given a simulation's *.conv.sh file and a unit (either 'l' or 't')
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


def find_pulsar_last(sourcedir, folder, filename): 
#Extract pulsars id, spin, B field and mass from the last time step of the initial.pulsars.dat.
	data = np.genfromtxt(sourcedir+'/'+folder+'/'+filename, usecols = (1, 2, 15, 16, 23))
	time = data[:,0]; starid = data[:,1]; spin = data[:,2]; field = data[:,3]; mass = data[:,4]
	t=[]; ID=[]; P=[]; B=[]; M_1=[]
	for k in range(len(time)):
		if time[k] == time[-1]:
        		t.append(time[k])
               		ID.append(starid[k])
			P.append(spin[k])
			B.append(field[k])
			M_1.append(mass[k])
	np.savetxt(folder+'_pulsarf.dat', np.c_[ID, P, B, M_1], fmt ='%d %f %e %f', delimiter = ' ', header = '1.ID, 2.P, 3.B, 4.M1', comments = '#')



def find_history(sourcedir, ids):    #Find star history in snapshots.		
	snaps = np.sort(glob(sourcedir+'/'+'initial.snap*.dat.gz'))
	t_conv = conv('t',sourcedir+'/'+'initial.conv.sh')
	#ids=gen_snap_bnslist(sourcedir, folder)
	#ids=[148267]
	for j in range(len(ids)):
		id_temp=ids[j]
		m0=[]; m1=[]; id0=[]; id1=[]; a=[]; e=[]; k0=[]; k1=[]; P0=[]; P1=[]; B0=[]; B1=[]; t=[] 
		for i in range(len(snaps)):
			time = get_time(snaps[i])*t_conv
			with gzip.open(snaps[i],'r') as fi:
					for _ in xrange(2):
						next(fi)
					for line in fi:
						data1=line.split()
						if int(data1[7])==1:
							if int(data1[10])==id_temp or int(data1[11])==id_temp:
								m0.append(float(data1[8])); m1.append(float(data1[9])); id0.append(int(data1[10])); id1.append(int(data1[11])); a.append(float(data1[12])); e.append(float(data1[13])); k0.append(int(data1[17])); k1.append(int(data1[18])); P0.append(float(data1[45])); P1.append(float(data1[46])); B0.append(float(data1[47])); B1.append(float(data1[48]))
								t.append(time)
						if int(data1[7])!=1:
							if int(data1[0])==id_temp:
								m0.append(float(data1[1])); m1.append(float(-100)); id0.append(int(data1[0])); id1.append(int(-100)); k0.append(int(data1[14])); k1.append(int(-100)); a.append(float(-100)); e.append(float(-100)); P0.append(float(-100)); P1.append(float(-100)); B0.append(float(-100)); B1.append(float(-100))
								t.append(time)

		idname=str(id_temp)
		np.savetxt(sourcedir+'/'+'NS'+'/'+idname+'_snap.dat', np.c_[t, id0, id1, m0, m1, k0, k1, a, e, P0, P1, B0, B1], fmt ='%f %d %d %f %f %d %d %f %f %f %f %f %f', delimiter= ' ', header = '1.time, 2.id0, 3.id1, 4.m0, 5.m1, 6.k0, 7.k1, 8.a, 9.e, 10.p0, 11.p1, 12.B0, 13.B1', comments = '#')



def find_ns_mt(sourcedir, folder): #Find mass transfering neutron stars, including collision stars.
	t_conv = conv('t',sourcedir+'/'+folder+'/'+'initial.conv.sh')
	snaps = np.sort(glob(sourcedir+'/'+folder+'/'+'initial.snap*.dat.gz'))
	#m0_temp=[]; m1_temp=[]; id0_temp=[]; id1_temp=[]; m_temp=[]; id_temp=[]; t_temp=[]; t0_temp=[]; t1_temp=[]
	nsmt, tnsmt = [], []
	nsmt.append(0)
	for i in range(len(snaps)-1):
		print i
		tnsmt.append(get_time(snaps[i])*t_conv/1000)    #Time unit: Gyr
		m_tempa, id_tempa = find_ns(snaps[i]); m_tempb, id_tempb = find_ns(snaps[i+1])
		Nmt = compare_mass(m_tempa, m_tempb, id_tempa, id_tempb)
		nsmt.append(Nmt)
        #data_ns_temp = np.genfromtxt(snaps[i],usecols=(0,1,7,8,9,10,11,14,17,18)) 
        #star_id=data_ns_temp[:,0]; star_m=data_ns_temp[:,1]; binflag=data_ns_temp[:,2]; mass0=data_ns_temp[:,3]; mass1=data_ns_temp[:,4]; id_0=data_ns_temp[:,5]; id_1=data_ns_temp[:,6]; st=data_ns_temp[:,7]; st1=data_ns_temp[:,8]; st2=data_ns_temp[:,9]
		#for j in range(len(binflag)):
        #    		if binflag[j] == 0 and st[j]  == 13: 
		#		id_temp.append(star_id[j]); m_temp.append(star_m[j]); t_temp.append(get_time(snaps[i])*t_conv/1000)
        #    		if binflag[j] == 1 and st1[j] == 13: 
		#		id0_temp.append(id_0[j]); m0_temp.append(mass0[j]); t0_temp.append(get_time(snaps[i])*t_conv/1000) 
        #    		if binflag[j] == 1 and st2[j] == 13: 
		#		id1_temp.append(id_1[j]); m1_temp.append(mass1[j]); t1_temp.append(get_time(snaps[i])*t_conv/1000)

   # np.savetxt(savename+'.txt', np.c_[tns, ns], delimiter = ' ', header = '1.t, 2.NS', comments = '#')

	tnsmt.append(get_time(snaps[-1])*t_conv/1000)    #Time unit: Gyr
	tnsmt = np.delete(tnsmt, 0)
	nsmt = np.delete(nsmt, 0)
	#print tnsmt, nsmt 
	return tnsmt, nsmt


def compare_mass(m_a, m_b, id_a, id_b, idcoll=544644, idbefore=313668):   #idm=544644(mm=0) id1=313668(m1=1.26071)
	mt=0
	for k in range(len(id_b)):
		x=id_b[k]; y=m_b[k]
		for m in range(len(id_a)):
			if id_a[m]==x and m_a[m]<y: mt+=1; print x
		
	return mt



def find_ns(snapshot): 
	m0_temp=[]; m1_temp=[]; id0_temp=[]; id1_temp=[]; m_temp=[]; id_temp=[]#; t_temp=[]; t0_temp=[]; t1_temp=[]
	data_ns_temp = np.genfromtxt(snapshot,usecols=(0,1,7,8,9,10,11,14,17,18))
	star_id=data_ns_temp[:,0]; star_m=data_ns_temp[:,1]; binflag=data_ns_temp[:,2]; mass0=data_ns_temp[:,3]; mass1=data_ns_temp[:,4]; id_0=data_ns_temp[:,5]; id_1=data_ns_temp[:,6]; st=data_ns_temp[:,7]; st1=data_ns_temp[:,8]; st2=data_ns_temp[:,9]
	for j in range(len(binflag)):
        	if st[j]  == 13:
			id_temp.append(int(star_id[j])); m_temp.append(star_m[j])#; t_temp.append(get_time(snaps[i])*t_conv/1000)
        	if st1[j] == 13:
			id0_temp.append(int(id_0[j])); m0_temp.append(mass0[j])#; t0_temp.append(get_time(snaps[i])*t_conv/1000)
        	if st2[j] == 13:
			id1_temp.append(int(id_1[j])); m1_temp.append(mass1[j])#; t1_temp.append(get_time(snaps[i])*t_conv/1000)

	M_temp = np.concatenate((m_temp, m0_temp, m1_temp), axis=0); ID_temp = np.concatenate((id_temp, id0_temp, id1_temp), axis=0)

	return ID_temp, M_temp


def find_ns2D(pathlist, start, end):
	sourcedir=np.genfromtxt(pathlist, dtype=str)
	#r2d_ns=[]; model=[]
	fi=open('/projects/b1011/syr904/projects/PULSAR/kickgrid_runs/paper/data/ns2Dradius_newmodel_new.dat', 'a+', 0)
	for i in range(start, end):
		pref='initial'
		filepath=sourcedir[i]
		filestr=filepath+'/'+pref
		snapproj=np.sort(glob(filestr+'.snap*.2Dproj.dat'))
		lastproj=snapproj[-1]
	
		dataproj=np.genfromtxt(lastproj)
		R2D=dataproj[:,0]; kstar=dataproj[:,3]; k0=dataproj[:,5]; k1=dataproj[:,6]
		for j in range(len(R2D)):
			if kstar[j]==13 or k0[j]==13 or k1[j]==13: 
				r2d_ns=R2D[j]
				fi.write('%d %f\n'%(i, r2d_ns))

	fi.close()




def plot_massdistri_hist(Ms):
	plt.figure()
	weights=np.ones_like(Ms)/float(len(Ms))
	plt.hist(Ms, 20, weights=weights, alpha=0.5)
	plt.yscale('log')
	plt.xlabel(r'$M_\odot$')
	plt.show()



def find_ns_position(sourcedir, folder, no):
	snap = sourcedir+'/'+folder+'/'+'initial.snap'+no+'.dat.gz' 
        rns=[]
        data_r = np.genfromtxt(snap,usecols=(2,7,14,17,18))
        r = data_r[:,0]; binflag=data_r[:,1]; st=data_r[:,2]; st1=data_r[:,3]; st2=data_r[:,4]
        for j in range(len(binflag)):
        	if binflag[j] == 0 and st[j]  == 13: rns.append(r[j])
        	if binflag[j] == 1 and st1[j] == 13: rns.append(r[j])
                if binflag[j] == 1 and st2[j] == 13: rns.append(r[j])	
	
	#print rns
	return rns


def find_nsbinary(sourcedir, folder):     #Find NS binary in all the snapshots
	yearsc=31557600
	twopi=6.283185307179586
	snaps = np.sort(glob(sourcedir+'/'+folder+'/'+'initial.snap*.dat.gz'))
	t_conv = conv('t',sourcedir+'/'+folder+'/'+'initial.conv.sh')
	t=[]; k0=[]; k1=[]; m0=[]; m1=[]; id0=[]; id1=[]; P0=[]; P1=[]; B0=[]; B1=[]
	for i in range(len(snaps)):
		time=t_conv*get_time(snaps[i])   #Myr
		data=np.genfromtxt(snaps[i], usecols=(7, 8, 9, 10, 11, 17, 18, 45, 46, 47, 48))
		binflag=data[:,0]; mass0=data[:,1]; mass1=data[:,2]; sid0=data[:,3]; sid1=data[:,4]; sk0=data[:,5]; sk1=data[:,6]; spin0=data[:,7]; spin1=data[:,8]; b0=data[:,9]; b1=data[:,10]
		for j in range(len(binflag)):
			if binflag[j]==1 and (sk0[j]==13 or sk1[j]==13):
				t.append(time)
				m0.append(mass0[j]); m1.append(mass1[j]); id0.append(sid0[j]); id1.append(sid1[j]); k0.append(sk0[j]); k1.append(sk1[j]); P0.append(twopi*yearsc/spin0[j]); P1.append(twopi*yearsc/spin1[j]); B0.append(b0[j]); B1.append(b1[j])

	np.savetxt('snap_bns.dat', np.c_[t, m0, m1, id0, id1, k0, k1, P0, P1, B0, B1], fmt ='%f %f %f %d %d %d %d %f %f %f %f', delimiter= ' ', header = '1.time, 2.m0, 3.m1, 4.id0, 5.id1, 6.k0, 7.k1, 8.P0(s), 9.P1(s), 10.B0, 11.B1', comments = '#')

	

def gen_snap_bnslist(sourcedir, folder):      #Find the list of ids of binary ns in snapshots
	data1=np.genfromtxt(sourcedir+'/'+folder+'/'+'snap_bns.dat')
	t=data1[:,0]; m0=data1[:,1]; m1=data1[:,2]; id0=data1[:,3]; id1=data1[:,4]; k0=data1[:,5]; k1=data1[:,6]; p0=data1[:,7]; p1=data1[:,8]
	idlist=[]
	for i in range(len(id0)):
		if k0[i]==13:idlist.append(int(id0[i]))
		if k1[i]==13:idlist.append(int(id1[i]))

	snap_idlist=Counter(idlist).keys()
	#print snap_idlist
	return snap_idlist

	#p_idlist=[]
	#pt=[]; pm0=[]; pm1=[]; 
	#for j in range(len(snap_idlist)):
	#	no=snap_idlist[j]
	#	with open(sourcedir+'/'+folder+'/'+'initial.morepulsars.dat') as f:
	#		next(f)
	#		for line in f:
	#			data2=line.split()
	#			if int(data2[8])==0:
	#				if float(data2[9])==no or float(data2[10])==no:



def gen_idhistory(sourcedir, folder):
	list=gen_snap_bnslist(sourcedir, folder)
	for i in range(len(list)):
		no=int(list[i])
		pt=[]; pm0=[]; pm1=[]; pid0=[]; pid1=[]; pk0=[]; pk1=[]; pp=[]; pB=[]; pa=[]; pe=[]
		with open(sourcedir+'/'+folder+'/'+'initial.morepulsars.dat') as f:
			next(f)
			for line in f:
				data2=line.split()
				if int(data2[8])==1:
					if int(data2[9])==no or int(data2[10])==no:
						pt.append(float(data2[1])); pid0.append(int(data2[9])); pid1.append(int(data2[10])); pm0.append(float(data2[11])); pm1.append(float(data2[12])); pk0.append(int(data2[19])); pk1.append(int(data2[20])); pp.append(float(data2[17])); pB.append(float(data2[15])); pa.append(float(data2[21])); pe.append(float(data2[22]))
					#if int(data2[10])==no:
						#pt.append(float(data2[1])); pm0.append(float(data2[11])); pm1.append(float(data2[12])); pk0.append(float(data2[19])); pk1.append(float(data2[20])); pp.append(float(data2[18])); pB.append(float(data2[16])); pa.append(float(data2[21])); pe.append(float(data2[22]))

		name=str(no)
		np.savetxt(sourcedir+'/'+folder+'/'+'history'+'/'+name+'_pulsar.dat', np.c_[pt, pm0, pm1, pid0, pid1, pk0, pk1, pp, pB, pa, pe], fmt ='%f %f %f %d %d %f %e %f %f', delimiter= ' ', header = '1.time, 2.m0, 3.m1, 4.id0, 5.id1, 6.k0, 7.k1, 8.P, 9.B, 10.a, 11.e', comments = '#')

	
def get_snap_Nns(snapshot, dynlim):
	nmtb=0; nmsp=0; npulsarsin=0; npulsarbin=0; nns=0; nnsdyn=0
	mspid=[]; mspcomid=[]

	#snaps=np.sort(glob(sourcedir+'/'+'initial.snap*.dat.gz'))
	#lastsnap=snaps[-1]
	#print lastsnap
	with gzip.open(snapshot, 'r') as fsnap:
		for _ in xrange(2):
			next(fsnap)
		for line in fsnap:
			datasnap=line.split()
			if int(datasnap[7])!=1:
				if int(datasnap[14])==13:
					nns+=1
					if int(datasnap[0])>dynlim: nnsdyn+=1
					spin=twopi*yearsc/float(datasnap[59])
					deathcut=(spin**2)*(0.17*10**12)
					if deathcut<float(datasnap[60]): npulsarsin+=1
					if spin<=0.03: 
						nmsp+=1; mspid.append(int(datasnap[0])); mspcomid.append(-100)
						print int(datasnap[0]), -100
			if int(datasnap[7])==1:
				if int(datasnap[17])==13:
					nns+=1
					if int(datasnap[10])>dynlim: nnsdyn+=1
					spin0=twopi*yearsc/float(datasnap[45])
					deathcut0=(spin0**2)*(0.17*10**12)
					if deathcut0<float(datasnap[47]): npulsarbin+=1
					if float(datasnap[44])>=1: nmtb+=1	
					if spin0<=0.03: 
						nmsp+=1; mspid.append(int(datasnap[10])); mspcomid.append(int(datasnap[11]))
						print int(datasnap[10]), int(datasnap[11])
				if int(datasnap[18])==13:
					nns+=1
					if int(datasnap[11])>dynlim: nnsdyn+=1
					spin1=twopi*yearsc/float(datasnap[46])
					deathcut1=(spin1**2)*(0.17*10**12)
					if deathcut1<float(datasnap[48]): npulsarbin+=1
					if float(datasnap[43])>=1: nmtb+=1
					if spin1<=0.03: 
						nmsp+=1; mspid.append(int(datasnap[11])); mspcomid.append(int(datasnap[10]))
						print int(datasnap[11]), int(datasnap[10])

	return npulsarsin, npulsarbin, nmsp, nmtb, mspid, mspcomid, nns, nnsdyn


def get_allsnap_Nns(sourcedir):
	snaps=np.sort(glob(sourcedir+'/'+'initial.snap*.dat.gz'))
	databh=np.genfromtxt(sourcedir+'/'+'initial.bh.dat')
	timetot=databh[:,1]; nbhtot=databh[:,2]
	NP=[]; NMSP=[]; NMT=[]; NBH=[]; T=[]
	for k in range(len(snaps)):
		Np, Nmsp, Nmt=get_snap_Nns(snaps[k])
		NP.append(Np); NMSP.append(Nmsp); NMT.append(Nmt)
		time=get_time(snaps[k])
		t_conv = conv('t',sourcedir+'/'+'initial.conv.sh')
		T.append(time*t_conv)
		for j in range(len(nbhtot)):
			if timetot[j]==time: 
				NBH.append(float(nbhtot[j]))

	
	return NP, NMSP, NMT, NBH, T


def get_allmodel_MSPID(pathlist, start, end):
	sourcedir=np.genfromtxt(pathlist, dtype=str)
	MSPID=[]; MSPCOMID=[]; model=[]
	for i in range(start, end):
		filepath=sourcedir[i]
		pref='initial'
		filestr=filepath+'/'+pref
		snaps=np.sort(glob(filestr+'.snap*.dat.gz'))
		lastsnap=snaps[-1]
		Npulsar, Nmsp, Nmtb, Mspid, Mspcomid, Nns, Nnsdyn=get_snap_Nns(lastsnap, 800000)
		MSPID=MSPID+Mspid; MSPCOMID=MSPCOMID+Mspcomid
		temp=list(np.full_like(Mspid, fill_value=i))
		model=model+temp
		print i

	print len(model), len(MSPID)

	np.savetxt('/projects/b1011/syr904/projects/PULSAR/kickgrid_runs/paper/data/kickgrid_mspid_newmodel.dat', np.c_[model, MSPID, MSPCOMID], fmt='%d %d %d', header='Model ID0 ID1', comments='#', delimiter= ' ')


def get_allmodel_MSPID_10to12Gyr(pathlist, start, end):
	dynno=[800000,800000,800000,800000,800000,800000,800000,800000,800000,800000,800000,800000,800000,800000,800000,800000,800000,800000,800000,800000,800000,800000,800000,800000,800000,3500000]
	sourcedir=np.genfromtxt(pathlist, dtype=str)
	MSPID=[]; MSPCOMID=[]; model=[]; time=[]
	for i in range(start, end):
		if end==1: sourcedir=[str(sourcedir)]; print len(sourcedir)
		filepath=sourcedir[i]
		pref='initial'
		filestr=filepath+'/'+pref
		t_conv=conv('t', filestr+'.conv.sh')
		snaps=np.sort(glob(filestr+'.snap*.dat.gz'))
		for j in range(len(snaps)):
			t=get_time(snaps[j])*t_conv/1000.
			print t, j
			if t>=10.:
				Npulsar, Nmsp, Nmtb, Mspid, Mspcomid, Nns, Nnsdyn=get_snap_Nns(snaps[j], dynno[i])
				MSPID=MSPID+Mspid; MSPCOMID=MSPCOMID+Mspcomid
				model_temp=list(np.full_like(Mspid, fill_value=i))
				model=model+model_temp
				time_temp=list(np.full_like(Mspid, fill_value=t, dtype=np.double))
				time=time+time_temp
				
				#break

		print i

	print len(model), len(MSPID)

	np.savetxt('/projects/b1011/syr904/projects/PULSAR/kickgrid_runs/paper/data/kickgrid_mspid_10to12Gyr_newmodel.dat', np.c_[time, model, MSPID, MSPCOMID], fmt='%f %d %d %d', header='Time Model ID0 ID1', comments='#', delimiter=' ')



def print_Nbh_Npulsar(start, end, pathlist):
	dynno=[800000,800000,800000,800000,800000,800000,800000,800000,800000,800000,800000,800000,800000,800000,800000,800000,800000,800000,800000,800000,800000,800000,800000,800000,800000,3500000]
	path=np.genfromtxt(pathlist, dtype='|S')
	NBH=[]; NP=[]; NMSP=[]; NMT=[]; NTOT=[]; NMSPSIN=[]; NMSPBIN=[]; NNS=[]; NNSDYN=[]; NPSIN=[]; NPBIN=[]
	for k in range(start, end):
		print path[k]
		snaps=np.sort(glob(path[k]+'/'+'initial.snap*.dat.gz'))
		lastsnap=snaps[-1]
		pref='initial'
		filestr=path[k]+'/'+pref
		Nbh, Ntot=dyn.find_NBH_NTOT_last(filestr)
		Npulsin, Npulbin, Nmsp, Nmt, MSPid, MSPcomid, Nns, Nnsdyn=get_snap_Nns(lastsnap,dynno[k])

		Nmspsin=0; Nmspbin=0
		for i in range(len(MSPid)):
			if int(MSPcomid[i])==-100: Nmspsin+=1
			else: Nmspbin+=1

		NBH.append(Nbh); NP.append(Npulsin+Npulbin); NMSP.append(Nmsp); NMT.append(Nmt); NTOT.append(Ntot); NPSIN.append(Npulsin); NPBIN.append(Npulbin)
		NNS.append(Nns); NMSPSIN.append(Nmspsin); NMSPBIN.append(Nmspbin); NNSDYN.append(Nnsdyn)

	
	NBH=np.array(NBH); NTOT=np.array(NTOT); NP=np.array(NP); NMSP=np.array(NMSP); NMT=np.array(NMT); NNS=np.array(NNS); NNSDYN=np.array(NNSDYN); NMSPSIN=np.array(NMSPSIN); NMSPBIN=np.array(NMSPBIN); NPSIN=np.array(NPSIN); NPBIN=np.array(NPBIN)

	print NBH, NTOT, NMT, NNS, NNSDYN, NP, NPSIN, NPBIN, NMSP, NMSPSIN, NMSPBIN 
	np.savetxt('/projects/b1011/syr904/projects/PULSAR/kickgrid_runs/paper/data/ns_number_newmodel.dat', np.c_[NBH, NTOT, NMT, NNS, NNSDYN, NP, NPSIN, NPBIN, NMSP, NMSPSIN, NMSPBIN], fmt ='%d %d %d %d %d %d %d %d %d %d %d', delimiter= ' ', header = '1.NBH 2.NTOT 3.NMT 4.NNS 5.NNSDYN 6.NP 7.sNP 8.bNP 9.NMSP 10.sNMSP 11.bNMSP', comments = '#')



def printout_Nbh_Npulsar_10to12Gyr(start, end, pathlist):
	dynno=[800000,800000,800000,800000,800000,800000,800000,800000,800000,800000,800000,800000,800000,800000,800000,800000,800000,800000,800000,800000,800000,800000,800000,800000,800000,3500000]
	path=np.genfromtxt(pathlist, dtype='|S')
	NBH=[]; NP=[]; NMSP=[]; NMT=[]; NTOT=[]; NMSPSIN=[]; NMSPBIN=[]; NNS=[]; NNSDYN=[]; M=[]; model=[]; T=[]
	f=open('/projects/b1011/syr904/projects/PULSAR/kickgrid_runs/paper/data/ns_number_9to12Gyr_newmodel.dat', 'a+', 0)
	f.write('#Model Age(Gyr) Mtot(Msun) NBH NTOT NMT NNS NNSDYN NP NMSP NMSPSIN NMSPBIN\n')
	for k in range(start, end):
		snaps=np.sort(glob(path[k]+'/'+'initial.snap*.dat.gz'))
		pref='initial'
		filestr=path[k]+'/'+pref
		t_conv=conv('t',filestr+'.conv.sh')
		m_conv=conv('m',filestr+'.conv.sh')
		if k==25: m_conv=2.166e+06
		for j in range(len(snaps)-1, 0, -1):
			t=get_time(snaps[j])
			Time=t_conv*t/1000.
			if Time>=9.0:
				Nbh, Ntot, ma=find_NBH_NTOT(filestr, t)
				mass=ma*m_conv
				print mass
				Npuls, Nmsp, Nmt, MSPid, MSPcomid, Nns, Nnsdyn=get_snap_Nns(snaps[j], dynno[k])
				NBH.append(Nbh); NP.append(Npuls); NMSP.append(Nmsp); NMT.append(Nmt); NTOT.append(Ntot); M.append(mass)
				model.append(k); T.append(Time)
				Nmspsin=0; Nmspbin=0
				for i in range(len(MSPid)):
					if int(MSPcomid[i])==-100: Nmspsin+=1
					else: Nmspbin+=1

				f.write('%d %f %e %d %d %d %d %d %d %d %d %d\n'%(k, Time, mass, Nbh, Ntot, Nmt, Nns, Nnsdyn, Npuls, Nmsp, Nmspsin, Nmspbin))

			if Time<9.0: break

		print k

	f.close()
	
	NBH=np.array(NBH); NTOT=np.array(NTOT); NP=np.array(NP); NMSP=np.array(NMSP); NMT=np.array(NMT)
	M=np.array(M); model=np.array(model); T=np.array(T)
	print NBH, NTOT, NMSP, NMT
	#np.savetxt('/projects/b1011/syr904/projects/PULSAR/kickgrid_runs/ns_number_10to12Gyr.dat', np.c_[model, T, NBH, NTOT, NMSP, NMT, M], fmt ='%d %f %d %d %d %d %f', delimiter= ' ', header = 'Model Age(Gyr) Nbh Ntot Nmsp Nnsmt M(msun)', comments = '#')



def plot_Nns_Nbh(pathlist, start, end):
        sourcedir=np.genfromtxt(pathlist, dtype='|S')
	NBH=[]; NMT=[]
	for k in range(start, end):
		filepath=sourcedir[k]
		pref='initial'
		filestr=filepath+'/'+pref
		snaps=np.sort(glob(filestr+'.snap*.dat.gz'))
                lastsnap=snaps[-1]

                databh=np.genfromtxt(filestr+'.bh.dat')
		Nbh=databh[-1,2]
		datans=np.genfromtxt(filestr+'.ns.dat')
		Nns_mtb=datans[-1, 4]
		NBH.append(float(Nbh)); NMT.append(float(Nns_mtb))
		print k
	
	plt.figure()
        plt.scatter(NBH, NMT, s=30, color='b')
        plt.xscale('log')
	plt.xlabel(r'$N_{BH}$')
	plt.ylabel(r'$N_{NS-MTB}$')
	#plt.xlim(xmin=0.1)
        plt.legend()
        #plt.show()
	plt.savefig('/projects/b1011/syr904/projects/nns-nbh.pdf', dpi=300)



def plot_t_Nns():
        sourcedir=['/projects/b1011/sourav/new_runs/kick_grid/rv1/kickscale_0.4', '/projects/b1011/sourav/new_runs/kick_grid/rv1/kickscale_0.04']
	titles=['Core-Collapse', 'Non Core-Collapse']
	f, axarr=plt.subplots(2)
	for i in range(2):
		filepath=sourcedir[i]
		pref='initial'
		filestr=filepath+'/'+pref
		datans=np.genfromtxt(filestr+'.ns.dat')
		Time=datans[:,0]; N_P=datans[:,5]; N_NS=datans[:,1]; N_MSP=datans[:,6]; N_NS_MTB=datans[:,4]
		t_conv=conv('t',filestr+'.conv.sh')
		Tns=[j * t_conv *0.001 for j in Time]

		databh=np.genfromtxt(filestr+'.bh.dat')
		Timebh=databh[:,1]; N_BH=databh[:,2]
		Tbh=[j * t_conv *0.001 for j in Timebh]

		axarr[i].scatter(Tns, N_NS_MTB, label=r'$N_{NS,MTB}$', color='b')
		axarr[i].scatter(Tbh, N_BH, label=r'$N_{BH}$', color='orange', s=5)
		axarr[i].set_yscale('log')
		#axarr[i].set_xscale('log')
		axarr[i].set_xlabel(r'$time(Gyr)$')
		axarr[i].set_xlim(-0.5, 13.)
		axarr[i].set_ylim(ymin=0.01)
		axarr[i].set_title(titles[i])
		axarr[i].legend(loc='best')

		print i

	plt.tight_layout()
	#plt.show()
	plt.savefig('/projects/b1011/syr904/projects/t-nmtb-nbh.pdf', dpi=300)


def plot_t_Nns_allmodel(pathlist, start, end):
	sourcedir=np.genfromtxt(pathlist, dtype='|S')
	plt.figure()
	for i in range(start, end):
		pref='initial'
		filestr=sourcedir[i]+'/'+pref
		databh=np.genfromtxt(filestr+'.bh.dat')
		Timebh=databh[:,1]; N_BH=databh[:,2]
		#Tbh=[j * t_conv *0.001 for j in Timebh]

		datans=np.genfromtxt(filestr+'.ns.dat')
		Time=datans[:,0]; N_P=datans[:,5]; N_NS=datans[:,1]; N_MSP=datans[:,6]; N_NS_MTB=datans[:,4]
		#print N_NS_MTB
		t_conv=conv('t',sourcedir[i]+'/'+'initial.conv.sh')
		Tns=[j * t_conv *0.001 for j in Time]

		if N_BH[-1]<=50:
			plt.scatter(Tns, N_NS_MTB, color='red')
		if N_BH[-1]>50:
			plt.scatter(Tns, N_NS_MTB, color='b')

	#blue_line = mlines.Line2D([], [], color='blue', label=r'$N_{BH} \gtrsim 200}$')
	#red_line = mlines.Line2D([], [], color='red', label = r'$N_{BH} \lesssim 10}$')
	red_dot, = plt.plot(z, "ro", markersize=8, label=r'$N_{BH} \lesssim 10}$')
	blue_dot, = plt.plot(z, "bo", markersize=8, label=r'$N_{BH} \gtrsim 200}$')
	plt.legend(handles=[red_line, blue_line])
	#plt.yscale('log')
	plt.xlabel(r'$time(Gyr)$')
	plt.xlim(-0.5, 13.)
	plt.legend()
	plt.show()




def print_Nns_snap(pathlist, start, end):
	sourcedir=np.genfromtxt(pathlist, dtype='|S')
	#sourcedir=['/projects/b1011/syr904/cmc/cmc-mpi-08/rundir_1-3/8e5rv1fb5kick1.0']
	
	for i in range(start, end):
		pref='initial'
		filestr=sourcedir[i]+'/'+pref
		snaps=np.sort(glob(filestr+'.snap*.dat.gz'))
		fhandle=open(filestr+'.ns.dat', 'w+', 0)
		fhandle.write('#1.Totaltime, 2.Nns,tot, 3.Nns,single, 4.Nns,binary, 5.Nns,mtb, 6.Npulsar, 7.Nmsp, 8.Nns-ns, 9.Nns-bh, 10.Nns-wd, 11.Nns-ms, 12.Nns-postms\n')
		for j in range(len(snaps)):
			N_NS=0; N_NS_SIN=0; N_NS_BIN=0; N_NS_MTB=0; N_PULS=0; N_MSP=0; N_NSNS=0; N_NSBH=0; N_NSWD=0; N_NSMS=0; N_NSPOSTMS=0
			T=get_time(snaps[j])
			print j
			with gzip.open(snaps[j], 'r') as fsnap:
				for _ in xrange(2):
					next(fsnap)
				for line in fsnap:
					datasnap=line.split()
					if int(datasnap[7])!=1:
						if int(datasnap[14])==13: 
							N_NS+=1; N_NS_SIN+=1
							spin=twopi*yearsc/float(datasnap[59])
							deathcut=(spin**2)*(0.17*10**12)
							if deathcut<float(datasnap[60]): N_PULS+=1
							if spin<=0.03: N_MSP+=1
					if int(datasnap[7])==1:
						if int(datasnap[17])==13:
							N_NS+=1; N_NS_BIN+=1
							spin0=twopi*yearsc/float(datasnap[45])
							deathcut0=(spin0**2)*(0.17*10**12)
							if deathcut0<float(datasnap[47]): N_PULS+=1
							if float(datasnap[44])>=1: N_NS_MTB+=1	
							if spin0<=0.03: N_MSP+=1

							if int(datasnap[18])<2: N_NSMS+=1
							elif int(datasnap[18])>=10 and int(datasnap[18])<=12: N_NSWD+=1
							elif int(datasnap[18])==13: N_NSNS+=1
							elif int(datasnap[18])==14: N_NSBH+=1
							else: N_NSPOSTMS+=1

						if int(datasnap[18])==13 and int(datasnap[17])!=13:
							N_NS+=1; N_NS_BIN+=1
							spin1=twopi*yearsc/float(datasnap[46])
							deathcut1=(spin1**2)*(0.17*10**12)
							if deathcut1<float(datasnap[48]): N_PULS+=1
							if float(datasnap[43])>=1: N_NS_MTB+=1
							if spin1<=0.03: N_MSP+=1

							if int(datasnap[17])<2: N_NSMS+=1
							elif int(datasnap[17])>=10 and int(datasnap[17])<=12: N_NSWD+=1
							elif int(datasnap[17])==13: N_NSNS+=1
							elif int(datasnap[17])==14: N_NSBH+=1
							else: N_NSPOSTMS+=1
			fhandle.write('%f %d %d %d %d %d %d %d %d %d %d %d\n'%(T, N_NS, N_NS_SIN, N_NS_BIN, N_NS_MTB, N_PULS, N_MSP, N_NSNS, N_NSBH, N_NSWD, N_NSMS, N_NSPOSTMS))
		fhandle.close()
		
		print i


def print_NSMTB_snap(pathlist, start, end):
	sourcedir=np.genfromtxt(pathlist, dtype='|S')
	fhandle=open('/projects/b1011/syr904/projects/PULSAR/kickgrid_runs/'+'NS_MTB_1e6_kickgrid.dat', 'w+', 0)
	fhandle.write('#1.Model, 2.Snapno, 3.Totaltime, 4.id0, 5.id1, 6.k0, 7.k1, 8.m0, 9.m1, 10.a, 11.ecc, 12.radrol0, 13.radrol1, 14.dmdt0, 15.dmdt1, 16.P0, 17.P1, 18.B0, 19.B1\n')
	for i in range(start, end):
		pref='initial'
		filestr=sourcedir[i]+'/'+pref
		snaps=np.sort(glob(filestr+'.snap*.dat.gz'))
		t_conv=conv('t',filestr+'.conv.sh')
		for j in range(len(snaps)):
			model=-100; snapno=-100; t=-100; id0=-100; id1=-100; k0=-100; k1=-100; m0=-100; m1=-100; a=-100; ecc=-100; radrol0=-100; radrol1=-100; dmdt0=-100; dmdt1=-100; p0=-100; p1=-100; b0=-100; b1=-100
			T=get_time(snaps[j])
			age=T*t_conv*0.001
			if age>=10.0:
				with gzip.open(snaps[j], 'r') as fsnap:
					for _ in xrange(2):
						next(fsnap)
					for line in fsnap:
						datasnap=line.split()
						if int(datasnap[7])==1:
							if (int(datasnap[17])==13 and float(datasnap[44])>=1) or (int(datasnap[18])==13 and float(datasnap[43])>=1):
								model=i; snapno=str(j).zfill(4); t=age
								id0=int(datasnap[10]); id1=int(datasnap[11]); k0=int(datasnap[17]); k1=int(datasnap[18]); m0=float(datasnap[8]); m1=float(datasnap[9]); a=float(datasnap[12]); ecc=float(datasnap[13]); radrol0=float(datasnap[43]); radrol1=float(datasnap[44]); dmdt0=float(datasnap[41]); dmdt1=float(datasnap[42]); p0=twopi*yearsc/float(datasnap[45]); p1=twopi*yearsc/float(datasnap[46]); b0=float(datasnap[47]); b1=float(datasnap[48]);
								fhandle.write('%d %s %f %d %d %d %d %f %f %f %f %f %f %f %f %f %f %e %e\n'%(model, snapno, t, id0, id1, k0, k1, m0, m1, a, ecc, radrol0, radrol1, dmdt0, dmdt1, p0, p1, b0, b1))

	fhandle.close()



def get_snap_BP(snapshot):
	Bs=[]; Bb=[]; Ps=[]; Pb=[]

	#snaps=np.sort(glob(sourcedir+'/'+'initial.snap*.dat.gz'))
	#lastsnap=snaps[-1]
	#print lastsnap
	with gzip.open(snapshot, 'r') as fsnap:
		for _ in xrange(2):
			next(fsnap)
		for line in fsnap:
			datasnap=line.split()
			if int(datasnap[7])!=1:
				if int(datasnap[14])==13:
					spin=twopi*yearsc/float(datasnap[59])
					Bs.append(float(datasnap[60])); Ps.append(spin)
			if int(datasnap[7])==1:
				if int(datasnap[17])==13:
					spin0=twopi*yearsc/float(datasnap[45])
					Bb.append(float(datasnap[47])); Pb.append(spin0)
				if int(datasnap[18])==13:
					spin1=twopi*yearsc/float(datasnap[46])
					Bb.append(float(datasnap[48])); Pb.append(spin1)

	return Bs, Bb, Ps, Pb


def sub_MSP(sourcedir, folder):   #Find sub-MSP in the model
	sunreal=[]; bunreal=[]; cpunreal=[]
	snaps=np.sort(glob(sourcedir+'/'+folder+'/'+'initial.snap*.dat.gz'))
	with open(sourcedir+'/'+folder+'/'+'end1e6.dat') as f:
		#next(f)
		for line in f:
			data=line.split()
			#print binflag
			if int(data[0])==834358:
				if int(data[8])==0:
					if float(data[6])<0.001:
						sunreal.append(int(data[2]))
				else:
					if int(data[19])==13:
						if float(data[17])<0.001:
							bunreal.append(int(data[9]))
							cpunreal.append(int(data[10]))
					if int(data[20])==13:
						if float(data[18])<0.001:
							bunreal.append(int(data[10]))
							cpunreal.append(int(data[9]))

	#print 'single=',sunreal, 'binary=', bunreal, 'companion=', cpunreal

	#snapno=int(len(snaps)-1)
	#ti=0.0
	#for k in range(len(bunreal)):
	#	id0=bunreal[k]; id1=cpunreal[k]
	#	history_maker_full4.history_maker(id0, id1, sourcedir+'/'+folder, snapno, ti)

	idlist=np.concatenate((sunreal,bunreal))
	for k in range(len(idlist)):
		sid=int(idlist[k])
		find_history(sid, sourcedir, folder)

	
		
def find_allNS(pathlist, start, end):
	sourcedir=np.genfromtxt(pathlist, dtype=str)
	#Nns=[]; Nnsesc=[]; Nnsprimd=[]; Nnsescprimd=[]
	fallns=open('/projects/b1011/syr904/projects/PULSAR/kickgrid_runs/paper/data/all_ns_new_26.dat', 'a+', 0)
	fallns.write('#1.Nns 2.Nns_primd 3.Nnsesc 4.Nnsesc_primd 5.Nnsccsn 6.Nnsecsn 7.Nnsothers\n')
	for i in range(start, end):
		Nns=0; Nnsesc=0; Nnsprimd=0; Nnsescprimd=0; Nnsccsn=0; Nnsecsn=0; Nnsothers=0
		filepath=sourcedir[i]
		filestr=filepath+'/initial'
		snaps=np.sort(glob(filestr+'.snap*.dat.gz'))
		nsid=[]; nsccsnid=[]; nsecsnid=[]; nsotherid=[]
		for j in range(len(snaps)):
			with gzip.open(snaps[j], 'r') as fsnap:
				for _ in xrange(2):
					next(fsnap)
				for line in fsnap:
					datasnap=line.split()
					if int(datasnap[14])==13: 
						nsid.append(int(datasnap[0]))
						if int(datasnap[61])==4: nsccsnid.append(int(datasnap[0]))
						elif int(datasnap[61])==5 or int(datasnap[61])==6 or int(datasnap[61])==7:
							nsecsnid.append(int(datasnap[0]))
						else: nsotherid.append(int(datasnap[0]))

					if int(datasnap[17])==13: 
						nsid.append(int(datasnap[10]))
						if int(datasnap[49])==4: nsccsnid.append(int(datasnap[10]))
						elif int(datasnap[49])==5 or int(datasnap[49])==6 or int(datasnap[49])==7:
							nsecsnid.append(int(datasnap[10]))
						else: nsotherid.append(int(datasnap[10]))

					if int(datasnap[18])==13: 
						nsid.append(int(datasnap[11]))
						if int(datasnap[50])==4: nsccsnid.append(int(datasnap[11]))
						elif int(datasnap[50])==5 or int(datasnap[50])==6 or int(datasnap[50])==7:
							nsecsnid.append(int(datasnap[11]))
						else: nsotherid.append(int(datasnap[11]))


		print len(nsid)

		#np.savetxt('/projects/b1011/syr904/projects/PULSAR/kickgrid_runs/primd_ns/new_intmd_kickgrid/primd_ns_'+str(i)+'.dat', np.c_[nsid], fmt='%d')

		NSID=Counter(nsid).keys()
		print len(NSID)
		n_primd=0; n_nonprimd=0
		for k in range(len(NSID)):
			if i<25:
				if NSID[k]<=800000: n_primd+=1
				else: n_nonprimd+=1
			else:
				if NSID[k]<=3500000: n_primd+=1
				else: n_nonprimd+=1

		NSCCSNID=Counter(nsccsnid).keys()
		nnsccsn=len(NSCCSNID)
		NSECSNID=Counter(nsecsnid).keys()
		nnsecsn=len(NSECSNID)
		NSOTHERID=Counter(nsotherid).keys()
		nnsothers=len(NSOTHERID)

		#print n_primd, n_nonprimd 

		nsescid=[]
		with open(filestr+'.esc.dat', 'r') as fesc:
			next(fesc)
			for line in fesc:
				dataesc=line.split()
				if int(dataesc[14])!=1:
					if int(dataesc[21])==13: nsescid.append(int(dataesc[13]))
				else:
					if int(dataesc[22])==13: nsescid.append(int(dataesc[17]))
					if int(dataesc[23])==13: nsescid.append(int(dataesc[18]))

		print len(nsescid)

		#np.savetxt('/projects/b1011/syr904/projects/PULSAR/kickgrid_runs/primd_ns/new_intmd_kickgrid/primd_nsesc_'+str(i)+'.dat', np.c_[nsescid], fmt='%d')

		NSESCID=Counter(nsescid).keys()
		print len(NSESCID)
		nesc_primd=0; nesc_nonprimd=0
		for o in range(len(NSESCID)):
			if i<25:
				if NSESCID[o]<=800000: nesc_primd+=1
				else: nesc_nonprimd+=1
			else:
				if NSESCID[o]<=3500000: nesc_primd+=1
				else: nesc_nonprimd+=1

		nesc_snap_primd=0; nesc_snap_nonprimd=0; nesc_ccsn=0; nesc_ecsn=0; nesc_other=0
		for p in range(len(NSESCID)):
			escid=int(NSESCID[p])
			for q in range(len(NSID)):
				if i<25:
					if int(NSID[q]==escid):
						if escid<=800000: nesc_snap_primd+=1
						else: nesc_snap_nonprimd+=1
				else:
					if int(NSID[q]==escid):
						if escid<=3500000: nesc_snap_primd+=1
						else: nesc_snap_nonprimd+=1

			for r in range(len(NSCCSNID)):
				if int(NSCCSNID[r]==escid):
					nesc_ccsn+=1

			for s in range(len(NSECSNID)):
				if int(NSECSNID[s]==escid):
					nesc_ecsn+=1

			for t in range(len(NSOTHERID)):
				if int(NSOTHERID[t]==escid):
					nesc_other+=1
			

		print nesc_snap_primd, nesc_snap_nonprimd


		#print nesc_primd, nesc_nonprimd

		Nns=n_primd+n_nonprimd-(nesc_snap_primd+nesc_snap_nonprimd); Nnsesc=nesc_primd+nesc_nonprimd
		Nnsprimd=n_primd-nesc_snap_primd; Nnsescprimd=nesc_primd
		Nnsccsn=nnsccsn-nesc_ccsn; Nnsecsn=nnsecsn-nesc_ecsn; Nnsothers=nnsothers-nesc_other

		fallns.write('%d %d %d %d %d %d %d\n'%(Nns, Nnsprimd, Nnsesc, Nnsescprimd, Nnsccsn, Nnsecsn, Nnsothers))

		print i

	#np.savetxt('/projects/b1011/syr904/projects/PULSAR/kickgrid_runs/paper/data/all_ns.dat', np.c_[Nns, Nnsprimd, Nnsesc, Nnsescprimd], fmt='%d %d %d %d', header='Nns Nns_primd Nnsesc Nnsesc_primd', comments='#')
	fallns.close()


def find_Nns_primd(sourcedir):
	nsfile=np.sort(glob(sourcedir+'/'+'primd_ns_*.dat'))
	nsescfile=np.sort(glob(sourcedir+'/'+'primd_nsesc_*.dat'))

	Nns=[]; Nnsesc=[]; Nnsprimd=[]; Nnsescprimd=[]
	for i in range(6):
		datans=np.genfromtxt(nsfile[i])
		datansesc=np.genfromtxt(nsescfile[i])

		NSID=Counter(datans).keys()
		#print len(NSID)
		n_primd=0; n_nonprimd=0
		for k in range(len(NSID)):
			if NSID[k]<=800000: n_primd+=1
			else: n_nonprimd+=1

		print n_primd, n_nonprimd 

		NSESCID=Counter(datansesc).keys()
		#print len(NSESCID)
		nesc_primd=0; nesc_nonprimd=0
		for o in range(len(NSESCID)):
			if NSESCID[o]<=800000: nesc_primd+=1
			else: nesc_nonprimd+=1

		print nesc_primd, nesc_nonprimd

		nskick=0; nskick_nonprimd=0
		for j in range(len(NSID)):
			if NSID[j] in NSESCID: 
				nskick+=1
				if NSID[j]>800000: nskick_nonprimd+=1


		print nskick, nskick_nonprimd

		Nns.append(n_primd+n_nonprimd-nskick); Nnsesc.append(nesc_primd+nesc_nonprimd)
		Nnsprimd.append(n_primd-(nskick-nskick_nonprimd)); Nnsescprimd.append(nesc_primd)

		#print i

	print Nns, Nnsprimd, Nnsesc, Nnsescprimd
	np.savetxt('/projects/b1011/syr904/projects/PULSAR/kickgrid_runs/primd_ns/new_intmd_kickgrid/primordial_ns.dat', np.c_[Nns, Nnsprimd, Nnsesc, Nnsescprimd], fmt='%d %d %d %d', header='#Nns Nnsprimd Nnsesc Nnsescprimd')


 
def find_primordialbin(pathlist, msplist):   #Find NS Primordial binaries
	datamsp=np.genfromtxt(msplist)
	id0msp=datamsp[:,1]; id1msp=datamsp[:,2]; modelno=datamsp[:,0]
	sourcedir=np.genfromtxt(pathlist, dtype=str)
	model_pri=[]; id0_pri=[]; id1_pri=[]
	for i in range(len(modelno)):
		no=int(modelno[i])
		filepath=sourcedir[no]
		filestr=filepath+'/'+'initial'
		t_conv = conv('t',filestr+'.conv.sh')

		data = np.genfromtxt(filestr+'.snap0000.dat.gz', usecols=(7,10,11))
		binflag=data[:,0]
		for j in range(len(binflag)):
			if binflag[j]==1:
				id0=data[:,1][j]; id1=data[:,2][j]
				if (id0msp[i]==id0 and id1msp[i]==id1) or (id1msp[i]==id0 and id0msp[i]==id1):
					model_pri.append(no); id0_pri.append(id0msp[i]); id1_pri.append(id1msp[i])
					break
		print i


	np.savetxt('/projects/b1011/syr904/projects/PULSAR/kickgrid_runs/paper/data/primordial_MSPbinary.dat', np.c_[model_pri, id0_pri, id1_pri], fmt='%d %d %d', delimiter=' ', comments = '#')
	print model_pri, id0_pri, id1_pri


	#snaps=np.sort(glob(sourcedir+'/'+folder+'/'+'initial.snap*.dat.gz'))
	#for j in range(len(id0)):
	#	x0=id0[j]; x1=id1[j]
	#	t=[]; m0=[]; m1=[]; i0=[]; i1=[]; a=[]; e=[]; k0=[]; k1=[]; P0=[]; P1=[]; B0=[]; B1=[]
	#	for k in range(len(snaps)):
	#		mass0=[]; mass1=[]; sid0=[]; sid1=[]; semima=[]; ecc=[]; sk0=[]; sk1=[]; spin0=[]; spin1=[]; #field0=[]; field1=[]
	#		with gzip.open(snaps[k],'r') as fi:
	#			for _ in xrange(2):
	#				next(fi)
	#			for line in fi:
	#				data1=line.split()
	#				if int(data1[7])==1:
	#					mass0.append(float(data1[8])); mass1.append(float(data1[9])); sid0.append(int(data1[#10])); sid1.append(int(data1[11])); semima.append(float(data1[12])); ecc.append(#float(data1[13])); sk0.append(int(data1[17])); sk1.append(int(data1[18])); spin0#.append(float(data1[46])); spin1.append(float(data1[47])); field0.append(float(#data1[48])); field1.append(float(data1[49]))
	#		
	#		time=t_conv*get_time(snaps[k])
	#		for l in range(len(mass0)):
	#			if (sid0[l]==x0 and sid1[l]==x1) or (sid0[l]==x1 and sid1[l]==x0):
	#				m0.append(mass0[l]); m1.append(mass1[l]); i0.append(sid0[l]); i1.append(sid1[l]); a.#append(semima[l]); e.append(ecc[l]); k0.append(sk0[l]); k1.append(sk1[l]); P0.append#(spin0[l]); P1.append(spin1[l]); B0.append(field0[l]); B1.append(field1[l]); t.#append(time)
	#				break
#
	#	y0=k0[-1]; y1=k1[-1]
	#	if y0==13 or y1==13:
	#		name=str(x0)+'_'+str(x1)
	#		np.savetxt(sourcedir+'/'+folder+'/'+'history'+'/'+name+'_primordial.dat', np.c_[t, m0, m1, i0, #i1, a, e, k0, k1, P0, P1, B0, B1], fmt ='%f %f %f %d %d %f %f %d %d %f %f %e %e', delimiter=# ' ', header = '1.time, 2.m0, 3.m1, 4.id0, 5.id1, 6.a, 7.e, 8.k0, 9.k1, 10.P0, 11.P1, #12.B0, 13.B1', comments = '#')
	#		#print x0, x1
#



def NSnumber(sourcedir, folder):
	Nret=[]; tret=[]; Nesc=[]; tesc=[]
	snaps=np.sort(glob(sourcedir+'/'+folder+'/'+'initial.snap*.dat.gz'))
	new_snaps=np.delete(snaps, 0)
	t_conv = conv('t',sourcedir+'/'+folder+'/'+'initial.conv.sh')
	for k in range(len(new_snaps)):
		time=t_conv*get_time(new_snaps[k])
		numret=0
		tret.append(time)
		with gzip.open(new_snaps[k],'r') as fi:
					for _ in xrange(2):
						next(fi)
					for line in fi:
						data=line.split()
						if int(data[7])==1:
							if int(data[17])==13:numret+=1	
							if int(data[18])==13:numret+=1
						if int(data[7])==0:
							if int(data[14])==13:numret+=1
								

		Nret.append(numret)

	#tret[0]+=1


	data1=np.genfromtxt(sourcedir+'/'+folder+'/'+'initial.esc.dat')
	tesc_temp=data1[:,1]; binflag=data1[:,14]; k=data1[:,21]; k0=data1[:,22]; k1=data1[:,23]
	for j in range(len(tesc_temp)):
		#numesc=0
		if binflag[j]==1:
			if k0[j]==13 or k1[j]==13: tesc.append(tesc_temp[j]*t_conv)#; numesc+=1 
		if binflag[j]==0:
			if k[j]==13: tesc.append(tesc_temp[j]*t_conv)#; numesc+=1

		#Nesc.append(numesc)

	tesc_key=Counter(tesc).keys()
	Nesc=Counter(tesc).values()
	tesc_key, Nesc=zip(*sorted(zip(tesc_key,Nesc)))

	#Nret=np.cumsum(Nret)
	#Nesc=np.cumsum(Nesc)
	#max=max(Nret)+10
	plt.figure()
	plt.xscale('log')
	plt.yscale('log')
	#plt.ylim((0.1, 6000))
	plt.scatter(tret, Nret, color='b', label='retained')
	plt.scatter(tesc_key, Nesc, color='r', alpha=0.7, label='escaped')
	plt.xlabel('t(Myr)')
	plt.ylabel('Number of NSs')
	plt.legend(loc='upper left')
	plt.savefig('NS'+'/'+'Nns_time.pdf')
		



def find_pulsar(sourcedir):
	sispin=[]; bispin=[]; Bs=[]; Bb=[]; idns=[]; idcom=[]; B=[]
	with open(sourcedir+'/'+'end1e6.dat') as f:
		#next(f)
		for line in f:
			data=line.split()
			#print binflag
			if int(data[0])==135549:
				if int(data[8])==0:
					if float(data[6])>0.05:
						#x1=float(data[6])*0.15*10**12
						#x2=float(data[6])*0.09*10**12
						y=float(data[6])*0.05*10**12   #High B
						#if float(data[5])>x2 and float(data[5])<x1:
						if float(data[5])>y:
					#B.append(float(data[5]))
							idns.append(int(data[2]))
							idcom.append(int(-1))
				else:
					if int(data[19])==13:
						if float(data[17])>0.05:
							#x1=float(data[17])*0.15*10**12
							#x2=float(data[17])*0.09*10**12
							y=float(data[17])*0.05*10**12
							#if float(data[15])>x2 and float(data[15])<x1:
							if float(data[15])>y:
						#B.append(float(data[15]))
								idns.append(int(data[9]))
								idcom.append(int(data[10]))
					if int(data[20])==13:
						if float(data[18])>0.05:
							#x1=float(data[18])*0.15*10**12
							#x2=float(data[18])*0.09*10**12
							y=float(data[18])*0.05*10**12
							#if float(data[16])>x2 and float(data[16])<x1:
							if float(data[16])>y:
						#B.append(float(data[16]))
								idns.append(int(data[10]))
								idcom.append(int(data[9]))

	
	#print len(idns)
	return idns, idcom
	#for i in range(len(idns)):	
		#history_maker_full5.history_maker(idns[i], idcom[i], sourcedir, 'initial', 149, 0.0)
	#find_history(sourcedir, idns)



def highB(sourcedir):
	idns, idcom=find_pulsar(sourcedir)
	indx=[]
	indx.append(idns.index(-6196740879348)); indx.append(idns.index(-6203282366687)); indx.append(idns.index(0))
	idns=np.delete(idns, indx); idcom=np.delete(idcom, indx)
	#print idns
	tis=[]; tib=[]; tfs=[]; tfb=[]; idis=[]; idib=[]; mis=[]; mib=[]; idfs=[]; idfb=[]; mfs=[]; mfb=[]
	for k in range(len(idns)):
		data0=np.genfromtxt(sourcedir+'/'+'history/snaphistory/'+str(idns[k])+'_snap.dat')
		first=data0[0]
		if idcom[k]==-1:
			if first[1]==idns[k]:
				tis.append(first[0]); mis.append(first[3])
			if first[2]==idns[k]:
				tis.append(first[0]); mis.append(first[4])
		else:
			if first[1]==idns[k]:
				tib.append(first[0]); mib.append(first[3])
			if first[2]==idns[k]:
				tib.append(first[0]); mib.append(first[4])

	tis=np.asarray(tis)
	tib=np.asarray(tib)
	

	Bfs=[]; Bfb=[]; Pfs=[]; Pfb=[]
	with open(sourcedir+'/'+'end1e6.dat') as f:
		for line in f:
			data=line.split()
			if float(data[0])==834358:
				for j in range(len(idns)):
					no=idns[j]
					if int(data[8])==0:
						if int(data[2])==no:
							idfs.append(no)
							Pfs.append(float(data[6]))
							Bfs.append(float(data[5]))
					if int(data[8])==1:
						if int(data[9])==no:
							idfs.append(no)
							Pfb.append(float(data[17]))
							Bfb.append(float(data[15]))
						if int(data[10])==no:
							idfs.append(no)
							Pfb.append(float(data[18]))
							Bfb.append(float(data[16]))

	#print idfs

	plt.figure()
	plt.xscale('log')
	plt.yscale('log')
	plt.scatter(Pfs, Bfs, s=35, marker='o', c=mis, label='single')
	plt.scatter(Pfb, Bfb, s=35, marker='^', c=mib, label='binary')
	plt.colorbar()
	#plt.xlim(10**-4, 100.)
	#plt.ylim(10**7, 10**15)
	plt.xlabel(r'$P(sec)$')
	plt.ylabel(r'$B(G)$')
	plt.legend(loc='upper left')
	plt.savefig('NS'+'/'+'highB_m.pdf')



def straightline(sourcedir):
	snaps=np.sort(glob(sourcedir+'/'+'initial.snap*.dat.gz'))
	sispin=[]; bispin=[]; Bs=[]; Bb=[]; idns=[]; idcom=[]; B=[]
	ns=0; nb=0
	with open(sourcedir+'/'+'end1e6.dat') as f:
		for line in f:
			data=line.split()
			if int(data[0])==834358:
				if int(data[8])==0:
					if float(data[6])>0.04:
						x1=float(data[6])*0.15*10**12
						x2=float(data[6])*0.09*10**12
						#y=float(data[6])*0.12*10**12
						if float(data[5])>x2 and float(data[5])<x1:
						#if float(data[5])>y:
					#B.append(float(data[5]))
							idns.append(int(data[2]))
							idcom.append(1.0)
							ns+=1
				else:
					if int(data[19])==13:
						if float(data[17])>0.04:
							x1=float(data[17])*0.15*10**12
							x2=float(data[17])*0.09*10**12
							#y=float(data[17])*0.12*10**12
							if float(data[15])>x2 and float(data[15])<x1:
							#if float(data[15])>y:
						#B.append(float(data[15]))
								idns.append(int(data[9]))
								idcom.append(int(data[10]))
								nb+=1
					if int(data[20])==13:
						if float(data[18])>0.04:
							x1=float(data[18])*0.15*10**12
							x2=float(data[18])*0.09*10**12
							#y=float(data[18])*0.12*10**12
							if float(data[16])>x2 and float(data[16])<x1:
							#if float(data[16])>y:
						#B.append(float(data[16]))
								idns.append(int(data[10]))
								idcom.append(int(data[9]))
								nb+=1


	datai=np.genfromtxt(snaps[1])
	binflagi=datai[:,7]
	dataf=np.genfromtxt(snaps[-1])
	binflagf=dataf[:,7]
	mt=[]; c=0
	for i in range(len(idns)):
		no=idns[i]
		for j in range(len(binflagi)):
			if binflagi[j]==0:
				if int(datai[:,0][j])==no: mi=datai[:,1][j]; c+=1
			if binflagi[j]==nap:
				if int(datai[:,10][j])==no: mi=datai[:,8][j]; c+=1
				if int(datai[:,11][j])==no: mi=datai[:,9][j]; c+=1
		if c==0:
			for h in range(len(snaps)-2):
				dataext=np.genfromtxt(snaps[h+2])
				binflagext=dataext[:,7]
				for g in range(len(binflagfext)):
					if binflagext[g]==0:
						if int(dataext[:,0][g])==no: mi=dataext[:,1][g]; c+=1
					if binflagext[g]==1:
						if int(dataext[:,10][g])==no: mi=dataext[:,8][g]; c+=1
						if int(dataext[:,11][g])==no: mi=dataext[:,9][g]; c+=1

				if c==1: break


		for k in range(len(binflagf)):
			if binflagf[k]==0:
				if int(dataf[:,0][k])==no: mf=dataf[:,1][k]; break
			if binflagf[k]==1:
				if int(dataf[:,10][k])==no: mf=dataf[:,8][k]; break
				if int(dataf[:,11][k])==no: mf=dataf[:,9][k]; break


		if mi==mf: mt.append(int(0))
		if mi<mf:
			a=0
			for l in range(len(snaps)-1):
				data1=np.genfromtxt(snaps[l])
				data2=np.genfromtxt(snaps[l+1])
				bf1=data1[:,7]; bf2=data2[:,7]
				for m in range(len(bf1)):
					if bf1[m]==0:
						if int(data1[:,0][m])==no: m1=data1[:,1][m]
					if bf1[m]==1:
						if int(data1[:,10][m])==no: m1=data1[:,8][m]
						if int(data1[:,11][m])==no: m1=data1[:,9][m]

				for n in range(len(bf2)):
					if bf2[n]==0:
						if int(data2[:,0][n])==no: m2=data2[:,1][n]
					if bf2[n]==1:
						if int(data2[:,10][n])==no: m2=data2[:,8][n]
						if int(data2[:,11][n])==no: m2=data1[:,9][n]

				if m2>m1: a+=1
			mt.append(a)


	print ns, nb
	print mt
	#return mt
	plt.figure()
	plt.hist(mt)
	plt.yscale('log')
	plt.title('Mass Transfering')
	plt.savefig('NS/'+'mt.pdf')



def find_merger_startype(ids ,sourcedir):
	typem=int(-1); type1=int(-1); type2=int(-1); type3=int(-1); type4=int(-1)
	idkm1=[]; idkm6=[]; idkm11=[]; idkm12=[]
	with open(sourcedir+'/'+'initial.collision.log') as fi:
		next(fi)
		for line in fi:
			idm=re.findall('idm=([\d]+)', line)[0]	
			if idm==str(ids):
				typem=re.findall('typem=([\d]+)', line)[0]
				type1=re.findall('type1=([\d]+)', line)[0]
				type2=re.findall('type2=([\d]+)', line)[0]
				if re.findall('type3', line):
					type3=re.findall('type3=([\d]+)', line)[0]
					if re.findall('type4', line):
						type4=re.findall('type4=([\d]+)', line)[0]

				if typem==1: idkm1.append(ids)
				if typem==6: idkm6.append(ids)
				if typem==11: idkm11.append(ids)
				if typem==12: idkm12.append(ids)

			 
		if typem==-1:
			print 'Cannot find ', ids

				
	return int(typem), int(type1), int(type2), int(type3), int(type4)
	print idkm1, idkm6, idkm11, idkm12


def get_nenc(pathlist, start, end):
	sourcedir=np.genfromtxt(pathlist, dtype=str)
	enc_bhpoor=[]; enc_bhrich=[]; enc_bhmid=[]
	if end==10: 
		f=open('/projects/b1011/syr904/projects/PULSAR/kickgrid_runs/paper/data/enc_bhrich_new.dat', 'a+', 0)
	elif end==19: f=open('/projects/b1011/syr904/projects/PULSAR/kickgrid_runs/paper/data/enc_bhmid_new.dat', 'a+', 0)
	else: f=open('/projects/b1011/syr904/projects/PULSAR/kickgrid_runs/paper/data/enc_bhpoor_new.dat', 'a+', 0)
	for i in range(start, end):
		pref='initial'
		filestr=sourcedir[i]+'/'+pref
		snaps=np.sort(glob(filestr+'.snap*.dat.gz'))
		lastsnap=snaps[-1]
		nsid, m=find_ns(lastsnap)
		print len(nsid)
    		
		for j in range(len(nsid)):
			#print MSPid[j]
			history=hic.history_maker([nsid[j]], [1], 'initial', sourcedir[i], 1.0)
			intact_num=len(history[nsid[j]]['binint']['binint'])
			if i<10: 
				enc_bhrich.append(intact_num)
				f.write('%d\n'%(intact_num))
			elif i<19: 
				enc_bhmid.append(intact_num)
				f.write('%d\n'%(intact_num))
			else: 
				enc_bhpoor.append(intact_num)
				f.write('%d\n'%(intact_num))

		print i


def get_nenc_msp(pathlist, mspfile):
	sourcedir=np.genfromtxt(pathlist, dtype=str)
	enc_msp=[]
	datamsp=np.genfromtxt(mspfile)
	modelno=datamsp[:,0]; id0=datamsp[:,1]; id1=datamsp[:,2]
	for i in range(len(modelno)):
		pref='initial'
		no=int(modelno[i])
		filestr=sourcedir[no]+'/'+pref
		snaps=np.sort(glob(filestr+'.snap*.dat.gz'))
		lastsnap=snaps[-1]
    		
		history=hic.history_maker([id0[i]], [1], 'initial', sourcedir[no], 1.0)
		intact_num=len(history[id0[i]]['binint']['binint'])
		enc_msp.append(intact_num)


		print i

	np.savetxt('/projects/b1011/syr904/projects/PULSAR/kickgrid_runs/paper/data/enc_msp.dat', np.c_[modelno, id0, id1, enc_msp], fmt='%d %d %d %d', header='Model ID0 ID1 Encounters', delimiter='', comments='#')

            


def get_mtb_numofencounter(sourcedir):
	data = np.genfromtxt(sourcedir)
	pathlist=np.genfromtxt('/projects/b1011/syr904/projects/PULSAR/kickgrid_runs/kickgrid_path.dat', dtype='|S')

	model=data[:,0]; totaltime=data[:,2]; id0=data[:,3]; id1=data[:,4]; k0=data[:,5]; k1=data[:,6]
	modelall_key=Counter(model).keys()
	modelall_value=Counter(model).values()
	model_key, model_value=zip(*sorted(zip(modelall_key,modelall_value)))

	model_value_cumsum=np.cumsum(model_value)
	model_key=np.array(model_key)


	interact_num_bhrich=[]
	interact_num_bhpoor=[]

	for i in range(len(model_key)):
		modelno=int(model_key[i])
		startlen=int(model_value_cumsum[i])-1
		modellen=int(model_value[i])
		lasttime=totaltime[startlen]
		print lasttime
		path=pathlist[modelno]

		j=startlen

		while totaltime[j]==lasttime:
			if k0[j]==13:
				print id0[j], data[:,1][j]	
				history=hic.history_maker([id0[j]], [1], 'initial', path, 1.0)
				intact_num=len(history[id0[j]]['binint']['binint'])
				if modelno>10: interact_num_bhpoor.append(intact_num)
				else: interact_num_bhrich.append(intact_num)
			if k1[j]==13:
				print id1[j], data[:,1][j]
				history=hic.history_maker([id1[j]], [1], 'initial', path, 1.0)
				intact_num=len(history[id1[j]]['binint']['binint'])
				if modelno>10: interact_num_bhpoor.append(intact_num)
				else: interact_num_bhrich.append(intact_num)

			j-=1

	return interact_num_bhrich, interact_num_bhpoor



def plot_mtb_encounter_cdf(sourcedir):
	enc_bhrich, enc_bhpoor=get_mtb_numofencounter(sourcedir)
	print enc_bhrich, enc_bhpoor
	#enc_bhrich=np.cumsum(enc_bhrich); enc_bhpoor=np.cumsum(enc_bhpoor)
	#print enc_bhrich, enc_bhpoor	

	plt.figure()
	plt.hist(enc_bhrich, 10, normed=1, histtype='step', alpha=0.7, label=r'$N_{BH} \gtrsim 200$')
	plt.hist(enc_bhpoor, 10, normed=1, histtype='step', alpha=0.7, label=r'$N_{BH} \lesssim 10$')
	plt.legend()
	plt.ylim(0, 0.4)
	plt.show()


def plot_encounter_cdf(pathlist, start, end):
	sourcedir=np.genfromtxt(pathlist, dtype=str)
	enc_bhpoor=[]; enc_bhrich=[]; enc_bhmid=[]
	for i in range(start, end):
		pref='initial'
		filestr=sourcedir[i]+'/'+pref
		snaps=np.sort(glob(filestr+'.snap*.dat.gz'))
		lastsnap=snaps[-1]
		Npulsar, Nmsp, Nmtb, MSPid=get_snap_Nns(lastsnap)
		for j in range(len(MSPid)):
			print MSPid[j]
			history=hic.history_maker([MSPid[j]], [1], 'initial', sourcedir[i], 1.0)
			intact_num=len(history[MSPid[j]]['binint']['binint'])
			if i>18: enc_bhpoor.append(intact_num)
			elif i<11: enc_bhrich.append(intact_num)
			else: enc_bhmid.append(intact_num)

		print i

	#print len(enc_bhrich), len(enc_bhpoor), len(enc_bhmid)
	enc_bhrich=np.pad(enc_bhrich, (0, 37), 'constant')
	enc_bhmid=np.pad(enc_bhmid, (0, 31), 'constant')
	#print len(enc_bhrich), len(enc_bhpoor), len(enc_bhmid)

	weights_bhrich= np.ones_like(enc_bhrich)/float(len(enc_bhrich))
	weights_bhpoor= np.ones_like(enc_bhpoor)/float(len(enc_bhpoor))
	weights_bhmid= np.ones_like(enc_bhmid)/float(len(enc_bhmid))


	plt.figure()
	plt.hist(enc_bhrich, bins=10, histtype='step', lw=1.8, weights=weights_bhrich, cumulative=True, range=(0.0, max(enc_bhpoor)), label=r'$N_{BH} \gtrsim 200$')
	plt.hist(enc_bhmid, bins=10, histtype='step', lw=1.8, weights=weights_bhmid, cumulative=True, range=(0.0, max(enc_bhpoor)), label=r'$ 10 \lesssim N_{BH} \lesssim 200$')
	plt.hist(enc_bhpoor, bins=10, histtype='step', lw=1.8, weights=weights_bhpoor, cumulative=True, range=(0.0, max(enc_bhpoor)), label=r'$N_{BH} \lesssim 10$')
	plt.ylim(0.0, 1.2)
	plt.xlim(-1., 29.)
	plt.ylabel(r'$Probability\ Density$')
	plt.xlabel(r'$N_{encounter}$')
	plt.legend(loc='lower right')
	#plt.show()
	plt.savefig('/projects/b1011/syr904/projects/PULSAR/kickgrid_runs/encounter_cdf.pdf', dpi=300)

	#plt.figure()
	#values_rich, base_rich = np.histogram(enc_bhrich)
	#cumulative_rich = np.cumsum(values_rich)
	#plt.plot(base_rich[:-1], cumulative_rich)
	#values_mid, base_mid = np.histogram(enc_bhmid)
	#cumulative_mid = np.cumsum(values_mid)
	#plt.plot(base_mid[:-1], cumulative_mid)
	#values_poor, base_poor = np.histogram(enc_bhpoor)
	#cumulative_poor = np.cumsum(values_poor)
	#plt.plot(base_poor[:-1], cumulative_poor)
	#plt.ylim(0.0, 1.2)
	#plt.ylabel(r'$Probability Density$')
	#plt.xlabel(r'$N_{encounter}$')
	#plt.legend(loc='upper right')
	##plt.show()
	#plt.savefig('/projects/b1011/syr904/projects/PULSAR/kickgrid_runs/encounter_cdf_alter.pdf', dpi=300)



##Find the fraction of NSs in different formation channels
def NS_formation_fraction(pathlist, start, end):
	sourcedir=np.genfromtxt(pathlist, dtype=str)
	f=open('/projects/b1011/syr904/projects/PULSAR/kickgrid_runs/paper/data/NS_formation_fraction.dat', 'a+', 0)
	f.write('#1.CC, 2.EIC, 3.AIC, 4.MIC, 5.Others\n')
	for i in range(len(sourcedir)):
		N0=0; N4=0; N5=0; N6=0; N7=0
		filestr=sourcedir[i]+'/initial'
		snaps=np.sort(glob(filestr+'.snap*.dat.gz'))
		lastsnap=snaps[-1]

		with gzip.open(lastsnap, 'r') as flast:
			for _ in xrange(2): next(flast)
			for line in flast:
				datalast=line.split()
				if int(datalast[7])!=1:
					if int(datalast[14])==13:
						if int(datalast[61])==0: N0+=1
						if int(datalast[61])==4: N4+=1
						if int(datalast[61])==5: N5+=1
						if int(datalast[61])==6: N6+=1
						if int(datalast[61])==7: N7+=1
				else:
					if int(datalast[17])==13:
						if int(datalast[49])==0: N0+=1
						if int(datalast[49])==4: N4+=1
						if int(datalast[49])==5: N5+=1
						if int(datalast[49])==6: N6+=1
						if int(datalast[49])==7: N7+=1
					if int(datalast[18])==13:
						if int(datalast[50])==0: N0+=1
						if int(datalast[50])==4: N4+=1
						if int(datalast[50])==5: N5+=1
						if int(datalast[50])==6: N6+=1
						if int(datalast[50])==7: N7+=1


		f.write('%d %d %d %d %d %d\n'%(i, N4, N5, N6, N7, N0))

		print i

	f.close()



def get_formationchannel(snapshot):
	fc_si=[]; fc_bi=[]
	with gzip.open(snapshot, 'r') as fsnap:
                for _ in xrange(2):
                	next(fsnap)
                for line in fsnap:
                        datasnap=line.split()
                        if int(datasnap[7])!=1:
                                if int(datasnap[14])==13:
                                        fc_si.append(int(datasnap[61]))
                        if int(datasnap[7])==1:
                                if int(datasnap[17])==13:
                                        fc_bi.append(int(datasnap[49]))
                                if int(datasnap[18])==13:
                                        fc_bi.append(int(datasnap[50]))

        return fc_si, fc_bi


def plot_formationchannel(sourcedir):
	pref='initial'
	filestr=sourcedir+'/'+pref
	snaps=np.sort(glob(filestr+'.snap*.dat.gz'))
	lastsnap=snaps[-1]
	FCsi, FCbi=get_formationchannel(lastsnap)
	FCbix=Counter(FCbi).keys(); FCbiy = Counter(FCbi).values()	
	FCsix=Counter(FCsi).keys(); FCsiy = Counter(FCsi).values()
	#print FCbix, FCsix	
		
	FCbiy.insert(1, 0); FCbiy.insert(4, 0); FCsiy.insert(4,0) 
	print FCbiy, FCsiy
	
	x=[1,2,3,4,5]
	xlabel=['Others', 'Type 2', 'ECSN', 'AIC', 'MIC']
	
	plt.figure()
	plt.scatter(x,FCsiy, c='purple', s=200, label='Single')
	plt.scatter(x, FCbiy, c='orange', s=200, marker='^', label='Binary')
	plt.xticks(x, xlabel)
	plt.yscale('log')
	plt.ylim(1, 1000)
	plt.legend(loc='best')
	plt.show()

	

def get_formation_BP(sourcedir):
	pref='initial'
	filestr=sourcedir+'/'+pref
	snaps=np.sort(glob(filestr+'.snap*.dat.gz'))
	ids=[]; idb=[]; Bs=[]; Bb=[]; Ps=[]; Pb=[]; FCs=[]; FCb=[]
	with gzip.open(snaps[2], 'r') as fsnap:
		for _ in xrange(2):
			next(fsnap)
		for line in fsnap:
			datasnap=line.split()
			#if int(datasnap[7])!=1:
				#if int(datasnap[14])==13:
					#ids.append(int(datasnap[0]))
					#Bs.append(float(datasnap[60])); Ps.append((twopi*yearsc)/float(datasnap[59])); FCs.append(int(datasnap[61]))
			if int(datasnap[7])==1:
				if int(datasnap[17])==13:
                                	idb.append(int(datasnap[10]))
                                        Bb.append(float(datasnap[47])); Pb.append((twopi*yearsc)/float(datasnap[45])); FCb.append(int(datasnap[49]))	
				if int(datasnap[18])==13:
                                	idb.append(int(datasnap[11]))
                                	Bb.append(float(datasnap[48])); Pb.append((twopi*yearsc)/float(datasnap[46])); FCb.append(int(datasnap[50]))
	
	for i in range(3, len(snaps)):
		idtot=ids+idb
		#print idtot
		with gzip.open(snaps[i], 'r') as fsnap:
                        for _ in xrange(2):
                                next(fsnap)
                        for line in fsnap:
				control=0
                                datasnap=line.split()
                                #if int(datasnap[7])!=1:
                                        #if int(datasnap[14])==13:
						#for j in range(len(idtot)):
							#if int(datasnap[0])==idtot[j]: control=1
						#if control==0:
                                                	#ids.append(int(datasnap[0]))
                                                	#Bs.append(float(datasnap[60])); Ps.append((twopi*yearsc)/float(datasnap[59])); FCs.append(int(datasnap[61]))
                                if int(datasnap[7])==1:
                                        if int(datasnap[17])==13:
						for j in range(len(idtot)):
                                                        if int(datasnap[10])==idtot[j]: control=1
						if control==0:
                                                	idb.append(int(datasnap[10]))
                                                	Bb.append(float(datasnap[47])); Pb.append((twopi*yearsc)/float(datasnap[45])); FCb.append(int(datasnap[49]))
                                        if int(datasnap[18])==13:
						for j in range(len(idtot)):
                                                       if int(datasnap[11])==idtot[j]: control=1
						if control==0: 
                                                	idb.append(int(datasnap[11]))
                                                	Bb.append(float(datasnap[48])); Pb.append((twopi*yearsc)/float(datasnap[46])); FCb.append(int(datasnap[50]))						
		
		print i, len(idtot)

	#print ids, idb, Bs, Bb, Ps, Pb, FCs, FCb
	idtot=ids+idb; B=Bs+Bb; P=Ps+Pb; FC=FCs+FCb
        np.savetxt('/projects/b1011/syr904/projects/PULSAR/bse_change/newformedNS_kconst.dat', np.c_[idtot, B, P, FC], fmt ='%ld %e %f %d', delimiter= ' ', header = 'id B P FC', comments = '#' )
	return ids, idb, Bs, Bb, Ps, Pb, FCs, FCb
	

def plot_formation_BP():
	#idsi, idbi, Bsi, Bbi, Psi, Pbi, FCsi, FCbi=get_formation_BP(sourcedir)
	dataBP=np.genfromtxt('/projects/b1011/syr904/projects/PULSAR/bse_change/newformedNS_kconst.dat')
	B=dataBP[:,1]; P=dataBP[:,2]; FC=dataBP[:,3]
	#print B, P
	
	FCmark=[]
	for i in range(len(FC)):
                if FC[i]==0: FCmark.append('D')
                if FC[i]==4: FCmark.append('o')
                if FC[i]==5: FCmark.append('s')
                if FC[i]==6: FCmark.append('^')
                if FC[i]==7: FCmark.append('*')
	
	B6=[]; P6=[]
	for j in range(len(B)):
		if FC[j]==6: B6.append(B[j]); P6.append(P[j])
				
	plt.figure()
	for k in range(len(P)):
		plt.scatter(P[k], B[k], marker=FCmark[k], alpha=0.5, s=15)
	#plt.scatter(P6, B6)
	plt.xlabel('P(sec)')
	plt.ylabel('B(G)')
	plt.xscale('log')
	plt.yscale('log')
	plt.xlim(0.0001, 1000)
	plt.ylim(10**7, 10**14)
	#plt.savefig('/projects/b1011/syr904/projects/PULSAR/bse_change/formation.pdf', dpi=300)	
	plt.show()



def get_id_BP(sourcedir, theid, modelno, savepath):
	pref='initial'
	filestr=sourcedir+'/'+pref
	snaps=np.sort(glob(filestr+'.snap*.dat.gz'))
	print len(snaps)
	t_conv = conv('t',sourcedir+'/'+'initial.conv.sh')
	Age=[]; B=[]; P=[]; FC=[]; m0=[]; m1=[]; k0=[]; k1=[]; a=[]; ecc=[]; id0=[]; id1=[]; radrol0=[]; radrol1=[]

	for i in range(len(snaps)):
		time=get_time(snaps[i])*t_conv
		with gzip.open(snaps[i], 'r') as fsnap:
			for _ in xrange(2):
				next(fsnap)
			for line in fsnap:
				datasnap=line.split()
				if int(datasnap[7])!=1:
					if int(datasnap[14])==13:
						if int(datasnap[0])==theid:
							Age.append(float(time)); B.append(float(datasnap[60])); P.append((twopi*yearsc)/float(datasnap[59])); FC.append(int(datasnap[61])); m0.append(float(datasnap[1])); m1.append(-100)
							k0.append(int(datasnap[14])); k1.append(-100); a.append(-100); ecc.append(-100)
							id0.append(int(datasnap[0])); id1.append(-100); radrol0.append(-100); radrol1.append(-100)
				if int(datasnap[7])==1:
					if int(datasnap[17])==13:
						if int(datasnap[10])==theid:
							Age.append(float(time)); B.append(float(datasnap[47])); P.append((twopi*yearsc)/float(datasnap[45])); FC.append(int(datasnap[49])); m0.append(float(datasnap[8])); m1.append(float(datasnap[9]))	
							k0.append(int(datasnap[17])); k1.append(int(datasnap[18])); a.append(float(datasnap[12])); ecc.append(float(datasnap[13]))
							id0.append(int(datasnap[10])); id1.append(int(datasnap[11])); radrol0.append(float(datasnap[43])); radrol1.append(float(datasnap[44]))
					if int(datasnap[18])==13:
						if int(datasnap[11])==theid:
							Age.append(float(time)); B.append(float(datasnap[48])); P.append((twopi*yearsc)/float(datasnap[46])); FC.append(int(datasnap[50])); m0.append(float(datasnap[9])); m1.append(float(datasnap[8]))
							k0.append(int(datasnap[18])); k1.append(int(datasnap[17])); a.append(float(datasnap[12])); ecc.append(float(datasnap[13]))
							id0.append(int(datasnap[11])); id1.append(int(datasnap[10])); radrol0.append(float(datasnap[44])); radrol1.append(float(datasnap[43]))
		print i
	np.savetxt(savepath+'/'+str(modelno)+'_'+str(theid)+'.dat', np.c_[Age, B, P, FC, id0, id1, m0, m1, k0, k1, radrol0, radrol1, a, ecc], fmt ='%f %e %f %d %d %d %f %f %d %d %f %f %f %f', delimiter= ' ', header = '1.Age(Myr) 2.B(G) 3.P(sec) 4.FC 5.id0 6.id1 7.m0 8.m1 9.k0 10.k1 11.radrol0 12.radrol1 13.a[AU] 14.ecc', comments = '#')
	#return B, P, FC, Age



def get_allid_BP(idfile, pathfile, savepath, start, end):
	dataid=np.genfromtxt(idfile, dtype=int)
	model=dataid[:,0]; id0=dataid[:,1]; id1=dataid[:,2]
	path=np.genfromtxt(pathfile, dtype=str)
	for k in range(0,len(model)):
		sourcedir=path[int(model[k])]
		print sourcedir
		if model[k]<end and model[k]>=start:
			get_id_BP(sourcedir, id0[k], model[k], savepath)
			print id0[k]

	print k
				


def get_binary_evolution(sourcedir, pathlist, savepath, num):
	dataid=np.genfromtxt(sourcedir)
	model=dataid[:,0]; id0=dataid[:,1]; id1=dataid[:,2]
	sourcedir=np.genfromtxt(pathlist, dtype=str)
	path=sourcedir[num]
	snaps=np.sort(glob(path+'/'+',snap*.dat.gz'))
	for i in range(len(id0)):
		if model[i]==float(num):
			if id1[i]!=-100:
				snapno=len(snaps)-1
				hi4.history_maker(int(id0[i]), int(id1[i]), path, snapno, 0.0, savepath)
				print id0[i]


def get_history_inpsrfile(sourcedir, theid):
	pref='initial'
	filestr=sourcedir+'/'+pref
	t_conv = conv('t',filestr+'.conv.sh')
	fname=filestr+'.morepulsars.dat'
	Age=[]; B=[]; P=[]; m0=[]; m1=[]; k0=[]; k1=[]; a=[]; ecc=[]
	id0=[]; id1=[]; radrol0=[]; radrol1=[]
	with open(fname, 'r') as fpsr:
		next(fpsr)
		for line in fpsr:
			datapsr=line.split()
			if int(datapsr[8])!=1:
				if int(datapsr[2])==theid:
					Age.append(float(datapsr[1])*t_conv); B.append(float(datapsr[5])); P.append(float(datapsr[6]))
					M.append(float(datapsr[3])); Mcom.append(-100); Kcom.append(-100)
					a.append(-100); ecc.append(-100); id0.append(int(datapsr[2]))
			if int(datapsr[8])==1:
				if int(datapsr[9])==theid:
					Age.append(float(datapsr[1])*t_conv); B.append(float(datapsr[15])); P.append(float(datapsr[17]))
					M.append(float(datapsr[11])); Mcom.append(float(datapsr[12])); Kcom.append(int(datapsr[20]))
					a.append(float(datapsr[21])); ecc.append(float(datapsr[22]))
				if int(datapsr[10])==theid:
					Age.append(float(datapsr[1])*t_conv); B.append(float(datapsr[16])); P.append(float(datapsr[18]))
					M.append(float(datapsr[12])); Mcom.append(float(datapsr[11])); Kcom.append(int(datapsr[19]))
					a.append(float(datapsr[21])); ecc.append(float(datapsr[22]))

	return Age, M, Mcom, B, P, Kcom, a, ecc


def get_NNS_inpsrfile(pathlist, start, end):   ##including retained and escaped pulsars
	sourcedir=np.genfromtxt(pathlist, dtype=str)
	NNS=[]; NPSR=[]; NMSP=[]; model=[]; time=[]
	fn=open('/projects/b1011/syr904/projects/PULSAR/kickgrid_runs/paper/data/nsnumber_newmodel_inpsr_9to12Gyr.dat', 'a+', 0)
	fn.write('#1.Model, 2.Time, 3.Nns, 4.Npsr, 5.Nmsp\n')
	for i in range(start, end):
		filestr=sourcedir[i]+'/'
		t_conv=conv('t',filestr+'initial.conv.sh')
		with open(filestr+'psrlast8e7.dat', 'r') as fpsr:
			for line in fpsr:
				datapsr=line.split()
				t=float(datapsr[1])
				tmyr=t*t_conv
				if tmyr>=9000.0:
					print tmyr, t
					tcheck=int(datapsr[0])
					print tcheck
					break

		nns=0; npsr=0; nmsp=0
		with open(filestr+'psrlast8e7.dat', 'r') as fpsr:
			for line in fpsr:
				datapsr=line.split()
				if str.isdigit(datapsr[0])==False: print 'no'; continue
				#print datapsr[0]
				if int(datapsr[0])==tcheck:
					#print int(datapsr[0]), float(datapsr[1])
					nns+=1
					if int(datapsr[19])==13:
						p=float(datapsr[17]); b=float(datapsr[15]) 
						if b>=(p**2)*(0.17*10**12): npsr+=1
						if p<=0.03: nmsp+=1
					if int(datapsr[20])==13:
						p=float(datapsr[18]); b=float(datapsr[16]) 
						if b>=(p**2)*(0.17*10**12): npsr+=1
						if p<=0.03: nmsp+=1
					if int(datapsr[7])==13:
						p=float(datapsr[6]); b=float(datapsr[5]) 
						if b>=(p**2)*(0.17*10**12): npsr+=1
						if p<=0.03: nmsp+=1

				if int(datapsr[0])>tcheck: 
					tcheck=int(datapsr[0])
					#print tcheck
					NNS.append(nns); NPSR.append(npsr); NMSP.append(nmsp)
					model.append(i); time.append(float(datapsr[1])*t_conv)
					tmyr=float(datapsr[1])*t_conv
					fn.write('%d %f %d %d %d\n'%(i, tmyr, nns, npsr, nmsp))
					nns=0; npsr=0; nmsp=0
					nns+=1
					if int(datapsr[19])==13:
						p=float(datapsr[17]); b=float(datapsr[15]) 
						if b>=(p**2)*(0.17*10**12): npsr+=1
						if p<=0.03: nmsp+=1
					if int(datapsr[20])==13:
						p=float(datapsr[18]); b=float(datapsr[16]) 
						if b>=(p**2)*(0.17*10**12): npsr+=1
						if p<=0.03: nmsp+=1
					if int(datapsr[7])==13:
						p=float(datapsr[6]); b=float(datapsr[5]) 
						if b>=(p**2)*(0.17*10**12): npsr+=1
						if p<=0.03: nmsp+=1

			NNS.append(nns); NPSR.append(npsr); NMSP.append(nmsp)
			model.append(i); time.append(float(datapsr[1])*t_conv)

		print NMSP
		print i

	fn.close()

	#np.savetxt('/projects/b1011/syr904/projects/PULSAR/kickgrid_runs/paper/data/nsnumber_newmodel_inpsr_9to12Gyr.dat', np.c_[model, time, NNS, NPSR, NMSP], fmt='%d %f %d %d %d', header='1.Model, 2.Time, 3.Nns, 4.Npsr, 5.Nmsp', delimiter='', comments='#')





def plot_BPhistory_inpsrfile(sourcedir, theid):
	Time, M0, M1, Bfield, Spin, K, A, ECC=get_history_inpsrfile(sourcedir, theid)
	plt.figure()
	plt.scatter(Spin, Bfield, 'o--')
	plt.yscale('log')
	plt.xlabel('P(sec)')
	plt.ylabel('B(G)')
	#plt.show()
	plt.savefig('projects/b1011/syr904/cmc/cmc-ns-bse/rundir_kconst/8e5rv1fb10kick1.0_kconst49/history/'+str(theid)+'_BP.pdf', dpi=300)


def plot_MBPchange(sourcedir):
	dataid=np.genfromtxt('/projects/b1011/syr904/projects/PULSAR/bse_change/MSPID_8e5kconst_fc.dat')
	id0=dataid[:,0]
	Mini=[]; Mfinl=[]; Pini=[]; Pfinl=[]; Bini=[]; Bfinl=[]
	for i in range(6,7):
		Time, M0, M1, Bfield, Spin, K, A, ECC=get_history_inpsrfile(sourcedir, id0[i])
		
		Mini.append(M0[0]); Mfinl.append(M0[-1])
		Pini.append(Spin[0]); Pfinl.append(Spin[-1]); Bini.append(Bfield[0]); Bfinl.append(Bfield[-1])

		
		fig, (ax0, ax1, ax2)=plt.subplots(nrows=3, sharex=True)
		ax0.plot(Time, M0, lw=1.5)
		#plt.xlabel(r'$Time(Myr)$')
		ax0.set_ylabel(r'$M_\odot$')
		ax0.set_title(str(id0[i]))
		ax0.set_xlim(10715,10725)
		##plt.show()

		ax1.plot(Time, Spin, lw=1.5)
		#plt.xlabel(r'$Time(Myr)$')
		ax1.set_ylabel(r'$P(sec)$')
		#plt.title(str(id0[i]))
		ax1.set_yscale('log')
		ax1.set_xlim(10715,10725)
		##plt.show()

		ax2.plot(Time, Bfield, lw=1.5)
		ax2.set_xlabel(r'$Time(Myr)$')
		ax2.set_ylabel(r'$B(G)$')
		ax2.set_xlim(10715,10725)
		##plt.title(str(id0[i]))
		ax2.set_yscale('log')

		plt.tight_layout()

		plt.savefig('/projects/b1011/syr904/projects/PULSAR/bse_change/history_kconst_fc/'+str(id0[i])+'_zoomin.pdf', dpi=300)
		##plt.show()

		print i

	#print Mini, Mfinl, Pini, Pfinl, Bini, Bfinl
	#Mini=np.array(Mini); Mfinl=np.array(Mfinl); Pini=np.array(Pini); Pfinl=np.array(Pfinl); Bini=np.array(Bini); Bfinl=np.array(Bfinl)

	#plt.figure()
	#plt.plot(Mfinl-Mini, (Bini**2-Bfinl**2)/(Pini-Pfinl)
	#plt.savefig('/projects/b1011/syr904/projects/PULSAR/bse_change/history_kconst_fc/M-P.pdf', dpi=300)


def get_id_position_BP(theid, snapshot):
	check=0
	with gzip.open(snapshot, 'r') as fsnap:
		for _ in xrange(2):
			next(fsnap)
		for line in fsnap:
			data=line.split()
			if int(data[7])==1:
				if int(data[10])==theid:
					r=float(data[2]); bfield=float(data[47]); spin=twopi*yearsc/float(data[45])
					dmdt0=np.float64(data[41]); dmdt1=np.float64(data[42]); radrol0=float(data[43]); radrol1=float(data[44]) 
					m0=float(data[8]); m1=float(data[9]); k0=int(data[17]); k1=int(data[18])
					a=float(data[12]); e=float(data[13]); fc=int(data[49])
					l=float(data[19])
					check=1
				if int(data[11])==theid:
					r=float(data[2]); bfield=float(data[48]); spin=twopi*yearsc/float(data[46])
					dmdt0=np.float64(data[42]); dmdt1=np.float64(data[41]); radrol0=float(data[44]); radrol1=float(data[43])
					m0=float(data[9]); m1=float(data[8]); k0=int(data[18]); k1=int(data[17])
					a=float(data[12]); e=float(data[13]); fc=int(data[50])
					l=float(data[20])
					check=1
			if int(data[7])!=1:
				if int(data[0])==theid:
					r=float(data[2]); bfield=float(data[60]); spin=twopi*yearsc/float(data[59])
					dmdt0=-100; dmdt1=-100; radrol0=-100; radrol1=-100
					m0=float(data[1]); m1=-100; k0=int(data[14]); k1=-100
					a=float(data[12]); e=float(data[13]); fc=int(data[61])
					l=float(data[15])
					check=1

			#if check==0: 

			if check==1: break		

	return r, bfield, spin, dmdt0, dmdt1, radrol0, radrol1, m0, m1, k0, k1, a, e, fc, l


def get_id_allmodel_position(idlist, pathlist):
	dataid=np.genfromtxt(idlist)
	model=dataid[:,0]; id0=dataid[:,1]; id1=dataid[:,2]
	sourcedir=np.genfromtxt(pathlist, dtype=str)
	rposition=[]; rc=[]; B=[]; P=[]; Nbh=[]; Ntot=[]
	f=open('/projects/b1011/syr904/projects/PULSAR/kickgrid_runs/paper/data/kickgrid_msp_newmodel_new.dat', 'a+', 0)
	f.write('#1.Model 2.ID0 3.ID1 4.r(pc) 5.B(G) 6.P(sec) 7.rc(pc) 8.Nbh 9.Ntot 10.dmdt0 11.dmdt1 12.rolrad0 13.rolrad1 14.m0 15.m1 16.k0 17.k1 18.a(AU) 19.ecc 20.Formation 21.L(mJy*kpc^2)\n')
	for i in range(len(model)):
		num=int(model[i])
		filepath=sourcedir[num]
		pref='initial'
		filestr=filepath+'/'+pref
		snaps=np.sort(glob(filestr+'.snap*.dat.gz'))
		lastsnap=snaps[-1]

		l_conv = conv('l', filestr+'.conv.sh')

		r_temp, b, s, dmt0, dmt1, rol0, rol1, m_0, m_1, k_0, k_1, A, E, F, L=get_id_position_BP(id0[i], lastsnap)
		rp=r_temp*l_conv
		L=L*Lsun
		rposition.append(rp); B.append(b); P.append(s)

		datadyn=np.genfromtxt(filestr+'.dyn.dat')
		rc_tempt=datadyn[-1][7]*l_conv
		rc.append(rc_tempt)

		nbh, ntot=dyn.find_NBH_NTOT_last(filestr)
		Nbh.append(nbh); Ntot.append(ntot)

		f.write('%d %d %d %f %e %f %f %d %d %e %e %f %f %f %f %d %d %f %f %d %f\n'%(num, id0[i], id1[i], rp, b, s, rc_tempt, nbh, ntot, dmt0, dmt1, rol0, rol1, m_0, m_1, k_0, k_1, A, E, F, L))

		print model[i]

	#np.savetxt('/projects/b1011/syr904/projects/PULSAR/kickgrid_runs/kickgrid_msp.dat', np.c_[model, id0, id1, rposition, B, P, rc, Nbh, Ntot], fmt='%d %d %d %f %e %f %f %d %d', header='Model ID0 ID1 r(pc) B(G) P(sec) rc(pc) Nbh Ntot', comments='#', delimiter= ' ')

	return rposition, rc, Nbh, Ntot


def get_id_allmodel_position_10to12Gyr(idlist, pathlist):
	dataid=np.genfromtxt(idlist)
	model=dataid[:,1]; id0=dataid[:,2]; id1=dataid[:,3]; snapt=dataid[:,0]
	#print snapt
	sourcedir=np.genfromtxt(pathlist, dtype=str)
	#sourcedir=[str(sourcedir)]   ##Use when there is only one path in the list
	#print sourcedir
	rposition=[]; rc=[]; B=[]; P=[]; Nbh=[]; Ntot=[]; time=[]
	m=25
	while (m <len(sourcedir)):
		f=open('/projects/b1011/syr904/projects/PULSAR/kickgrid_runs/paper/data/kickgrid_msp_10to12Gyr/msp_10to12Gyr_'+str(m)+'.dat', 'a+', 0)
		f.write('#Time Model ID0 ID1 r(pc) B(G) P(sec) rc(pc) Nbh Ntot dmdt0(codeunit) dmdt1(codeunit) radrol0 radrol1 m0 m1 k0 k1 a(AU) e Formation\n')
		filepath=sourcedir[m]
		#print filepath
		pref='initial'
		filestr=filepath+'/'+pref
		snaps=np.sort(glob(filestr+'.snap*.dat.gz'))
		#print snaps
		l_conv = conv('l', filestr+'.conv.sh')
		t_conv = conv('t', filestr+'.conv.sh')
		datadyn=np.genfromtxt(filestr+'.dyn.dat')
		databh=np.genfromtxt(filestr+'.bh.dat')
		print 'model=', m
		for i in range(len(model)):
			num=int(model[i])
			#if num>=20: continue
			if num==m:
			
				for j in range(len(snaps)-1, 0, -1):
					t_real=round(get_time(snaps[j])*t_conv/1000., 6)
					#print t_real
					if t_real==snapt[i]:
						print t_real
						t=get_time(snaps[j])
						#print t
						r_temp, b, s, dmt0, dmt1, rol0, rol1, m_0, m_1, k_0, k_1, A, E, F=get_id_position_BP(id0[i], snaps[j])
						r_temp=r_temp*l_conv
						rposition.append(r_temp); B.append(b); P.append(s); time.append(t)

						for k in range(len(datadyn[:,7])):
							tdyn=datadyn[:,0][k]
							if tdyn==t:
								rc_tempt=datadyn[:,7][k]*l_conv
								rc.append(rc_tempt)
								n_tot=datadyn[:,3][k]

						for h in range(len(databh[:,2])):
							tbh=databh[:,1][h]
							if tbh==t:
								n_bh=databh[:,2][h]

						Nbh.append(n_bh); Ntot.append(n_tot)
						#print "yes"
						f.write('%f %d %d %d %f %e %f %f %d %d %f %f %f %f %f %f %d %d %f %f %d\n'%(t_real, num, id0[i], id1[i], r_temp, b, s, rc_tempt, n_bh, n_tot, dmt0, dmt1, rol0, rol1, m_0, m_1, k_0, k_1, A, E, F))

				print model[i]

			if num>m: 
				m+=1; f.close()
				print 'break'; break

	#np.savetxt('/projects/b1011/syr904/projects/PULSAR/kickgrid_runs/kickgrid_msp10Gyr.dat', np.c_[model, id0, id1, rposition, B, P, rc, Nbh, Ntot], fmt='%d %d %d %f %e %f %f %d %d', header='Model ID0 ID1 r(pc) B(G) P(sec) rc(pc) Nbh Ntot', comments='#', delimiter= ' ')

	return rposition, rc, Nbh, Ntot


def get_interact_t_type(theid, pref, filepath):
	filestr=filepath+'/'+pref
	t_conv=conv('t', filestr+'.conv.sh')
	hdict=hic.history_maker([theid], [1], pref, filepath, 1.0)
	ts=[]; types=[]
	for i in hdict[theid]['binint']['binint'].keys():
		nopar=hdict[theid]['binint']['binint'][i]['interaction']['type']['nopars']
		time=hdict[theid]['binint']['binint'][i]['interaction']['type']['time']
		tp=hdict[theid]['binint']['binint'][i]['interaction']['type']['type']
		ts.append(time*t_conv); types.append(tp)

	return ts, types


def find_normalpsr(snapshot, modelno, t):
	time=[]; id0=[]; id1=[]; m0=[]; m1=[]; k0=[]; k1=[]; B=[]; P=[]; a=[]; ecc=[]; FC=[]; model=[]
	with gzip.open(snapshot, 'r') as fsnap:
		for _ in xrange(2):
			next(fsnap)
		for line in fsnap:
			datasnap=line.split()
			if int(datasnap[7])!=1:
				if int(datasnap[14])==13:
					spin=twopi*yearsc/float(datasnap[59])
					deathcut=(spin**2)*(0.17*10**12)
					if spin>0.03 and float(datasnap[60])>=deathcut:
						time.append(t); id0.append(int(datasnap[0])); id1.append(-100); m0.append(float(datasnap[1])); m1.append(-100)
						k0.append(int(datasnap[14])); k1.append(-100); FC.append(int(datasnap[61])); B.append(float(datasnap[60])); P.append(spin)
						a.append(-100); ecc.append(-100)
						model.append(modelno)

			if int(datasnap[7])==1:
				if int(datasnap[17])==13:
					spin0=twopi*yearsc/float(datasnap[45])
					deathcut0=(spin0**2)*(0.17*10**12)
					if spin0>0.03 and float(datasnap[47])>=deathcut0:
						time.append(t); id0.append(int(datasnap[10])); id1.append(int(datasnap[11])); m0.append(float(datasnap[8])); m1.append(float(datasnap[9]))
						k0.append(int(datasnap[17])); k1.append(int(datasnap[18])); FC.append(int(datasnap[49])); B.append(float(datasnap[47])); P.append(spin0)
						a.append(float(datasnap[12])); ecc.append(float(datasnap[13]))
						model.append(modelno)

				if int(datasnap[18])==13:
					spin1=twopi*yearsc/float(datasnap[46])
					deathcut1=(spin1**2)*(0.17*10**12)
					if spin1>0.03 and float(datasnap[48])>=deathcut1:
						time.append(t); id0.append(int(datasnap[11])); id1.append(int(datasnap[10])); m0.append(float(datasnap[9])); m1.append(float(datasnap[8]))
						k0.append(int(datasnap[18])); k1.append(int(datasnap[17])); FC.append(int(datasnap[50])); B.append(float(datasnap[48])); P.append(spin1)
						a.append(float(datasnap[12])); ecc.append(float(datasnap[13]))
						model.append(modelno)

	return model, time, B, P, FC, id0, id1, m0, m1, k0, k1, a, ecc


def get_normalpsr_last(pathlist):   ##using /projects/b1011/syr904/projects/PULSAR/kickgrid_runs/path_normalpsr_paper.dat for pathlist
	Md=[]; T=[]; BF=[]; S=[]; F=[]; ID0=[]; ID1=[]; M_0=[]; M_1=[]; K_0=[]; K_1=[]; Aaxis=[]; E=[]
	sourcedir=np.genfromtxt(pathlist, dtype=str)
	for i in range(len(sourcedir)):
		pref='initial'
		filepath=sourcedir[i]
		filestr=filepath+'/'+pref
		t_conv=conv('t', filestr+'.conv.sh')
		snaps=np.sort(glob(filestr+'.snap*.dat.gz'))
		lastsnap=snaps[-1]
		t=get_time(lastsnap)*t_conv

		Model, Time, Bf, Spin, Fc, Id0, Id1, M0, M1, K0, K1, A, Ecc=find_normalpsr(lastsnap, i, t)
		Md=Md+Model; T=T+Time; BF=BF+Bf; S=S+Spin; F=F+Fc; ID0=ID0+Id0; ID1=ID1+Id1; M_0=M_0+M0; M_1=M_1+M1; K_0=K_0+K0; K_1=K_1+K1; Aaxis=Aaxis+A; E=E+Ecc

		print i

	np.savetxt('/projects/b1011/syr904/projects/PULSAR/kickgrid_runs/paper/data/kickgrid_normalpsr_newmodel.dat', np.c_[Md, ID0, ID1, BF, S, F, M_0, M_1, K_0, K_1, Aaxis, E], fmt ='%d %f %e %f %d %d %d %f %f %d %d %f %f', delimiter= ' ', header = 'Model id0 id1 B(G) P(sec) FC m0 m1 k0 k1 a[AU] ecc', comments = '#')


def get_normalpsr_10Gyr(filepath, savepath, modelno):
	pref='initial'
	filestr=filepath+'/'+pref
	t_conv=conv('t', filestr+'.conv.sh')
	snaps=np.sort(glob(filestr+'.snap*.dat.gz'))
	for i in range(len(snaps)-1, 0, -1):
		ti=get_time(snaps[i])*t_conv/1000.
		print i
		if ti<10.:
			Model, Time, Bf, Spin, Fc, Id0, Id1, M0, M1, K0, K1, A, Ecc=find_normalpsr(snaps[i+1], modelno, ti)
			break


	np.savetxt(savepath+'/kickgrid_normalpsr_10Gyr'+str(modelno)+'.dat', np.c_[Model, Time, Bf, Spin, Fc, Id0, Id1, M0, M1, K0, K1, A, Ecc], fmt ='%d %f %e %f %d %d %d %f %f %d %d %f %f', delimiter= ' ', header = 'Model Age(Myr) B(G) P(sec) FC id0 id1 m0 m1 k0 k1 a[AU] ecc', comments = '#')



def get_normalpsr_ini(pathlist, savepath):
	ffinl=np.sort(glob(savepath+'/'+'finl*.dat'))
	sourcedir=np.genfromtxt(pathlist, dtype=str)
	for i in range(len(ffinl)):
		datafinl=np.genfromtxt(ffinl[i])
		ids=datafinl[:,4]

		time=[]; id0=[]; id1=[]; m0=[]; m1=[]; k0=[]; k1=[]; B=[]; P=[]; a=[]; ecc=[]; FC=[]
		filepath=sourcedir[i]
		pref='initial'
		filestr=filepath+'/'+pref
		t_conv=conv('t', filestr+'.conv.sh')
		snaps=np.sort(glob(filestr+'.snap*.dat.gz'))

		fini=open(savepath+'/ini'+str(i)+'.dat', 'a+', 0)
		fini.write('#Age(Myr) B(G) P(sec) FC id0 id1 m0 m1 k0 k1 a[AU] ecc\n')

		for j in range(len(ids)):
			theid=int(ids[j])
			print 'theid=',theid
			check=0
			for k in range(len(snaps)):
				print k
				if check==0:
					with gzip.open(snaps[k], 'r') as fsnap:
						for _ in xrange(2):
							next(fsnap)
						for line in fsnap:
							datasnap=line.split()
							if int(datasnap[0])==theid:
								if int(datasnap[14])==13:
									check=1
									t=get_time(snaps[k])*t_conv
									spin=twopi*yearsc/float(datasnap[59])
									time=t; id0=int(datasnap[0]); id1=-100; m0=float(datasnap[1]); m1=-100
									k0=int(datasnap[14]); k1=-100; FC=-100; B=float(datasnap[60]); P=spin
									a=-100; ecc=-100
									print theid
									fini.write('%f %e %f %d %d %d %f %f %d %d %f %f\n'%(time, B, P, FC, id0, id1, m0, m1, k0, k1, a, ecc))
								break

							if int(datasnap[10])==theid:
								if int(datasnap[17])==13:
									check=1
									t=get_time(snaps[k])*t_conv
									spin0=twopi*yearsc/float(datasnap[45])
									time=t; id0=int(datasnap[10]); id1=int(datasnap[11]); m0=float(datasnap[8]); m1=float(datasnap[9])
									k0=int(datasnap[17]); k1=int(datasnap[18]); FC=int(datasnap[49]); B=float(datasnap[47]); P=spin0
									a=float(datasnap[12]); ecc=float(datasnap[13])
									print theid
									fini.write('%f %e %f %d %d %d %f %f %d %d %f %f\n'%(time, B, P, FC, id0, id1, m0, m1, k0, k1, a, ecc))
								break

							if int(datasnap[11])==theid:
								if int(datasnap[18])==13:
									check=1
									t=get_time(snaps[k])*t_conv
									spin1=twopi*yearsc/float(datasnap[46])
									time=t; id0=int(datasnap[11]); id1=int(datasnap[10]); m0=float(datasnap[9]); m1=float(datasnap[8])
									k0=int(datasnap[18]); k1=int(datasnap[17]); FC=int(datasnap[50]); B=float(datasnap[48]); P=spin1
									a=float(datasnap[12]); ecc=float(datasnap[13])
									print theid
									fini.write('%f %e %f %d %d %d %f %f %d %d %f %f\n'%(time, B, P, FC, id0, id1, m0, m1, k0, k1, a, ecc))
								break

				else: break

		fini.close()
		print i
		#np.savetxt(savepath+'/ini'+str(i)+'.dat', np.c_[time, B, P, FC, id0, id1, m0, m1, k0, k1, a, ecc], fmt ='%f %e %f %d %d %d %f %f %d %d %f %f', delimiter= ' ', header = 'Age(Myr) B(G) P(sec) FC id0 id1 m0 m1 k0 k1 a[AU] ecc', comments = '#') 


def get_2Dradius(projfile, obsfile, mspids):
	rbh=[]; rmsp=[]; rns=[]; mbh=[]; mns=[]; mmsp=[]
	dataobs=np.genfromtxt(obsfile)
	rc=dataobs[0, 7]; rhl=dataobs[0, 8]

	with open(projfile, 'r') as fproj:
		for _ in xrange(2):
			next(fproj)
		for line in fproj:
			dataproj=line.split()
			if int(dataproj[3])==13:
				rns.append(float(dataproj[0]))
				mns.append(float(dataproj[9]))
				try:
					ids=mspids.index(int(dataproj[12]))
				except ValueError:
					continue
				else:
					#print int(dataproj[12])
					rmsp.append(float(dataproj[0]))
					mmsp.append(float(dataproj[9]))
			if int(dataproj[3])==14:
				if float(dataproj[0])<=rc:
					rbh.append(float(dataproj[0]))
					mbh.append(float(dataproj[9]))

			if int(dataproj[5])==13 or int(dataproj[6])==13:
				rns.append(float(dataproj[0]))
				mns.append(float(dataproj[9]))
				if int(dataproj[5])==13:
					try:
						ids=mspids.index(int(dataproj[13]))
					except ValueError:
						continue
					else:
						#print int(dataproj[13])
						rmsp.append(float(dataproj[0]))
						mmsp.append(float(dataproj[9]))
				if int(dataproj[6])==13:
					try:
						ids=mspids.index(int(dataproj[14]))
					except ValueError:
						continue
					else:
						#print int(dataproj[14])
						rmsp.append(float(dataproj[0]))
						mmsp.append(float(dataproj[9]))

			if int(dataproj[5])==14 or int(dataproj[6])==14:
				if float(dataproj[0])<=rc:
					rbh.append(float(dataproj[0]))
					mbh.append(float(dataproj[9]))

	#print rmsp
	#print mbh, mns, mmsp

	#print rc, rhl

	return rbh, rns, rmsp, mbh, mns, mmsp, rc, rhl


def get_mean2Dradius_allmodels(pathlist, msplist, start, end):
	sourcedir=np.genfromtxt(pathlist, dtype=str)
	#print sourcedir
	datamsp=np.genfromtxt(msplist)
	model=datamsp[:,0]; nsid=datamsp[:,1]

	rbh_mean=[]; rns_mean=[]; rmsp_mean=[]
	dbh=[]; dbh_mean=[]
	RC=[]; RHL=[]


	for i in range(start, end):
		#print i
		pref='initial'
		filepath=sourcedir[i]
		filestr=filepath+'/'+pref
		projs=np.sort(glob(filestr+'.snap*.2Dproj.dat'))
		observs=np.sort(glob(filestr+'.snap*.obs_params.dat'))
		#print projs, observs
		lastproj=projs[-1]
		lastobs=observs[-1]

		mspids=[]
		for j in range(len(model)):
			if int(model[j])==i:
				mspids.append(int(nsid[j]))

		#print mspids

		Rbh, Rns, Rmsp, Mbh, Mns, Mmsp, Rc, Rhl=get_2Dradius(lastproj, lastobs, mspids)
		RC.append(Rc); RHL.append(Rhl)

		rbh_mean.append(np.mean(Rbh)); rns_mean.append(np.mean(Rns)); rmsp_mean.append(np.mean(Rmsp))
		rbh_max=np.max(Rbh); mbh_sum=np.sum(Mbh); mbh_mean=np.mean(Mbh)
		dbh.append(mbh_sum/(np.pi*rbh_max**2))
		dbh_mean.append(mbh_mean/(np.pi*rbh_mean[i]**2))

	#print rbh_mean, rns_mean, rmsp_mean, dbh, dbh_mean
	return rbh_mean, rns_mean, rmsp_mean, dbh, dbh_mean


def get_3Dradius(snapfile, rc3D):
	rbh3D=[]; mbh3D=[]	
	with gzip.open(snapfile, 'r') as fsnap:
		for _ in xrange(2):
			next(fsnap)
		for line in fsnap:
			datasnap=line.split()
			#print datasnap
			if int(datasnap[14])==14:
				if float(datasnap[2])<=rc3D:
					rbh3D.append(float(datasnap[2]))
					mbh3D.append(float(datasnap[1]))

			if int(datasnap[17])==14 or int(datasnap[18])==14:
				if float(datasnap[2])<=rc3D:
					rbh3D.append(float(datasnap[2]))
					mbh3D.append(float(datasnap[1]))

	return rbh3D, mbh3D


def get_mean3Dradius_allmodels(pathlist, start, end):
	sourcedir=np.genfromtxt(pathlist, dtype=str)
	rbh3D_mean=[]; RC3D=[]
	dbh3D=[]; dbh3D_mean=[]
	for i in range(start, end):
		#print i
		pref='initial'
		filepath=sourcedir[i]
		filestr=filepath+'/'+pref
		snaps=np.sort(glob(filestr+'.snap*.dat.gz'))
		lastsnap=snaps[-1]
		#print lastsnap

		l_conv=conv('l', filestr+'.conv.sh')

		datadyn=np.genfromtxt(filestr+'.dyn.dat')
		Rc3D=datadyn[-1][7]

		Rbh3D, Mbh3D=get_3Dradius(lastsnap, Rc3D)
		Rbh3D=np.array(Rbh3D)*l_conv

		RC3D.append(Rc3D*l_conv)

		rbh3D_mean.append(np.mean(Rbh3D))
		rbh3D_max=np.max(Rbh3D); mbh3D_sum=np.sum(Mbh3D); mbh3D_mean=np.mean(Mbh3D)
		dbh3D.append(mbh3D_sum/(np.pi*rbh3D_max**2))
		dbh3D_mean.append(mbh3D_mean/(np.pi*rbh3D_mean[i]**2))

	return rbh3D_mean, dbh3D, dbh3D_mean


def get_id_formationchannel(pathlist, msplist):
	fc=[]
	sourcedir=np.genfromtxt(pathlist, dtype=str)
	datamsp=np.genfromtxt(msplist)
	model=datamsp[:,0]; ids=datamsp[:,1]
	for i in range(len(model)):
		no=int(model[i])
		if no==24:
			theid=int(ids[i])
			filepath=sourcedir[no]
			filestr=filepath+'/initial'
			snaps=np.sort(glob(filestr+'.snap*.dat.gz'))
			check=0
			for j in range(len(snaps)-1, 0, -2):
				with gzip.open(snaps[j], 'r') as fsnap:
					for _ in xrange(2):
						next(fsnap)
					for line in fsnap:
						datasnap=line.split()
						if int(datasnap[10])==theid: fc.append(int(datasnap[49])); check=1; break
						if int(datasnap[11])==theid: fc.append(int(datasnap[50])); check=1; break

				if check==1: break

			if check==0: fc.append(-100)

			print i

	print fc
	return fc


####################
###Comparison to Observations###
####################

def XrayLum(M_A, M_D, A, K1):  ##Equation (29) in Rappaport et al 1982
##K=0, alpha=q, beta=1 for WD donor
	Const_grav=6.67408*10**-8
	c_lightspeed=3*10**10
	solarmass=2*10**33 ##g
	M_A=M_A*solarmass
	M_D=M_D*solarmass
	R_A=1.4e-05*6.96e10
	A=A*1.5e13
	Porb=2*np.pi*math.sqrt(A**3/(Const_grav*(M_D+M_A)))
	q=M_A/M_D
	Jterm=-(32./5.)*(4*np.pi**2)**(4./3.)*Const_grav**(5./3.)*c_lightspeed**-5*(M_A+M_D)**(-1./3.)*M_A*M_D*Porb**(-8./3.)
	M_Ddot=np.abs(M_D*Jterm/(2./3.-1/q))
	print M_Ddot/solarmass*365*24*3600

	eta_bol=0.55; epsilon=1.0  ##Kremer et al 2018
	Lx=eta_bol*epsilon*Const_grav*M_A*M_Ddot/R_A   ##in erg/s

	if K1<10: f=1
	else: f=6

	Mcrit=5.3e-11*f*(M_A/solarmass)**0.3*(Porb/3600.)**1.4

	return Lx, M_Ddot/solarmass*365*24*3600, Mcrit



def get_r_v(idlist, pathlist):
	dataid=np.genfromtxt(idlist, dtype=int)
	model=dataid[:,0]; id0=dataid[:,1]; id1=dataid[:,2]
	sourcedir=np.genfromtxt(pathlist, dtype=str)
	r2D=[]; vz=[]; r3D=[]; vt=[]; Mr2D=[]
	for i in range(len(id0)):
		pref='initial'
		num=int(model[i])
		filepath=sourcedir[num]
		filestr=filepath+'/'+pref
		snapproj=np.sort(glob(filestr+'.snap*.2Dproj.dat'))
		snaps=np.sort(glob(filestr+'.snap*.dat.gz'))
		lastproj=snapproj[-1]
		lastsnap=snaps[-1]
		Mtot=0
		with open(lastproj, 'r') as fproj:
			for _ in xrange(2):
				next(fproj)
			for line in fproj:
				dataproj=line.split()
				Mtot+=float(dataproj[9])
				if int(float(dataproj[12]))==id0[i] or int(float(dataproj[13]))==id0[i] or int(float(dataproj[14]))==id0[i]:
					Mtot=Mtot-float(dataproj[9])
					r2D.append(float(dataproj[0])); vz.append(float(dataproj[20])); Mr2D.append(Mtot)
					print i
				
		with gzip.open(lastsnap, 'r') as fsnap:
			for _ in xrange(2):
				next(fsnap)
			for line in fsnap:
				datasnap=line.split()
				if int(datasnap[0])==id0[i] or int(datasnap[10])==id0[i] or int(datasnap[11])==id0[i]:
					r3D.append(float(datasnap[2])); vt.append(float(datasnap[4]))
					print i


	print len(model), len(r2D), len(vz), len(r3D), len(vt), len(Mr2D)


	np.savetxt('/projects/b1011/syr904/projects/PULSAR/kickgrid_runs/paper/data/kickgrid_accl_normalpsr_newmodel.dat', np.c_[model, id0, id1, r2D, vz, r3D, vt, Mr2D], fmt ='%d %d %d %f %f %f %f %f', delimiter = ' ', header = '1.model, 2.id0, 3.id1, 4.r2D(pc), 5.vl(km/s), 6.r3D(codeunit), 7.vt(codeunit), 8. Mtot(<r2D)(Msun)', comments = '#')

						
def add_acceleration(r_2d, v_l, pdot0, p0, mr2d, a, e, m0, m1):
	if m1==-100:
		##for single pulsar acceleration is from GC potential (Phinney 1992 2.4)
		maxal=1.1*Gconst*(mr2d*Msun)/(np.pi*(r_2d*PC)**2)
		pmax=(1+abs(v_l)*10**5/clight)*p0
		pmin=(1-abs(v_l)*10**5/clight)*p0
		#print abs(v_l)*10**5/clight
		pdotmax=pmax*(pdot0/p0+maxal/clight)
		pdotmin=pmin*(pdot0/p0-maxal/clight)


	else:
		##for binary pulsar need to add doppler shift from binary motions
		rp=(a*AU)*(1-e)  ##periapsis distance in cm
		miu=Gconst*(m0+m1)*Msun  ##standard gravitational parameter
		#print miu*(2./rp-1./(a*AU))
		vrel=math.sqrt(miu*(2./rp-1./(a*AU)))   ##in cm/s
		v0=m1*vrel/(m0+m1)
		v1=m0*vrel/(m0+m1)
		a0=Gconst*m1*Msun/rp**2; a1=Gconst*m0*Msun/rp**2

		maxal=1.1*Gconst*(mr2d*Msun)/(np.pi*(r_2d*PC)**2)+a0
		v_lmax=abs(v_l)*10**5+v0    ##in cm/s
		#print v_lmax/clight
		pmax=(1+v_lmax/clight)*p0
		pmin=(1-v_lmax/clight)*p0
		pdotmax=pmax*(pdot0/p0+maxal/clight)
		pdotmin=pmin*(pdot0/p0-maxal/clight)

	#print maxal, m1

	return  pmin, pmax, pdotmin, pdotmax



def add_acceleration_GCpotential(r_2d, v_l, pdot0, p0, mr2d, a, e, m0, m1):
	maxal=1.1*Gconst*(mr2d*Msun)/(np.pi*(r_2d*PC)**2)
	p=(1+v_l*10**5/clight)*p0
	#p=(1-abs(v_l)*10**5/clight)*p0
	#print p, maxal
	pdotmax=p*(pdot0/p0+maxal/clight)
	pdotmin=p*(pdot0/p0-maxal/clight)

	return  pdotmin, pdotmax, p


def get_acceleration(acclfile, psrfile):
	dataacc=np.genfromtxt(acclfile)
	datapsr=np.genfromtxt(psrfile)

	model=dataacc[:,0]; id0=dataacc[:,1]; id1=dataacc[:,2]; r2d=dataacc[:,3]; vl=dataacc[:,4]; Mr2d=dataacc[:,7]
	
	A=datapsr[:,17]; ecc=datapsr[:,18]; M0=datapsr[:,13]; M1=datapsr[:,14]; P=datapsr[:,5]; B=datapsr[:,4]

	Pdot=Kconst*yearsc*np.array(B)*np.array(B)/np.array(P)
	#print Pdot/np.array(P)

	pminus=[]; pplus=[]; pdotminus=[]; pdotplus=[]; pdotl=[]; pdotu=[]; paccl=[]
	for i in range(len(id0)):
		#Pmin, Pmax, Pdotmin, Pdotmax=add_acceleration(r2d[i], vl[i], Pdot[i], P[i], Mr2d[i], A[i], ecc[i], M0[i], M1[i])
		Pdotmin, Pdotmax, Paccl=add_acceleration_GCpotential(r2d[i], vl[i], Pdot[i], P[i], Mr2d[i], A[i], ecc[i], M0[i], M1[i])
		#print Pmin, Pmax, Pdotmin, Pdotmax
		#print Pdotmin, Pdotmax
		#pminus.append(P[i]-Pmin); pplus.append(Pmax-P[i])
		pdotminus.append(Pdot[i]-Pdotmin); pdotplus.append(Pdotmax-Pdot[i])
		pdotl.append(Pdotmin); pdotu.append(Pdotmax)
		paccl.append(Paccl)

		#print i

	#return pminus, pplus, pdotminus, pdotplus
	return pdotminus, pdotplus, pdotl, pdotu, paccl


##Gammas
def get_EncounterRates_lastsnap(pathlist):
	sourcedir=np.genfromtxt(pathlist, dtype=str)
	Rc=[]; Rho0=[]; V0=[]; GAMMA=[]; model=[]
	for i in range(len(sourcedir)):
		filepath=sourcedir[i]
		pref='initial'
		filestr=filepath+'/'+pref
		snaps=np.sort(glob(filestr+'.snap*.dat.gz'))
		obs=np.sort(glob(filestr+'.snap*.obs_params.dat'))
		projs=np.sort(glob(filestr+'.snap*.2Dproj.dat'))
		#snapno_max=str(int(len(snaps)-1)).zfill(4)
		#projfile=filestr+'.snap'+snapno_max+'.2Dproj.dat'
		#obsfile=filestr+'.snap'+snapno_max+'.obs_params.dat'
		projfile=projs[-1]
		obsfile=obs[-1]
		print obsfile
		
		dataobs=np.genfromtxt(obsfile)
		rc=dataobs[0,7]; rhl=dataobs[0,8]     ###in pc
		print rc, rhl
		Rc.append(rc)
		
		dataproj=np.genfromtxt(projfile)
		r2d=dataproj[:,0]; Ltot=dataproj[:,1]; binflag=dataproj[:,2]; k=dataproj[:,3]
		k0=dataproj[:,5]; k1=dataproj[:,6]; vx=dataproj[:,18]; vy=dataproj[:,19]; vz=dataproj[:,20]

		Ltotc=0; v1=[]; v2=[]; v3=[]
		for j in range(len(r2d)):
			if r2d[j]<=rc:
				if binflag[j]!=1 and k[j]<13:
					if Ltot[j]<=20:
						Ltotc+=Ltot[j]
						v1.append(vx[j]); v2.append(vy[j]); v3.append(vz[j])
				if binflag[j]==1 and k0[j]<13:
					if Ltot[j]<=20:
						Ltotc+=Ltot[j]
						v1.append(vx[j]); v2.append(vy[j]); v3.append(vz[j])
				if binflag[j]==1 and k1[j]<13:
					if Ltot[j]<=20:
						Ltotc+=Ltot[j]
						v1.append(vx[j]); v2.append(vy[j]); v3.append(vz[j])

		Rho0.append(Ltotc/rc**3)

		v1, v2, v3=np.array(v1), np.array(v2), np.array(v3)
		vsigma=[np.std(v1), np.std(v2), np.std(v3)]
		vc=np.mean(vsigma)

		V0.append(vc)

		model.append(i)

	Rc, Rho0, V0=np.array(Rc), np.array(Rho0), np.array(V0)

	Gamma=Rho0**2*Rc**3/V0
	sgamma=Rho0/V0
	#Nest=np.exp(-1.1+1.5*np.log10(Gamma))

	print Rc, Rho0, V0, Gamma, sgamma

	np.savetxt('/projects/b1011/syr904/projects/PULSAR/kickgrid_runs/paper/data/GammaModel_newmodel.dat', np.c_[model, Rho0, Rc, V0, Gamma, sgamma], fmt ='%d %f %f %f %f %f', delimiter = ' ', header = '1.model, 2. rho_c[L_sun/pc^3], 3.r_c[pc], 4.v_c[km/s], 5.Gamma, 6.gamma', comments = '#')

	return Gamma, sgamma



def get_EncounterRates_10to12Gyr(pathlist):
	sourcedir=np.genfromtxt(pathlist, dtype=str)
	fout=open('/projects/b1011/syr904/projects/PULSAR/kickgrid_runs/paper/data/GammaModel_10to12Gyr_newmodel.dat', 'a+', 0)
	fout.write('#1.Model, 2.t, 3.Gamma\n')
	for i in range(len(sourcedir)):
		filepath=sourcedir[i]
		pref='initial'
		filestr=filepath+'/'+pref
		snaps=np.sort(glob(filestr+'.snap*.dat.gz'))
		snapobs=np.sort(glob(filestr+'.snap*.obs_params.dat'))
		snapproj=np.sort(glob(filestr+'.snap*.2Dproj.dat'))
		snapno_max=str(int(len(snaps)-1)).zfill(4)
		Rc=[]; Rho0=[]; V0=[]; t=[]
		for k in range(len(snapobs)-1, 0, -1):	
			projfile=snapproj[k]
			obsfile=snapobs[k]
			dataobs=np.genfromtxt(obsfile)
			time=dataobs[0,10]
			if time>=10000.0:
				t.append(time)
				rc=dataobs[0,7]; rhl=dataobs[0,8]     ###in pc
				#print rc, rhl
				Rc.append(rc)
		
				dataproj=np.genfromtxt(projfile)
				r2d=dataproj[:,0]; Ltot=dataproj[:,1]; binflag=dataproj[:,2]; k=dataproj[:,3]
				k0=dataproj[:,5]; k1=dataproj[:,6]; vx=dataproj[:,18]; vy=dataproj[:,19]; vz=dataproj[:,20]

				Ltotc=0; v1=[]; v2=[]; v3=[]
				for j in range(len(r2d)):
					if r2d[j]<=rc:
						if binflag[j]!=1 and k[j]<13:
							if Ltot[j]<=20:
								Ltotc+=Ltot[j]
								v1.append(vx[j]); v2.append(vy[j]); v3.append(vz[j])
						if binflag[j]==1 and k0[j]<13:
							if Ltot[j]<=20:
								Ltotc+=Ltot[j]
								v1.append(vx[j]); v2.append(vy[j]); v3.append(vz[j])
						if binflag[j]==1 and k1[j]<13:
							if Ltot[j]<=20:
								Ltotc+=Ltot[j]
								v1.append(vx[j]); v2.append(vy[j]); v3.append(vz[j])

				Rho0.append(Ltotc/rc**3)

				##Calculate the velocity dispersion at the center of the model
				v1, v2, v3=np.array(v1), np.array(v2), np.array(v3)
				vsigma=[np.std(v1), np.std(v2), np.std(v3)]
				vc=np.mean(vsigma)
				V0.append(vc)

				gamma=(Ltotc/rc**3)**2*rc**3/vc

				fout.write('%d %f %f\n'%(i, time, gamma))
				print k

			if time<10000.0: break

		Rc, Rho0, V0=np.array(Rc), np.array(Rho0), np.array(V0)

		Gamma=Rho0**2*Rc**3/V0   ##this is an array
		print i, t, Gamma

		#np.savetxt(fout, np.c_[t, Gamma], fmt ='%f %f', delimiter = ' ', header = '1.t, 2.Gamma', comments = '#')

	fout.close()

		

def DoubleNS(path):
	filestr=path+'/initial'
	snaps=np.sort(glob(filestr+'.snap*.dat.gz'))
	lastsnap=snaps[-1]

	B0=[]; B1=[]; P0=[]; P1=[]; M0=[]; M1=[]; A=[]; E=[]; ID0=[]; ID1=[]; FC0=[]; FC1=[]
	with gzip.open(lastsnap, 'r') as fsnap:
		for _ in xrange(2):
			next(fsnap)
		for line in fsnap:
			data=line.split()
			if int(data[17])==13 and int(data[18])==13:
				B0.append(float(data[47])); B1.append(float(data[48])); P0.append(twopi*yearsc/float(data[45])); P1.append(twopi*yearsc/float(data[46])); M0.append(float(data[8])); M1.append(float(data[9])); A.append(float(data[12])); E.append(float(data[13])); ID0.append(int(data[10])); ID1.append(int(data[11])); FC0.append(int(data[49])); FC1.append(int(data[50]))


	np.savetxt('/projects/b1011/syr904/projects/PULSAR/kickgrid_runs/paper/data/doubelNS.dat', np.c_[ID0, ID1, M0, M1, B0, B1, P0, P1, FC0, FC1, A, E], fmt='%d %d %f %f %e %e %f %f %d %d %f %f', header='1.ID0, 2.ID1, 3.M0(Msun), 4.M1(Msun), 5.B0(G), 6.B1(G), 7.P0(sec), 8.P1(sec), 9.FC0, 10.FC1, 11.a(AU), 12.ecc', delimiter='', comments='#')


def find_DNS_NSBH(pathlist, start, end):
	sourcedir=np.genfromtxt(pathlist, dtype=str)
	pref='initial'
	Model=[]; Time=[]; ID0=[]; ID1=[]; P0=[]; P1=[]; B0=[]; B1=[]; M0=[]; M1=[]; A=[]; E=[]; Type=[]
	for i in range(start, end):
		filestr=sourcedir[i]+'/'+pref
		t_conv=conv('t', filestr+'.conv.sh')
		snaps=np.sort(glob(filestr+'.snap*.dat.gz'))
		with open(filestr+'.ns.dat', 'r') as fns:
			next(fns)
			for line in fns:
				datans=line.split()
				if float(datans[0])*t_conv>=9000.0:
					if int(datans[7])>0:
						model=i; t=float(datans[0])
						print model, t, 'DNS'
						for j in range(len(snaps)-1, 0, -1):
							time=get_time(snaps[j])
							if round(time, 6)==round(t, 6):
								print model, t, 'DNS'
								with gzip.open(snaps[j], 'r') as fsnap:
									next(fsnap)
									next(fsnap)
									for line in fsnap:
										datasnap=line.split()
										if int(datasnap[17])==13 and int(datasnap[18])==13:
											Model.append(i); Time.append(time); ID0.append(int(datasnap[10])); ID1.append(int(datasnap[11])); P0.append(twopi*yearsc/float(datasnap[45])); P1.append(twopi*yearsc/float(datasnap[46])); B0.append(float(datasnap[47])); B1.append(float(datasnap[48])); M0.append(float(datasnap[8])); M1.append(float(datasnap[9])); A.append(float(datasnap[12])); E.append(float(datasnap[13])); Type.append(1313) 



					if int(datans[8])>0: 
						model=i; t=float(datans[0])
						print model, t, 'NSBH'
						for j in range(len(snaps)-1, 0, -1):
							time=get_time(snaps[j])
							if round(time, 6)==round(t, 6):
								print model, t, 'NSBH'
								with gzip.open(snaps[j], 'r') as fsnap:
									next(fsnap)
									next(fsnap)
									for line in fsnap:
										datasnap=line.split()
										if int(datasnap[17])==13 and int(datasnap[18])==14:
											Model.append(i); Time.append(time); ID0.append(int(datasnap[10])); ID1.append(int(datasnap[11])); P0.append(twopi*yearsc/float(datasnap[45])); P1.append(-100); B0.append(float(datasnap[47])); B1.append(-100); M0.append(float(datasnap[8])); M1.append(float(datasnap[9])); A.append(float(datasnap[12])); E.append(float(datasnap[13])); Type.append(1314)
										if int(datasnap[17])==14 and int(datasnap[18])==13:
											Model.append(i); Time.append(time); ID0.append(int(datasnap[10])); ID1.append(int(datasnap[11])); P0.append(-100); P1.append(twopi*yearsc/float(datasnap[46])); B0.append(-100); B1.append(float(datasnap[48])); M0.append(float(datasnap[8])); M1.append(float(datasnap[9])); A.append(float(datasnap[12])); E.append(float(datasnap[13])); Type.append(1413)

		#print Model, ID0, ID1, P0, P1, B0, B1, M0, M1, A, E, Type

	np.savetxt('/projects/b1011/syr904/projects/PULSAR/kickgrid_runs/paper/data/DNS-NSBH.dat', np.c_[Model, Time, ID0, ID1, P0, P1, B0, B1, M0, M1, A, E, Type], fmt='%d %f %d %d %f %f %e %e %f %f %f %f %d', header='1.Model, 2.Time, 3.ID0, 4.ID1, 5.P0, 6.P1, 7.B0, 8.B1, 9.M0, 10.M1, 11.a, 12.ecc, 13.Type', delimiter='', comments='#')


def DNS_NSBH_Unique(datfile, pathlist):
	data=np.genfromtxt(datfile)
	model=data[:,0]; id0=data[:,2]; id1=data[:,3]
	print np.array(data[0,:])
	allmodel=Counter(model).keys()

	lines=[]; lines_needed=[]
	for i in range(len(allmodel)):
		modelno=int(allmodel[i])
		idmodel=[]	
		for j in range(len(model)):
			if model[j]==modelno:
				idmodel.append(str(id0[j])+str(id1[j]))

		idkey=Counter(idmodel).keys()
		print idkey


		for k in range(len(idkey)):
			for l in range(len(model)):
				if model[l]==modelno:
					ids=str(id0[l])+str(id1[l])
					if ids==idkey[k]:
						lines.append(l)
			lines_needed.append(lines[-1])

		print modelno

	print lines_needed

	with open('/projects/b1011/syr904/projects/PULSAR/kickgrid_runs/paper/data/DNS_NSBH_Unique.dat', 'a+') as out_file:
		for m in range(len(lines_needed)):
			listm=' '.join(str(e) for e in data[lines_needed[m],:])
			out_file.write(listm+'\n')



def BNS_NSBH_Appearperiod(datfile,pathlist):
	data=np.genfromtxt(datfile, usecols=(0, 1, 2, 3, 12))
	model=data[:,0]; time=data[:,1]; id0=data[:,2]; id1=data[:,3]; tp=data[:,4]
	allmodel=Counter(model).keys()
	sourcedir=np.genfromtxt(pathlist, dtype=str)

	ids=[]; periods=[]; types=[]
	for i in range(len(allmodel)):
		modelno=int(allmodel[i])
		t_conv=conv('t', sourcedir[modelno]+'/initial.conv.sh')

		idmodel=[]; pmodel=[]; tmodel=[]
		for j in range(len(model)):
			if model[j]==modelno:
				idmodel.append(str(id0[j])+str(id1[j]))
				pmodel.append(time[j])
				tmodel.append(tp[j])

		idkey=Counter(idmodel).keys()
		ids=ids+idkey

		for k in range(len(idkey)):
			p_temp=[]; t_temp=[]
			for l in range(len(idmodel)):
				if idmodel[l]==idkey[k]:
					p_temp.append(pmodel[l])
					t_temp.append(tmodel[l])

			#print p_temp

			if len(p_temp)==1: periods.append(30); types.append(t_temp[-1])
			if len(p_temp)>1: periods.append((p_temp[-1]-p_temp[0])*t_conv); types.append(t_temp[-1])


	return ids, periods, types


def Ave_stellardensity(pathlist):
	sourcedir=np.genfromtxt(pathlist, dtype=str)
	median=[2.86, 0.88, 0.25]
	mode=[1.29,0.22,0.034]
	Npoor_median=[]; Nmid_median=[]; Nrich_median=[]
	Npoor_mode=[]; Nmid_mode=[]; Nrich_mode=[]
	for i in range(len(sourcedir)):
		n_poor_median=0; n_mid_median=0; n_rich_median=0
		n_poor_mode=0; n_mid_mode=0; n_rich_mode=0
		path=sourcedir[i]
		filestr=path+'/'+'initial'
		projs=np.sort(glob(filestr+'.snap*.2Dproj.dat'))
		lastproj=projs[-1]
		with open(lastproj, 'r') as fproj:
			next(fproj)
			next(fproj)
			for line in fproj:
				dataproj=line.split()
				if i<10:
					if float(dataproj[0])<=median[0]: n_rich_median+=1
					if float(dataproj[0])<=mode[0]: n_rich_mode+=1
				elif i<19:
					if float(dataproj[0])<=median[1]: n_mid_median+=1
					if float(dataproj[0])<=mode[1]: n_mid_mode+=1
				elif i<25:
					if float(dataproj[0])<=median[2]: n_poor_median+=1
					if float(dataproj[0])<=mode[2]: n_poor_mode+=1


		if i<10: Nrich_median.append(n_rich_median); Nrich_mode.append(n_rich_mode)
		elif i<19: Nmid_median.append(n_mid_median); Nmid_mode.append(n_mid_mode)
		elif i<25: Npoor_median.append(n_poor_median); Npoor_mode.append(n_poor_mode)

		print i
		
	Nrich_median=np.array(Nrich_median); Nmid_median=np.array(Nmid_median); Npoor_median=np.array(Npoor_median)

	print np.mean(Nrich_median)/(median[0])**3, np.mean(Nmid_median)/(median[1])**3, np.mean(Npoor_median)/(median[2])**3


	return Nrich_median, Nmid_median, Npoor_median, Nrich_mode, Nmid_mode, Npoor_mode



def find_DWD_MS_Triple(pathlist, start, end):
	sourcedir=np.genfromtxt(pathlist, dtype=str)
	Ntriple=[]
	for i in range(len(sourcedir)):
		filestr=sourcedir[i]+'/initial'
		h=scripts3.read_binint(filestr+'.binint.log')
		ntri=0
		for j in range(len(h)):
			loutoput=len(h[j]['output'])
			for k in range(loutoput):
				if h[j]['output'][k]['type']=='triple':
					#print 'yes'
					mtot=float(h[j]['output'][k]['m'][0])+float(h[j]['output'][k]['m'][1])
					kin0=int(re.findall(r'\d+', h[j]['output'][k]['ktype'][0])[0])
					kin1=int(re.findall(r'\d+', h[j]['output'][k]['ktype'][1])[0])
					kout=int(re.findall(r'\d+', h[j]['output'][k]['ktype'][2])[0])
					if (kin0>=10 and kin0<=12) and (kin1>=10 and kin1<=12):
						print 'DWD'
					 	if mtot>=1.4:
							print kin0, kin1, kout
							ntri+=1


		Ntriple.append(ntri)
		print i
	
	print Ntriple
				

















		
