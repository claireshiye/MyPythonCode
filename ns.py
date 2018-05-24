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
import seaborn as sns
import gzip
import math
import re
import history_maker_full4 as hi4
import history_maker_full5 as hi5
import history_cmc as hic
import dynamics as dyn


#path = '/projects/b1011/kyle/cmc/branches/cmc-mpi/m22_project/rundir/m10'
#fo = 'm10_400000e5_4.5_1.0_0.05_FULL'
#fname = 'initial.pulsars.dat'
yearsc=31557600
twopi=6.283185307179586


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
	Bs=K*np.sqrt(np.abs(Pdots)*np.array(Ps)*0.001)
	Bb=K*np.sqrt(np.abs(Pdotb)*np.array(Pb)*0.001)
	Ps=np.array(Ps); Pb=np.array(Pb); Bs=np.array(Bs); Bb=np.array(Bb)
	#print Bs, Bb
	return Ps, Pb, Bs, Bb


def conv_dict(): return {'l':15, 't':19}    #?


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
        	if binflag[j] != 1 and st[j]  == 13:
			id_temp.append(int(star_id[j])); m_temp.append(star_m[j])#; t_temp.append(get_time(snaps[i])*t_conv/1000)
        	if binflag[j] == 1 and st1[j] == 13:
			id0_temp.append(int(id_0[j])); m0_temp.append(mass0[j])#; t0_temp.append(get_time(snaps[i])*t_conv/1000)
        	if binflag[j] == 1 and st2[j] == 13:
			id1_temp.append(int(id_1[j])); m1_temp.append(mass1[j])#; t1_temp.append(get_time(snaps[i])*t_conv/1000)

	M_temp = np.concatenate((m_temp, m0_temp, m1_temp), axis=0); ID_temp = np.concatenate((id_temp, id0_temp, id1_temp), axis=0)

	return ID_temp, M_temp


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

	
def get_snap_Nns(snapshot):
	nmtb=0; nmsp=0; npulsar=0
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
					spin=twopi*yearsc/float(datasnap[59])
					deathcut=(spin**2)*(0.17*10**12)
					if deathcut<float(datasnap[60]): npulsar+=1
					if spin<=0.02: 
						nmsp+=1; mspid.append(int(datasnap[0])); mspcomid.append(-100)
						print int(datasnap[0]), -100
			if int(datasnap[7])==1:
				if int(datasnap[17])==13:
					spin0=twopi*yearsc/float(datasnap[45])
					deathcut0=(spin0**2)*(0.17*10**12)
					if deathcut0<float(datasnap[47]): npulsar+=1
					if float(datasnap[44])>=1: nmtb+=1	
					if spin0<=0.02: 
						nmsp+=1; mspid.append(int(datasnap[10])); mspcomid.append(int(datasnap[11]))
						print int(datasnap[10]), int(datasnap[11])
				if int(datasnap[18])==13:
					spin1=twopi*yearsc/float(datasnap[46])
					deathcut1=(spin1**2)*(0.17*10**12)
					if deathcut1<float(datasnap[48]): npulsar+=1
					if float(datasnap[43])>=1: nmtb+=1
					if spin1<=0.02: 
						nmsp+=1; mspid.append(int(datasnap[11])); mspcomid.append(int(datasnap[10]))
						print int(datasnap[11]), int(datasnap[10])

	return npulsar, nmsp, nmtb, mspid, mspcomid


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
		Npulsar, Nmsp, Nmtb, Mspid, Mspcomid=get_snap_Nns(lastsnap)
		MSPID=MSPID+Mspid; MSPCOMID=MSPCOMID+Mspcomid
		temp=list(np.full_like(Mspid, fill_value=i))
		model=model+temp
		print i

	print len(model), len(MSPID)

	#np.savetxt('/projects/b1011/syr904/projects/PULSAR/kickgrid_runs/kickgrid_mspid.dat', np.c_[model, MSPID, MSPCOMID], fmt='%d %d %d', header='Model ID0 ID1', comments='#', delimiter= ' ')


def get_allmodel_MSPID_10Gyr(pathlist, start, end):
	sourcedir=np.genfromtxt(pathlist, dtype=str)
	MSPID=[]; MSPCOMID=[]; model=[]
	for i in range(start, end):
		filepath=sourcedir[i]
		pref='initial'
		filestr=filepath+'/'+pref
		t_conv=conv('t', filestr+'.conv.sh')
		snaps=np.sort(glob(filestr+'.snap*.dat.gz'))
		for j in range(len(snaps)-1, 0, -1):
			t=get_time(snaps[j])*t_conv/1000.
			print j
			if t<10.:
				Npulsar, Nmsp, Nmtb, Mspid, Mspcomid=get_snap_Nns(snaps[j+1])
				MSPID=MSPID+Mspid; MSPCOMID=MSPCOMID+Mspcomid
				temp=list(np.full_like(Mspid, fill_value=i))
				model=model+temp
				break

		print i

	print len(model), len(MSPID)

	np.savetxt('/projects/b1011/syr904/projects/PULSAR/kickgrid_runs/kickgrid_mspid_10Gyr.dat', np.c_[model, MSPID, MSPCOMID], fmt='%d %d %d', header='Model ID0 ID1', comments='#', delimiter= ' ')



def plot_Nbh_Npulsar(start, end, pathlist):
	path=np.genfromtxt(pathlist, dtype='|S')
	NBH=[]; NP=[]; NMSP=[]; NMT=[]; NTOT=[]
	for k in range(start, end):
		print path[k]
		snaps=np.sort(glob(path[k]+'/'+'initial.snap*.dat.gz'))
		lastsnap=snaps[-1]
		pref='initial'
		filestr=path[k]+'/'+pref
		Nbh, Ntot=dyn.find_NBH_NTOT_last(filestr)
		Npuls, Nmsp, Nmt, MSPid=get_snap_Nns(lastsnap)
		NBH.append(Nbh); NP.append(Npuls); NMSP.append(Nmsp); NMT.append(Nmt); NTOT.append(Ntot)
	
	NBH=np.array(NBH); NTOT=np.array(NTOT); NP=np.array(NP); NMSP=np.array(NMSP); NMT=np.array(NMT)	
	print NBH, NTOT, NMSP, NMT
	np.savetxt('/projects/b1011/syr904/projects/PULSAR/kickgrid_runs/rvgrid/ns_number.dat', np.c_[NBH, NTOT, NMSP, NMT], fmt ='%d %d %d %d', delimiter= ' ', header = 'Nbh Ntot Nmsp Nnsmt', comments = '#')

	plt.figure()
	#plt.scatter(NBH, NP, color='b', label=r'$N_{pulsars}$')
	#plt.scatter(NBH, NMSP, color='orange', label=r'$N_{MSP}$')
	plt.scatter(NBH, NMSP, color='green', s=100, alpha=0.7)
	plt.xscale('log')
	#plt.yscale('symlog')
	plt.ylabel(r'$N_{MSP}$')
	plt.xlabel(r'$N_{BH}$')
	plt.legend()
	plt.show()


def plot_Nbh_Npulsar_10to12Gyr(start, end, pathlist):
	path=np.genfromtxt(pathlist, dtype='|S')
	NBH=[]; NP=[]; NMSP=[]; NMT=[]; NTOT=[]
	for k in range(start, end):
		print path[k]
		snaps=np.sort(glob(path[k]+'/'+'initial.snap*.dat.gz'))
		pref='initial'
		filestr=path[k]+'/'+pref
		t_conv=conv('t',filestr+'.conv.sh')
		for j in range(len(snaps)-1, 0, -1):
			t=get_time(snaps[j])
			Time=t_conv*t/1000.
			if Time>=10.0:
				Nbh, Ntot=dyn.find_NBH_NTOT(filestr, t)
				Npuls, Nmsp, Nmt, MSPid=get_snap_Nns(snaps[j])
				NBH.append(Nbh); NP.append(Npuls); NMSP.append(Nmsp); NMT.append(Nmt); NTOT.append(Ntot)

			if Time<10.0: break
	
	NBH=np.array(NBH); NTOT=np.array(NTOT); NP=np.array(NP); NMSP=np.array(NMSP); NMT=np.array(NMT)	
	print NBH, NTOT, NMSP, NMT
	np.savetxt('/projects/b1011/syr904/projects/PULSAR/kickgrid_runs/ns_number_10to12Gyr.dat', np.c_[NBH, NTOT, NMSP, NMT], fmt ='%d %d %d %d', delimiter= ' ', header = 'Nbh Ntot Nmsp Nnsmt', comments = '#')

	plt.figure()
	plt.scatter(np.log(NBH/NTOT), np.log(NMSP/NTOT), color='b')
	#plt.xscale('log')
	#plt.yscale('log')
	plt.xlabel(r'$ln\ N_{BH}/N_{TOT}$')
	plt.ylabel(r'$ln\ N_{MSP}/N_{TOT}$')
	plt.savefig('/projects/b1011/syr904/projects/PULSAR/kickgrid_runs/nbhnmsp_10to12Gyr.pdf', dpi=300)



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
	for i in range(start, end):
		pref='initial'
		filestr=sourcedir[i]+'/'+pref
		snaps=np.sort(glob(filestr+'.snap*.dat.gz'))
		fhandle=open(filestr+'.ns.dat', 'w+', 0)
		fhandle.write('#1.Totaltime, 2.Nns,tot, 3.Nns,single, 4.Nns,binary, 5.Nns,mtb, 6.Npulsar, 7.Nmsp, 8.Nns-ns, 9.Nns-bh, 10.Nns-wd, 11.Nns-ms, 12.Nns-postms\n')
		for j in range(len(snaps)):
			N_NS=0; N_NS_SIN=0; N_NS_BIN=0; N_NS_MTB=0; N_PULS=0; N_MSP=0; N_NSNS=0; N_NSBH=0; N_NSWD=0; N_NSMS=0; N_NSPOSTMS=0
			T=get_time(snaps[j])
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
							if spin<=0.02: N_MSP+=1
					if int(datasnap[7])==1:
						if int(datasnap[17])==13:
							N_NS+=1; N_NS_BIN+=1
							spin0=twopi*yearsc/float(datasnap[45])
							deathcut0=(spin0**2)*(0.17*10**12)
							if deathcut0<float(datasnap[47]): N_PULS+=1
							if float(datasnap[44])>=1: N_NS_MTB+=1	
							if spin0<=0.02: N_MSP+=1

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
							if spin1<=0.02: N_MSP+=1

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

	
		

 
def find_primordialbin(sourcedir, folder):   #Find NS Primordial binaries
	yearsc=31557600
	twopi=6.283185307179586
	t_conv = conv('t',sourcedir+'/'+folder+'/'+'initial.conv.sh')

	data = np.genfromtxt(sourcedir+'/'+folder+'/'+'initial.snap0000.dat.gz', usecols=(7,10,11))
	binflag=data[:,0]
	id0=[]; id1=[]
	for i in range(len(binflag)):
		if binflag[i]==1:
			id0.append(data[:,1][i]); id1.append(data[:,2][i])

	snaps=np.sort(glob(sourcedir+'/'+folder+'/'+'initial.snap*.dat.gz'))
	for j in range(len(id0)):
		x0=id0[j]; x1=id1[j]
		t=[]; m0=[]; m1=[]; i0=[]; i1=[]; a=[]; e=[]; k0=[]; k1=[]; P0=[]; P1=[]; B0=[]; B1=[]
		for k in range(len(snaps)):
			mass0=[]; mass1=[]; sid0=[]; sid1=[]; semima=[]; ecc=[]; sk0=[]; sk1=[]; spin0=[]; spin1=[]; field0=[]; field1=[]
			with gzip.open(snaps[k],'r') as fi:
				for _ in xrange(2):
					next(fi)
				for line in fi:
					data1=line.split()
					if int(data1[7])==1:
						mass0.append(float(data1[8])); mass1.append(float(data1[9])); sid0.append(int(data1[10])); sid1.append(int(data1[11])); semima.append(float(data1[12])); ecc.append(float(data1[13])); sk0.append(int(data1[17])); sk1.append(int(data1[18])); spin0.append(float(data1[46])); spin1.append(float(data1[47])); field0.append(float(data1[48])); field1.append(float(data1[49]))
			
			time=t_conv*get_time(snaps[k])
			for l in range(len(mass0)):
				if (sid0[l]==x0 and sid1[l]==x1) or (sid0[l]==x1 and sid1[l]==x0):
					m0.append(mass0[l]); m1.append(mass1[l]); i0.append(sid0[l]); i1.append(sid1[l]); a.append(semima[l]); e.append(ecc[l]); k0.append(sk0[l]); k1.append(sk1[l]); P0.append(spin0[l]); P1.append(spin1[l]); B0.append(field0[l]); B1.append(field1[l]); t.append(time)
					break

		y0=k0[-1]; y1=k1[-1]
		if y0==13 or y1==13:
			name=str(x0)+'_'+str(x1)
			np.savetxt(sourcedir+'/'+folder+'/'+'history'+'/'+name+'_primordial.dat', np.c_[t, m0, m1, i0, i1, a, e, k0, k1, P0, P1, B0, B1], fmt ='%f %f %f %d %d %f %f %d %d %f %f %e %e', delimiter= ' ', header = '1.time, 2.m0, 3.m1, 4.id0, 5.id1, 6.a, 7.e, 8.k0, 9.k1, 10.P0, 11.P1, 12.B0, 13.B1', comments = '#')
			#print x0, x1



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



def get_id_BP(sourcedir, theid, savepath):
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
                                        		Age.append(float(time)); B.append(float(datasnap[60])); P.append((twopi*yearsc)/float(datasnap[59])); FC.append(-100); m0.append(float(datasnap[1])); m1.append(-100)
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
	np.savetxt(savepath+'/'+str(theid)+'.dat', np.c_[Age, B, P, FC, id0, id1, m0, m1, k0, k1, radrol0, radrol1, a, ecc], fmt ='%f %e %f %d %d %d %f %f %d %d %f %f %f %f', delimiter= ' ', header = 'Age(Myr) B(G) P(sec) FC id0 id1 m0 m1 k0 k1 radrol0 radrol1 a[AU] ecc', comments = '#')
	#return B, P, FC, Age



def get_allid_BP(sourcedir, idfile, savepath):
	dataid=np.genfromtxt(idfile, dtype=int)
	model=dataid[:,0]; id0=dataid[:,1]; id1=dataid[:,2]
	for i in range(0,len(model)):
		if model[i]==23:
			get_id_BP(sourcedir, id0[i], savepath)
			print id0[i]
				


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
					check=1
				if int(data[11])==theid:
					r=float(data[2]); bfield=float(data[48]); spin=twopi*yearsc/float(data[46])
					check=1
			if int(data[7])!=1:
				if int(data[0])==theid:
					r=float(data[2]); bfield=float(data[60]); spin=twopi*yearsc/float(data[59])
					check=1

			if check==1: break

	return r, bfield, spin


def get_id_allmodel_position(idlist, pathlist):
	dataid=np.genfromtxt(idlist)
	model=dataid[:,0]; id0=dataid[:,1]; id1=dataid[:,2]
	sourcedir=np.genfromtxt(pathlist, dtype=str)
	rposition=[]; rc=[]; B=[]; P=[]; Nbh=[]; Ntot=[]
	for i in range(len(model)):
		num=int(model[i])
		filepath=sourcedir[num]
		pref='initial'
		filestr=filepath+'/'+pref
		snaps=np.sort(glob(filestr+'.snap*.dat.gz'))
		lastsnap=snaps[-1]

		l_conv = conv('l', filestr+'.conv.sh')

		r_temp, b, s=get_id_position_BP(id0[i], lastsnap)
		rposition.append(r_temp*l_conv); B.append(b); P.append(s)

		datadyn=np.genfromtxt(filestr+'.dyn.dat')
		rc_tempt=datadyn[-1][7]
		rc.append(rc_tempt*l_conv)

		nbh, ntot=dyn.find_NBH_NTOT_last(filestr)
		Nbh.append(nbh); Ntot.append(ntot)

		print model[i]

	np.savetxt('/projects/b1011/syr904/projects/PULSAR/kickgrid_runs/kickgrid_msp1.dat', np.c_[model, id0, id1, rposition, B, P, rc, Nbh, Ntot], fmt='%d %d %d %f %e %f %f %d %d', header='Model ID0 ID1 r(pc) B(G) P(sec) rc(pc) Nbh Ntot', comments='#', delimiter= ' ')

	return rposition, rc, Nbh, Ntot


def get_id_allmodel_position_10Gyr(idlist, pathlist):
	dataid=np.genfromtxt(idlist)
	model=dataid[:,0]; id0=dataid[:,1]; id1=dataid[:,2]
	sourcedir=np.genfromtxt(pathlist, dtype=str)
	rposition=[]; rc=[]; B=[]; P=[]; Nbh=[]; Ntot=[]
	for i in range(len(model)):
		num=int(model[i])
		filepath=sourcedir[num]
		pref='initial'
		filestr=filepath+'/'+pref
		snaps=np.sort(glob(filestr+'.snap*.dat.gz'))
		l_conv = conv('l', filestr+'.conv.sh')
		t_conv = conv('t', filestr+'.conv.sh')
		datadyn=np.genfromtxt(filestr+'.dyn.dat')
		databh=np.genfromtxt(filestr+'.bh.dat')
		for j in range(len(snaps)-1, 0, -1):
			t=get_time(snaps[j])*t_conv/1000.
			print j
			if t<10.:
				t10=get_time(snaps[j+1])
				print t10
				r_temp, b, s=get_id_position_BP(id0[i], snaps[j+1])
				rposition.append(r_temp*l_conv); B.append(b); P.append(s)

				for k in range(len(datadyn[:,7])):
					tdyn=datadyn[:,0][k]
					if tdyn==t10:
						rc_tempt=datadyn[:,7][k]
						rc.append(rc_tempt*l_conv)
						n_tot=datadyn[:,3][k]

				for h in range(len(databh[:,2])):
					tbh=databh[:,1][h]
					if tbh==t10:
						n_bh=databh[:,2][h]

				Nbh.append(n_bh); Ntot.append(n_tot)
				print "yes"
				break

		print model[i]

	np.savetxt('/projects/b1011/syr904/projects/PULSAR/kickgrid_runs/kickgrid_msp10Gyr.dat', np.c_[model, id0, id1, rposition, B, P, rc, Nbh, Ntot], fmt='%d %d %d %f %e %f %f %d %d', header='Model ID0 ID1 r(pc) B(G) P(sec) rc(pc) Nbh Ntot', comments='#', delimiter= ' ')

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
					if spin>0.03 and float(datasnap[60])>deathcut:
						time.append(t); id0.append(int(datasnap[0])); id1.append(-100); m0.append(float(datasnap[1])); m1.append(-100)
						k0.append(int(datasnap[14])); k1.append(-100); FC.append(-100); B.append(float(datasnap[60])); P.append(spin)
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


def get_normalpsr_finl(filepath, savepath, modelno):
	time=[]; id0=[]; id1=[]; m0=[]; m1=[]; k0=[]; k1=[]; B=[]; P=[]; a=[]; ecc=[]; FC=[]; model=[]
	pref='initial'
	filestr=filepath+'/'+pref
	t_conv=conv('t', filestr+'.conv.sh')
	snaps=np.sort(glob(filestr+'.snap*.dat.gz'))
	lastsnap=snaps[-1]
	t=get_time(lastsnap)*t_conv

	Model, Time, Bf, Spin, Fc, Id0, Id1, M0, M1, K0, K1, A, Ecc=find_normalpsr(lastsnap, modelno, t)


	np.savetxt(savepath+'/kickgrid_normalpsr'+str(modelno)+'.dat', np.c_[Model, Time, Bf, Spin, Fc, Id0, Id1, M0, M1, K0, K1, A, Ecc], fmt ='%d %f %e %f %d %d %d %f %f %d %d %f %f', delimiter= ' ', header = 'Model Age(Myr) B(G) P(sec) FC id0 id1 m0 m1 k0 k1 a[AU] ecc', comments = '#')


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




		






		
