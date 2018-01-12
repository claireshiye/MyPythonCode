import numpy as np
from glob import glob
import collections
from collections import Counter
import os,sys
import subprocess
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
import gzip
import math
import re
import history_maker_full5


#path = '/projects/b1011/kyle/cmc/branches/cmc-mpi/m22_project/rundir/m10'
#fo = 'm10_400000e5_4.5_1.0_0.05_FULL'
#fname = 'initial.pulsars.dat'

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
						if int(data1[7])==0:
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
        	if binflag[j] == 0 and st[j]  == 13:
			id_temp.append(int(star_id[j])); m_temp.append(star_m[j])#; t_temp.append(get_time(snaps[i])*t_conv/1000)
        	if binflag[j] == 1 and st1[j] == 13:
			id0_temp.append(int(id_0[j])); m0_temp.append(mass0[j])#; t0_temp.append(get_time(snaps[i])*t_conv/1000)
        	if binflag[j] == 1 and st2[j] == 13:
			id1_temp.append(int(id_1[j])); m1_temp.append(mass1[j])#; t1_temp.append(get_time(snaps[i])*t_conv/1000)

	M_temp = np.concatenate((m_temp, m0_temp, m1_temp), axis=0); ID_temp = np.concatenate((id_temp, id0_temp, id1_temp), axis=0)

	return ID_temp, #M_temp




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

	



def spin_PPplot(sourcedir):   #Generate spin distribution and PP plot.
	allspin=[]; sispin=[]; bispin=[]; unreal=[]; Bs=[]; Bb=[]; Bsi=[]; Bbi=[]; Psi=[]; Pbi=[]
	with open(sourcedir+'/'+'end1e6.dat') as f:
		#next(f)
		for line in f:
			data=line.split()
			#print binflag
			if int(data[0])==135549:
				if int(data[8])==0:
					allspin.append(float(data[6]))
					sispin.append(float(data[6]))
					Bs.append(float(data[5]))
					#print float(data[6])
					#if float(data[6])<=0.001:
						#unreal.append(float(data[6]))
				else:
					if int(data[19])==13:
						allspin.append(float(data[17]))
						bispin.append(float(data[17]))
						Bb.append(float(data[15]))
						#print float(data[17])
						#if float(data[17])<=0.001:
							#unreal.append(float(data[17]))
					if int(data[20])==13:
						allspin.append(float(data[18]))
						bispin.append(float(data[18]))
						Bb.append(float(data[16]))
						#print float(data[18])
						#if float(data[18])<=0.001:
							#unreal.append(float(data[18]))

	#with open(sourcedir+'/'+folder+'/'+'head1e6.dat') as fi:
	#	next(fi)
	#	for line in fi:
	#		datai=line.split()
	#		if int(datai[0])==50:
	#			if int(datai[8])==0:
	#				Psi.append(float(datai[6]))
	#				Bsi.append(float(datai[5]))
	#			else:
	#				if int(datai[19])==13:
	#					Pbi.append(float(datai[17]))
	#					Bbi.append(float(datai[15]))
	#				if int(data[20])==13:
	#					Pbi.append(float(datai[18]))
	#					Bbi.append(float(datai[16]))


	#allspin=np.asarray(allspin)
	#plt.figure()
	#plt.xscale('log')
	#plt.xlabel(r'$SPIN(s)$')
	#plt.hist(allspin, bins=10**np.linspace(np.log10(min(allspin)), np.log10(max(allspin)), 25))
	#plt.savefig('6.5838288_spin.pdf')

	#sispin=np.asarray(sispin); bispin=np.asarray(bispin); unreal=np.asarray(unreal)
	#plt.figure()
	#plt.xscale('log')
	#plt.xlabel(r'$SPIN(s)$')
	#plt.hist(sispin, bins=10**np.linspace(np.log10(min(sispin)), np.log10(max(sispin)), 25), histtype='step', color='b', label='single')
	#plt.hist(bispin, bins=10**np.linspace(np.log10(min(bispin)), np.log10(max(bispin)), 25), histtype='step', alpha=0.5, color='r', label='binary')
	##plt.hist(unreal, bins=10**np.linspace(np.log10(min(unreal)), np.log10(max(unreal)), 25), histtype='step', alpha=0.5, color='k')
	#plt.legend()
	#plt.savefig('NS'+'/'+'7.5538408_spin.pdf')

	#Death Line
	x=np.logspace(-4.0, 2.0, num=50)

	Ps=np.asarray(sispin); Pb=np.asarray(bispin); Bs=np.asarray(Bs); Bb=np.asarray(Bb); Bsi=np.asarray(Bsi); Bbi=np.asarray(Bbi); Psi=np.asarray(Psi); Pbi=np.asarray(Pbi)
	plt.figure()
	plt.xscale('log')
	plt.yscale('log')
	plt.plot(x, (x**2)*(0.17*10**12), 'k--')    #Deadline
	#plt.plot(x, (4*10**8)*(x/(8.6*10**-4))**(7./6.), 'r--')   #Spin-up line
	#plt.plot(x, x*0.12*10**12, 'b--')  #Artificial line
	#plt.plot(x, x*0.15*10**12, 'g--')  #Artificial line
	plt.plot(x, x*0.04*10**12, 'r--')  #Artificial line
	#plt.plot(x, 10**(1.81*np.log10(x)+6.9), 'k')   #Eq.(18)
	#plt.plot(x, 10**(1.75*np.log10(x)-11.64), 'k.')   #Eq.(19)
	#plt.scatter(Psi, Bsi, marker='.', color='k', label='15 Myr')
	#plt.scatter(Pbi, Bbi, marker='.', color='k')
	plt.scatter(Ps, Bs, color='purple', label='single')
	plt.scatter(Pb, Bb, color='orange', label='binary', alpha=0.7)
	plt.xlim(10**-4, 100.)
	plt.ylim(10**7, 10**15)
	plt.xlabel(r'$P(sec)$')
	plt.ylabel(r'$B(G)$')
	plt.legend(loc='upper left')
	plt.savefig(sourcedir+'/'+'NS'+'/'+'PP_9.3158381.pdf')
	
	#return max(Bs), max(Bb)


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
			if binflagi[j]==1:
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






							
		


	
				
