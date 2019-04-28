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
import random
import history_maker_full4 as hi4
import history_maker_full5 as hi5
import history_cmc as hic
import dynamics as dyn
import scripts1
import scripts2
import scripts3
import LISA_calculations as lisa
import ecc_calc as gwcalc

savepath='/projects/b1095/syr904/projects/SGRB'

##Colllecting useful models from the set of new models
def find_models():
	allmodels=np.genfromtxt('/projects/b1095/syr904/projects/SGRB/newruns/path_final.dat', dtype=str)
	status=allmodels[:,1]; paths=allmodels[:,0]

	sourcedir=[]
	for i in range(len(paths)):
		if status[i]=='2' or status[i]=='3': sourcedir.append(paths[i])

	np.savetxt(savepath+'/newruns/path_newruns.dat', np.c_[sourcedir], fmt='%s')



##Find any mergers in cluster
def find_merger_nsns(pathlist, start, end):
	sourcedir=np.genfromtxt(pathlist, dtype=str)
	model_coll=[]; model_merge=[]

	idcoll=[]; id0coll=[]; id1coll=[]; id2coll=[]; id3coll=[]
	idsemerge=[]; id0merge=[]; id1merge=[]
	typecoll=[]; k0=[]; k1=[]; k2=[]; k3=[]
	timecoll=[]; timesemerge=[]
	timecoll_myr=[]; timesemerge_myr=[]
	mf_merge=[]; m0_merge=[]; m1_merge=[]
	mf_coll=[]; m0_coll=[]; m1_coll=[]; m2_coll=[]; m3_coll=[]
	r_merge=[]; r_coll=[]


	for i in range(start, end):
		filestr=sourcedir[i]+'/initial'
		collfile=filestr+'.collision.log'
		binintfile=filestr+'.binint.log'
		semergefile=filestr+'.semergedisrupt.log'

		t_conv=dyn.conv('t', filestr+'.conv.sh')

		ncoll2=0; ncoll3=0; ncoll4=0
		nsemerge=0

		##Check in-cluster merger in the collision file
		colldata=scripts1.readcollfile(collfile)
		for j in range(len(colldata)):
			line=colldata[j].split()
			if int(line[2])==2:  ##Single-single star collision
				if int(line[11])==13 and int(line[12])==13:
					model_coll.append(i)
					ncoll2+=1
					idcoll.append(int(line[3])); id0coll.append(int(line[5])); id1coll.append(int(line[7])); id2coll.append(-100); id3coll.append(-100)
					typecoll.append(int(line[2]))
					timecoll.append(float(line[0])); timecoll_myr.append(t_conv*float(line[0]))
					mf_coll.append(float(line[4])); m0_coll.append(float(line[6])); m1_coll.append(float(line[8])); m2_coll.append(-100); m3_coll.append(-100)
					r_coll.append(float(line[9]))
					k0.append(int(line[11])); k1.append(int(line[12])); k2.append(-100); k3.append(-100)
					

			if int(line[2])==3:  ##Binary-single star collision
				if len(line)==16:  ##Three stars collision
					if (int(line[13])==13 and int(line[14])==13) or (int(line[14])==13 and int(line[15])==13) or (int(line[13])==13 and int(line[15])==13):
						model_coll.append(i)
						ncoll3+=1
						idcoll.append(int(line[3])); id0coll.append(int(line[5])); id1coll.append(int(line[7])); id2coll.append(int(line[9])); id3coll.append(-100)
						typecoll.append(int(line[2]))
						timecoll.append(float(line[0])); timecoll_myr.append(t_conv*float(line[0]))
						mf_coll.append(float(line[4])); m0_coll.append(float(line[6])); m1_coll.append(float(line[8])); m2_coll.append(float(line[10])); m3_coll.append(-100)
						r_coll.append(float(line[11]))
						k0.append(int(line[13])); k1.append(int(line[14])); k2.append(int(line[15])); k3.append(-100)
						
						
				if len(line)==13: ##Two stars collision
					if int(line[11])==13 and int(line[12])==13:
						model_coll.append(i)
						ncoll3+=1
						idcoll.append(int(line[3])); id0coll.append(int(line[5])); id1coll.append(int(line[7])); id2coll.append(-100); id3coll.append(-100)
						typecoll.append(int(line[2]))
						timecoll.append(float(line[0])); timecoll_myr.append(t_conv*float(line[0]))
						mf_coll.append(float(line[4])); m0_coll.append(float(line[6])); m1_coll.append(float(line[8])); m2_coll.append(-100); m3_coll.append(-100)
						r_coll.append(float(line[9]))
						k0.append(int(line[11])); k1.append(int(line[12])); k2.append(-100); k3.append(-100)


			if int(line[2])==4:  ##Binary-binary star collision
				if len(line)==19: ##Four stars collision
					if (int(line[15])==13 and int(line[16])==13) or (int(line[16])==13 and int(line[17])==13) or (int(line[17])==13 and int(line[18])==13) or (int(line[15])==13 and int(line[17])==13) or (int(line[15])==13 and int(line[18])==13) or (int(line[16])==13 and int(line[18])==13):
						model_coll.append(i)
						ncoll4+=1
						idcoll.append(int(line[3])); id0coll.append(int(line[5])); id1coll.append(int(line[7])); id2coll.append(int(line[9])); id3coll.append(int(line[11]))
						typecoll.append(int(line[2]))
						timecoll.append(float(line[0])); timecoll_myr.append(t_conv*float(line[0]))
						mf_coll.append(float(line[4])); m0_coll.append(float(line[6])); m1_coll.append(float(line[8])); m2_coll.append(float(line[10])); m3_coll.append(float(line[12]))
						r_coll.append(float(line[13]))
						k0.append(int(line[15])); k1.append(int(line[16])); k2.append(int(line[17])); k3.append(int(line[18]))


				if len(line)==16:  ##Three stars collision
					if (int(line[13])==13 and int(line[14])==13) or (int(line[14])==13 and int(line[15])==13) or (int(line[13])==13 and int(line[15])==13):
						model_coll.append(i)
						ncoll4+=1
						idcoll.append(int(line[3])); id0coll.append(int(line[5])); id1coll.append(int(line[7])); id2coll.append(int(line[9])); id3coll.append(-100)
						typecoll.append(int(line[2]))
						timecoll.append(float(line[0])); timecoll_myr.append(t_conv*float(line[0]))
						mf_coll.append(float(line[4])); m0_coll.append(float(line[6])); m1_coll.append(float(line[8])); m2_coll.append(float(line[10])); m3_coll.append(-100)
						r_coll.append(float(line[11]))
						k0.append(int(line[13])); k1.append(int(line[14])); k2.append(int(line[15])); k3.append(-100)
						
						
				if len(line)==13: ##Two stars collision
					if int(line[11])==13 and int(line[12])==13:
						model_coll.append(i)
						ncoll4+=1
						idcoll.append(int(line[3])); id0coll.append(int(line[5])); id1coll.append(int(line[7])); id2coll.append(-100); id3coll.append(-100)
						typecoll.append(int(line[2]))
						timecoll.append(float(line[0])); timecoll_myr.append(t_conv*float(line[0]))
						mf_coll.append(float(line[4])); m0_coll.append(float(line[6])); m1_coll.append(float(line[8])); m2_coll.append(-100); m3_coll.append(-100)
						r_coll.append(float(line[9]))
						k0.append(int(line[11])); k1.append(int(line[12])); k2.append(-100); k3.append(-100)


		print 'collfile:', ncoll2, ncoll3, ncoll4#, typecoll, timecoll
	
		##Check in-cluster merger in the semerge file
		semergedata=scripts2.readmergefile(semergefile)
		for k in range(len(semergedata)):
			line=semergedata[k].split()
			if int(line[1])<3:
				if int(line[-1])==13 and int(line[-2])==13:
					model_merge.append(i)
					nsemerge+=1
					timesemerge.append(float(line[0])); timesemerge_myr.append(t_conv*float(line[0]))
					idsemerge.append(int(line[2])); id0merge.append(int(line[4])); id1merge.append(int(line[6]))
					mf_merge.append(float(line[3])); m0_merge.append(float(line[5])); m1_merge.append(float(line[7]))
					r_merge.append(float(line[8]))


		print 'semerge:', nsemerge#, idsemerge, timesemerge

	np.savetxt(savepath+'/extrememodels/GWcap.dat', np.c_[model_coll, timecoll, timecoll_myr, typecoll, idcoll, id0coll, id1coll, id2coll, id3coll, mf_coll, m0_coll, m1_coll, m2_coll, m3_coll, r_coll, k0, k1, k2, k3], fmt='%d %f %f %d %d %d %d %d %d %f %f %f %f %f %f %d %d %d %d', delimiter='', header='1.Model 2.Time(code) 3.Time(Myr) 4.Type 5.IDM 6.ID0 7.ID1 8.ID2 9.ID3 10.MM 11.M0 12.M1 13.M2 14.M3 15.R 16.K0 17.K1 18.K2 19.K3', comments='#')

	np.savetxt(savepath+'/extrememodels/Incluster.dat', np.c_[model_merge, timesemerge, timesemerge_myr, idsemerge, id0merge, id1merge, mf_merge, m0_merge, m1_merge, r_merge], fmt='%d %f %f %d %d %d %f %f %f %f', delimiter='', header='1.Model 2.Time(code) 3.Time(Myr) 4.IDM 5.ID0 6.ID1 7.MM 8.M0 9.M1 10.R', comments='#')



def find_escmerger_nsns(pathlist, start, end):
	sourcedir=np.genfromtxt(pathlist, dtype=str)

	model_esc=[]; timeesc=[]; timeesc_myr=[]; tins=[]; m0=[]; m1=[]; id0=[]; id1=[]; a=[]; ecc=[]
	for i in range(start, end):
		filestr=sourcedir[i]+'/initial'
		escfile=filestr+'.esc.dat'

		t_conv=dyn.conv('t', filestr+'.conv.sh')

		nesc=0; nescmerger=0

		with open(escfile, 'r') as fesc:
			next(fesc)
			for line in fesc:
				dataesc=line.split()
				if int(dataesc[14])==1:
					if int(dataesc[22])==13 and int(dataesc[23])==13:
						nesc+=1
						t_inspiral=lisa.inspiral_time_peters(float(dataesc[19]), float(dataesc[20]), float(dataesc[15]), float(dataesc[16]))
						tesc=float(dataesc[1])*t_conv/1000.   ##In Gyr
						if t_inspiral+tesc<=14.:
							nescmerger+=1
							#print int(dataesc[17]), int(dataesc[18]), tesc, t_inspiral
							model_esc.append(i)
							timeesc_myr.append(float(dataesc[1])*t_conv); tins.append(t_inspiral*1000.)
							timeesc.append(float(dataesc[1])); m0.append(float(dataesc[15])); m1.append(float(dataesc[16])); id0.append(int(dataesc[17])); id1.append(int(dataesc[18])); a.append(float(dataesc[19])); ecc.append(float(dataesc[20]))


		print nesc, nescmerger		
		#print i

	np.savetxt(savepath+'/extrememodels/Escmerger.dat', np.c_[model_esc, timeesc, timeesc_myr, tins, m0, m1, id0, id1, a, ecc], fmt='%d %f %f %f %f %f %d %d %f %f', header='1.Model 2.Time_esc(code) 3.Time_esc(Myr) 4.T_insp(Myr) 5.M0 6.M1 7.ID0 8.ID1 9.A 10.ECC', delimiter='', comments='#')



def find_formationtime(ids, filepath): ##Find formation time of merging BNSs
	snaps=np.sort(glob(filepath+'*.snap*.dat.gz'))
	#t_conv=dyn.conv('t', filestring+'.conv.sh')
	tform=-100; snapno=-100
	for m in range(len(snaps)-1, 0, -1):
		check1=0; check2=0
		with gzip.open(snaps[m], 'r') as fsnap:
			next(fsnap)
			next(fsnap)
			for line in fsnap:
				datasnap=line.split()
				if (int(datasnap[10])==ids[0] and int(datasnap[11])==ids[1]) or (int(datasnap[10])==ids[1] and int(datasnap[11])==ids[0]) and (int(datasnap[17])==13 and int(datasnap[18])==13):
					tform=dyn.get_time(snaps[m])
					#tform=tform*t_conv  ##in Myr
					snapno=m
					check2=1

		with gzip.open(snaps[m-1], 'r') as fsnap:
			next(fsnap)
			next(fsnap)
			for line in fsnap:
				datasnap=line.split()
				if (int(datasnap[10])==ids[0] and int(datasnap[11])==ids[1]) or (int(datasnap[10])==ids[1] and int(datasnap[11])==ids[0]) and (int(datasnap[17])==13 and int(datasnap[18])==13):
					tform=dyn.get_time(snaps[m-1])
					#tform=tform*t_conv  ##in Myr
					snapno=m-1
					check1=0


		if check1==0 and check2==1:
			print 'snapno=', m 
			break

		print m
 
	return tform, snapno  ##in code unit



def find_all_mergers(pathlist, start, end):   ##Both in the cluster and escaped
	sourcedir=np.genfromtxt(pathlist, dtype=str)

	##Mergers in the clusters
	model_coll=[]; model_merge=[]

	idcoll=[]; id0coll=[]; id1coll=[]; id2coll=[]; id3coll=[]
	idsemerge=[]; id0merge=[]; id1merge=[]
	typecoll=[]; k0=[]; k1=[]; k2=[]; k3=[]
	timecoll=[]; timesemerge=[]
	timecoll_myr=[]; timesemerge_myr=[]
	mf_merge=[]; m0_merge=[]; m1_merge=[]
	mf_coll=[]; m0_coll=[]; m1_coll=[]; m2_coll=[]; m3_coll=[]
	r_merge=[]; r_coll=[]
	tform_semerge=[]; snapno_semerge=[]

	##Mergers outside of the clusters
	#model_esc=[]; timeesc=[]; timeesc_myr=[]; tins=[]; m0=[]; m1=[]; id0=[]; id1=[]; a=[]; ecc=[]
	#tform_esc=[]; snapno_esc=[]


	##Numbers
	Ncoll2=[]; Ncoll3=[]; Ncoll4=[]; Ncoll=[]
	Nsemerge=[]
	Nesc=[]; Nescmerge=[]
	Models=[]


	####For mergers that happen in the clusters####
	for i in range(start, end):
		filestr=sourcedir[i]+'/initial'
		collfile=filestr+'.collision.log'
		binintfile=filestr+'.binint.log'
		semergefile=filestr+'.semergedisrupt.log'

		t_conv=dyn.conv('t', filestr+'.conv.sh')

		ncoll2=0; ncoll3=0; ncoll4=0
		nsemerge=0
		nesc=0; nescmerger=0

		##Check in-cluster merger in the collision file
		colldata=scripts1.readcollfile(collfile)
		for j in range(len(colldata)):
			line=colldata[j].split()
			if int(line[2])==2:  ##Single-single star collision
				if int(line[11])==13 and int(line[12])==13:
					model_coll.append(i)
					ncoll2+=1
					idcoll.append(int(line[3])); id0coll.append(int(line[5])); id1coll.append(int(line[7])); id2coll.append(-100); id3coll.append(-100)
					typecoll.append(int(line[2]))
					timecoll.append(float(line[0])); timecoll_myr.append(t_conv*float(line[0]))
					mf_coll.append(float(line[4])); m0_coll.append(float(line[6])); m1_coll.append(float(line[8])); m2_coll.append(-100); m3_coll.append(-100)
					r_coll.append(float(line[9]))
					k0.append(int(line[11])); k1.append(int(line[12])); k2.append(-100); k3.append(-100)
					

			if int(line[2])==3:  ##Binary-single star collision
				if len(line)==16:  ##Three stars collision
					if (int(line[13])==13 and int(line[14])==13) or (int(line[14])==13 and int(line[15])==13) or (int(line[13])==13 and int(line[15])==13):
						model_coll.append(i)
						ncoll3+=1
						idcoll.append(int(line[3])); id0coll.append(int(line[5])); id1coll.append(int(line[7])); id2coll.append(int(line[9])); id3coll.append(-100)
						typecoll.append(int(line[2]))
						timecoll.append(float(line[0])); timecoll_myr.append(t_conv*float(line[0]))
						mf_coll.append(float(line[4])); m0_coll.append(float(line[6])); m1_coll.append(float(line[8])); m2_coll.append(float(line[10])); m3_coll.append(-100)
						r_coll.append(float(line[11]))
						k0.append(int(line[13])); k1.append(int(line[14])); k2.append(int(line[15])); k3.append(-100)
						
						
				if len(line)==13: ##Two stars collision
					if int(line[11])==13 and int(line[12])==13:
						model_coll.append(i)
						ncoll3+=1
						idcoll.append(int(line[3])); id0coll.append(int(line[5])); id1coll.append(int(line[7])); id2coll.append(-100); id3coll.append(-100)
						typecoll.append(int(line[2]))
						timecoll.append(float(line[0])); timecoll_myr.append(t_conv*float(line[0]))
						mf_coll.append(float(line[4])); m0_coll.append(float(line[6])); m1_coll.append(float(line[8])); m2_coll.append(-100); m3_coll.append(-100)
						r_coll.append(float(line[9]))
						k0.append(int(line[11])); k1.append(int(line[12])); k2.append(-100); k3.append(-100)


			if int(line[2])==4:  ##Binary-binary star collision
				if len(line)==19: ##Four stars collision
					if (int(line[15])==13 and int(line[16])==13) or (int(line[16])==13 and int(line[17])==13) or (int(line[17])==13 and int(line[18])==13) or (int(line[15])==13 and int(line[17])==13) or (int(line[15])==13 and int(line[18])==13) or (int(line[16])==13 and int(line[18])==13):
						model_coll.append(i)
						ncoll4+=1
						idcoll.append(int(line[3])); id0coll.append(int(line[5])); id1coll.append(int(line[7])); id2coll.append(int(line[9])); id3coll.append(int(line[11]))
						typecoll.append(int(line[2]))
						timecoll.append(float(line[0])); timecoll_myr.append(t_conv*float(line[0]))
						mf_coll.append(float(line[4])); m0_coll.append(float(line[6])); m1_coll.append(float(line[8])); m2_coll.append(float(line[10])); m3_coll.append(float(line[12]))
						r_coll.append(float(line[13]))
						k0.append(int(line[15])); k1.append(int(line[16])); k2.append(int(line[17])); k3.append(int(line[18]))


				if len(line)==16:  ##Three stars collision
					if (int(line[13])==13 and int(line[14])==13) or (int(line[14])==13 and int(line[15])==13) or (int(line[13])==13 and int(line[15])==13):
						model_coll.append(i)
						ncoll4+=1
						idcoll.append(int(line[3])); id0coll.append(int(line[5])); id1coll.append(int(line[7])); id2coll.append(int(line[9])); id3coll.append(-100)
						typecoll.append(int(line[2]))
						timecoll.append(float(line[0])); timecoll_myr.append(t_conv*float(line[0]))
						mf_coll.append(float(line[4])); m0_coll.append(float(line[6])); m1_coll.append(float(line[8])); m2_coll.append(float(line[10])); m3_coll.append(-100)
						r_coll.append(float(line[11]))
						k0.append(int(line[13])); k1.append(int(line[14])); k2.append(int(line[15])); k3.append(-100)
						
						
				if len(line)==13: ##Two stars collision
					if int(line[11])==13 and int(line[12])==13:
						model_coll.append(i)
						ncoll4+=1
						idcoll.append(int(line[3])); id0coll.append(int(line[5])); id1coll.append(int(line[7])); id2coll.append(-100); id3coll.append(-100)
						typecoll.append(int(line[2]))
						timecoll.append(float(line[0])); timecoll_myr.append(t_conv*float(line[0]))
						mf_coll.append(float(line[4])); m0_coll.append(float(line[6])); m1_coll.append(float(line[8])); m2_coll.append(-100); m3_coll.append(-100)
						r_coll.append(float(line[9]))
						k0.append(int(line[11])); k1.append(int(line[12])); k2.append(-100); k3.append(-100)


		print 'collfile:', ncoll2, ncoll3, ncoll4#, typecoll, timecoll
		Ncoll2.append(ncoll2); Ncoll3.append(ncoll3); Ncoll4.append(ncoll4)
		Ncoll.append(ncoll2+ncoll3+ncoll4)

	
		##Check in-cluster merger in the semerge file
		semergedata=scripts2.readmergefile(semergefile)
		for k in range(len(semergedata)):
			line=semergedata[k].split()
			if int(line[1])<3:
				if int(line[-1])==13 and int(line[-2])==13:
					model_merge.append(i)
					nsemerge+=1
					timesemerge.append(float(line[0])); timesemerge_myr.append(t_conv*float(line[0]))
					idsemerge.append(int(line[2])); id0merge.append(int(line[4])); id1merge.append(int(line[6]))
					mf_merge.append(float(line[3])); m0_merge.append(float(line[5])); m1_merge.append(float(line[7]))
					r_merge.append(float(line[8]))
					#tf_semerge, sno_semerge=find_formationtime([int(line[4]), int(line[6])], filestr)
					#if tf_semerge!=-100:
					#	tform_semerge.append(tf_semerge*t_conv); snapno_semerge.append(sno_semerge)
					#else:
					#	tform_semerge.append(t_conv*float(line[0])); snapno_semerge.append(-100)


		print 'semerge:', nsemerge#, idsemerge, timesemerge
		Nsemerge.append(nsemerge)

		#fescbns=open(savepath+'/extrememodels/Escmerger.dat', 'a+', 0)
		#fescbns.write('#1.Model 2.Time_esc(code) 3.Time_esc(Myr) 4.T_insp(Myr) 5.M0 6.M1 7.ID0 8.ID1 9.A 10.ECC 11.Tform(Myr) 12.Snapno\n')
		####For mergers that happen outside of the clusters####
		escfile=filestr+'.esc.dat'
		with open(escfile, 'r') as fesc:
			next(fesc)
			for line in fesc:
				dataesc=line.split()
				if int(dataesc[14])==1:
					if int(dataesc[22])==13 and int(dataesc[23])==13:
						nesc+=1
						#print float(dataesc[19]), float(dataesc[20]), float(dataesc[15]), float(dataesc[16])
						t_inspiral=gwcalc.t_inspiral_2(float(dataesc[19]), float(dataesc[20]), float(dataesc[15]), float(dataesc[16]), 0, 0, 0, 1100)/10**6 ##in Myr
						tesc=float(dataesc[1])*t_conv  ##In Myr
						#model_esc=i
						#timeesc_myr=float(dataesc[1])*t_conv; tins=t_inspiral
						#timeesc=float(dataesc[1]); m0=float(dataesc[15]); m1=float(dataesc[16]); id0=int(dataesc[17]); id1=int(dataesc[18]); a=float(dataesc[19]); ecc=float(dataesc[20])
						#tf_esc, sno_esc=find_formationtime([int(dataesc[17]), int(dataesc[18])], filestr)
						#if tf_esc!=-100:
						#	tform_esc=tf_esc*t_conv; snapno_esc=sno_esc
						#else:
						#	tform_esc=float(dataesc[1])*t_conv; snapno_esc=-100

						#fescbns.write('%d %f %f %f %f %f %d %d %f %f %f %d\n'%(model_esc, timeesc, timeesc_myr, tins, m0, m1, id0, id1, a, ecc, tform_esc, snapno_esc))


						if t_inspiral+tesc<=14000.:
							nescmerger+=1
							#print int(dataesc[17]), int(dataesc[18]), tesc, t_inspiral
							
							#model_esc=i
							#timeesc_myr=float(dataesc[1])*t_conv; tins=t_inspiral
							#timeesc=float(dataesc[1]); m0=float(dataesc[15]); m1=float(dataesc[16]); id0=int(dataesc[17]); id1=int(dataesc[18]); a=float(dataesc[19]); ecc=float(dataesc[20])
							#tf_esc, sno_esc=find_formationtime([int(dataesc[17]), int(dataesc[18])], filestr)
							#if tf_esc!=-100:
							#	tform_esc=tf_esc*t_conv; snapno_esc=sno_esc
							#else:
							#	tform_esc=float(dataesc[1])*t_conv; snapno_esc=-100

							#fescbns.write('%d %f %f %f %f %f %d %d %f %f %f %d\n'%(model_esc, timeesc, timeesc_myr, tins, m0, m1, id0, id1, a, ecc, tform_esc, snapno_esc))

		#fescbns.close()
		print 'escaped:', nesc, nescmerger
		Nesc.append(nesc); Nescmerge.append(nescmerger)


		Models.append(i)


	##Output files
	#np.savetxt('/projects/b1011/syr904/projects/PULSAR/kickgrid_runs/MSPBHinGC/data'+'/GWcap.dat', np.c_[model_coll, timecoll, timecoll_myr, typecoll, idcoll, id0coll, id1coll, id2coll, id3coll, mf_coll, m0_coll, m1_coll, m2_coll, m3_coll, r_coll, k0, k1, k2, k3], fmt='%d %f %f %d %d %d %d %d %d %f %f %f %f %f %f %d %d %d %d', delimiter='', header='1.Model 2.Time(code) 3.Time(Myr) 4.Type 5.IDM 6.ID0 7.ID1 8.ID2 9.ID3 10.MM 11.M0 12.M1 13.M2 14.M3 15.R 16.K0 17.K1 18.K2 19.K3', comments='#')

	#np.savetxt('/projects/b1011/syr904/projects/PULSAR/kickgrid_runs/MSPBHinGC/data'+'/Incluster.dat', np.c_[model_merge, timesemerge, timesemerge_myr, idsemerge, id0merge, id1merge, mf_merge, m0_merge, m1_merge, r_merge, tform_semerge, snapno_semerge], fmt='%d %f %f %d %d %d %f %f %f %f %f %d', delimiter='', header='1.Model 2.Time(code) 3.Time(Myr) 4.IDM 5.ID0 6.ID1 7.MM 8.M0 9.M1 10.R 11.Tform(Myr) 12.Snapno', comments='#')

	#np.savetxt(savepath+'/extrememodels/Escmerger.dat', np.c_[model_esc, timeesc, timeesc_myr, tins, m0, m1, id0, id1, a, ecc, tform_esc, snapno_esc], fmt='%d %f %f %f %f %f %d %d %f %f %f %d', header='1.Model 2.Time_esc(code) 3.Time_esc(Myr) 4.T_insp(Myr) 5.M0 6.M1 7.ID0 8.ID1 9.A 10.ECC 11.Tform(Myr) 12.Snapno', delimiter='', comments='#')

	np.savetxt('/projects/b1011/syr904/projects/PULSAR/kickgrid_runs/MSPBHinGC/data'+'/num_bnsmerger.dat', np.c_[Models, Ncoll, Ncoll2, Ncoll3, Ncoll4, Nsemerge, Nesc, Nescmerge], fmt='%d %d %d %d %d %d %d %d', header='1.Model 2.Ncoll 3.Ncoll2 4.Ncoll3 5.Ncoll4 6.Nsemerge 7.Nesc 8.Nescmerge', delimiter='', comments='#')

	#np.savetxt(savepath+'/extrememodels/Escbns.dat', np.c_[model_esc, timeesc, timeesc_myr, tins, m0, m1, id0, id1, a, ecc], fmt='%d %f %f %f %f %f %d %d %f %f', header='1.Model 2.Time_esc(code) 3.Time_esc(Myr) 4.T_insp(Myr) 5.M0 6.M1 7.ID0 8.ID1 9.A 10.ECC', delimiter='', comments='#')



def find_all_mergers_pnmodels(pathlist, start, end):   
##Both in the cluster and escaped. Code adjusted for Carl's PN models since their are extension files from reruns.
	sourcedir=np.genfromtxt(pathlist, dtype=str)
	pref=sourcedir[:,1]
	sourcedir=sourcedir[:,0]

	##Mergers in the clusters
	model_coll=[]; model_merge=[]

	idcoll=[]; id0coll=[]; id1coll=[]; id2coll=[]; id3coll=[]
	idsemerge=[]; id0merge=[]; id1merge=[]
	typecoll=[]; k0=[]; k1=[]; k2=[]; k3=[]
	timecoll=[]; timesemerge=[]
	timecoll_myr=[]; timesemerge_myr=[]
	mf_merge=[]; m0_merge=[]; m1_merge=[]
	mf_coll=[]; m0_coll=[]; m1_coll=[]; m2_coll=[]; m3_coll=[]
	r_merge=[]; r_coll=[]
	tform_semerge=[]; snapno_semerge=[]

	##Mergers outside of the clusters
	model_esc=[]; timeesc=[]; timeesc_myr=[]; tins=[]; m0=[]; m1=[]; id0=[]; id1=[]; a=[]; ecc=[]
	tform_esc=[]; snapno_esc=[]


	##Numbers
	Ncoll2=[]; Ncoll3=[]; Ncoll4=[]; Ncoll=[]
	Nsemerge=[]
	Nesc=[]; Nescmerge=[]
	Models=[]


	####For mergers that happen in the clusters####
	for i in range(start, end):
		filestr=sourcedir[i]+pref[i]
		collfile=filestr+'.collision.log'
		collfileext=filestr+'-extension.collision.log'
		binintfile=filestr+'.binint.log'
		binintfileext=filestr+'-extension.binint.log'
		semergefile=filestr+'.semergedisrupt.log'
		semergefileext=filestr+'-extension.semergedisrupt.log'

		t_conv=dyn.conv('t', filestr+'.conv.sh')

		ncoll2=0; ncoll3=0; ncoll4=0
		nsemerge=0
		nesc=0; nescmerger=0

		##Check in-cluster merger in the collision file
		colldata=scripts1.readcollfile(collfile)
		if os.path.isfile(collfileext) and os.path.getsize(collfileext) > 0:
			colldataext=scripts1.readcollfile(collfileext)
			colldata=colldata+colldataext

		for j in range(len(colldata)):
			line=colldata[j].split()
			if int(line[2])==2:  ##Single-single star collision
				if int(line[11])==13 and int(line[12])==13:
					model_coll.append(i)
					ncoll2+=1
					idcoll.append(int(line[3])); id0coll.append(int(line[5])); id1coll.append(int(line[7])); id2coll.append(-100); id3coll.append(-100)
					typecoll.append(int(line[2]))
					timecoll.append(float(line[0])); timecoll_myr.append(t_conv*float(line[0]))
					mf_coll.append(float(line[4])); m0_coll.append(float(line[6])); m1_coll.append(float(line[8])); m2_coll.append(-100); m3_coll.append(-100)
					r_coll.append(float(line[9]))
					k0.append(int(line[11])); k1.append(int(line[12])); k2.append(-100); k3.append(-100)
					

			if int(line[2])==3:  ##Binary-single star collision
				if len(line)==16:  ##Three stars collision
					if (int(line[13])==13 and int(line[14])==13) or (int(line[14])==13 and int(line[15])==13) or (int(line[13])==13 and int(line[15])==13):
						model_coll.append(i)
						ncoll3+=1
						idcoll.append(int(line[3])); id0coll.append(int(line[5])); id1coll.append(int(line[7])); id2coll.append(int(line[9])); id3coll.append(-100)
						typecoll.append(int(line[2]))
						timecoll.append(float(line[0])); timecoll_myr.append(t_conv*float(line[0]))
						mf_coll.append(float(line[4])); m0_coll.append(float(line[6])); m1_coll.append(float(line[8])); m2_coll.append(float(line[10])); m3_coll.append(-100)
						r_coll.append(float(line[11]))
						k0.append(int(line[13])); k1.append(int(line[14])); k2.append(int(line[15])); k3.append(-100)
						
						
				if len(line)==13: ##Two stars collision
					if int(line[11])==13 and int(line[12])==13:
						model_coll.append(i)
						ncoll3+=1
						idcoll.append(int(line[3])); id0coll.append(int(line[5])); id1coll.append(int(line[7])); id2coll.append(-100); id3coll.append(-100)
						typecoll.append(int(line[2]))
						timecoll.append(float(line[0])); timecoll_myr.append(t_conv*float(line[0]))
						mf_coll.append(float(line[4])); m0_coll.append(float(line[6])); m1_coll.append(float(line[8])); m2_coll.append(-100); m3_coll.append(-100)
						r_coll.append(float(line[9]))
						k0.append(int(line[11])); k1.append(int(line[12])); k2.append(-100); k3.append(-100)


			if int(line[2])==4:  ##Binary-binary star collision
				if len(line)==19: ##Four stars collision
					if (int(line[15])==13 and int(line[16])==13) or (int(line[16])==13 and int(line[17])==13) or (int(line[17])==13 and int(line[18])==13) or (int(line[15])==13 and int(line[17])==13) or (int(line[15])==13 and int(line[18])==13) or (int(line[16])==13 and int(line[18])==13):
						model_coll.append(i)
						ncoll4+=1
						idcoll.append(int(line[3])); id0coll.append(int(line[5])); id1coll.append(int(line[7])); id2coll.append(int(line[9])); id3coll.append(int(line[11]))
						typecoll.append(int(line[2]))
						timecoll.append(float(line[0])); timecoll_myr.append(t_conv*float(line[0]))
						mf_coll.append(float(line[4])); m0_coll.append(float(line[6])); m1_coll.append(float(line[8])); m2_coll.append(float(line[10])); m3_coll.append(float(line[12]))
						r_coll.append(float(line[13]))
						k0.append(int(line[15])); k1.append(int(line[16])); k2.append(int(line[17])); k3.append(int(line[18]))


				if len(line)==16:  ##Three stars collision
					if (int(line[13])==13 and int(line[14])==13) or (int(line[14])==13 and int(line[15])==13) or (int(line[13])==13 and int(line[15])==13):
						model_coll.append(i)
						ncoll4+=1
						idcoll.append(int(line[3])); id0coll.append(int(line[5])); id1coll.append(int(line[7])); id2coll.append(int(line[9])); id3coll.append(-100)
						typecoll.append(int(line[2]))
						timecoll.append(float(line[0])); timecoll_myr.append(t_conv*float(line[0]))
						mf_coll.append(float(line[4])); m0_coll.append(float(line[6])); m1_coll.append(float(line[8])); m2_coll.append(float(line[10])); m3_coll.append(-100)
						r_coll.append(float(line[11]))
						k0.append(int(line[13])); k1.append(int(line[14])); k2.append(int(line[15])); k3.append(-100)
						
						
				if len(line)==13: ##Two stars collision
					if int(line[11])==13 and int(line[12])==13:
						model_coll.append(i)
						ncoll4+=1
						idcoll.append(int(line[3])); id0coll.append(int(line[5])); id1coll.append(int(line[7])); id2coll.append(-100); id3coll.append(-100)
						typecoll.append(int(line[2]))
						timecoll.append(float(line[0])); timecoll_myr.append(t_conv*float(line[0]))
						mf_coll.append(float(line[4])); m0_coll.append(float(line[6])); m1_coll.append(float(line[8])); m2_coll.append(-100); m3_coll.append(-100)
						r_coll.append(float(line[9]))
						k0.append(int(line[11])); k1.append(int(line[12])); k2.append(-100); k3.append(-100)


		print 'collfile:', ncoll2, ncoll3, ncoll4#, typecoll, timecoll
		Ncoll2.append(ncoll2); Ncoll3.append(ncoll3); Ncoll4.append(ncoll4)
		Ncoll.append(ncoll2+ncoll3+ncoll4)

	
		##Check in-cluster merger in the semerge file
		semergedata=scripts2.readmergefile(semergefile)
		if os.path.isfile(semergefileext) and os.path.getsize(semergefileext)>0:
			semergedataext=scripts2.readmergefile(semergefileext)
			semergedata=semergedata+semergedataext

		for k in range(len(semergedata)):
			line=semergedata[k].split()
			if int(line[1])<3:
				if int(line[-1])==13 and int(line[-2])==13:
					model_merge.append(i)
					nsemerge+=1
					timesemerge.append(float(line[0])); timesemerge_myr.append(t_conv*float(line[0]))
					idsemerge.append(int(line[2])); id0merge.append(int(line[4])); id1merge.append(int(line[6]))
					mf_merge.append(float(line[3])); m0_merge.append(float(line[5])); m1_merge.append(float(line[7]))
					r_merge.append(float(line[8]))
					tf_semerge, sno_semerge=find_formationtime([int(line[4]), int(line[6])], sourcedir[i])
					if tf_semerge!=-100:
						tform_semerge.append(tf_semerge*t_conv); snapno_semerge.append(sno_semerge)
					else:
						tform_semerge.append(t_conv*float(line[0])); snapno_semerge.append(-100)


		print 'semerge:', nsemerge#, idsemerge, timesemerge
		Nsemerge.append(nsemerge)


		####For mergers that happen outside of the clusters####
		escfile=filestr+'.esc.dat'
		with open(escfile, 'r') as fesc:
			next(fesc)
			for line in fesc:
				dataesc=line.split()
				if int(dataesc[14])==1:
					if int(dataesc[22])==13 and int(dataesc[23])==13:
						nesc+=1
						#print float(dataesc[19]), float(dataesc[20]), float(dataesc[15]), float(dataesc[16])
						t_inspiral=gwcalc.t_inspiral_2(float(dataesc[19]), float(dataesc[20]), float(dataesc[15]), float(dataesc[16]), 0, 0, 0, 1100)/10**6 ##in Myr
						tesc=float(dataesc[1])*t_conv  ##In Myr
						model_esc.append(i)
						timeesc_myr.append(float(dataesc[1])*t_conv); tins.append(t_inspiral)
						timeesc.append(float(dataesc[1])); m0.append(float(dataesc[15])); m1.append(float(dataesc[16])); id0.append(int(dataesc[17])); id1.append(int(dataesc[18])); a.append(float(dataesc[19])); ecc.append(float(dataesc[20]))
						tf_esc, sno_esc=find_formationtime([int(dataesc[17]), int(dataesc[18])], filestr)
						if tf_esc!=-100:
							tform_esc.append(tf_esc*t_conv); snapno_esc.append(sno_esc)
						else:
							tform_esc.append(float(dataesc[1])*t_conv); snapno_esc.append(-100)

						if t_inspiral+tesc<=14000.:
							nescmerger+=1
							#print int(dataesc[17]), int(dataesc[18]), tesc, t_inspiral
							#model_esc.append(i)
							#timeesc_myr.append(float(dataesc[1])*t_conv); tins.append(t_inspiral)
							#timeesc.append(float(dataesc[1])); m0.append(float(dataesc[15])); m1.append(float(dataesc[16])); id0.append(int(dataesc[17])); id1.append(int(dataesc[18])); a.append(float(dataesc[19])); ecc.append(float(dataesc[20]))
							#tf_esc, sno_esc=find_formationtime([int(dataesc[17]), int(dataesc[18])], filestr)
							#if tf_esc!=-100:
							#	tform_esc.append(tf_esc*t_conv); snapno_esc.append(sno_esc)
							#else:
							#	tform_esc.append(float(dataesc[1])*t_conv); snapno_esc.append(-100)

		escfileext=filestr+'-extension.esc.dat'
		if os.path.isfile(escfileext) and os.path.getsize(escfileext) > 0:
			print 'model=', i
			with open(escfileext, 'r') as fesc:
				next(fesc)
				for line in fesc:
					dataesc=line.split()
					if int(dataesc[14])==1:
						if int(dataesc[22])==13 and int(dataesc[23])==13:
							nesc+=1
							#print float(dataesc[19]), float(dataesc[20]), float(dataesc[15]), float(dataesc[16])
							t_inspiral=gwcalc.t_inspiral_2(float(dataesc[19]), float(dataesc[20]), float(dataesc[15]), float(dataesc[16]), 0, 0, 0, 1100)/10**6 ##in Myr
							tesc=float(dataesc[1])*t_conv  ##In Myr
							model_esc.append(i)
							timeesc_myr.append(float(dataesc[1])*t_conv); tins.append(t_inspiral)
							timeesc.append(float(dataesc[1])); m0.append(float(dataesc[15])); m1.append(float(dataesc[16])); id0.append(int(dataesc[17])); id1.append(int(dataesc[18])); a.append(float(dataesc[19])); ecc.append(float(dataesc[20]))
							tf_esc, sno_esc=find_formationtime([int(dataesc[17]), int(dataesc[18])], filestr)
							if tf_esc!=-100:
								tform_esc.append(tf_esc*t_conv); snapno_esc.append(sno_esc)
							else:
								tform_esc.append(float(dataesc[1])*t_conv); snapno_esc.append(-100)

							if t_inspiral+tesc<=14000.:
								nescmerger+=1
								#print int(dataesc[17]), int(dataesc[18]), tesc, t_inspiral
								#model_esc.append(i)
								#timeesc_myr.append(float(dataesc[1])*t_conv); tins.append(t_inspiral)
								#timeesc.append(float(dataesc[1])); m0.append(float(dataesc[15])); m1.append(float(dataesc[16])); id0.append(int(dataesc[17])); id1.append(int(dataesc[18])); a.append(float(dataesc[19])); ecc.append(float(dataesc[20]))
								#tf_esc, sno_esc=find_formationtime([int(dataesc[17]), int(dataesc[18])], sourcedir[i])
								#if tf_esc!=-100:
								#	tform_esc.append(tf_esc*t_conv); snapno_esc.append(sno_esc)
								#else:
								#	tform_esc.append(float(dataesc[1])*t_conv); snapno_esc.append(-100)	


		print 'escaped:', nesc, nescmerger
		Nesc.append(nesc); Nescmerge.append(nescmerger)


		Models.append(i)


	##Output files
	np.savetxt(savepath+'/pnmodels/GWcap.dat', np.c_[model_coll, timecoll, timecoll_myr, typecoll, idcoll, id0coll, id1coll, id2coll, id3coll, mf_coll, m0_coll, m1_coll, m2_coll, m3_coll, r_coll, k0, k1, k2, k3], fmt='%d %f %f %d %d %d %d %d %d %f %f %f %f %f %f %d %d %d %d', delimiter='', header='1.Model 2.Time(code) 3.Time(Myr) 4.Type 5.IDM 6.ID0 7.ID1 8.ID2 9.ID3 10.MM 11.M0 12.M1 13.M2 14.M3 15.R 16.K0 17.K1 18.K2 19.K3', comments='#')

	np.savetxt(savepath+'/pnmodels/Incluster.dat', np.c_[model_merge, timesemerge, timesemerge_myr, idsemerge, id0merge, id1merge, mf_merge, m0_merge, m1_merge, r_merge, tform_semerge, snapno_semerge], fmt='%d %f %f %d %d %d %f %f %f %f %f %d', delimiter='', header='1.Model 2.Time(code) 3.Time(Myr) 4.IDM 5.ID0 6.ID1 7.MM 8.M0 9.M1 10.R 11.Tform(Myr) 12.Snapno', comments='#')

	#np.savetxt(savepath+'/pnmodels/Escmerger.dat', np.c_[model_esc, timeesc, timeesc_myr, tins, m0, m1, id0, id1, a, ecc, tform_esc, snapno_esc], fmt='%d %f %f %f %f %f %d %d %f %f %f %d', header='1.Model 2.Time_esc(code) 3.Time_esc(Myr) 4.T_insp(Myr) 5.M0 6.M1 7.ID0 8.ID1 9.A 10.ECC 11.Tform(Myr) 12.Snapno', delimiter='', comments='#')

	np.savetxt(savepath+'/pnmodels/num_bnsmerger.dat', np.c_[Models, Ncoll, Ncoll2, Ncoll3, Ncoll4, Nsemerge, Nesc, Nescmerge], fmt='%d %d %d %d %d %d %d %d', header='1.Model 2.Ncoll 3.Ncoll2 4.Ncoll3 5.Ncoll4 6.Nsemerge 7.Nesc 8.Nescmerge', delimiter='', comments='#')

	np.savetxt(savepath+'/pnmodels/Escbns.dat', np.c_[model_esc, timeesc, timeesc_myr, tins, m0, m1, id0, id1, a, ecc, tform_esc, snapno_esc], fmt='%d %f %f %f %f %f %d %d %f %f %f %d', header='1.Model 2.Time_esc(code) 3.Time_esc(Myr) 4.T_insp(Myr) 5.M0 6.M1 7.ID0 8.ID1 9.A 10.ECC 11.Tform(Myr) 12.Snapno', delimiter='', comments='#')



def find_formationchannel_escbns(pathlist, escbnsfile, start, end):
	databns=np.genfromtxt(escbnsfile)
	id0=databns[:,6]; id1=databns[:,7]; model=databns[:,0]
	sourcedir=np.genfromtxt(pathlist, dtype=str)
	ffm=open('/projects/b1011/syr904/projects/SGRB/extrememodels/fmchannel.dat', 'a+', 0)
	ffm.write('#1.Model 2.ID0 3.ID1 4.Primordial 5.Dynamics\n')
	##6.ExcNS 7.ExcNonNS 8.CollNS 9.EvolNS 10.DynEject 11.SNEject
	for i in range(start, end):
		filestr=sourcedir[i]+'/initial'
		modelno=i
		for j in range(len(model)):
			if model[j]==modelno:
				firstsnap=filestr+'.snap0000.dat.gz'
				dyn=1; pri=0
				with gzip.open(firstsnap, 'r') as ffirst:
					next(ffirst)
					next(ffirst)
					for line in ffirst:
						datafirst=line.split()
						if (int(datafirst[10])==id0[j] and int(datafirst[11])==id1[j]) or (int(datafirst[11])==id0[j] and int(datafirst[10])==id1[j]):
							pri=1; dyn=0
							break
				
				#if dyn==1:	
						#hdict=hic.history_maker([long(id0[j])], [1], 'initial', sourcedir[i], 1)
						##Can add extra checking of histories here


				ffm.write('%d %d %d %d %d\n'%(modelno, id0[j], id1[j], pri, dyn))

		print i

	ffm.close()


##Calculate the maximum event rate from the "best" model
def find_eventrate_maximum(pathlist, totaldraw, deltat):  ##deltat in Myr
	sourcedir=np.genfromtxt(pathlist, dtype=str)
	gwfile=np.genfromtxt('/projects/b1011/syr904/projects/SGRB/extrememodels/GWcap.dat')
	inclusterfile=np.genfromtxt('/projects/b1011/syr904/projects/SGRB/extrememodels/Incluster.dat')
	escbnsfile=np.genfromtxt('/projects/b1011/syr904/projects/SGRB/extrememodels/Escbns_new.dat')
	modelgw=gwfile[:,0]; tmggw=gwfile[:,2] ##All times in Myr
	modelin=inclusterfile[:,0]; tmgin=inclusterfile[:,2]
	modelesc=escbnsfile[:,0]; tesc=escbnsfile[:,2]; tinsp=escbnsfile[:,3]
	tmgesc=np.array(tesc)+np.array(tinsp)

	agedist=np.genfromtxt('/projects/b1011/syr904/projects/SGRB/AgeDistribution/metallicity_0.001.dat')
	agelist=agedist[:,0]

	##Maximum event rate
	model_maxbns=sourcedir[2]
	filestr_max=model_maxbns+'/initial'
	collmax=filestr_max+'.collision.log'
	semergemax=filestr_max+'.semergedisrupt.log'
	escmax=filestr_max+'.esc.dat'

	draw_coll=[]; draw_semerge=[]; draw_esc=[]
	ndraw=0
	i=0
	while i < totaldraw:
		ncoll=0; nsemerge=0; nesc=0
		ageno=random.randint(0, len(agelist)-1)
		agemin=agelist[ageno]*1000. ##In Myr
		agemax=agemin+deltat  ##in Myr
		#print i
		if agemax>12000.: continue

		i+=1
		ndraw+=1
		#print agemin, agemax
		for j in range(len(modelgw)):
			if modelgw[j]==2:
				if tmggw[j]>=agemin and tmggw[j]<=agemax:
					ncoll+=1
		for k in range(len(modelin)):
			if modelin[k]==2:
				if tmgin[k]>=agemin and tmgin[k]<=agemax:
					nsemerge+=1
		for h in range(len(modelesc)):
			if modelesc[h]==2:
				if tmgesc[h]>=agemin and tmgesc[h]<=agemax:
					nesc+=1

		draw_coll.append(ncoll); draw_semerge.append(nsemerge); draw_esc.append(nesc)

	print np.mean(draw_coll), np.mean(draw_semerge), np.mean(draw_esc)
	print ndraw

	mean_coll=np.mean(draw_coll); mean_semerge=np.mean(draw_semerge); mean_esc=np.mean(draw_esc)

	ntot=mean_coll+mean_semerge+mean_esc

	est1=ntot*0.3*10**9/(deltat*10**6)  ##Gpc^-3 yr^-1
	est2=ntot*0.7*10**9/(deltat*10**6)
	est3=ntot*2.1*10**9/(deltat*10**6)

	print est1, est2, est3



##Calculate the average event rate from a small set of models with no weighting on the models
def find_eventrate_average(pathlist, totaldraw, deltat, modeldraw):   ##deltat in Myr
	sourcedir=np.genfromtxt(pathlist, dtype=str)
	pref=sourcedir[:,1]; sourcedir=sourcedir[:,0]
	gwfile=np.genfromtxt('/projects/b1011/syr904/projects/SGRB/extrememodels/GWcap.dat')
	inclusterfile=np.genfromtxt('/projects/b1011/syr904/projects/SGRB/extrememodels/Incluster.dat')
	escbnsfile=np.genfromtxt('/projects/b1011/syr904/projects/SGRB/extrememodels/Escbns.dat')
	modelgw=gwfile[:,0]; tmggw=gwfile[:,2] ##All times in Myr
	modelin=inclusterfile[:,0]; tmgin=inclusterfile[:,2]
	modelesc=escbnsfile[:,0]; tesc=escbnsfile[:,2]; tinsp=escbnsfile[:,3]
	tmgesc=np.array(tesc)+np.array(tinsp)

	agedist=np.genfromtxt('/projects/b1011/syr904/projects/SGRB/AgeDistribution/metallicity_0.001.dat')
	agelist=agedist[:,0]

	mean_coll=[]; mean_semerge=[]; mean_esc=[]

	n=0; nn=0
	while n < modeldraw:
		modelno=random.randint(0, len(sourcedir)-1)
		model=sourcedir[modelno]
		filestr=model+pref[modelno]

		draw_coll=[]; draw_semerge=[]; draw_esc=[]
		ni=0

		n+=1; nn+=1
		i=0
		while i < totaldraw:
			ncoll=0; nsemerge=0; nesc=0
			ageno=random.randint(0, len(agelist)-1)
			agemin=agelist[ageno]*1000. ##In Myr
			agemax=agemin+deltat  ##in Myr
			#print i
			if agemax>12000.: continue

			i+=1
			ni+=1
			#print agemin, agemax
			if len(modelgw)>0:
				for j in range(len(modelgw)):
					if modelgw[j]==modelno:
						if tmggw[j]>=agemin and tmggw[j]<=agemax:
							ncoll+=1.

			if len(modelin)>0:			
				for k in range(len(modelin)):
					if modelin[k]==modelno:
						if tmgin[k]>=agemin and tmgin[k]<=agemax:
							nsemerge+=1.

			if len(modelesc)>0:	
				for h in range(len(modelesc)):
					if modelesc[h]==modelno:
						if tmgesc[h]>=agemin and tmgesc[h]<=agemax:
							nesc+=1.

			draw_coll.append(ncoll); draw_semerge.append(nsemerge); draw_esc.append(nesc)

		#print np.mean(draw_coll), np.mean(draw_semerge), np.mean(draw_esc)
		#print ni

		mean_coll.append(np.mean(draw_coll)); mean_semerge.append(np.mean(draw_semerge)); mean_esc.append(np.mean(draw_esc))

		print n

	#print nn


	totmean_coll=np.mean(mean_coll); totmean_semerge=np.mean(mean_semerge); totmean_esc=np.mean(mean_esc)
	print totmean_coll, totmean_semerge, totmean_esc
	ntot=totmean_coll+totmean_semerge+totmean_esc

	est1=ntot*0.3*10**9/(deltat*10**6)  ##Gpc^-3 yr^-1
	est2=ntot*0.7*10**9/(deltat*10**6)  ##Gpc^-3 yr^-1
	est3=ntot*2.1*10**9/(deltat*10**6)  ##Gpc^-3 yr^-1

	print est1, est2, est3
	





















			




