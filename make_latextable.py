import numpy as np
import re
import collections
from collections import Counter


##Read in the relevant data files
sourcedir=np.genfromtxt('/projects/b1095/syr904/projects/SGRB/newruns/finaldata/path_allfinished_newruns_maingrid.dat', dtype=str)
paths=sourcedir[:,0]; status=sourcedir[:,1]

cluster_properties=np.genfromtxt('/projects/b1095/syr904/projects/SGRB/newruns/finaldata/clusterproperty_12Gyr_maingrid.dat')
DNS_NSBH=np.genfromtxt('/projects/b1095/syr904/projects/SGRB/newruns/finaldata/DNS_NSBH_Unique_9to12Gyr_maingrid.dat')
numbbh=np.genfromtxt('/projects/b1095/syr904/projects/SGRB/newruns/finaldata/num_merger_BBH_maingrid.dat')
numdns=np.genfromtxt('/projects/b1095/syr904/projects/SGRB/newruns/finaldata/num_merger_DNS_maingrid.dat')
numnsbh=np.genfromtxt('/projects/b1095/syr904/projects/SGRB/newruns/finaldata/num_merger_NSBH_maingrid.dat')


N=[]; Rv=[]; Z=[]; Rg=[]
Npsr=[]; Nmsp=[]; NSBH=[]; DNS=[]
model_name=[]
#rc=[]; rhl=[]; Mtot=[]; Nbh=[]; Nns=[]
#Ngwcap=[]; Nincluster=[]; Nesc=[]; Nescmerge=[]


##Properties at 12 Gyr
Nbh=cluster_properties[:,1]; Mtot=np.array(cluster_properties[:,2])/100000.; rc=cluster_properties[:,5]; rhl=cluster_properties[:,6]; Nns=cluster_properties[:,7]#; status=cluster_properties[:,10]
Npsr=cluster_properties[:,10]; Nmsp=cluster_properties[:,11]


Ngwcap_dns=numdns[:,1]; Nincluster_dns=numdns[:,5]; Nesc_dns=numdns[:,6]; Nescmerge_dns=numdns[:,7]
Ngwcap_nsbh=numnsbh[:,1]; Nincluster_nsbh=numnsbh[:,5]; Nesc_nsbh=numnsbh[:,6]; Nescmerge_nsbh=numnsbh[:,7]

modelno=DNS_NSBH[:,0]; types=DNS_NSBH[:,12]

for i in range(len(paths)):
	##Initial Conditions
	s=paths[i].split('/')
	n_star=float(s[-2])
	z=float(s[-3][1:])
	rg=int(s[-4][2:])
	rv=float(s[-5][2:])

	N.append(n_star/100000.); Rv.append(rv); Z.append(z); Rg.append(rg) 
	model_name.append('N'+str(int(n_star/100000.))+'-RV'+str(rv)+'-RG'+str(int(rg))+'-Z'+str(z/0.02))
	print model_name[i]


	#if status[i]==1:
	#	model.append(str(i+1)+'*')
	#else: model.append(str(i+1))
	#model.append(i)


	###NS numbers
	#filestr=paths[i]+'initial'
	#nsfile=filestr+'.ns.dat'
	#with open(nsfile,'r') as fns:
	#	for line in fns:pass
	#	lastns=line
	#datans=lastns.split()
	#Npsr.append(int(datans[5])); Nmsp.append(int(datans[6]))


	##DNS and NSBH
	if int(status[i])==1:
		ndns=0; nnsbh=0
		for j in range(len(modelno)):
			if int(modelno[j])==i:
				if types[j]==1313.: ndns+=1
				else: nnsbh+=1
	else:
		ndns=-100; nnsbh=-100

	NSBH.append(nnsbh); DNS.append(ndns)

#print newmodel


#np.savetxt('/projects/b1095/syr904/projects/SGRB/newruns/table_allnums_BBH.txt', np.c_[model, N, Rv, Z, Rg, rc, rhl, Mtot, Nbh, Nns, Npsr, Nmsp, Ngwcap, Nincluster, Nesc, Nescmerge, DNS, NSBH, status], header='1.Model 2.N 3.rv(pc) 4.z 5.rg(kpc) 6.rc(pc) 7.rhl(pc) 8.mtot(Msun) 9.Nbh 10.Nns 11.Npsr 12.Nmsp 13.Ngwcap 14.Nincluster 15.Nesc 16.Nescmerge 17.DNS 18.NSBH 19.status', fmt='%d %d %g %g %d %.2f %.2f %.2f %d %d %d %d %d %d %d %d %d %d %d', delimiter=' ', comments='#')

f=open('/projects/b1095/syr904/projects/SGRB/newruns/finaldata/latextable_maingrid.txt', 'w+')
for i in range(len(paths)):
	if int(status[i])==1:
		f.write('%s & %.2f & %.2f & %.2f & %d & %d & %d & %d & %d & %d & %d & %d & %d & %d & %d & %d & %d & %d\\\\\n'%(model_name[i], rc[i], rhl[i], Mtot[i], Nbh[i], Nns[i], Npsr[i], Nmsp[i], Ngwcap_dns[i], Nincluster_dns[i], Nesc_dns[i], Nescmerge_dns[i], Ngwcap_nsbh[i], Nincluster_nsbh[i], Nesc_nsbh[i], Nescmerge_nsbh[i], DNS[i], NSBH[i]))
	else:
		f.write('%s & %.2f & %.2f & %.2f & %d & %d & %d & %d & %d & %d & %d & %d & %d & %d & %d & %d & %s\\\\\n'%(model_name[i], rc[i], rhl[i], Mtot[i], Nbh[i], Nns[i], Npsr[i], Nmsp[i], Ngwcap_dns[i], Nincluster_dns[i], Nesc_dns[i], Nescmerge_dns[i], Ngwcap_nsbh[i], Nincluster_nsbh[i], Nesc_nsbh[i], Nescmerge_nsbh[i], "\multicolumn{2}{c}{dissolved}"))


f.close()

		









