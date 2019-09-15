import numpy as np
import os, sys
import dynamics as dyn

modelfiles = np.genfromtxt('/projects/b1091/CMC_Grid_March2019/rundir/path.dat', dtype=str, delimiter='')
sourcedir=modelfiles[:,0]

status=[]
for i in range(len(sourcedir)):
	filestr1=sourcedir[i]+'initial'
	dynfile1=filestr1+'.dyn.dat'

	filestr2=sourcedir[i]+'initial2'
	dynfile2=filestr2+'.dyn.dat'


	if os.path.isfile(dynfile2) and os.path.getsize(dynfile2)>0:	
		t_conv=dyn.conv('t', filestr1+'.conv.sh')

		with open(dynfile2, 'r') as fdyn:
			next(fdyn); next(fdyn)
			for line in fdyn:
				datadyn=line.split()
				totaln=float(datadyn[3])
				break

		with open(dynfile2, 'r') as fdyn:
			for line in fdyn: pass
			lastdyn=line.split()
	
		time=float(lastdyn[0])*t_conv
		finaln=float(lastdyn[3])

		if time>=11999.9 and time<14005.: 
			status.append(1)
		if time<11999.9 and finaln<=0.1*totaln:
			status.append(2)
		if time<11999.9 and finaln>0.1*totaln: 
			status.append(3)
		if time>=14005.: 
			status.append(4)


	elif os.path.isfile(dynfile1) and os.path.getsize(dynfile1)>0:
		t_conv=dyn.conv('t', filestr1+'.conv.sh')

		with open(dynfile1, 'r') as fdyn:
			next(fdyn); next(fdyn)
			for line in fdyn:
				datadyn=line.split()
				totaln=float(datadyn[3])
				break

		with open(dynfile1, 'r') as fdyn:
			for line in fdyn: pass
			lastdyn=line.split()
	
		time=float(lastdyn[0])*t_conv
		finaln=float(lastdyn[3])

		if time>=11999.9 and time<14005.: 
			status.append(1)
		if time<11999.9 and finaln<=0.1*totaln:
			status.append(2)
		if time<11999.9 and finaln>0.1*totaln: 
			status.append(3)
		if time>=14005.: 
			status.append(4)
	else:
		status.append(0)

#arrays=np.c_[sourcedir, status].astype(str)
#import pdb; pdb.set_trace()

f=open('/projects/b1091/CMC_Grid_March2019/rundir/path_final.dat', 'a+', 0)
f.write('#1 means done; 2 means dissolved; 3 means dissolved but with final number of stars > 0.1 initial number of stars; 4 means model runs past 14 Gyr; 0 means missing dyn file\n')
for j in range(len(sourcedir)):
	f.write('%s %d\n'%(sourcedir[j], status[j]))

f.close()








