import numpy as np
import os, sys
import dynamics as dyn

modelfiles = np.genfromtxt('/projects/b1091/CMC_Grid_March2019/rundir/path_copy.dat', dtype=str, delimiter='')
sourcedir=modelfiles[:,0]

status=[]
for i in range(len(sourcedir)):
	filestr=sourcedir[i]+'initial'
	dynfile=filestr+'.dyn.dat'

	if os.path.isfile(dynfile) and os.path.getsize(dynfile)>0:
		t_conv=dyn.conv('t', filestr+'.conv.sh')

		with open(dynfile, 'r') as fdyn:
			next(fdyn); next(fdyn)
			for line in fdyn:
				datadyn=line.split()
				totaln=float(datadyn[3])
				break

		with open(dynfile, 'r') as fdyn:
			for line in fdyn: pass
			lastdyn=line.split()
	
		time=float(lastdyn[0])*t_conv
		finaln=float(lastdyn[3])

		if time>=13999.9: 
			status.append(2)
		if time<13999.9 and finaln<=0.1*totaln:
			status.append(3)
		if time<13999.9 and finaln>0.1*totaln: 
			status.append(1)
	else:
		status.append(0)

#arrays=np.c_[sourcedir, status].astype(str)
#import pdb; pdb.set_trace()

f=open('/projects/b1091/CMC_Grid_March2019/rundir/path_final.dat', 'a+', 0)
for j in range(len(sourcedir)):
	f.write('%s %d\n'%(sourcedir[j], status[j]))

f.close()








