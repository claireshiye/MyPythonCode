import numpy as np
import os, sys
import dynamics as dyn

modelfiles = np.genfromtxt('/projects/b1091/CMC_Grid_March2019/rundir/path_copy.dat', dtype=str)
sourcedir=modelfiles[:,0]
#print sourcedir

status=[]; paths=[]
for i in range(len(sourcedir)):
	filestr=sourcedir[i]+'initial'
	dynfile=filestr+'dyn.dat'

	if os.path.isfile(dynfile) and os.path.getsize(dynfile)>0:
		t_conv=dyn.conv('t', filestr+'.conv.sh')
		with open(dynfile, 'r') as fdyn:
			next(fdyn)
			for line in fdyn:pass
        		lastdyn=line
        	datalast=lastdyn.split()
        	time=float(datalast[0])*t_conv
        	if time>=14000.: 
        		status.append(2)
        	else: status.append(1)

	else:
		status.append(0)

	#paths.append(str(sourcedir[i]))


np.savetxt('/projects/b1091/CMC_Grid_March2019/rundir/'+'path_final.dat', np.c_[sourcedir, status],fmt='%s %d', delimiter='')





