import numpy as np
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt 
from glob import glob

path=np.genfromtxt('/projects/b1011/sourav/new_runs/kick_grid/kick_grid_path.dat', dtype='|S')
dm=np.genfromtxt('/projects/b1011/syr904/projects/NGC3201/MWGC.dat')
dl=np.genfromtxt('/projects/b1011/syr904/projects/NGC3201/LMC.dat')

##Filter super young LMC clusters
agel=[]; rcl=[]; rhl=[]
for i in range(len(dl[:,0])):
	if float(dl[:,0][i])>0.1*float(dl[:,3][i]):
		agel.append(float(dl[:,0][i])); rcl.append(float(dl[:,1][i])); rhl.append(float(dl[:,2][i]))


##Find rc, rh and age in models
def find_rcrh(no):
	filepath=path[no]
	snapobs=np.sort(glob(filepath+'/'+'initial.snap*.obs_params.dat'))
	x=[]; y=[]; z=[]
	for k in range(len(snapobs)):
		dataobs=np.genfromtxt(snapobs)
		x.append(float(dataobs[0, 7])); y.append(float(dataobs[0, 8])); z.append(float(dataobs[0, 10]))

	return x, y, z


rc=[]; rh=[]; age=[]
n=0
for _ in xrange(16):
	a, b, c=find_rcrh(n)
	rc.append(a); rh.append(b); age.append(c)
	n+=1


plt.figure(1)
for k in range(16):
	plt.plot(age[k], rc[k], color='grey', linewidth=0.5)

plt.scatter(dm[:,0], dm[:,1], color='b', marker='*', label='MW')
plt.scatter(agel, rcl, color='red', marker='^', label='LMC/SMC/Fornax')
plt.xscale('log')
plt.yscale('log')
plt.legend(loc='best')
plt.savefig('/projects/b1011/syr904/projects/NGC3201/rc_new.png', dpi=300)

plt.figure(2)
for j in range(16):
	plt.plot(age[j], rh[j], color='grey', linewidth=0.5)

plt.scatter(dm[:,0], dm[:,2], color='b', marker='*', label='MW')
plt.scatter(agel, rhl, color='red', marker='^', label='LMC/SMC/Fornax')
plt.xscale('log')
plt.yscale('log')
plt.legend(loc='best')
plt.savefig('/projects/b1011/syr904/projects/NGC3201/rh_new.png', dpi=300)
