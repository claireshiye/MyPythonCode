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
		dataobs=np.genfromtxt(snapobs[k])
		x.append(float(dataobs[0, 7])); y.append(float(dataobs[0, 8])); z.append(float(dataobs[0, 10])/1000.0)

	return x, y, z


rc=[]; rh=[]; age=[]
n=0
for _ in xrange(16):
	a, b, c=find_rcrh(n)
	rc.append(a); rh.append(b); age.append(c)
	n+=1


plt.figure(1)
for k in range(11):
	plt.plot(age[k], rc[k], color='k', linewidth=0.5)

for k1 in range(11, 16):
	plt.plot(age[k1], rc[k1], color='b', linewidth=0.5)

plt.scatter(dm[:,0], dm[:,1], color='cyan', marker='*', label='MW')
plt.scatter(agel, rcl, color='red', marker='^', label='LMC/SMC/Fornax')
plt.xlim(0.02, 20.)
plt.ylim(0.02, 12.)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('$Age(Gyr)$')
plt.ylabel('$r_c(pc)$')
plt.legend(loc='best')
#plt.savefig('/projects/b1011/syr904/projects/NGC3201/rc_new.png', dpi=300)

plt.figure(2)
for j in range(11):
	plt.plot(age[j], rh[j], color='k', linewidth=0.5)

for j1 in range(11, 16):
	plt.plot(age[j1], rh[j1], color='b', linewidth=0.5)

plt.scatter(dm[:,0], dm[:,2], color='cyan', marker='*', label='MW')
plt.scatter(agel, rhl, color='red', marker='^', label='LMC/SMC/Fornax')
plt.xscale('log')
plt.yscale('log')
plt.xlim(0.02, 20.)
plt.ylim(0.2, 25.)
plt.xlabel('$Age(Gyr)$')
plt.ylabel('$r_h(pc)$')
plt.legend(loc='upper left')
#plt.savefig('/projects/b1011/syr904/projects/NGC3201/rh_new.png', dpi=300)

f, axarr=plt.subplots(2, sharex=True, sharey=False)
for k in range(11):
	axarr[0].plot(age[k], rc[k], color='k', linewidth=1, zorder=1)

for k1 in range(11,16):
	axarr[0].plot(age[k1], rc[k1], color='b', linewidth=1, zorder=1)

axarr[0].scatter(dm[:,0], dm[:,1], color='orange', marker='^', s=15, label='$MW$', zorder=2)
axarr[0].plot(dm[6,0], dm[6,1], color='yellow', marker='*', markersize=14, label='$NGC\ 3201$', zorder=3)
axarr[0].plot([1000, 10000], color='k', label='$Models\ 1-11$')
axarr[0].plot([1000, 10000], color='b', label='$Models\ 12-16$')
axarr[0].set_xscale('log')
axarr[0].set_yscale('log')
axarr[0].set_ylabel(r'$\rm{r_c(pc)}$', fontsize=18)
axarr[0].set_ylim(0.02, 12)

for j in range(11):
	axarr[1].plot(age[j], rh[j], color='k', linewidth=1, zorder=1)

for j1 in range(11,16):
	axarr[1].plot(age[j1], rh[j1], color='b', linewidth=1, zorder=1)

axarr[1].scatter(dm[:,0], dm[:,2], color='orange', marker='^', s=15, label='$MW$', zorder=2)
axarr[1].plot(dm[6,0], dm[6,2], color='yellow', marker='*', markersize=14, label='$NGC\ 3201$', zorder=3)
axarr[1].plot([1000, 10000], color='k', label='$Models\ 1-11$')
axarr[1].plot([1000, 10000], color='b', label='$Models\ 12-16$')
axarr[1].set_xscale('log')
axarr[1].set_yscale('log')
axarr[1].set_xlabel(r'$\rm{time(Gyr)}$',fontsize=18)
axarr[1].set_ylabel(r'$\rm{r_h(pc)}$', fontsize=18)
axarr[1].set_xlim(0.08,20)
axarr[1].set_ylim(0.2, 25)

f.subplots_adjust(hspace=0)
for ax in axarr.flat:
    ax.label_outer()	

plt.legend(loc='best', prop={'size': 10}, numpoints=1)
#plt.legend([yellow_star, red_triangle], ["$Models 1-11$", "Attr A+B"])

plt.savefig('/projects/b1011/syr904/projects/NGC3201/rcrh_new.png', dpi=300)
