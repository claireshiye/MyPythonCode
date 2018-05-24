import numpy as np
import matplotlib
#matplotlib.use('PDF')
import matplotlib.pyplot as plt 
from glob import glob
import pandas as pd
import math
import conversions as conv

path=np.genfromtxt('/projects/b1011/sourav/new_runs/kick_grid/kick_grid_path.dat', dtype='|S')
dm=np.genfromtxt('/projects/b1011/syr904/projects/NGC3201/MWGC.dat')
dl=np.genfromtxt('/projects/b1011/syr904/projects/NGC3201/LMC.dat')
path1=np.append(path, ['/projects/b1011/syr904/cmc/cmc-mpi-04/rundir/NGC3201/8e5w5rv1fb50'])
df1=pd.read_fwf('/projects/b1011/syr904/projects/NGC3201/harris_table1.dat', header=None)
df2=pd.read_fwf('/projects/b1011/syr904/projects/NGC3201/harris_table2.dat', header=None)
df3=pd.read_fwf('/projects/b1011/syr904/projects/NGC3201/harris_table3.dat', header=None)
names=df1[0].values; magha=df2[7].values; rcha=df3[8].values; rhha=df3[9].values; dsunha=df1[10].values; dgcha=df1[11].values
#print dsunha, dgcha
#print magha, rcha, rhha
#print np.isnan(magharris)
#print path[32]


##Filter super young LMC clusters
agel=[]; rcl=[]; rhl=[]
for i in range(len(dl[:,0])):
	if float(dl[:,0][i])>0.1*float(dl[:,3][i]):
		agel.append(float(dl[:,0][i])); rcl.append(float(dl[:,1][i])); rhl.append(float(dl[:,2][i]))


##Find rc, rh and age in models
def find_rcrh(no):
	filepath=path1[no]
	#print filepath
	snapobs=np.sort(glob(filepath+'/'+'initial.snap*.obs_params.dat'))
	x=[]; y=[]; z=[]
	for k in range(len(snapobs)):
		dataobs=np.genfromtxt(snapobs[k])
		x.append(float(dataobs[0, 7])); y.append(float(dataobs[0, 8])); z.append(float(dataobs[0, 10])/1000.0)

	return x, y, z

def select_mwgc_mass(dm):
	agem=dm[:,0]; rcm=dm[:,1]; rhm=dm[:,2]; mm=dm[:,3]; dgcm=dm[:,4]
	Agem=[]; Rcm=[]; Rhm=[]
	for l in range(len(agem)):
		if dgcm[l]>=5. and dgcm[l]<=13.:
			if mm[l]>=80000. and mm[l]<=750000.:
				Agem.append(agem[l]); Rcm.append(rcm[l]); Rhm.append(rhm[l])		
	return Agem, Rcm, Rhm

def select_mwgc_mag(magharris, rcharris, rhharris, dsunharris, dgcharris, nameharris):
	#magharris=harris2[7].values; rcharris=harris3[8].values; rhharris=harris3[9].values
	gcname=[]; maghr=[]; rchr=[]; rhhr=[]; dsunhr=[]; dgchr=[]
	Agem=[]; Rcm=[]; Rhm=[]
	for l in range(157):
		if math.isnan(float(magharris[l]))==False and math.isnan(float(rcharris[l]))==False and math.isnan(float(rhharris[l]))==False:
			gcname.append(str(nameharris[l])); maghr.append(float(magharris[l])); rchr.append(float(rcharris[l])); rhhr.append(float(rhharris[l])); dsunhr.append(float(dsunharris[l])); dgchr.append(float(dgcharris[l]))
	
	plt.figure()
	plt.hist(maghr, bins=20)
	plt.show()

	plt.figure()
	plt.hist(dgchr, bins=50)
	plt.xlim(0, 20.)
	plt.show()
	
	#print maghr, rchr, rhhr, dsunhr, len(maghr), len(rchr), len(rhhr), len(dsunhr)
	##Conversion from arcmin to pc
	rchrpc=[]; rhhrpc=[]
	for h in range(len(rchr)):
		rchrpc.append(conv.arcmin_to_pc(rchr[h], dsunhr[h]))
		rhhrpc.append(conv.arcmin_to_pc(rhhr[h], dsunhr[h]))

	#print rchrpc, rhhrpc, maghrpc, dsunhrpc
	##Selection
	xsel=0; ysel=0	
	for m in range(len(maghr)):
		#if maghr[m]<=-6.45 and maghr[m]>=-8.45:
		if maghr[m]<=-5.95 and maghr[m]>=-8.95:
			xsel+=1
			if dgchr[m]>=6. and dgchr[m]<=12.0:
				ysel+=1
				print gcname[m]
				#print maghr[m]
				Agem.append(12.5); Rcm.append(rchrpc[m]); Rhm.append(rhhrpc[m]) 
	#print Agem, Rcm, Rhm		
	print xsel, ysel	
	return Agem, Rcm, Rhm
##Find rc, rh and age is the models
rc=[]; rh=[]; age=[]
n=0
for _ in xrange(16):
	a, b, c=find_rcrh(n)
	rc.append(a); rh.append(b); age.append(c)
	n+=1

rc50, rh50, age50=find_rcrh(32)


def plot_separatefig():
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

#AGEM, RCM, RHM=select_mwgc_mass(dm)
AGEM, RCM, RHM=select_mwgc_mag(magha, rcha, rhha, dsunha, dgcha, names)
print len(AGEM)
f, axarr=plt.subplots(2, sharex=True, sharey=False)
for k in range(11):
	if k!=4:
		axarr[0].plot(age[k], rc[k], color='k', linewidth=1, zorder=1)
axarr[0].plot(age[4], rc[4], color='red', linewidth=1.5, zorder=1)

for k1 in range(11,16):
	axarr[0].plot(age[k1], rc[k1], color='b', linewidth=1, zorder=1)

#axarr[0].plot(age50, rc50, color='red', linewidth=1, zorder=1)
#axarr[0].errorbar(12.5,RCM[15], xerr=0.5, capsize=0, color='green', fmt='.', ms=1, lw=1.5, label=r'$\rm{NGC 6397}$', zorder=2)
axarr[0].errorbar(AGEM, RCM, xerr=0.5, capsize=0, color='orange', fmt='.', ms=1, lw=1.5, label=r'$\rm{Other\ MW\ GCs}$', zorder=2)
axarr[0].errorbar(12.5,RCM[15], xerr=0.5, capsize=0, color='limegreen', fmt='.', ms=1, lw=1.5, label=r'$\rm{NGC 6397}$', zorder=2)
#axarr[0].scatter(AGEM, RCM, color='orange', marker='^', s=27, label='$MW$', zorder=2)
#axarr[0].scatter(dm[:,0], dm[:,1], color='orange', marker='^', s=15, label='$MW$', zorder=2)
axarr[0].plot(11.2, dm[6,1], color='red', marker='*', markersize=22, label=r'$\rm{NGC}\ 3201$', zorder=3)
axarr[0].plot([1000, 10000], color='k', label=r'$\rm{Models\ with\ N_{BH} \gtrsim 200}$')
axarr[0].plot([1000, 10000], color='b', label=r'$\rm{Models\ with N_{BH} \lesssim 10}$')
#axarr[0].set_xscale('log')
#axarr[0].set_yscale('log')
axarr[0].set_ylabel(r'$\rm{r_c(pc)}$', fontsize=18)
axarr[0].set_ylim(-0.4, 4.0)

for j in range(11):
	if j!=4:
		axarr[1].plot(age[j], rh[j], color='k', linewidth=1, zorder=1)
axarr[1].plot(age[4], rh[4], color='red', linewidth=1.5, zorder=1)

for j1 in range(11,16):
	axarr[1].plot(age[j1], rh[j1], color='b', linewidth=1, zorder=1)

#axarr[1].plot(age50, rh50, color='red', linewidth=1, zorder=1)
#axarr[1].errorbar(12.5,RHM[15], xerr=0.5, capsize=0, color='green', fmt='.', ms=1, lw=1.5, label=r'$\rm{NGC 6397}$', zorder=2)
axarr[1].errorbar(AGEM, RHM, xerr=0.5, capsize=0, color='orange', fmt='.', ms=1, lw=1.5, label=r'$\rm{Other\ MW\ GCs}$', zorder=2)
axarr[1].errorbar(12.5,RHM[15], xerr=0.5, capsize=0, color='limegreen', fmt='.', ms=1, lw=1.5, label=r'$\rm{NGC\ 6397}$', zorder=2)
#plt.axhline(RHM[15], xmin=12., xmax=13., color='limegreen', lw=1.5, ls='--', label=r'$\rm{NGC\ 6397}$', zorder=2)
#axarr[1].scatter(AGEM, RHM, color='orange', marker='^', s=27, label='$MW$', zorder=2)
#axarr[1].scatter(dm[:,0], dm[:,2], color='orange', marker='^', s=15, label='$MW$', zorder=2)
axarr[1].plot(11.2, dm[6,2], color='red', marker='*', markersize=22, markeredgewidth=1, label=r'$\rm{NGC}\ 3201$', zorder=3)
axarr[1].plot([1000, 10000], color='k', label=r'$\rm{Models\ with}\ \it{N}_{BH} \gtrsim 200$')
axarr[1].plot([1000, 10000], color='b', label=r'$\rm{Models\ with}\ \it{N}_{BH} \lesssim 10$')
#axarr[1].set_xscale('log')
#axarr[1].set_yscale('log')
axarr[1].set_xlabel(r'$\rm{time(Gyr)}$',fontsize=18)
axarr[1].set_ylabel(r'$\rm{r_{hl}(pc)}$', fontsize=18)
#axarr[1].set_xlim(0.08,15)
axarr[1].set_xlim(-0.5,13.2)
axarr[1].set_ylim(0.0, 5.8)

f.subplots_adjust(hspace=0)
for ax in axarr.flat:
    ax.label_outer()	

axarr[0].xaxis.set_tick_params(width=1,length=6,which='major')
axarr[0].xaxis.set_tick_params(width=1,length=2,which='minor')
axarr[0].yaxis.set_tick_params(width=1,length=6, which='major')

for tick in axarr[0].xaxis.get_major_ticks():
        tick.label.set_fontsize(14)
for tick in axarr[0].yaxis.get_major_ticks():
        tick.label.set_fontsize(14)

axarr[1].xaxis.set_tick_params(width=1,length=6,which='major')
axarr[1].xaxis.set_tick_params(width=1,length=2,which='minor')
axarr[1].yaxis.set_tick_params(width=1,length=6, which='major')

for tick in axarr[1].xaxis.get_major_ticks():
        tick.label.set_fontsize(14)
for tick in axarr[1].yaxis.get_major_ticks():
        tick.label.set_fontsize(14)

plt.legend(loc='upper left', prop={'size': 12}, numpoints=1, frameon=False)
#plt.tight_layout()
f.set_size_inches(7,9)
plt.show()
#plt.savefig('/projects/b1011/syr904/projects/NGC3201/rcrh_new.pdf', dpi=350)
