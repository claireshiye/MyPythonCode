import numpy as np
from glob import glob
import collections
from collections import Counter
import os,sys
import subprocess
import matplotlib
import matplotlib.lines as mlines
matplotlib.use('PDF')
import matplotlib.pyplot as plt
import seaborn as sns
import gzip
import math
import re
import history_cmc as hic
import dynamics as dyn    

yearsc=31557600
twopi=6.283185307179586

biout=np.sort(glob('/projects/b1011/syr904/cmc/cmc-mpi/bse_wrap/bse/history/alloutput/binary*.dat'))

B=[]; P=[]; FC=[]
Nns=0
for i in range(len(biout)):
	databi=np.genfromtxt(biout[i])
	lastline=databi[-1]
	k0=lastline[1]; k1=lastline[15]; m0=lastline[3]; m1=lastline[17]
	if k0==13:
		B.append(lastline[32]); P.append(twopi*yearsc/lastline[12]); FC.append(lastline[34])
		Nns+=1	
	if k1==13:
		B.append(lastline[33]); P.append(twopi*yearsc/lastline[26]); FC.append(lastline[35])
		Nns+=1

for j in range(len(P)):
	if P[j]<=0.02: print j; print P[j], B[j]

print B, P, FC, Nns

#Death Line
x=np.logspace(-4.0, 2.0, num=50)

plt.figure()
plt.xscale('log')
plt.yscale('log')
plt.plot(x, (x**2)*(0.17*10**12), 'k--')    #Deadline

plt.scatter(P, B, alpha=0.6, s=50, c=FC, cmap='bwr')
#plt.scatter(Pb, Bb, color='orange', label='binary', alpha=0.7)
plt.xlim(10**-4, 100.)
plt.ylim(10**7, 10**15)
plt.xlabel(r'$P(sec)$')
plt.ylabel(r'$B(G)$')
plt.legend(loc='upper left')

#plt.show()
plt.savefig('/projects/b1011/syr904/projects/PULSAR/bse_change/bse_ns_BP.pdf', dpi=300)
