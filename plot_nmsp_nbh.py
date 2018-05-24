import numpy as np
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt

data=np.genfromtxt('/projects/b1011/syr904/projects/PULSAR/kickgrid_runs/ns_number_10to12Gyr.dat')
nbh=np.array(data[:,0]); ntot=np.array(data[:,1]); nmsp=np.array(data[:,2])
#print np.log(nbh/ntot), np.log(nmsp/ntot)

plt.figure()
plt.scatter(np.log(nbh/ntot), np.log(nmsp/ntot), color='b')
#plt.xscale('log')
#plt.yscale('log')
plt.xlabel(r'$ln\ N_{BH}/N_{TOT}$')
plt.ylabel(r'$ln\ N_{MSP}/N_{TOT}$')
plt.savefig('/projects/b1011/syr904/projects/PULSAR/kickgrid_runs/nbhnmsp_10to12Gyr.pdf', dpi=300)


