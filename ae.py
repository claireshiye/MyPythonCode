import numpy as np
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math

f = open('/projects/b1011/syr904/projects/PULSAR/bse_change/history_kconst/IDae.dat','r')

lines = f.readlines()
print lines
colors = cm.rainbow(np.linspace(0, 1, len(lines)))

for i in range(len(lines)):
	line = lines[i]
	line = line.split('\n')
	line = line[0]	
	print line
	data = np.loadtxt('/projects/b1011/syr904/cmc/cmc-ns-bse/rundir_micaic/8e5rv1kick_1.0_kconst49/history/'+line)	
	t_form = data[-1,0]/1000.
	RLOflag = 0
	#print line, len(data)
	a_1 = []
	e_1 = []
	a_2 = []
        e_2 = []
	a = []
	e = []
	minmax = 0
	NUM = 0
	NUM_triples = 0
	NUM_triples_stable = 0
	MTflag = 0
	plt.figure()
	ax1 = plt.subplot(1,1,1)
	for j in range(len(data)-1):
		a = []
		e = []
		a.append(data[j,4])
		a.append(data[j+1,4])
		e.append(data[j,5])
                e.append(data[j+1,5])

		if data[j,6] == 0:
			if minmax == 0:
				a_max = data[j-1,4]
				e_max = data[j-1,5]
				minmax = 1

			if data[j,1] == 0:
                                plt.plot(a,e,linewidth=2.0,ls='-',c='black',label='BSE')
				a.append(data[j,4])
				e.append(data[j,5])
			if data[j,1] == 1:
				plt.plot(a,e,linewidth=2.0,ls='--',c='blue',label='Binary-single')
				a_1.append(data[j,4])
				e_1.append(data[j,5])
				a.append(data[j,4])
				e.append(data[j,5])
				NUM = NUM + 1
			if data[j,1] == 2:
				plt.plot(a,e,linewidth=2.0,ls='--',c='red',label='Binary-binary')
				a_2.append(data[j,4])
                                e_2.append(data[j,5])
				a.append(data[j,4])
                                e.append(data[j,5])
				NUM = NUM + 1
			if data[j,1] == 3:
				plt.plot(a,e,linewidth=2.0,ls='--',c='orange',label='Triple')
				NUM_triples = NUM_triples + 1
				if data[j,7] == 1:
					NUM_triples_stable = NUM_triples_stable + 1
		if data[j,6] == 0:
                        if RLOflag == 0:
                                t_RLO = data[j-1,0]/1000.
                                a = []
				e = []
				a.append(data[j-1,4])
				a.append(data[j,4])
				e.append(data[j-1,5])
				e.append(data[j,5])
				RLOflag = 1
				#### Plot the first RLO ###
				if data[j-1,1] == 0:
					plt.plot(a,e,linewidth=2.0,ls='-',c='black')
				if data[j-1,1] == 1:
					plt.plot(a,e,linewidth=2.0,ls='--',c='blue')
				if data[j-1,1] == 2:
					plt.plot(a,e,linewidth=2.0,ls='--',c='red')
				if data[j-1,1] == 3:
					plt.plot(a,e,linewidth=2.0,ls='--',c='orange')				


	print 'TRIPLES:', NUM_triples_stable, 'out of', NUM_triples, 'are dynamically stable.'
	size = 25
	a_min = data[-1,4]
        e_min = data[-1,5]
	plt.scatter(a_1,e_1,marker='o',s=size,color='black',label='Dynamical interaction')
	plt.scatter(a_2,e_2,marker='o',s=size,color='black')
	#plt.quiver(a_full[-1], e_full[-1], a_full[1]-a_full[-1], e_full[1]-e_full[-1], scale_units='xy', angles='xy', scale=1,color=colors[i])

        plt.scatter(a_min,e_min,marker='*',s=300,color='black',label='Formation of binary')
        plt.scatter(a_max,e_max,marker='o',s=250,color='black',label='Onset of RLO')
	plt.xscale('log')
	plt.title(r'$\rm{Binary\,ID}$: $%i$, $t_{\rm{form}} = %.2f$ $\rm{Gyr}$, $t_{\rm{RLO}} = %.2f$ $\rm{Gyr}$, $N_{\rm{enc}}=%i$' %(i+1, t_form,t_RLO,NUM),fontsize=18,y=1.02)
	plt.xlim(0.0008,400)
	plt.ylim(-0.04,1.03)
	plt.xlabel(r'$a\, (\rm{AU})$',fontsize=24)
	plt.ylabel(r'$\rm{Eccentricity}$',fontsize=24)
	[i.set_linewidth(2.0) for i in ax1.spines.itervalues()]
	ax1.xaxis.set_tick_params(width=2)
	ax1.yaxis.set_tick_params(width=2)
	for tick in ax1.xaxis.get_major_ticks():
	        tick.label.set_fontsize(14)
	for tick in ax1.yaxis.get_major_ticks():
	        tick.label.set_fontsize(14)
	#plt.legend(loc=1, numpoints=1, frameon=False)

	plt.savefig('/projects/b1011/syr904/projects/PULSAR/bse_change/history_kconst/'+line+'.pdf')
	#plt.show()


