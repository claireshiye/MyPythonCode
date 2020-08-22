import numpy as np
import matplotlib.pyplot as plt
import scripts
import glob
import subprocess
import constants
import scipy.integrate
import json
import gzip

def plot_t_vs_rc(filestrings, plotfilename='rc_comparison.pdf', YSCALE='linear', XSCALE='linear', XLIM=[0.,13.], YLIM=[0.,4.], IND1=1, IND2=1):
    import scripts
    #plotfilename = 'rc_comparison.pdf'
    fig = plt.figure()
    fig.subplots_adjust(hspace=0.,wspace=0.)
    #XLIM = [0., 12.]
    #YLIM = [0., 4.]
    colors = ['black', 'red', 'blue', 'gold', 'green', 'magenta']
    #IND1, IND2 = MULTIPLE, len(filestrings)/MULTIPLE
    for i in range(len(filestrings)):
        print 'reading from', filestrings[i]
        units = scripts.read_units(filestrings[i])
        filename = filestrings[i]+'.dyn.pruned.dat'
        data = np.loadtxt(filename, usecols=(0,7,20))
        rc = data[:,1]*units[0]['l_pc']
	rh = data[:,2]*units[0]['l_pc']
        t = data[:,0]*units[0]['t_myr']*1e-3
    
        ax = fig.add_subplot(IND1, IND2, i+1)
        ax.plot(t, rc, ls='solid', color='black')
	ax.plot(t, rh, ls='dashed', color='blue')
	
	ax.set_xlim(XLIM)
	ax.set_ylim(YLIM)
	ax.set_xscale(XSCALE)
    	ax.set_yscale(YSCALE)
	if i+1>=IND1:
		ax.set_xlabel(r'$t\ ({\rm Gyr})$', size=18)
	else:
		ax.set_xticklabels([])
	if i%IND2==0:
                ax.set_ylabel(r'$r_c,\ r_h\,\rm{(pc)}$', size=18)
    		#ax.set_xlabel(r'$t\ ({\rm Gyr})$', size=16)
	else:
		ax.set_yticklabels([])
    	
    #ax.set_xlim([0., 12000.])
    #ax.set_ylim([0, 4.])

    plt.savefig(plotfilename)
    plt.show()


def plot_t_vs_Nbh(filestrings, plotfilename='t_vs_Nbh_comparison.pdf', YSCALE='log', XSCALE='log', XLIM=[1,30000.], YLIM=[1., 5e3]):
    import scripts
    #plotfilename = 't_vs_Nbh_comparison.pdf'
    fig = plt.figure()
    fig.subplots_adjust(hspace=0.,wspace=0.)
    colors = ['black', 'red', 'blue', 'green', 'magenta']
    for i in range(len(filestrings)):
        units = scripts.read_units(filestrings[i])
        filename = filestrings[i]+'.bh.dat'
        print 'reading from', filename
        #data = np.loadtxt(filename, usecols=(1,2,3,4,5,6,))
	data = np.loadtxt(filename)
        Nbh = data[:,2]
	Nbh_sin = data[:,3]
	Nbh_bin = data[:,4]
        Nbh_bin_bhbh = data[:,5]
        Nbh_bin_bhnbh = data[:,6]
	if YSCALE=='log':
		Nbh += 1e-15
		Nbh_sin += 1e-15
		Nbh_bin += 1e-15
		Nbh_bin_bhbh += 1e-15
		Nbh_bin_bhnbh += 1e-15
        t = data[:,1]*units[0]['t_myr']/1e3
   	print np.min(Nbh_bin_bhnbh), np.max(Nbh_bin_bhnbh) 
        #ax = fig.add_subplot(1, len(filestrings), i+1)
	ax = fig.add_subplot(1,1,1)
        #ax.plot(t, Nbh, ls='solid', color='black')
	#ax.plot(t, Nbh_sin, ls='dotted', color='blue', lw=2)
	#ax.plot(t, Nbh_bin, ls='dashed', color='red', lw=2)
        #ax.plot(t, Nbh_bin_bhbh, ls='dashed', dashes=(10, 10, 2, 10, 10, 10, 2, 10), color='green')
	#ax.plot(t, Nbh_bin_bhnbh, ls='dashed', dashes=(10,5,10,5), color='magenta')
	ax.plot(t, Nbh_bin_bhbh, ls='solid', color=colors[i])
	ax.plot(t, Nbh_bin_bhnbh, ls='dashed', color=colors[i])
  
        ax.set_xlim(XLIM)
	ax.set_ylim(YLIM)
	ax.set_xscale(XSCALE)
    	ax.set_yscale(YSCALE)
	ax.set_xlabel(r'$t\ ({\rm Gyr})$', size=16)

	if i==0:
                #ax.set_ylabel(r'$N_{bh}$, $N_{bh,b}$', size=16)
		ax.set_ylabel(r'$N_{\rm{BH}}$', size=18)
    		#ax.set_xlabel(r'$t\ ({\rm Gyr})$', size=16)
	else:
		ax.set_yticklabels([])

    	ax.set_xlabel(r'$t\ ({\rm Gyr})$', size=18)
    ax.set_ylabel(r'$N_{\rm{BH}}$', size=16)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(XLIM)
    ax.set_ylim(YLIM)
    plt.tight_layout()
     

    plt.savefig(plotfilename)
    plt.show()


def plot_t_vs_Nbh_esc(filestrings, plotfilename='t_vs_Nbh_esc_comparison.pdf', YSCALE='log', XSCALE='log', XLIM=[1,30000.], YLIM=[1., 5e3]):
    import scripts
    #plotfilename = 't_vs_Nbh_comparison.pdf'
    fig = plt.figure()
    fig.subplots_adjust(hspace=0.,wspace=0.)
    colors = ['black', 'red', 'blue', 'gold', 'green', 'magenta']
    for i in range(len(filestrings)):
        units = scripts.read_units(filestrings[i])
        filename = filestrings[i]+'.esc.bh.dat'
        print 'reading from', filename
        #data = np.loadtxt(filename, usecols=(1,2,3,4,5,6,))
	data = np.loadtxt(filename)
        Nbh = data[:,2]
	Nbh_sin = data[:,3]
	Nbh_bin = data[:,4]
        Nbh_bin_bhbh = data[:,5]
        Nbh_bin_bhnbh = data[:,6]
	if YSCALE=='log':
		Nbh += 1e-15
		Nbh_sin += 1e-15
		Nbh_bin += 1e-15
		Nbh_bin_bhbh += 1e-15
		Nbh_bin_bhnbh += 1e-15
        t = data[:,1]*units[0]['t_myr']
   	print np.min(Nbh_bin_bhnbh), np.max(Nbh_bin_bhnbh) 
        ax = fig.add_subplot(1, len(filestrings), i+1)
	ax.set_xscale(XSCALE)
    	ax.set_yscale(YSCALE)
    	ax.set_xlim(XLIM)
    	ax.set_ylim(YLIM)
	ax.set_xlabel(r'$t\ ({\rm Myr})$', size=18)
	if i==0:
		ax.set_ylabel(r'$N_{\rm{BH}}$', size=18)
	else:
		ax.set_yticklabels([])
	
        ax.plot(t, Nbh, ls='solid', color='black')
	ax.plot(t, Nbh_sin, ls='dotted', color='blue', lw=2)
	ax.plot(t, Nbh_bin, ls='dashed', color='red', lw=2)
        ax.plot(t, Nbh_bin_bhbh, ls='dashed', dashes=(10, 10, 2, 10, 10, 10, 2, 10), color='green')
	ax.plot(t, Nbh_bin_bhnbh, ls='dashed', dashes=(10,5,10,5), color='magenta')
	
    plt.savefig(plotfilename)
    plt.show()


def plot_t_vs_Nbh_esc_bound(filestrings, plotfilename='t_vs_Nbh_esc_bound.pdf', YSCALE='log', XSCALE='log', XLIM=[1,30000.], YLIM=[1., 5e3]):
    import scripts
    #plotfilename = 't_vs_Nbh_comparison.pdf'
    fig = plt.figure()
    fig.subplots_adjust(hspace=0.,wspace=0.)
    colors = ['black', 'red', 'blue', 'gold', 'green', 'magenta']
    for i in range(len(filestrings)):
        units = scripts.read_units(filestrings[i])
        filename = filestrings[i]+'.esc.bh.dat'
	filename1 = filestrings[i]+'.bh.dat'
	#BOUND BHs
        print 'reading from', filename
        #data = np.loadtxt(filename, usecols=(1,2,3,4,5,6,))
	data = np.loadtxt(filename1)
        Nbh = data[:,2]
	Nbh_sin = data[:,3]
	Nbh_bin = data[:,4]
        Nbh_bin_bhbh = data[:,5]
        Nbh_bin_bhnbh = data[:,6]
	if YSCALE=='log':
		Nbh += 1e-15
		Nbh_sin += 1e-15
		Nbh_bin += 1e-15
		Nbh_bin_bhbh += 1e-15
		Nbh_bin_bhnbh += 1e-15
        t = data[:,1]*units[0]['t_myr']
	#ESCAPED BHs
        print 'reading from', filename1
        #data = np.loadtxt(filename, usecols=(1,2,3,4,5,6,))
	data1 = np.loadtxt(filename)
        Nbh1 = data1[:,2]
	Nbh_sin1 = data1[:,3]
	Nbh_bin1 = data1[:,4]
        Nbh_bin_bhbh1 = data1[:,5]
        Nbh_bin_bhnbh1 = data1[:,6]
	if YSCALE=='log':
		Nbh1 += 1e-15
		Nbh_sin1 += 1e-15
		Nbh_bin1 += 1e-15
		Nbh_bin_bhbh1 += 1e-15
		Nbh_bin_bhnbh1 += 1e-15
        t1 = data1[:,1]*units[0]['t_myr']
   	print np.min(Nbh_bin_bhnbh), np.max(Nbh_bin_bhnbh) 
        ax = fig.add_subplot(2, len(filestrings), i+1)
	ax.set_xscale(XSCALE)
    	ax.set_yscale(YSCALE)
    	ax.set_xlim(XLIM)
    	ax.set_ylim(YLIM)
	
	#ax.set_xlabel(r'$t\ ({\rm Myr})$', size=18)
	if i==0:
		ax.set_ylabel(r'$N_{\rm{BH}}$', size=14)
	else:
		ax.set_yticklabels([])
	
	ax.set_xticklabels([])	
        ax.plot(t, Nbh, ls='solid', color='black')
	ax.plot(t, Nbh_sin, ls='dotted', color='blue', lw=2)
	ax.plot(t, Nbh_bin, ls='dashed', color='red', lw=2)
        ax.plot(t, Nbh_bin_bhbh, ls='dashed', dashes=(10, 10, 2, 10, 10, 10, 2, 10), color='green')
	ax.plot(t, Nbh_bin_bhnbh, ls='dashed', dashes=(10,5,10,5), color='magenta')
	
	ax1 = fig.add_subplot(2, len(filestrings), len(filestrings)+i+1)
	ax1.set_xscale(XSCALE)
    	ax1.set_yscale(YSCALE)
    	ax1.set_xlim(XLIM)
    	ax1.set_ylim(YLIM)
	
	#ax.set_xlabel(r'$t\ ({\rm Myr})$', size=18)
	if i==0:
		ax1.set_ylabel(r'$N_{\rm{BH, esc}}$', size=14)
	else:
		ax1.set_yticklabels([])
	
        ax1.plot(t1, Nbh1, ls='solid', color='black')
	ax1.plot(t1, Nbh_sin1, ls='dotted', color='blue', lw=2)
	ax1.plot(t1, Nbh_bin1, ls='dashed', color='red', lw=2)
        ax1.plot(t1, Nbh_bin_bhbh1, ls='dashed', dashes=(10, 10, 2, 10, 10, 10, 2, 10), color='green')
	ax1.plot(t1, Nbh_bin_bhnbh1, ls='dashed', dashes=(10,5,10,5), color='magenta')

	ax1.set_xlabel(r'$t\ ({\rm Myr})$', size=14)
	


    plt.savefig(plotfilename)
    plt.show()


def plot_t_vs_M(filestrings, plotfilename='M_comparison.pdf', YSCALE='log', XSCALE='log', XLIM=[1e-3,2e4], YLIM=[1e-2,1.2], IND1=1, IND2=1):
    import scripts
    #plotfilename = 'rc_comparison.pdf'
    fig = plt.figure()
    fig.subplots_adjust(hspace=0.,wspace=0.)
    #XLIM = [0., 12.]
    #YLIM = [0., 4.]
    colors = ['black', 'red', 'blue', 'gold', 'green', 'magenta']
    IND1, IND2 = MULTIPLE, len(filestrings)/MULTIPLE
    for i in range(len(filestrings)):
        print 'reading from', filestrings[i]
        units = scripts.read_units(filestrings[i])
        filename = filestrings[i]+'.dyn.dat'
        data = np.loadtxt(filename, usecols=(0,4))
        M_msun = data[:,1]*units[0]['m_msun']
	M_code = data[:,1]
        t = data[:,0]*units[0]['t_myr']
	print t[-1], M_code[-1]
    
        ax = fig.add_subplot(IND1, IND2, i+1)
        ax.plot(t, M_code, ls='solid', color='black')
	#ax.axvline(tdiss[i], color='red', ls='dotted', lw=1)
	
	ax.set_xlim(XLIM)
	ax.set_ylim(YLIM)
	ax.set_xscale(XSCALE)
    	ax.set_yscale(YSCALE)
	if i>=IND2:
		ax.set_xlabel(r'$t\ ({\rm Myr})$', size=16)
	if i==0:
                ax.set_ylabel(r'$M/M_i$', size=16)
    		#ax.set_xlabel(r'$t\ ({\rm Gyr})$', size=16)
	else:
		ax.set_yticklabels([])
    	
    #ax.set_xlim([0., 12000.])
    #ax.set_ylim([0, 4.])

    plt.savefig(plotfilename)
    plt.show()


def plot_t_vs_fbbh(filestrings):
    import scripts
    plotfilename = 't_vs_fbbh_comparison.pdf'
    fig = plt.figure()
    colors = ['black', 'red', 'blue', 'gold', 'green', 'magenta']
    for i in range(len(filestrings)):
        units = scripts.read_units(filestrings[i])
        filename = filestrings[i]+'.bh.dat'
        print 'reading from', filename
        data = np.loadtxt(filename, usecols=(1,12))
        Nbh = data[:,1]
        t = data[:,0]*units[0]['t_myr']
    
        ax = fig.add_subplot(111)
        ax.plot(t, Nbh, ls='solid', color=colors[i])

    ax.set_xlabel(r'$t$ (Myr)')
    ax.set_ylabel(r'$f_{b,\rm{BH}}$')
    ax.set_xscale('log')
    ax.set_yscale('linear')

    plt.savefig(plotfilename)
    plt.show()


def get_rc_sim(file_location, filestring, TWINDOW=1.):
    file=file_location+'/'+filestring+'.dyn.dat'
    data=np.loadtxt(file, usecols=(0,7,))
    convfile=file_location+'/'+filestring
    units=scripts.read_units(convfile)
    rcfarr = []
    t_final = data[-1,0]*units[0]['t_myr']
    for i in range(len(data)):
    	if data[i,0]*units[0]['t_myr']>=t_final-TWINDOW:
       		rcfarr.append(data[i,1])
    rcfarr = np.array(rcfarr)

    #rc_i, rc_f = data[0,1]*units[0]['l_pc'], mean(data[-6:-1,1])*units[0]['l_pc']
    rc_i, rc_f, rc_f_std = data[0,1]*units[0]['l_pc'], np.mean(rcfarr)*units[0]['l_pc'], np.std(rcfarr)*units[0]['l_pc']
    return rc_i, rc_f, rc_f_std
	

def get_rh_sim(file_location, filestring, TWINDOW=1.):
    file=file_location+'/'+filestring+'.dyn.dat'
    data=np.loadtxt(file, usecols=(0,20,))
    convfile=file_location+'/'+filestring
    units=scripts.read_units(convfile)
    rhfarr = []
    t_final = data[-1,0]*units[0]['t_myr']
    for i in range(len(data)):
    	if data[i,0]*units[0]['t_myr']>=t_final-TWINDOW:
       		rhfarr.append(data[i,1])
    rhfarr = np.array(rhfarr)

    #rc_i, rc_f = data[0,1]*units[0]['l_pc'], mean(data[-6:-1,1])*units[0]['l_pc']
    rh_i, rh_f, rh_f_std = data[0,1]*units[0]['l_pc'], np.mean(rhfarr)*units[0]['l_pc'], np.std(rhfarr)*units[0]['l_pc']
    return rh_i, rh_f, rh_f_std


def get_BH_sim(file_location, filestring, TWINDOW=1.):
    filename=file_location+'/'+filestring+'.bh.dat'
    data=np.loadtxt(filename)
    convfile=file_location+'/'+filestring
    units=scripts.read_units(convfile)
    nbharr = []
    fbharr = []
    nbhbharr = []
    nbhnbharr = []
    t_final = data[-1,1]*units[0]['t_myr']
    for i in range(len(data)):
    	if data[i,1]*units[0]['t_myr']>=t_final-TWINDOW:
       		nbharr.append(data[i,2])
		fbharr.append(data[i,12])
     		nbhbharr.append(data[i,5])
  		nbhnbharr.append(data[i,6])
    nbharr = np.array(nbharr)
    fbharr = np.array(fbharr)
    nbhbharr = np.array(nbhbharr)
    nbhnbharr = np.array(nbhnbharr)

    nbh, nbherr, fbh, fbherr, nbhbh, nbhbherr, nbhnbh, nbhnbherr, nbhbhrange, nbhnbhrange = np.mean(nbharr), np.std(nbharr), np.mean(fbharr), np.std(fbharr), np.mean(nbhbharr), np.std(nbhbharr), np.mean(nbhnbharr), np.std(nbhnbharr), [np.min(nbhbharr), np.max(nbhbharr)], [np.min(nbhnbharr), np.max(nbhnbharr)]
    return  nbh, nbherr, fbh, fbherr, nbhbh, nbhbherr, nbhnbh, nbhnbherr, nbhbhrange, nbhnbhrange

def get_M_sim(file_location, filestring, TWINDOW=1.):
    filename = file_location+'/'+filestring+'.dyn.dat'
    data = np.loadtxt(filename, usecols=(0,4))
    convfile=file_location+'/'+filestring
    units=scripts.read_units(convfile)
    #Marr = []
    #t_final = data[-1,0]*units[0]['t_myr']
    #for i in range(len(data)):
    #	if data[i,1]*units[0]['t_myr']>=t_final-TWINDOW:
    #		Marr.append(data[i,1])
    #Marr = np.array(Marr)
    #M_i, M_f, M_f_err = data[0,1]*units[0]['m_msun'], np.mean(Marr)*units[0]['m_msun'], np.std(Marr)*units[0]['m_msun']
    M_i, M_f = data[0,1]*units[0]['m_msun'], data[-1,1]*units[0]['m_msun']
    #return M_i, M_f, M_f_err
    return M_i, M_f


def get_N_sim(file_location, filestring, TWINDOW=1.):
    filename = file_location+'/'+filestring+'.dyn.dat'
    data = np.loadtxt(filename, usecols=(0,3))
    convfile=file_location+'/'+filestring
    units=scripts.read_units(convfile)
    Narr = []
    t_final = data[-1,0]*units[0]['t_myr']
    for i in range(len(data)):
    	if data[i,1]*units[0]['t_myr']>=t_final-TWINDOW:
		Narr.append(data[i,1])
    Narr = np.array(Narr)
    N_i, N_f, N_f_err = data[0,1], np.mean(Narr), np.std(Narr)
    return N_i, N_f, N_f_err


def extract_birthkick(filestring, readfilename, writefilename):
    import scripts
    units = scripts.read_units(filestring)
    f = open(readfilename, 'r')
    f1 = open(writefilename, 'w')
    f1.write('#1.startype 2.startype1 3.v1(km/s) 4.v2(km/s) 5.v3(km/s) 6.v(km/s) 7.kickno 8.binary?\n')
    kms = units[0]['l_cgs']/units[0]['nbt_cgs']/1e5
    print kms
    for line in f:
        if line.rfind('birth kick(iso)')>-1: #For single objects
            binary = 0
            stype1 = -1
            v1 = float(line.split('vs[1]=')[1].split(',')[0])
            v1 = v1*kms
            v2 = float(line.split('vs[2]=')[1].split(',')[0])
            v2 = v2*kms
            v3 = float(line.split('vs[3]=')[1].split(',')[0])
            v3 = v3*kms
            v = (v1*v1 + v2*v2 + v3*v3)**0.5
            stype = int(line.split('type=')[1].split()[0])
            kickno = int(line.split('vs[0]=')[1].split(',')[0])
            f1.write("%d %d %g %g %g %g %d %d\n" %(stype, stype1, v1, v2, v3, v, kickno, binary))
        elif line.rfind('birth kick(bin)')>-1: #For binary objects
            binary = 1
            v1 = float(line.split('vs[1]=')[1].split(',')[0])
            v1 = v1*kms
            v2 = float(line.split('vs[2]=')[1].split(',')[0])
            v2 = v2*kms
            v3 = float(line.split('vs[3]=')[1].split(',')[0])
            v3 = v3*kms
            v = (v1*v1 + v2*v2 + v3*v3)**0.5
            stype = int(line.split('type1=')[1].split()[0])
            stype1 = int(line.split('type2=')[1].split()[0])
            kickno = int(line.split('vs[0]=')[1].split(',')[0])
            f1.write("%d %d %g %g %g %g %d %d\n" %(stype, stype1, v1, v2, v3, v, kickno, binary))

    f.close()
    f1.close()

def plot_kicks(filename, plotfilename):
    data = np.loadtxt(filename)
    NSv = []
    BHv = []
    for i in range(len(data)):
        if data[i,7]==0.: #single system
            if data[i,0]==13.: #NS
                NSv.append(data[i,5])
            elif data[i,0]==14.: #BH
                BHv.append(data[i,5])
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(111)
    print BHv
    ax.hist(NSv, range=[0.1, 1e5], bins=5000, histtype='step', color='black', label='NS')
    ax.hist(BHv, range=[0.1, 1e5], bins=5000, histtype='step', color='blue', label='BH')
    ax.set_xlabel('SN kick velocity')
    ax.set_ylabel('dn/dv')
    ax.set_xscale('log')
    ax.set_xlim([1,1e5])

    ax.legend(loc='best')
    plt.savefig(plotfilename)
    plt.show()

    return NSv, BHv


def plot_SBP_L_comparison(filestrings, plotfilename='SBP.pdf'):
    colors = ['black', 'red', 'blue', 'gold', 'green', 'magenta']
    labels = ['0.01', '0.1', '1', '10']
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(211)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$r_{\rm{2D}}\ (\rm{pc})$', size=16)
    ax.set_ylabel(r'${\rm SBP}\ (L_\odot\ {\rm pc}^{-2})$', size=16)
    ax.set_xlim([1e-3, 1e2])
    ax.set_ylim([1e-3, 1e6])
    for i in range(len(filestrings)):
        filename = filestrings[i]+'_sbp.dat'
        data = np.loadtxt(filename)
        #ax.errorbar(data[:,0], data[:,1], data[:,2], ecolor='grey', barsabove=True, elinewidth=1, marker='o', mfc='none', mec=colors[i], ms=4, label=labels[i], ls='None') 
	ax.errorbar(data[:,0], data[:,1], data[:,2], ecolor='grey', barsabove=True, elinewidth=1, marker='o', mfc='red', mec='black', ms=6, ls='None')
    ax.legend(loc='best')
    plt.savefig(plotfilename)
    plt.show()


def plot_SBP_L_comparison_same_run_multisnap(filestring, snapnos, plotfilename):
    import scripts
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$r_{\rm{2D}}\,\rm{pc}$')
    ax.set_ylabel(r'$L/L_\odot$')
    for i in range(len(snapnos)):
        filename = filestring+'.snap'+snapnos[i]+'2D_sbp.dat'
	t_myr = scripts.find_t_myr(filestring, snapnos[i])
        data = np.loadtxt(filename)
        ax.errorbar(data[:,0], data[:,1], data[:,2], ecolor='grey', barsabove=True, elinewidth=1, marker='o', mfc='none', mec='black', ms=4, label=str(t_myr), ls='solid') 
    ax.legend(loc='best')
    plt.savefig(plotfilename)
    plt.show()
	

def plot_SBP_n_comparison(filestrings, plotfilename):
    colors = ['black', 'red', 'blue', 'gold', 'green', 'magenta']
    labels = ['0.01', '0.1', '1', '10']
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$r_{\rm{2D}}\,\rm{pc}$')
    ax.set_ylabel(r'$n/\rm{pc}^2$')
    for i in range(len(filestrings)):
        filename = filestrings[i]+'_sbp.dat'
        data = np.loadtxt(filename)
        ax.errorbar(data[:,0], data[:,3], data[:,4], ecolor='grey', barsabove=True, elinewidth=1, marker='o', mfc='none', mec=colors[i], ms=4, label=labels[i], ls='None') 
    ax.legend(loc='best')
    plt.savefig(plotfilename)
    plt.show()
       

def find_t_myr(fileloc, filestring, typestring, snapno):
	"""goes in the given snapshot and finds the physical time corresponding to that snap"""
	import gzip
	import scripts
	convfilestring = fileloc+'/'+filestring
	snapfile = fileloc+'/'+filestring+'.'+typestring+snapno+'.dat.gz'

	f=gzip.open(snapfile,'rb')
	line=f.readline()
	a=line.split()
	b=a[1]
	c=b.split('=')
	t=float(c[1])
	d=scripts.read_units(convfilestring)
	t_myr=t*d['t_myr'][0]
	f.close()
	return t_myr



def plot_evolution_of_BH_massdist(filestring, plotfilestring, typestring='bhinfo'):
	import glob
	import scripts
	import constants
	
	globstring = filestring+'.'+typestring+'*.dat.gz'
	files = glob.glob(globstring)
	print files
	import matplotlib.pyplot as plt
	

	for i in range(len(files)):
		print 'Working on', files[1]
		snapno = files[i].split(typestring)[1].split('.dat.gz')[0]
		t_myr = find_t_myr(filestring, typestring, snapno)	
		data = np.loadtxt(files[i], usecols=[1,7,8,9,14,17,18])
		mbh, mbhs, mbhb = [], [], []
		for j in range(len(data)):
			if data[j,1]==1.:
				if data[j,5]==14.:
					mbh.append(data[j,2])
					mbhb.append(data[j,2])
				if data[j,6]==14.:
					mbh.append(data[j,3])
					mbhb.append(data[j,3])
			elif data[j,1]==0.:
				if data[j,4]==14.:
					mbh.append(data[j,0])
					mbhs.append(data[j,0])
		plotfilename = filestring+'.'+typestring+snapno+'.mdist.pdf'
		print 'Now Plot'
		fig = plt.figure()
		ax = fig.add_subplot(111)
		ax.hist(mbh, bins=50, range=[0,100], normed=True, histtype='step', color='black', ls='solid', lw=1, label='all')
		ax.hist(mbhs, bins=50, range=[0,100], normed=True, histtype='step', color='red', ls='dashed', lw=1, label='single')
		ax.hist(mbhb, bins=50, range=[0,100], normed=True, histtype='step', color='black', ls='dotted', lw=2, label='binary')
		ax.set_title(r'$t = %.3f\ \rm{Myr}$' %(t_myr))
		ax.set_xlabel(r'$M_{\rm{BH}}\ (M_\odot)$')
		ax.set_ylabel(r'$dp/dM_{\rm{BH}}$')
		ax.set_xlim([0,50])
		ax.legend(loc='best', frameon=False, numpoints=1)
		plt.savefig(plotfilename)
		plt.close()
		
			
def MoverL(snapfilename, plotfilename='MoverL.pdf'):
	import scripts
	filestring = snapfilename.split('.snap')[0]
	units = scripts.read_units(filestring)
	data = np.loadtxt(snapfilename)
	Marray, Larray, rarray, startypearray = [], [], [], []
	for i in range(len(data)):
		if data[i,7] == 0.:
			Marray.append(data[i,1])
			Larray.append(data[i,15])
			rarray.append(data[i,2]*units[0]['l_pc'])
			startypearray.append(data[i,14])
		else:
			Marray.append(data[i,8]+data[i,9])
			Larray.append(data[i,19]+data[i,20])
			rarray.append(data[i,2]*units[0]['l_pc'])
			if data[i,17]>=data[i,18]:
				temptype = data[i,18]
			else:
				temptype = data[i,17]
			startypearray.append(temptype)
	#moverl, moverlerr, moverllow, moverlhigh = [], [], [], []
	#tempm, templ = Marray[0], Larray[0]
	#moverl.append(tempm/templ)
	#moverlerr.append(tempm/templ)
	#moverllow.append(0.)
	#moverlhigh.append(2.*tempm/templ)
	
	#for i in range(1,len(rarray)):
	#	tempm += Marray[i]
	#	templ += Larray[i]
	#	moverl.append(tempm/templ)
	#	temperr = tempm/templ/(i+1)**0.5
	#	moverlerr.append(temperr)
	#	moverllow.append(tempm/templ - temperr)
	#	moverlhigh.append(tempm/templ + temperr)
		
	#print len(moverllow), len(moverlhigh), len(rarray)
	
	#under sample
	r = [0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 20., 30., 40., 50., 60., 70.]
	
	print "calculating stuff"
	moverl, moverlerr, moverllow, moverlhigh, centraldist, countarr, brightcountarr = [], [], [], [], [], [], []
	for i in range(len(r)):
		tempm, templ, count, countbright = 0., 0., 0, 0
		for j in range(len(rarray)):
			if rarray[j]<=r[i]:
				tempm += Marray[j]
				templ += Larray[j]
				count += 1
				if startypearray[j]<13.:
					countbright += 1
		if count>1 and countbright>1:
			tempmoverl = tempm/templ
			moverl.append(tempmoverl)
			tempmoverlerr = tempmoverl/count**0.5
			moverlerr.append(tempmoverlerr)
			moverllow.append(tempmoverl-tempmoverlerr)
			moverlhigh.append(tempmoverl+tempmoverlerr)
			centraldist.append(r[i])
			countarr.append(count)
			brightcountarr.append(countbright)
				
	print len(centraldist), len(moverl), len(moverlerr)	
	print "plotting stuff"
	
	fig = plt.figure()
	ax = fig.add_subplot(311)
	ax.set_xlabel(r'$r\ ({\rmpc})$', size=16)
	ax.set_ylabel(r'$M/L\ (M_\odot/L_\odot)$', size=16)
        ax.plot(centraldist, moverl, ls='solid', color='black', lw=2)
	ax.fill_between(centraldist, moverllow, moverlhigh, color='red', alpha=0.3)
	
    	ax.set_xscale('log')
    	ax.set_yscale('log')
	ax.set_xlim([1e-2, 1e2])
	ax.set_ylim([0.5, 1e2])

    	plt.savefig(plotfilename)
        plt.show()
	return moverl, moverlerr, centraldist, countarr, brightcountarr
					

def VSIGMA(snapfilename, plotfilename='vsigma.pdf'):
	import scripts
	filestring = snapfilename.split('.snap')[0]
	units = scripts.read_units(filestring)
	data = np.loadtxt(snapfilename)
	Marray, Larray, rarray, startypearray, vrarray, vtarray = [], [], [], [], [], []
	for i in range(len(data)):
		if data[i,7] == 0.:
			Marray.append(data[i,1])
			Larray.append(data[i,15])
			rarray.append(data[i,2]*units[0]['l_pc'])
			startypearray.append(data[i,14])
			tempvr = data[i,3]*1e-5*units[0]['l_cgs']/units[0]['nbt_cgs']
			vrarray.append(tempvr)
			tempvt = data[i,4]*1e-5*units[0]['l_cgs']/units[0]['nbt_cgs']
			vtarray.append(tempvt)
			
		else:
			Marray.append(data[i,8]+data[i,9])
			Larray.append(data[i,19]+data[i,20])
			rarray.append(data[i,2]*units[0]['l_pc'])
			tempvr = data[i,3]*1e-5*units[0]['l_cgs']/units[0]['nbt_cgs']
			vrarray.append(tempvr)
			tempvt = data[i,4]*1e-5*units[0]['l_cgs']/units[0]['nbt_cgs']
			vtarray.append(tempvt)
			if data[i,17]>=data[i,18]:
				temptype = data[i,18]
			else:
				temptype = data[i,17]
			startypearray.append(temptype)
	#under sample
	r = [0., 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 20., 30., 40., 50., 60., 70.]
	
	print "calculating stuff"
	
	POS, POSBRIGHT, DISPVR, DISPVT, DISP, DISPERR, DISPVRBRIGHT, DISPVTBRIGHT, DISPBRIGHT, DISPBRIGHTERR = [], [], [], [], [], [], [], [], [], []
	for i in range(len(r)-1):
		vrs, vts, vrbrights, vtbrights = [], [], [], []
		count, countbright = 0, 0
		tempr = (r[i]+r[i+1])/2.
		for j in range(len(rarray)):
			if rarray[j]<r[i+1] and rarray[j]>=r[i]:
				vrs.append(vrarray[j])
				vts.append(vtarray[j])
				count += 1
				if startypearray[j]<11.: #Only brights
					countbright += 1
					vrbrights.append(vrarray[j])
					vtbrights.append(vtarray[j])
		if count>2:
			dispvr = np.std(vrs)
			dispvt = np.std(vts)
			POS.append(tempr)
			DISPVR.append(dispvr)
			DISPVT.append(dispvt)
			temp = (dispvr**2. + dispvt**2.)**0.5
			DISP.append(temp)
			temperr = temp/count**0.5
			DISPERR.append(temperr)
		if countbright>2:
			dispvrbright = np.std(vrbrights)
			dispvtbright = np.std(vtbrights)
			DISPVRBRIGHT.append(dispvrbright)
			DISPVTBRIGHT.append(dispvtbright)
			POSBRIGHT.append(tempr)
			temp = (dispvrbright**2. + dispvtbright**2.)**0.5
			DISPBRIGHT.append(temp)
			temperr = temp/countbright**0.5
			DISPBRIGHTERR.append(temperr)
		
		
	print "plotting stuff"
	print len(POS), len(DISPVR), len(DISPVT), len(POSBRIGHT), len(DISPVRBRIGHT), len(DISPVTBRIGHT)

	DISPLOW, DISPHIGH, DISPBRIGHTLOW, DISPBRIGHTHIGH = [], [], [], []
	for i in range(len(DISP)):
		DISPLOW.append(DISP[i]-DISPERR[i])
		DISPHIGH.append(DISP[i]+DISPERR[i])
	for i in range(len(DISPBRIGHT)):
		DISPBRIGHTLOW.append(DISPBRIGHT[i]-DISPBRIGHTERR[i])
		DISPBRIGHTHIGH.append(DISPBRIGHT[i]+DISPBRIGHTERR[i])
	
	fig = plt.figure()
	ax = fig.add_subplot(311)
	ax.set_xlabel(r'$r\ ({\rm pc})$', size=16)
	ax.set_ylabel(r'$v_\sigma\ ({\rm kms}^{-1})$', size=16)
        ax.plot(POSBRIGHT, DISPBRIGHT, ls='solid', color='black', lw=2)
	ax.fill_between(POSBRIGHT, DISPBRIGHTLOW, DISPBRIGHTHIGH, color='red', alpha=0.3)
	ax.plot(POS, DISP, ls='dashed', color='black', lw=2)
	ax.fill_between(POS, DISPLOW, DISPHIGH, color='green', alpha=0.3)
	
    	ax.set_xscale('log')
    	ax.set_yscale('log')
	ax.set_xlim([1e-2, 1e2])
	ax.set_ylim([1e-1, 1e1])

    	plt.savefig(plotfilename)
        plt.show()
	return DISPVR 


def run_make_2D_projection_allsnaps(filestring, seedy=100, proj=(0,1), snapgap=10):
	np.random.seed(seedy)
	units = scripts.read_units(filestring)
	lpc = units[0]['l_pc']
	globstring = filestring+'.snap*.dat.gz'
	snapfilelist = glob.glob(globstring)
	snapfilelist_sorted = np.sort(snapfilelist)
	for i in range(124, len(snapfilelist_sorted)):
		snapno = snapfilelist_sorted[i].split('snap')[1].split('.dat.gz')[0]
		print 'looking at snapno:', snapno
		if int(snapno)%snapgap==0:
			make_2D_projection(filestring, snapno, seedy=seedy, proj=proj)
	

def make_2D_projection(filestring, snapno, seedy=100, proj=(0,1)):
	np.random.seed(seedy)
	units = scripts.read_units(filestring)
	lpc = units[0]['l_pc']
	kms = 1e-5 * units[0]['l_cgs']/units[0]['nbt_cgs']
	t_myr = scripts.find_t_myr(filestring, snapno) 

	writefilename=filestring+'.snap'+snapno+'.2Dproj_old.dat'
	writefile=open(writefilename, 'w')
	writefile.write("#t=%g\n#1.r2D(pc) 2.Ltot(Lsun) 3.binflag 4.startype 5.L(Lsun) 6.startype0 7.startype1 8.L0(Lsun) 9.L1(Lsun) 10.Mtot(Msun) 11.M0(Msun) 12.M1(Msun) 13.vr(km/s) 14.vt(km/s) 15.id 16.id0 17.id1\n" %(t_myr))

	#read the snapfile
	snapfile = filestring+'.snap'+snapno+'.dat.gz'
	colnos = (2, 7, 14, 15, 17, 18, 19, 20, 1, 8, 9, 3, 4, 0, 10, 11)
	#0-r, 1-binflag 2-startype 3-L 4-startype0 5-startype1 6-L0 7-L1 8-Mtot 9-M0 10-M1
	data = np.loadtxt(snapfile, usecols=colnos)
	r2darr, binflagarr, startypearr, Larr, startype0arr, startype1arr, L0arr, L1arr, mtotarr, m0arr, m1arr = [], [], [], [], [], [], [], [], [], [], []
	vrarr, vtarr, idarr, id0arr, id1arr = [], [], [], [], []
	valid_line = 1
	for i in range(len(data)):
		try:
			for j in range(len(data[i])):
				if str(data[i,j])=='nan' or str(data[i,j])=='inf':
					valid_line = 0
					raise StopIteration()
				else:
					valid_line = 1
		except StopIteration:
			print "there was a nan or inf"
			pass
		
		if valid_line==1:
			rs = scripts.make_3D(data[i,0])
			r2d = 0
			for j in proj:
				r2d += rs[j]**2. 
			r2d = r2d**0.5
			r2darr.append(r2d)
			binflagarr.append(data[i,1])
			startypearr.append(data[i,2])
			Larr.append(data[i,3])
			startype0arr.append(data[i,4])
			startype1arr.append(data[i,5])
			L0arr.append(data[i,6])
			L1arr.append(data[i,7])
			mtotarr.append(data[i,8])
			m0arr.append(data[i,9])	
			m1arr.append(data[i,10])
			vrarr.append(data[i,11])
			vtarr.append(data[i,12])
			idarr.append(data[i,13])
			id0arr.append(data[i,14])
			id1arr.append(data[i,15])
			
		

	r2d_sorted = np.argsort(r2darr)
	for i in range(len(r2d_sorted)):
		if binflagarr[r2d_sorted[i]]==1.:
			tempL = L0arr[r2d_sorted[i]] + L1arr[r2d_sorted[i]]
		elif binflagarr[r2d_sorted[i]]==0.:
			tempL = Larr[r2d_sorted[i]]
		writefile.write("%g %g %d %d %g %d %d %g %g %g %g %g %g %g %d %d %d\n" %( r2darr[r2d_sorted[i]]*units[0]['l_pc'], tempL, binflagarr[r2d_sorted[i]], startypearr[r2d_sorted[i]], Larr[r2d_sorted[i]], startype0arr[r2d_sorted[i]], startype1arr[r2d_sorted[i]], L0arr[r2d_sorted[i]], L1arr[r2d_sorted[i]], mtotarr[r2d_sorted[i]], m0arr[r2d_sorted[i]], m1arr[r2d_sorted[i]], vrarr[r2d_sorted[i]]*kms, vtarr[r2d_sorted[i]]*kms, idarr[r2d_sorted[i]], id0arr[r2d_sorted[i]], id1arr[r2d_sorted[i]] ))
	writefile.close()


def get_sbp_from_2D_projection(filestring, snapno, BINNO=50, LCUT=15, NCUT=0.6):
	filename = filestring+'.snap'+snapno+'.2Dproj.dat.gz'
	print filename
	projfile = gzip.open(filename, 'r')
	projfile.seek(0)
	line = projfile.readline()
	print line.split('t=')[1].split()
	t_myr = float(line.split('t=')[1].split()[0])
	print projfile
	data = np.loadtxt(projfile)
	writefilename=filestring+'.snap'+snapno+'.2D_SBP_NCUT.dat'
	writefile=open(writefilename, 'w')
	writefile.write("#t=%g\n#1.r2Dlow(pc) 2.r2Dmid(pc) 3.r2Dhigh 4.Sigma(L/pc^2) 5.Sigmaerr(L/pc^2) 6.Sigma_n(1/pc^2) 7.Sigma_nerr(1/pc^2) 8.N_density(1/pc^2) 9.N_densityerr(1/pc^2)\n" %(t_myr))
	writefilename1=filestring+'.snap'+snapno+'.2D_SBPLcut'+str(LCUT)+'_NCUT.dat'
	writefile1=open(writefilename1, 'w')
	writefile1.write("#t=%g\n#1.r2Dlow(pc) 2.r2Dmid(pc) 3.r2Dhigh 4.Sigma(L/pc^2) 5.Sigmaerr(L/pc^2) 6.Sigma_n(1/pc^2) 7.Sigma_nerr(1/pc^2) 8.N_density(1/pc^2) 9.N_densityerr(1/pc^2)\n" %(t_myr))

	lr2d = np.log10(data[:,0])
	lbinsize = (lr2d[-1]-lr2d[0])/float(BINNO)
	print lbinsize
        mtot=data[:,9]; m0=data[:,10]; m1=data[:,11]; biflag=data[:,2]; ksin=data[:,3]; k0=data[:,5]; k1=data[:,6]
	n2d_prev = 0
	for i in range(1, BINNO+1):
		lsum, lsumerr, n2d, n2derr = 0., 0., 0., 0.
		lsumcut, lsumcuterr, n2dcut, n2dcuterr = 0., 0., 0., 0.
		lr_high, lr_low = lr2d[0]+i*lbinsize, lr2d[0]+(i-1)*lbinsize
		lr_mid = (lr_low+lr_high)/2.
		area = np.pi * ((10.**lr_high)**2. - (10.**lr_low)**2.)
                nd, ndcut = 0., 0.
		try:
			for j in range(int(n2d_prev), len(lr2d)):
				#print n2d_prev, lr_low, lr_high, lr2d[j]
				if lr2d[j]<lr_high and lr2d[j]>=lr_low:
					lsum = lsum + data[j,1]
					n2d = n2d + 1
                                        if NCUT!=-1 and biflag[j]!=1 and ksin[j]<10 and mtot[j]>=NCUT:
                                                nd += 1
                                        if NCUT!=-1 and biflag[j]==1 and ((k0[j]<10 and m0[j]>=NCUT) or (k1[j]<10 and m1[j]>=NCUT)):
                                                nd += 1
					if data[j,1]<LCUT:
						lsumcut += data[j,1]
						n2dcut += 1
                                                if NCUT!=-1 and biflag[j]!=1 and ksin[j]<10 and mtot[j]>=NCUT:
                                                        ndcut += 1 
                                                if NCUT!=-1 and biflag[j]==1 and ((k0[j]<10 and m0[j]>=NCUT) or (k1[j]<10 and m1[j]>=NCUT)):
                                                        ndcut += 1
                                                                                                                                          
				else:
					raise StopIteration()
		except StopIteration:
			print 'got value:n2d=', n2d, 'n2d_prev=', n2d_prev, 'nd=', nd
		n2d_prev += n2d
		if n2d>2:
			sbp, sbperr = lsum/area, lsum/float(n2d)**0.5/area	
			snp, snperr = n2d/area, float(n2d)**0.5/area
                        sdp, sdperr = nd/area, float(nd)**0.5/area
			writefile.write('%g %g %g %g %g %g %g %g %g\n' %(10**lr_low, 10**lr_mid, 10**lr_high, sbp, sbperr, snp, snperr, sdp, sdperr))
			if n2dcut>2:
				############3	
				sbpcut, sbpcuterr = lsumcut/area, lsumcut/float(n2d)**0.5/area	
				snpcut, snpcuterr = n2dcut/area, float(n2dcut)**0.5/area
                                sdpcut, sdpcuterr = ndcut/area, float(ndcut)**0.5/area                              
				################
				writefile1.write('%g %g %g %g %g %g %g %g %g\n' %(10**lr_low, 10**lr_mid, 10**lr_high, sbpcut, sbpcuterr, snpcut, snpcuterr, sdpcut, sdpcuterr))

	projfile.close()
	writefile.close()
	writefile1.close()
	return t_myr
	

def run_get_sbp_from_2D_projection(fileloc, runstring='initial', BINNO=50, LCUT=20, WHICHSNAP='last'):
	filestring = fileloc+'/'+runstring
	projfilestring = filestring+'.snap*.2Dproj.dat'
	projfilelist = glob.glob(projfilestring)
	projfilelist_sorted = np.sort(projfilelist)
	print projfilelist_sorted
	t = []
	if WHICHSNAP=='last':
		snapno = projfilelist_sorted[-1].split('/')[-1].split('snap')[1].split('.2Dproj')[0]
		t_myr = get_sbp_from_2D_projection(filestring, snapno, BINNO=BINNO, LCUT=LCUT)
		t.append(t_myr)
	elif WHICHSNAP=='all':
		for i in range(len(projfilelist_sorted)):
			snapno = projfilelist_sorted[i].split('/')[-1].split('snap')[1].split('.2Dproj')[0]
			print snapno
			t_myr = get_sbp_from_2D_projection(filestring, snapno, BINNO=BINNO, LCUT=LCUT)
			t.append(t_myr)
	return t


def read_time(filename):
	f = open(filename, 'r')
	line = f.readline()
	t = float(line.split('=')[1])
	f.close()
	return t

def make_SBP_many_figures(fileloc, runstring='initial', plotloc=''):
	filestring = fileloc+'/'+runstring
	sbpfilestring = filestring+'.snap*.2D_SBPLcut20.dat'
	sbpfiles = glob.glob(sbpfilestring)
	sbpfiles_sorted = np.sort(sbpfiles)
	for i in range(len(sbpfiles_sorted)):
		snapno = sbpfiles_sorted[i].split('/')[-1].split('snap')[1].split('.2D')[0]
		print snapno
		plotfilename = plotloc+'/test_'+str(i)+'.png'
		make_SBP_figures(fileloc, snapno, runstring=runstring, XLIM=[1e-2, 1e2], YLIM=[1e-3, 1e6], XSCALE='log', YSCALE='log', plotfilename=plotfilename)


def make_SBP_figures(fileloc, snapno, runstring='initial', XLIM=[1e-4, 1e2], YLIM=[1e-3, 1e9], XSCALE='log', YSCALE='log', plotfilename='test.png'):
	filestring = fileloc+'/'+runstring
	sbpfilename = filestring+'.snap'+snapno+'.2D_SBPLcut20.dat'
	dynfile = filestring+'.dyn.dat'
	dyndata = np.loadtxt(dynfile, usecols=(0,7))
	units = scripts.read_units(filestring)
	rc = dyndata[:,1]*units[0]['l_pc']
	tdyn = dyndata[:,0]*units[0]['t_myr']

	t = read_time(sbpfilename)
	fig = plt.figure()
	ax = fig.add_subplot(2, 1, 1)
	data = np.loadtxt(sbpfilename)
	ax.errorbar(data[:,1], data[:,3], xerr=(data[:,2]-data[:,0])/2., yerr=data[:,4], ecolor='grey', marker='o', mfc='blue', mec='blue', ls='None')
	ax.set_xlabel(r'$r_{\rm{2D}}$ (pc)')
	ax.set_ylabel(r'$\Sigma$ ($L_\odot\ \rm{pc}^{-2}$)')
	ax.set_xscale(XSCALE)
	ax.set_yscale(YSCALE)
	ax.set_xlim(XLIM)
	ax.set_ylim(YLIM)
	ax.xaxis.set_label_position('top')
	ax1 = fig.add_subplot(2, 1, 2)
	ax1.plot(tdyn, rc, ls='solid', color='black', lw=1)
	ax1.axvline(t, ls='dotted', lw=1, color='red')
	ax1.set_xlabel(r'$t$ (Myr)')
	ax1.set_ylabel(r'$r_c$ (pc)')
	ax1.set_xscale('log')
	ax1.set_xlim([10., 2e4])
	ax1.set_ylim([1e-1, 5.])
	plt.savefig(plotfilename)
	plt.close()
	#plt.show()	
	
	

def run_make_2D_projection(filestringlist, runstring='initial'):
	for i in range(len(filestringlist)):
		print 'working at file:', filestringlist[i]
		filestring = filestringlist[i]+'/'+runstring
		#print filestring
		globstring = filestring+'.snap*.dat.gz'
		snapfilelist = glob.glob(globstring)
		#print snapfilelist
		snapfilelist_sorted = np.sort(snapfilelist)
		snapno = snapfilelist_sorted[-1].split('snap')[1].split('.dat.gz')[0]
		make_2D_projection(filestring, snapno, seedy=100, proj=(0,1))
		 
		
def get_half_light_obs(filename2dproj, Lcut=20.):
	file2dproj = open(filename2dproj, 'r')
	file2dproj.seek(0)
	t_myr = float(file2dproj.readline().split('=')[1])
	data = np.loadtxt(file2dproj)
	Ltot = np.sum(data[:,1])
	#print Ltot
	Ltotcut = 0
	for i in range(len(data)):
		if data[i,1]<20.:
			Ltotcut += data[i,1]
	L, Lcut = 0., 0.
	Llowcut, Lhighcut, Llow, Lhigh = 0, 0, 0, 0
	success = 0
	i = 0
	while success<2 and i<len(data):
		Llow = L
		L += data[i,1]
		Lhigh = L
		#print L, Llow, Lhigh
		if data[i,1]<20.:
			Llowcut = Lcut
			Lcut += data[i,1]
			Lhighcut = Lcut
			#print L, Llowcut, Lhighcut
		if Llow<=0.5*Ltot and Lhigh>=0.5*Ltot:
			print 'came here', Llow, Lhigh, 0.5*Ltot
			hllow = data[i-1,0]
			hlhigh = data[i,0]
			success += 1
		if Llowcut<=0.5*Ltotcut and Lhighcut>=0.5*Ltotcut:
			print 'came here', Llowcut, Lhighcut, 0.5*Ltotcut
			hllowcut = data[i-1, 0]
			hlhighcut = data[i,0]
			success += 1
		i+=1
	file2dproj.close()

	return hllow, hlhigh, hllowcut, hlhighcut, t_myr
			
		
def run_get_half_light_obs(filestringlist, writefilename='half_light.dat', runstring='initial', runstringlist=[], Lcut=20.):
	writefile = open(writefilename, 'w')
	writefile.write("#1.hllow 2.hlhigh 3.hllowcut 4.hlhighcut 5.t_myr 6.runname\n")
	for i in range(len(filestringlist)):
		try:
			print 'working at file:', filestringlist[i]
			if len(runstringlist)==0:
				filestring = filestringlist[i]+'/'+runstring
			else:
				filestring = filestringlist[i]+'/'+runstringlist[i]
			globstring = filestring+'.*2Dproj.dat'
			print globstring
			globfilelist = glob.glob(globstring)
			globfilelist_sorted = np.sort(globfilelist)
			datafilename = globfilelist_sorted[-1]
			hllow, hlhigh, hllowcut, hlhighcut, t_myr = get_half_light_obs(datafilename, Lcut=Lcut)
			writefile.write("%g %g %g %g %g %s\n" %(hllow, hlhigh, hllowcut, hlhighcut, t_myr, filestringlist[i]))
		except IndexError:
			writefile.write("#%s\n" %(filestringlist[i]))
			pass
	
	writefile.close()	


def fit_king_cum_func(x, p0, p1):
	Lcum = np.pi * p0 * p1 * p1 * np.log( np.abs(1 + (x/p1)**2.))
	return Lcum

def fit_king_cum_func1(x,p0,p1,p2):
	"""p0 is the central density
	 p1 are p2 are rc^2 and rt^2, respectively. 
	x is actually square of r2d. 
	"""
	#print 'test:', 1.+(x/p1)
	#Now making sure that curve-fit does not accidentally passes some values that will generate a negative in log
	#x, p1, p2 = x**2., p1**2., p2**2.
	#x, p1, p2 = x**0.5, p1**0.5, p2**0.5
	p1, p2, p0 = np.abs(p1), np.abs(p2), np.abs(p0)
	Lcum = np.pi * p0 * p1 * ( np.log(1.+(x/p1)) - 4.*( (1.+(x/p1))**0.5 - 1.)/(1.+ (p2/p1) )**0.5 + (x/p1)/(1.+(p2/p1) ) )
	return Lcum

def fit_king_cum_curvefit(r2d, Lcum, p0guess=[1e5, 1.]):
	import scipy.optimize as opt
	p_opt, p_cov = opt.curve_fit(fit_king_cum_func, r2d, Lcum, p0=p0guess)
	return p_opt, p_cov

def fit_king_cum_curvefit1(r2dsquare, Lcum, p0guess=[1e5, 1., 30.]):
	import scipy.optimize as opt
	p_opt, p_cov = opt.curve_fit(fit_king_cum_func1, r2dsquare, Lcum, p0=p0guess)
	return p_opt, p_cov

def get_kingfit_params(filename2dproj, rhl, p0guess=[1e5, 1.]):
	file2dproj = open(filename2dproj, 'r')
	file2dproj.seek(0)
	t_myr = float(file2dproj.readline().split('=')[1])
	data = np.loadtxt(file2dproj)
	r2d = []
	L = []
	Lcum = []
	tempLcum = 0.
	i, success = 0, 0
	while i<len(data) and data[i,0]<=rhl:
		r2d.append(data[i,0])
		L.append(data[i,1])
		tempLcum += data[i,1]
		Lcum.append(tempLcum)
		i+=1
	print len(r2d), len(L), len(Lcum)
	r2d = np.array(r2d)
	L = np.array(L)
	Lcum = np.array(Lcum)
	p_opt, p_cov = fit_king_cum_curvefit(r2d, Lcum, p0guess=p0guess)

	file2dproj.close()	
	return r2d, Lcum, p_opt, p_cov, t_myr


def get_kingfit_params1(filename2dproj, rhl, p0guess=[1e5, 1., 900.]):
	file2dproj = open(filename2dproj, 'r')
	file2dproj.seek(0)
	t_myr = float(file2dproj.readline().split('=')[1])
	data = np.loadtxt(file2dproj)
	r2d = []
	L = []
	Lcum = []
	tempLcum = 0.
	i, success = 0, 0
	while i<len(data) and data[i,0]<=rhl*rhl:
		r2d.append(data[i,0])
		L.append(data[i,1])
		tempLcum += data[i,1]
		Lcum.append(tempLcum)
		i+=1
	print len(r2d), len(L), len(Lcum)
	r2d = np.array(r2d)
	r2dsquare = r2d*r2d
	L = np.array(L)
	Lcum = np.array(Lcum)
	p_opt, p_cov = fit_king_cum_curvefit1(r2dsquare, Lcum, p0guess=p0guess)

	file2dproj.close()	
	return r2d, Lcum, p_opt, p_cov, t_myr



def run_get_kingfit_params(filestringlist, writefilename='rcobs_rhocobs.dat', runstring='initial', runstringlist=[], p0guess=[1e5,1.], RHLFILENAME='half_light.dat', factor=0.5):
	writefile = open(writefilename, 'w')
	writefile.write("#1.sigmacobs(Msun/pc^2) 2.sigmacobserr 3.rcobs(pc) 4.rcobserr(pc) 5.t_myr 6.runname\n")
	halflightfile = open(RHLFILENAME, 'r')
	halflightfile.seek(0)
	for i in range(len(filestringlist)):
		try:
			print 'working at file:', filestringlist[i]
			#extract estimated half light radius
			halflightfile.seek(0)
			for line in halflightfile:
				if line.rfind(filestringlist[i])>-1 and line.rfind('#')==-1:
					rhl = float(line.split()[1])
			if len(runstringlist)==0:
				filestring = filestringlist[i]+'/'+runstring
			else:
				filestring = filestringlist[i]+'/'+runstringlist[i]
			globstring = filestring+'.*2Dproj.dat'
			globfilelist = glob.glob(globstring)
			globfilelist_sorted = np.sort(globfilelist)
			datafilename = globfilelist_sorted[-1]
			#Now find a good initial guess
			if len(p0guess)==0:
				lastsnapno = datafilename.split('/')[-1].split('snap')[1].split('.')[0]
				globstring1 = filestring+'.snap'+lastsnapno+'*SBPLcut20.dat'
				globfilelist1 = glob.glob(globstring1)
				if len(globfilelist1)==0:
					run_get_sbp_from_2D_projection(filestringlist[i], runstring='initial', BINNO=50, LCUT=20, WHICHSNAP='last')
					sbpfilename = filestring+'.snap'+lastsnapno+'.2D_SBPLcut20.dat'
				else:
					sbpfilename = np.sort(globfilelist1)[-1]
					print 'sbp file exists:', sbpfilename
				sbpdata = np.loadtxt(sbpfilename)
				rhocobsguess = np.max(sbpdata[:,3]) 
				rcobsguess = 0.5*rhl
				tmp_p0guess = [rhocobsguess, rcobsguess]
				print 'p0guess', tmp_p0guess
			else:
				tmp_p0guess = p0guess
							
			r2d, Lcum, p_opt, p_cov, t_myr = get_kingfit_params(datafilename, rhl, p0guess=tmp_p0guess)
			p_err = np.sqrt(np.diag(p_cov))
			sigmac, sigmacerr, rc, rcerr = np.abs(p_opt[0]), np.abs(p_err[0]), np.abs(p_opt[1]), np.abs(p_err[1])
			writefile.write("%g %g %g %g %g %s\n" %(sigmac, sigmacerr, rc, rcerr, t_myr, filestringlist[i]))
		except ():
			writefile.write("#%s\n" %(filestringlist[i]))
			pass
	
	writefile.close()
	halflightfile.close()	


def run_get_kingfit_params1(filestringlist, writefilename='rcobs_rhocobs.dat', runstring='initial', p0guess=[1e5,1., 900.], RHLFILENAME='half_light.dat'):
	writefile = open(writefilename, 'w')
	writefile.write("#1.sigmacobs(Msun/pc^2) 2.sigmacobserr 3.rcobs(pc) 4.rcobserr(pc) 5.rtobs(pc) 6.rtobserr(pc) 7.log(rtoverrc) 8.log(rtoverrc)_err 9.t_myr 10.runname\n")
	halflightfile = open(RHLFILENAME, 'r')
	halflightfile.seek(0)
	for i in range(len(filestringlist)):
		try:
			print 'working at file:', filestringlist[i]
			#extract estimated half light radius
			halflightfile.seek(0)
			for line in halflightfile:
				if line.rfind(filestringlist[i])>-1 and line.rfind('#')==-1:
					rhl = float(line.split()[1])
			filestring = filestringlist[i]+'/'+runstring
			globstring = filestring+'.*2Dproj.dat'
			globfilelist = glob.glob(globstring)
			globfilelist_sorted = np.sort(globfilelist)
			datafilename = globfilelist_sorted[-1]
			#print 'test:', datafilename, rhl, p0guess
			r2d, Lcum, p_opt, p_cov, t_myr = get_kingfit_params1(datafilename, rhl, p0guess=p0guess)
			p_err = np.sqrt(np.diag(p_cov))
			sigmac, sigmacerr = np.abs(p_opt[0]), np.abs(p_err[0])
			rc = p_opt[1]**0.5
			rcerr = 0.5 * p_err[1] / rc
			rt = p_opt[2]**0.5
			rterr = 0.5 * p_err[2] / rt
			logrtoverrc = np.log10(rt/rc)
			rtoverrcerr = (rt/rc) * (rcerr/rc + rterr/rt)
			logrtoverrcerr = (np.log10(rt/rc + rtoverrcerr) - np.log10(rt/rc - rtoverrcerr))/2.
				
			#sigmac, sigmacerr, rc, rcerr = p_opt[0], p_err[0], p_opt[1]**0.5, p_err[1]
			writefile.write("%g %g %g %g %g %g %g %g %g %s\n" %(sigmac, sigmacerr, rc, rcerr, rt, rterr, logrtoverrc, logrtoverrcerr, t_myr, filestringlist[i]))
		except (IndexError, RuntimeError):
			writefile.write("#%s\n" %(filestringlist[i]))
			pass
	
	writefile.close()
	halflightfile.close()



def get_t_disrupt(fileloc, runstring='initial'):
	outputfilename = fileloc+'/output.out'
	convfilename = fileloc+'/initial'
	units = scripts.read_units(convfilename)
	outputfile = open(outputfilename, 'r')
	outputfile.seek(0)
	tarr, marr, trharr = [], [], []
	for line in outputfile:
		if line.rfind('TotalTime')>-1 and line.rfind('T_MAX_PHYS')==-1:
			print line
			t = line.split('TotalTime=')[1].split()[0]
			tarr.append(float(t))
		if line.rfind('Mtotal')>-1:
			m = line.split('Mtotal=')[1].split()[0]
			marr.append(float(m))
		if line.rfind('trh=')>-1:
			trh = line.split('trh=')[1].split()[0]
			trharr.append(float(trh))
	tarr = np.array(tarr)
	marr = np.array(marr)
	trharr = np.array(trharr)
	tarr = tarr * units[0]['t_myr']
	marr = marr * units[0]['m_msun']
	trharr = trharr * units[0]['t_myr']
	
	i=0
	success = 0	
	while success==0 and i<len(tarr):
		mdot = (marr[i]-marr[i+1])/(tarr[i+1]-tarr[i])
		#mdotlim = marr[i+1]/trharr[i+1]
		mdotlim = 10.*marr[i+1]/units[0]['t_myr']
		print success, i, mdot, mdotlim, tarr[i], trharr[i], marr[i], marr[i+1]
		if mdot>=mdotlim:
			tdiss = tarr[i+1]
			success = 1
		i += 1
	outputfile.close()
	return tdiss, tarr, trharr
	


def compare_trh_with_massloss_timescale(fileloc, runstring='initial'):
	writefilename = fileloc+'/trh_massloss_timescale.dat'
	outputfilename = fileloc+'/output.out'
	convfilename = fileloc+'/initial'
	units = scripts.read_units(convfilename)
	outputfile = open(outputfilename, 'r')
	outputfile.seek(0)
	tarr, marr, trharr = [], [], []
	for line in outputfile:
		if line.rfind('TotalTime')>-1 and line.rfind('T_MAX_PHYS')==-1:
			#print line
			t = line.split('TotalTime=')[1].split()[0]
			tarr.append(float(t))
		if line.rfind('Mtotal')>-1:
			m = line.split('Mtotal=')[1].split()[0]
			marr.append(float(m))
		if line.rfind('trh=')>-1:
			trh = line.split('trh=')[1].split()[0]
			trharr.append(float(trh))
	tarr = np.array(tarr)
	marr = np.array(marr)
	trharr = np.array(trharr)
	tarr = tarr * units[0]['t_myr']
	marr = marr * units[0]['m_msun']
	trharr = trharr * units[0]['t_myr']
	print len(tarr), len(trharr), len(marr)
	count = np.min(np.array([len(tarr), len(trharr), len(marr)]))
	print count
	#mdotarr = np.zeros(len(tarr))
	mdotarr = np.zeros(count)
	#for i in range(1, len(tarr)):
	for i in range(1, count):
		#print i, marr[i], tarr[i], trharr[i]
		mdotarr[i] = (marr[i-1] - marr[i])/(tarr[i]-tarr[i-1])
	mdotarr[0] = mdotarr[1]

	writefile = open(writefilename, 'w')
	writefile.write("#1.t(Myr) 2.trh(Myr) 3.Mbound(Msun) 4.Mdot(Msun) 5.Mbound/Mdot(Myr)\n")
	for i in range(count):
		if mdotarr[i]>0.:
			writefile.write("%g %g %g %g %g\n" %(tarr[i], trharr[i], marr[i], mdotarr[i], marr[i]/mdotarr[i]))
	
	outputfile.close()
	writefile.close()
	return tarr, trharr, marr


def running_average(dataarray, averagewidth=100):
	dataavearr = []
	errarr = []
	nbins = int(len(dataarray)/averagewidth)
	if nbins<=10:
		print 'decrease average width. datalength:', len(dataarray), 'averagewidth:', averagewidth
	else:
		for i in range(nbins):
			startindex = i*averagewidth
			if (i+1)*averagewidth>len(dataarray):
				endindex = -1
			else:
				endindex = (i+1)*averagewidth
			
			dataave = np.mean(dataarray[startindex:endindex])
			dataavearr.append(dataave)
			err = np.std(dataarray[startindex:endindex])
			errarr.append(err)
	dataavearr = np.array(dataavearr)
	errarr = np.array(errarr)
	return dataavearr, errarr



def running_average_percentiles(dataarray, averagewidth=100, PERCENTILE_LIST=[2.3, 97.7,]):
	datamedarr = []
	datalowarr = []
	datahigharr = []
	nbins = int(len(dataarray)/averagewidth)
	if nbins<=10:
		print 'decrease average width. datalength:', len(dataarray), 'averagewidth:', averagewidth
	else:
		for i in range(nbins):
			startindex = i*averagewidth
			if (i+1)*averagewidth>len(dataarray):
				endindex = -1
			else:
				endindex = (i+1)*averagewidth
			
			data_med = np.median(dataarray[startindex:endindex])
			datamedarr.append(data_med)
			data_low = np.percentile(dataarray[startindex:endindex], PERCENTILE_LIST[0])
			data_high = np.percentile(dataarray[startindex:endindex], PERCENTILE_LIST[1])
			datalowarr.append(data_low)
			datahigharr.append(data_high)
	datamedarr = np.array(datamedarr)
	datalowarr = np.array(datalowarr)
	datahigharr = np.array(datahigharr)
	return datamedarr, datalowarr, datahigharr 



def find_t_disrupt(fileloc, runstring='initial', N=100):
	snapstring = fileloc+'/'+runstring+'.snap*.dat.gz'
	snapfiles = glob.glob(snapstring)
	snapfiles_sorted = np.sort(snapfiles)
	snaptimes, snapnos = [], []
	filestring = fileloc+'/'+runstring
	typestring = 'snap'
	for i in range(len(snapfiles_sorted)):
		snapno = snapfiles_sorted[i].split('snap')[1].split('.dat')[0]
		timesnap = find_t_myr(filestring, typestring, snapno)
		snaptimes.append(timesnap)
		snapnos.append(snapno)
	#last_snap = snapfiles_sorted[-1]
	#snapno = last_snap.split('snap')[1].split('.dat')[0]
	#filestring = fileloc+'/'+runstring
	#typestring = 'snap'
	#last_t_myr = find_t_myr(filestring, typestring, snapno)
	#dt = (np.log10(last_t_myr))/100
	dt = np.log10(snaptimes[-1]) / N
	timeboundaries = []
	for i in range(101):
   		timeboundaries.append(10.**(i*dt))
	timescaledatafile = fileloc+'/trh_massloss_timescale.dat'
	data = np.loadtxt(timescaledatafile)
	counter, avetimes, avetrh, avemdottimescale = running_average_at_fixed_times(data[:,0], data[:,1], data[:,4], timeboundaries)
	#print len(counter), len(avetimes), len(avetrh), len(avemdottimescale)

	writefilename = fileloc+'/trh_massloss_timescale_average_fixed_timewindow.dat'
	writefile = open(writefilename, 'w')
	writefile.write("#1.t(Myr) 2.trh(Myr) 3.mcl/mdot(Myr)\n")
	for i in range(len(avetimes)):
		writefile.write("%g %g %g\n" %(avetimes[i], avetrh[i], avemdottimescale[i]))
	writefile.close()

	#Find t_disrupt
	if avemdottimescale[-1]>avetrh[-1]: #Cluster model did not disrupt at integration stopping time
		dissolved = 0
		snapofinterest = snapnos[-1]
		t_myr_of_interest = snaptimes[-1]
		t_disrupt, tbounds = 12000., [12000., 12000.]
		print 'test', snapofinterest, t_myr_of_interest, t_disrupt
	else:
		dissolved = 1
		try:
			for i in range(len(avetimes)-1, 0, -1):
				if avemdottimescale[i]>avetrh[i]:
					t_disrupt = avetimes[i]
					raise StopIteration()
		except StopIteration:
			print 'found t_disrupt', t_disrupt
			pass
		for i in range(len(timeboundaries)-1):
			if timeboundaries[i+1]>=t_disrupt>timeboundaries[i]:
				tbounds = [timeboundaries[i], timeboundaries[i+1]]		
		
		#Find the closest snapshot number
		try:
			for i in range(len(snaptimes)-1):
				if snaptimes[i]<=t_disrupt<=snaptimes[i+1]:
					snapofinterest = snapnos[i]
					t_myr_of_interest = snaptimes[i]
					raise StopIteration()
		except StopIteration:
			print 'found snapno and snaptime', snaptimes[i], snaptimes[i+1], snapnos[i], snapnos[i+1]



	return t_disrupt, tbounds, snapofinterest, t_myr_of_interest, dissolved 

def running_average_at_fixed_times(tarray, dataarray, timeboundaries):
	#times = np.zeros(len(timeboundaries)-1)
	#data = np.zeros(len(timeboundaries)-1)
	times = []
	data = []
	#data1 = []
	#counter = np.zeros(len(timeboundaries)-1)
	counter = [0]
	j = -1
	for i in range(len(tarray)-1):
		for j in range(len(timeboundaries)):
			if tarray[i+1]>=timeboundaries[j] and tarray[i]<timeboundaries[j]:
				#counter[j] = int(i)
				counter.append(i)
	counter.append(-1)
	counter_cleaned = [0]
	for i in range(1, len(counter)):
		if counter[i]!=counter[i-1]:
			counter_cleaned.append(counter[i])

	for i in range(len(counter_cleaned)-1):
		times.append( np.mean( tarray[ counter_cleaned[i]:counter_cleaned[i+1] ] ) )
		data.append(np.mean(dataarray[ counter_cleaned[i]:counter_cleaned[i+1] ] ) )
		#data1.append(np.mean(data1array[ counter_cleaned[i]:counter_cleaned[i+1] ] ) )
		
	#return counter_cleaned, times, data, data1
	return counter_cleaned, times, data



def running_median_at_fixed_times(tarray_unsorted, dataarray_unsorted, timeboundaries, PERCENTILES=[2.3, 16., 50., 84., 97.7]):
	#times = np.zeros(len(timeboundaries)-1)
	#data = np.zeros(len(timeboundaries)-1)
	#times = []
	#data = []
	#data1 = []
	#percentiles = 
	#counter = np.zeros(len(timeboundaries)-1)
	ind = np.argsort(tarray_unsorted)
	tarray, dataarray = np.zeros(len(ind)), np.zeros(len(ind))
	for i in range(len(ind)):
		tarray[i] = tarray_unsorted[ind[i]]
		dataarray[i] = dataarray_unsorted[ind[i]]


	counter = [0]
	j = -1
	for i in range(len(tarray)-1):
		for j in range(len(timeboundaries)):
			if tarray[i+1]>=timeboundaries[j] and tarray[i]<timeboundaries[j]:
				#counter[j] = int(i)
				counter.append(i)
	counter.append(-1)
	counter_cleaned = [0]
	for i in range(1, len(counter)):
		if counter[i]!=counter[i-1]:
			counter_cleaned.append(counter[i])

	times = np.zeros(len(counter_cleaned)-1)
	data = np.zeros((len(counter_cleaned)-1, len(PERCENTILES)))
	#data1 = np.zeros((len(counter_cleaned)-1, len(PERCENTILES)))

	for i in range(len(counter_cleaned)-1):
		#times.append( np.mean( tarray[ counter_cleaned[i]:counter_cleaned[i+1] ] ) )
		#data.append(np.median(dataarray[ counter_cleaned[i]:counter_cleaned[i+1] ] ) )
		#data1.append(np.median(data1array[ counter_cleaned[i]:counter_cleaned[i+1] ] ) )

		times[i] = np.mean( tarray[ counter_cleaned[i]:counter_cleaned[i+1] ] )
		
		for j in range(len(PERCENTILES)):
			data[i,j] = np.percentile( dataarray[ counter_cleaned[i]:counter_cleaned[i+1] ], PERCENTILES[j] )
			#data1[i,j] = np.percentile( data1array[ counter_cleaned[i]:counter_cleaned[i+1] ], PERCENTILES[j] )


		
	return counter_cleaned, times, data
		

def get_t_disrupt_and_plot_t_vs_M(fileloclist, runstring='initial'):
	filestringlist = []
	tdiss = []
	for i in range(len(fileloclist)):
		print 'finding tdiss at:', fileloclist[i]
		temptdiss, tarr, trharr = get_t_disrupt(fileloclist[i], runstring=runstring)
		filestringlist.append(fileloclist[i]+'/initial')
		tdiss.append(temptdiss)
		plotfilename = fileloclist[i]+'/M_comparison.pdf'
		plot_t_vs_M(filestringlist, tdiss, plotfilename=plotfilename, YSCALE='linear', XSCALE='linear', XLIM=[0,1e4], YLIM=[0,1.2], MULTIPLE=2)
	
	return tdiss	

def get_mdot(fileloc, runstring='initial', NPOINTS=1000):
	filestring = fileloc+'/'+runstring
	dynfilename = filestring+'.dyn.dat'
	units = scripts.read_units(filestring)
	data = np.loadtxt(dynfilename, usecols=(0,4))
	dt = data[-1,0]/NPOINTS

	j = 1
	counter = [0,]
	for i in range(len(data)-1):
		if data[i,0]<j*dt and data[i+1,0]>=j*dt:
			counter.append(i)
			j += 1
	counter.append(len(data)-1)
	t, mdot = [], []
	for i in range(len(counter)-1):
		t.append( (data[counter[i], 0] + data[counter[i+1], 0])/2. )
		mdot.append( (data[counter[i], 1]-data[counter[i+1], 1])/dt  )
	t = np.array(t)
	mdot = np.array(mdot)
	t = t * units[0]['t_myr']
	mdot = mdot * units[0]['m_msun']/units[0]['t_myr']
	
	return t, mdot


def get_tidal_mdot(fileloc, runstring='initial', N=1000):
	filestring = fileloc+'/'+runstring
	escfilename = filestring+'.esc.dat'
	units = scripts.read_units(filestring)
	data = np.loadtxt(escfilename, usecols=(0,1,2))
	t, mesc, mdotesc = [], [], []
	tempm = 0.
	for i in range(len(data)-1):
		if data[i,0]==data[i+1,0]: #same time
			tempm += data[i,2]
		else: #different time
			t.append((data[i,1]+data[i+1,1])/2.)
			mesc.append(tempm)
			tempm = 0

	print np.shape(t)
	dt = t[-1]/N
	j = 1
	counter = []
	for i in range(len(t)-1):
		if t[i]<j*dt and t[i+1]>=j*dt:
			counter.append(i)
			j += 1
	tnew, mdot = [], []
	for i in range(len(counter)-1):
		tnew.append( (t[counter[i]] + t[counter[i+1]])/2. )
		mdot.append( mesc[i]/dt  )
	return t, mesc, tnew, mdot




def plot_t_vs_M_flattop(filename0, filename1, plotfilename='M_comparison.pdf', YSCALE='linear', XSCALE='log', XLIM=[1e-3,20], YLIM=[0,1.2]):
	import scripts
	fig = plt.figure()
    	fig.subplots_adjust(hspace=0.,wspace=0.)
    	colors = ['black', 'red', 'blue', 'gold', 'green', 'magenta']

	dynfile0 = filename0+'.dyn.dat'
	dynfile1 = filename1+'.dyn.dat'
    	data0 = np.loadtxt(dynfile0, usecols=(0,4))
	data1 = np.loadtxt(dynfile1, usecols=(0,4))
	units0 = scripts.read_units(filename0)
	units1 = scripts.read_units(filename1)
    
        ax = fig.add_subplot(1, 1, 1)
	ax.plot(data0[:,0]*units0[0]['t_myr']*1e-3, data0[:,1], lw=2, ls='solid', color='black')
	ax.plot(data1[:,0]*units1[0]['t_myr']*1e-3, data1[:,1], lw=2, ls='dashed', color='blue')

	ax.set_xlim(XLIM)
	ax.set_ylim(YLIM)
	ax.set_xscale(XSCALE)
    	ax.set_yscale(YSCALE)
		
	ax.set_xlabel(r'$t\ ({\rm Gyr})$', size=16)
        ax.set_ylabel(r'$M/M_i$', size=16)

    	plt.savefig(plotfilename)
    	plt.show()



def plot_t_vs_rc_slide(filestrings, plotfilename='slide_rc_comparison.pdf', YSCALE='linear', XSCALE='linear', XLIM=[0.,13.], YLIM=[0.,4.]):
    import scripts
    #plotfilename = 'rc_comparison.pdf'
    fig = plt.figure()
    fig.subplots_adjust(hspace=0.,wspace=0.)
    #XLIM = [0., 12.]
    #YLIM = [0., 4.]
    colors = ['black', 'red', 'blue', 'green', 'magenta']
    for i in range(len(filestrings)):
        print 'reading from', filestrings[i]
        units = scripts.read_units(filestrings[i])
        filename = filestrings[i]+'.dyn.pruned.dat'
        data = np.loadtxt(filename, usecols=(0,7,20))
        rc = data[:,1]*units[0]['l_pc']
	rh = data[:,2]*units[0]['l_pc']
        t = data[:,0]*units[0]['t_myr']*1e-3
    
        ax = fig.add_subplot(1, 1, 1)
	#ax.yaxis.tick_right()
	#ax.yaxis.set_label_position("right")

        ax.plot(t, rc, ls='solid', color=colors[i])
	ax.plot(t, rh, ls='dashed', color=colors[i])
	
	ax.set_xlim(XLIM)
	ax.set_ylim(YLIM)
	ax.set_xscale(XSCALE)
    	ax.set_yscale(YSCALE)
	
	ax.set_xlabel(r'$t\ ({\rm Gyr})$', size=18)
        ax.set_ylabel(r'$r_c,r_h\,\rm{(pc)}$', size=18)
    	

    plt.savefig(plotfilename)
    plt.show()


def plot_t_vs_rcoverrh_slide(filestrings, plotfilename='slide_rcoverrh_comparison.pdf', YSCALE='linear', XSCALE='linear', XLIM=[0.,13.], YLIM=[0.,1.]):
    import scripts
    #plotfilename = 'rc_comparison.pdf'
    fig = plt.figure()
    fig.subplots_adjust(hspace=0.,wspace=0.)
    #XLIM = [0., 12.]
    #YLIM = [0., 4.]
    colors = ['black', 'red', 'blue', 'green', 'magenta']
    for i in range(len(filestrings)):
        print 'reading from', filestrings[i]
        units = scripts.read_units(filestrings[i])
        filename = filestrings[i]+'.dyn.pruned.dat'
        data = np.loadtxt(filename, usecols=(0,7,20))
        rc = data[:,1]*units[0]['l_pc']
	rh = data[:,2]*units[0]['l_pc']
        t = data[:,0]*units[0]['t_myr']*1e-3
    
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(t, rc/rh, ls='solid', color=colors[i])
	#ax.plot(t, rc/rh, ls='dashed', color=colors[i])
	
	ax.set_xlim(XLIM)
	ax.set_ylim(YLIM)
	ax.set_xscale(XSCALE)
    	ax.set_yscale(YSCALE)
	
	ax.set_xlabel(r'$t\ ({\rm Gyr})$', size=18)
        ax.set_ylabel(r'$r_c/r_h$', size=18)
    	

    plt.savefig(plotfilename)
    plt.show()


def plot_t_vs_rcoverrh_averaged(filestrings, avewidths, plotfilename='slide_rcoverrh_averaged.pdf', YSCALE='linear', XSCALE='linear', XLIM=[0.,13.], YLIM=[0.,1.]):
    import scripts
    #plotfilename = 'rc_comparison.pdf'
    fig = plt.figure()
    fig.subplots_adjust(hspace=0.,wspace=0.)
    #XLIM = [0., 12.]
    #YLIM = [0., 4.]
    colors = ['black', 'red', 'blue', 'green', 'magenta']
    for i in range(len(filestrings)):
        print 'reading from', filestrings[i]
        units = scripts.read_units(filestrings[i])
        filename = filestrings[i]+'.dyn.dat'
        data = np.loadtxt(filename, usecols=(0,7,20))
        rc = data[:,1]*units[0]['l_pc']
	rh = data[:,2]*units[0]['l_pc']
        t = data[:,0]*units[0]['t_myr']*1e-3
	
	tave, terr = running_average(t, averagewidth=avewidths[i])
	rcave, rcerr = running_average(rc, averagewidth=avewidths[i])
	rhave, rherr = running_average(rh, averagewidth=avewidths[i])

	ratio = rcave/rhave
	ratiofracerr = rcerr/rcave + rherr/rhave
	ratioerr = ratio * ratiofracerr	
    
        ax = fig.add_subplot(1, 1, 1)
	#ax.yaxis.tick_right()
	#ax.yaxis.set_label_position("right")

        ax.plot(tave, ratio, ls='solid', lw=1.5, color=colors[i])
	ax.fill_between(tave, ratio-ratioerr, ratio+ratioerr, facecolor=colors[i], alpha=0.2)

	
	ax.set_xlim(XLIM)
	ax.set_ylim(YLIM)
	ax.set_xscale(XSCALE)
    	ax.set_yscale(YSCALE)
	
	ax.set_xlabel(r'$t\ ({\rm Gyr})$', size=18)
        ax.set_ylabel(r'$r_c/r_h$', size=18)
    	

    plt.savefig(plotfilename)
    plt.show()



def plot_t_vs_Nbh_slide(filestrings, plotfilename='slide_t_vs_Nbh_comparison.pdf', YSCALE='log', XSCALE='log', XLIM=[1,30000.], YLIM=[1., 5e3]):
    import scripts
    fig = plt.figure()
    fig.subplots_adjust(hspace=0.,wspace=0.)
    colors = ['black', 'red', 'blue', 'green', 'magenta']
    for i in range(len(filestrings)):
        units = scripts.read_units(filestrings[i])
        filename = filestrings[i]+'.bh.dat'
        print 'reading from', filename
	data = np.loadtxt(filename)
        Nbh = data[:,2]
	Nbh_sin = data[:,3]
	Nbh_bin = data[:,4]
        Nbh_bin_bhbh = data[:,5]
        Nbh_bin_bhnbh = data[:,6]
	if YSCALE=='log':
		Nbh += 1e-15
		Nbh_sin += 1e-15
		Nbh_bin += 1e-15
		Nbh_bin_bhbh += 1e-15
		Nbh_bin_bhnbh += 1e-15
        t = data[:,1]*units[0]['t_myr']*1e-3
   	print np.min(Nbh_bin_bhnbh), np.max(Nbh_bin_bhnbh) 
        ax = fig.add_subplot(111)
        ax.plot(t, Nbh, ls='solid', color=colors[i])
  
        ax.set_xlim(XLIM)
	ax.set_ylim(YLIM)
	ax.set_xscale(XSCALE)
    	ax.set_yscale(YSCALE)
	ax.set_xlabel(r'$t\ ({\rm Gyr})$', size=16)

	ax.set_ylabel(r'$N_{\rm{BH}}$', size=18)

    	ax.set_xlabel(r'$t\ ({\rm Gyr})$', size=18)

    plt.savefig(plotfilename)
    plt.show()


def plot_SBP_slide(filestrings, SIGMA0=[1532.65,114898.,1616.99], RCOBS=[3.93941,0.235572,3.20052], RHOBS=[], plotfilename='slide_SBP_comparison.pdf', YSCALE='log', XSCALE='log', XLIM=[.1,1e2], YLIM=[1., 1e6], N=[10,15], RCPOS0=[0.4, 0.75, 0.4], RCPOS1=[0.55, 0.85, 0.55], RHPOS0=[0.4, 0.75, 0.4], RHPOS1=[0.55, 0.85, 0.55]):
    import scripts
    fig = plt.figure()
    fig.subplots_adjust(hspace=0.,wspace=0., bottom=0.15)
    #colors = ['black', 'blue', 'red', 'green', 'magenta', 'cyan']
    colors = ['black', 'red', 'blue', 'green', 'magenta']
    markers = ['o', '^', 's', '*', 'd', 'p']
    linestyles = ['solid', 'dashed', 'dashdot', 'dotted',]
    for i in range(len(filestrings)):
        globstring = filestrings[i]+'.snap*.2D_SBPLcut20.dat'
	globfiles = np.sort(glob.glob(globstring))
	print globfiles[-1]
	data = np.loadtxt(globfiles[-1])
        
	ax = fig.add_subplot(1, 1, 1)
	#ax.errorbar(data[:,1], data[:,3], xerr=(data[:,2]-data[:,0])/2., yerr=data[:,4], ecolor='grey', marker=markers[i], markersize=4, mfc=colors[i], mec=colors[i], ls=linestyles[i], color=colors[i])
	ax.errorbar(data[:,1], data[:,3], xerr=(data[:,2]-data[:,0])/2., yerr=data[:,4], ecolor='grey', marker=markers[i], markersize=4, mfc=colors[i], mec=colors[i], ls='None')
	#ax.plot(data[:,1], data[:,3], ls=linestyles[i], color=colors[i])
	#xmin, xmax = data[N[0],1], data[N[1],1]
	xmin, xmax = data[N[i,0],1], data[N[i,1],1]
	print SIGMA0[i], xmin, xmax
	ax.axhline(SIGMA0[i], xmin=xmin, xmax=xmax, ls=linestyles[i], color=colors[i], lw=3, alpha=0.6)

	#for j in range(len(data)-1):
	#	if data[j,3]<RCOBS[i] and data[j+1,3]>RCOBS[i]:
	#		ymin, ymax = data[j,1], data[j+1,1]

	
	ax.axvline(RCOBS[i], ymin=RCPOS0[i], ymax=RCPOS1[i], ls=linestyles[i], color=colors[i], lw=1, alpha=0.8)
	print RCOBS[i], RCPOS0[i], RCPOS1[i]
	
	ax.axvline(RHOBS[i], ymin=RHPOS0[i], ymax=RHPOS1[i], ls=linestyles[i], color=colors[i], lw=1, alpha=0.8)
	#ax.axvline(RHOBS[i], ls=linestyles[i], color=colors[i], lw=4, alpha=0.4)
	print RHOBS[i], RCPOS0[i], RCPOS1[i]
	#ax.vline(RCOBS[i], )
	
	ax.set_xlabel(r'$r_{\rm{2D}}\ \rm{(pc)}$', size=18)
	ax.set_ylabel(r'$\Sigma$ ($L_\odot\ \rm{pc}^{-2}$)', size=18)
	ax.set_xscale(XSCALE)
	ax.set_yscale(YSCALE)
	ax.set_xlim(XLIM)
	ax.set_ylim(YLIM)

    plt.savefig(plotfilename)
    plt.show()


def plot_BH_bound_formed_slide(filestrings, plotfilename='slide_bh_form_bound.pdf', YSCALE='log', XSCALE='log', XLIM=[4e-3,20.], YLIM=[1,1e5]):
	import scripts
	fig = plt.figure()
	fig.subplots_adjust(hspace=0., wspace=0., bottom=0.15, top=0.99)
	colors = ['black', 'red', 'blue', 'green', 'magenta']
	for i in range(len(filestrings)):
		print 'working at:', filestrings[i]
		units = scripts.read_units(filestrings[i])
		bhformfile = filestrings[i]+'.bhformation.dat'
		bhboundfile = filestrings[i]+'.bh.dat'
		bhformdata = np.loadtxt(bhformfile, usecols=(0,))
		bhformdata = bhformdata*units[0]['t_myr']*1e-3
		bhbounddata = np.loadtxt(bhboundfile, usecols=(1,2))
		t = bhbounddata[:,0]*units[0]['t_myr']*1e-3
		
		ax = fig.add_subplot(1, 1, 1)
		hist, edges = np.histogram(bhformdata, bins=10000, range=[1e-3,20.])
		binmids = np.zeros(len(hist)+1)
		for j in range(len(edges)-1):
			binmids[j] = (edges[j+1]+edges[i])/2.
		binmids[-1] = t[-1]
		#hist = hist+1e-15
		cumhist1 = np.cumsum(hist)
		cumhist = np.zeros(len(cumhist1)+1)
		for j in range(len(cumhist1)):
			cumhist[j] = cumhist1[j]
		cumhist[-1] = cumhist1[-1]
		#boundbh = bhbounddata[:,1]+1e-15
		boundbh = bhbounddata[:,1]
		#ax.hist(bhformdata, bins=10000, range=[1e-3,20.], cumulative=True, ls='solid', color=colors[i], lw=2, histtype='step')
		ax.plot(binmids, cumhist, color=colors[i], ls='solid', lw=2)
		ax.plot(t, boundbh, color=colors[i], ls='solid', lw=1)


	ax.set_xscale(XSCALE)
	ax.set_yscale(YSCALE)
	ax.set_xlim(XLIM)
	ax.set_ylim(YLIM)
	ax.set_xlabel(r'$t\ \rm{(Gyr)}$', size=18)
	ax.set_ylabel(r'$N_{\rm{BH}}$', size=18)
	plt.savefig(plotfilename)
	plt.show()

	return t, boundbh
				


def plot_t_vs_M_slide(filestrings, plotfilename='slide_M_comp.pdf', YSCALE='linear', XSCALE='log', XLIM=[4e-3,20], YLIM=[1e-2,1.2]):
    import scripts
    #plotfilename = 'rc_comparison.pdf'
    fig = plt.figure()
    fig.subplots_adjust(hspace=0.,wspace=0., bottom=0.15, top=0.99)
    ax = fig.add_subplot(1, 1, 1)
    #colors = ['black', 'red', 'blue', 'gold', 'green', 'magenta']
    colors = ['black', 'red', 'blue', 'green', 'black', 'red', 'green', 'green', 'magenta']
    linestyles = ['solid', 'solid', 'solid', 'dashed', 'dashed', 'dashed']
    for i in range(len(filestrings)):
        print 'reading from', filestrings[i]
        units = scripts.read_units(filestrings[i])
        filename = filestrings[i]+'.dyn.dat'
        data = np.loadtxt(filename, usecols=(0,4))
        #M_msun = data[:,1]*units[0]['m_msun']
	M_code = data[:,1]
        t = data[:,0]*units[0]['t_myr']*1e-3
	#print t[-1], M_code[-1]
    
        ax.plot(t, M_code, ls=linestyles[i], color=colors[i], lw=1)
	
    ax.set_xlabel(r'$t\ ({\rm Gyr})$', size=18)
    ax.set_ylabel(r'$M/M_i$', size=18)
    ax.set_xscale(XSCALE)
    ax.set_yscale(YSCALE)
    ax.set_xlim(XLIM)
    ax.set_ylim(YLIM)
    plt.savefig(plotfilename)
    plt.show()


def plot_Mbh_hist_slide(filestrings, plotfilename='slide_Mbh_hist.pdf', YSCALE='linear', XSCALE='linear', BINNO=50, RANGE=[1,40]):
	import scripts
	fig = plt.figure()
	fig.subplots_adjust(hspace=0., wspace=0.)
	ax = fig.add_subplot(1,1,1)	
	colors = ['black', 'red', 'blue', 'gold', 'green', 'magenta']
	for i in range(len(filestrings)):
		print 'reading from:', filestrings[i]
		units = scripts.read_units(filestrings[i])
		filename = filestrings[i]+'.bhformation.dat'
		data = np.loadtxt(filename, usecols=(5,))
		ax.hist(data, bins=BINNO, range=RANGE, histtype='step', ls='solid', lw=1.5, normed=True, color=colors[i], cumulative=True)
	plt.savefig(plotfilename)
	plt.show()


def extract_semimajor_initial(filestrings):
	fig = plt.figure()
	fig.subplots_adjust(hspace=0., wspace=0.)
	ax = fig.add_subplot(1,1,1)
	colors = ['black', 'red', 'blue', 'gold', 'green', 'magenta']
	for i in range(len(filestrings)):
		snapfile = filestrings[i]+'.snap0000.dat.gz'
		data = np.loadtxt(snapfile, usecols=(7,8,9,12,1))
		lowa, higha = [], []
		for j in range(len(data)):
			if data[j,0]==1:
				
				if data[j,1]>15. or data[j,2]>15.:
					
					higha.append(data[j,2])
				else:
					#print 'came here'
					lowa.append(data[j,2])
		#ax.hist(lowa, bins=50, normed=1, histtype='step')
		#ax.hist(higha, bins=50, normed=1, histtype='step')
	#plt.show()
	return data, higha, lowa



def estimate_detectable_nbincore(filestring, snapno, rcobs, mlow_ms=0.4):
	"""
		Find only binary MS stars between turn-off and mlow_ms mass 
		with a companion with q=0.5. Double this number to get Nb_detect. Count 
		Total number of MS single or binary stars Nms. fb_detect = Nb_detect/Nms. 
	"""
	#First get core radius and half light radius
	filename2dproj = filestring+'/initial.snap'+str(snapno)+'.2Dproj.dat'
	data = np.loadtxt(filename2dproj)
	totalmscore, nbindetect, starcountcore, nbinactualcore = 0, 0, 0, 0
	i=0
	while data[i,0]<=rcobs:
		starcountcore += 1
		if data[i,2]==0. and data[i,5]<2. and data[i,9]>=mlow_ms:
			totalmscore += 1
		if data[i,2]==1.:
			nbinactualcore += 1
			if data[i,5]<2. and data[i,6]<2.:
				if data[i,10]>data[i,11]:
					m1 = data[i,10]
					m2 = data[i,11]
				else:
					m1 = data[i,11]
					m2 = data[i,10]
				if m1>=mlow_ms:
					totalmscore += 1
					if m1/m2>=0.5:
						nbindetect += 1
		i+=1
	nbindetect = nbindetect*2
	fbdetect = float(nbindetect)/float(totalmscore)
	return fbdetect, nbindetect, totalmscore, starcountcore, nbinactualcore
				
				
def find_bss_core(filestring, snapno, rcobs, z=0.001, mcut=1.1):
	#snapfile=filestring+'.snap'+snapno+'.dat.gz'
	#convfile=filestring+'.conv.sh'
	snapfilestring = filestring+'/initial'
	t_Myr = scripts.find_t_myr(snapfilestring, snapno)
	mguess = scripts.find_MS_turnoff(float(t_Myr))
	m = scripts.find_MS_TO(t_Myr, z, mguess)
	m_cut = mcut*m
	
	filename2dproj = filestring+'/initial.snap'+str(snapno)+'.2Dproj.dat'
	data = np.loadtxt(filename2dproj)
	i=0
	total_bss, sing_bss, bin_bss = 0, 0, 0
	while data[i,0]<=rcobs:
		if data[i,2]==0. and data[i,3]<2. and data[i,9]>=m_cut:
			total_bss += 1
			sing_bss += 1
		if data[i,2]==1:
			if data[i,5]<2. and data[i,10]>=m_cut:
				bin_bss += 1
				total_bss += 1
			if data[i,6]<2. and data[i,11]>=m_cut:
				bin_bss += 1
				total_bss += 1
		i+=1
	return total_bss, sing_bss, bin_bss, t_Myr
			


def find_NBH(filestring, snapno, rcobs):
	filename2dproj = filestring+'/initial.snap'+str(snapno)+'.2Dproj.dat'
	data = np.loadtxt(filename2dproj)
	NBH, NBHcore = 0, 0
	NBHsing, NBHsingcore = 0, 0
	NBHbin, NBHbincore = 0, 0
	for i in range(len(data)):
		if data[i,2]==0. and data[i,3]==14.:
			NBH += 1
			NBHsing += 1
			if data[i,0]<=rcobs:
				NBHcore += 1
				NBHsingcore += 1
		if data[i,2]==1:
			if data[i,5]==14.:
				NBH += 1
				NBHbin += 1
				if data[i,0]<=rcobs:
					NBHcore += 1
					NBHbincore += 1
			if data[i,6]==14.:
				NBH+=1
				NBHbin += 1
				if data[i,0]<=rcobs:
					NBHcore += 1
					NBHbincore += 1
			if data[i,6]==14. and data[i,7]==14.:
				NBHbin -=1
				if data[i,0]<=rcobs:
					NBHbincore -= 1

	return NBH, NBHcore, NBHsing, NBHsingcore, NBHbin, NBHbincore


def find_bsscore_nbincore(clusterpropfilename, writefilename, mlow_ms=0.4, z=0.001, mcut=1.1):
	f = open(clusterpropfilename, 'r')
	filestrings, rcobslist, t_myrlist, rhocobslist = [], [], [], [] 
	for line in f:
		if line.rfind('#')==-1: #uncommented lines
			a = line.split()
			filestrings.append(a[-1])
			rcobslist.append(np.abs(float(a[2])))
			rhocobslist.append(np.abs(float(a[0])))
			t_myrlist.append(float(a[4]))
	f.close()
	writefile = open(writefilename, 'w')
	writefile.write("#1.rcobs(pc) 2.rhocobs(Lsun/pc^2) 3.fb_detectable 4.Nbin_detectable 5.N_MS_core 6.N_core 7.Nbin_actual_core 8.N_BSScore 9.N_BSScore_sing 10.N_BSScore_bin 11.t(Myr) 12.NBH 13.NBHcore 14.NBHsing 15.NBHsingcore 16.NBHbin 17.NBHbincore 18.snapno 19.runname\n")
	writefile.write("#lower limit in main sequence mass is: %f\n" %(mlow_ms) )	
	for i in range(len(filestrings)):
		snapstring = filestrings[i]+'/initial.snap*.dat.gz'
		snapfiles_sorted = np.sort(glob.glob(snapstring))
		snapfile = snapfiles_sorted[-1]
		snapno = snapfile.split('snap')[1].split('.dat')[0]
		print 'looking at file:', filestrings[i], snapno
		#Get binary info
		
		fbdetect, nbindetect, totalmscore, starcountcore, nbinactualcore = estimate_detectable_nbincore(filestrings[i], snapno, rcobslist[i], mlow_ms=mlow_ms)
		print 'found binary info:', fbdetect, nbindetect, totalmscore, starcountcore
		#Get BSS info
		total_bss, sing_bss, bin_bss, t_Myr = find_bss_core(filestrings[i], snapno, rcobslist[i], z=z, mcut=mcut)
		print 'found BSS info:', total_bss, sing_bss, bin_bss, t_Myr
		NBH, NBHcore, NBHsing, NBHsingcore, NBHbin, NBHbincore = find_NBH(filestrings[i], snapno, rcobslist[i])
		print 'found BH info:', NBH, NBHcore, NBHsing, NBHsingcore, NBHbin, NBHbincore
		writefile.write("%f %f %f %d %d %d %d %d %d %d %f %d %d %d %d %d %d %s %s\n" %(rcobslist[i], rhocobslist[i], fbdetect, nbindetect, totalmscore, starcountcore, nbinactualcore, total_bss, sing_bss, bin_bss, t_Myr, NBH, NBHcore, NBHsing, NBHsingcore, NBHbin, NBHbincore, snapno, filestrings[i]))
	writefile.close()
		


def observed_mcore_vs_nbsscore(obsfilename='Harris_Leigh.txt', filename='core_bss_bh_binary.dat', nbhcorecut=10):
	f = open(obsfilename, 'r')
	f.seek(0)
	mcore, nbsscore, name = [], [], []
	for line in f:
		if line.rfind('#')==-1:
			a = line.split()
			name.append(a[0])
			temprho = float(a[-6])
			temprc = float(a[-12])
			tempr_rsun = float(a[9])
			rcpc = temprc*tempr_rsun*1e3*np.pi/180./60.
			rhoc = 10.**temprho * 2.
			mcore.append(rhoc*4.*np.pi*rcpc*rcpc*rcpc/3)
			nbsscore.append(float(a[-5]))
	mcoreselect, nbsscoreselect = [], []
	for i in range(len(name)):
		if name[i] == 'NGC104' or name[i]=='NGC6254' or name[i]=='NGC6266':
			mcoreselect.append(mcore[i])
			nbsscoreselect.append(nbsscore[i])
	f.close()
	data = np.loadtxt(filename, usecols=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17])
	mc, nbssc, mchigh, nbsschigh = [], [], [], []
	for i in range(len(data)):
		mc.append(np.pi*data[i,0]*data[i,0]*2.*data[i,1])
		nbssc.append(data[i,7])
		if data[i,12]>nbhcorecut:
			mchigh.append(np.pi*data[i,0]*data[i,0]*2.*data[i,1])
			nbsschigh.append(data[i,7])
			
			
	
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.plot(np.log10(mcore), np.log10(nbsscore), ms=8, marker='o', mfc='None', mec='red', ls='None', lw=2, label='all GC')
	ax.plot(np.log10(mcoreselect), np.log10(nbsscoreselect), ms=14, marker='^', mfc='red', mec='black', ls='None', lw=2, label='GC w. BH cand.')
	ax.plot(np.log10(mc), np.log10(nbssc), ms=8, marker='o', mfc='None', mec='blue', ls='None', lw=2, label='all model')
	ax.plot(np.log10(mchigh), np.log10(nbsschigh), ms=14, marker='^', mfc='blue', mec='black', ls='None', lw=2, label=r'model w. $N_{BH}>%d$' %(nbhcorecut))
	
	ax.set_xlabel(r'$\log (M_c/M_\odot)$')
	ax.set_ylabel(r'$\log N_{\rm{BSS},c}$')
	ax.legend(loc='best', numpoints=1, frameon=0)
	plt.savefig('Nbsscore_vs_Mcore_obs.pdf')
	plt.show()
	

def observed_fbincore_vs_nbsscore(obsfilename='Harris_Leigh.txt', filename='core_bss_bh_binary.dat', nbhcorecut=10):
	f = open(obsfilename, 'r')
	f.seek(0)
	Nbincore, nbsscore, name = [], [], []
	for line in f:
		if line.rfind('#')==-1:
			a = line.split()
			name.append(a[0])
			temprho = float(a[-6])
			temprc = float(a[-12])
			tempr_rsun = float(a[9])
			rcpc = temprc*tempr_rsun*1e3*np.pi/180./60.
			rhoc = 10.**temprho * 2.
			mcore = rhoc*4.*np.pi*rcpc*rcpc*rcpc/3
			Nbin = float(a[-2])*mcore/0.5
			Nbincore.append(Nbin)
			nbsscore.append(float(a[-5]))
	Nbincoreselect, nbsscoreselect = [], []
	for i in range(len(name)):
		if name[i] == 'NGC104' or name[i]=='NGC6254' or name[i]=='NGC6266':
			Nbincoreselect.append(Nbincore[i])
			nbsscoreselect.append(nbsscore[i])
	f.close()
	data = np.loadtxt(filename, usecols=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17])
	Nbc, nbssc, Nbchigh, nbsschigh = [], [], [], []
	for i in range(len(data)):
		Nbc.append(data[i,3])
		nbssc.append(data[i,7])
		if data[i,12]>nbhcorecut:
			Nbchigh.append(data[i,3])
			nbsschigh.append(data[i,7])
			
			
	
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.plot(np.log10(Nbincore), np.log10(nbsscore), ms=8, marker='o', mfc='None', mec='red', ls='None', lw=2, label='all GC')
	ax.plot(np.log10(Nbincoreselect), np.log10(nbsscoreselect), ms=14, marker='^', mfc='red', mec='black', ls='None', lw=2, label='GC w. BH cand.')
	ax.plot(np.log10(Nbc), np.log10(nbssc), ms=8, marker='o', mfc='None', mec='blue', ls='None', lw=2, label='all model')
	ax.plot(np.log10(Nbchigh), np.log10(nbsschigh), ms=14, marker='^', mfc='blue', mec='black', ls='None', lw=2, label=r'model w. $N_{BH}>%d$' %(nbhcorecut))
	
	ax.set_xlabel(r'$\log (N_b)$')
	ax.set_ylabel(r'$\log N_{\rm{BSS},c}$')
	ax.legend(loc='best', numpoints=1, frameon=0)
	plt.savefig('Nbsscore_vs_Nbc_obs.pdf')
	plt.show()	


def check_binary_prop(filestring):
	datafile = filestring+'.snap0000.dat.gz'
	data = np.loadtxt(datafile)
	high, high_b, low, low_b = 0, 0, 0, 0
	high_mratio, low_mratio = [], []
	high_a, low_a = [], []
	for i in range(len(data)):
		if data[i,7]==1.: #binary
			if data[i,8]>data[i,9]:
				m0 = data[i,8]
				m1 = data[i,9]
			else:
				m0 = data[i,9]
				m1 = data[i,8]
			if m0>15.:
				high += 1	
				high_b += 1
				high_mratio.append(m1/m0)
				high_a.append(data[i,12])
			elif m0<=15:
				low += 1
				low_b += 1
				low_mratio.append(m1/m0)
				low_a.append(data[i,12])
		else:
			if data[i,1]>15.:
				high += 1
			elif data[i,1]<=15:
				low += 1
	return float(high), float(high_b), float(low), float(low_b), high_mratio, low_mratio, high_a, low_a
				
			
		
		
def convert_to_3d(r, vr, vt):
	costheta = np.random.uniform(-1, 1)
	sintheta = (1-costheta**2.)**0.5
	phi = np.random.uniform(0, 4*np.pi)
	rz = r*costheta
	rx = r*sintheta*np.cos(phi)
	ry = r*sintheta*np.sin(phi)
	
	anglev = np.random.uniform(0., 4.*np.pi)
	magv = (vr*vr + vt*vt)**0.5	
	thetadot = np.cos(anglev) * vt/r
	phidot = np.sin(anglev)*vt/r
	
	vx = vr * np.sin(np.arccos(costheta)) * np.cos(phi) + r * thetadot * costheta * np.cos(phi) - r * phidot * np.sin(np.arccos(costheta)) * np.sin(phi)
	vy = vr * np.sin(np.arccos(costheta)) * np.sin(phi) + r * thetadot * costheta * np.sin(phi) + r * phidot * np.sin(np.arccos(costheta)) * np.cos(phi)
	vz = vr * costheta - r * thetadot * np.sin(np.arccos(costheta))

	r3d = np.array([rx, ry, rz])
	v3d = np.array([vx, vy, vz])

	return r3d, v3d


def chirp_mass(m1, m2):
	return (m1*m2)**(0.6) / (m1+m2)**(0.2)


def read_lastline(filename):
	f = open(filename, 'r')
	lastline = f.readlines()[-1]
	lastarr = lastline.split()
	f.close()
	return lastarr, lastline

def read_last_few_lines(filename, N=10):
	f = open(filename, 'r')
	lines = f.readlines()
	lastline = lines[-1]
	lastarr = lastline.split()
	prop_arr = np.zeros((N, len(lastarr)))
	for i in range(-N, 0, 1):
		temp = lines[i]
		temparr = temp.split()
		for j in range(len(temparr)):
			prop_arr[i+N,j] = float(temparr[j])
	f.close()
	return prop_arr


def collect_nbh_final(filestring):
	#relevant header
	##1:tcount  #2:TotalTime  #3:Nbh,tot  #4:Nbh,single  #5:Nbinarybh  #6:Nbh-bh  #7:Nbh-nonbh  #8:Nbh-ns  #9:N_bh-wd  #10:N_bh-star  #11:Nbh-ms  #12:Nbh-postms #13:fb_bh [(# binaries containing a bh)/(total # systems containing a bh)
	filename = filestring+'.bh.dat'
	lastarr, lastline = read_lastline(filename)
	Nbh = float(lastarr[2])
	Nbh_sin = float(lastarr[3])
	Nbh_bin = float(lastarr[4])
        Nbh_bin_bhbh = float(lastarr[5])
        Nbh_bin_bhnbh = float(lastarr[6])

	return Nbh, Nbh_sin, Nbh_bin, Nbh_bin_bhbh, Nbh_bin_bhnbh


def collect_nbh_final_stat(filestring, twindow=1000.):
	#relevant header
	##1:tcount  #2:TotalTime  #3:Nbh,tot  #4:Nbh,single  #5:Nbinarybh  #6:Nbh-bh  #7:Nbh-nonbh  #8:Nbh-ns  #9:N_bh-wd  #10:N_bh-star  #11:Nbh-ms  #12:Nbh-postms #13:fb_bh [(# binaries containing a bh)/(total # systems containing a bh)
	filename = filestring+'.bh.dat'
	units = scripts.read_units(filestring)
	data = np.genfromtxt(filename)
	tmyr = data[:,1]*units[0]['t_myr']
	tempnbh, tempnbhbh, tempnbhnbh, tempnbhb, tempnbhs = [], [], [], [], []
	for i in range(len(data)):
		if tmyr[i]>tmyr[-1]-twindow:
			tempnbh.append(data[i,2])
			tempnbhbh.append(data[:,5])
			tempnbhnbh.append(data[:,6])
			tempnbhb.append(data[:,4])
			tempnbhs.append(data[:,3])
	tempnbh = np.array(tempnbh)
	tempnbhbh = np.array(tempnbhbh)
	tempnbhnbh = np.array(tempnbhnbh)
	tempnbhb = np.array(tempnbhb)
	tempnbhs = np.array(tempnbhs)

	Nbh, Nbherr = np.mean(tempnbh), np.std(tempnbh)
	Nbh_sin, Nbh_sinerr = np.mean(tempnbhs), np.std(tempnbhs)
	Nbh_bin, Nbh_binerr = np.mean(tempnbhb), np.std(tempnbhb)
	Nbh_bin_bhbh, Nbh_bin_bhbherr = np.mean(tempnbhbh), np.std(tempnbhbh)
	Nbh_bin_bhnbh, Nbh_bin_bhnbherr = np.mean(tempnbhnbh), np.std(tempnbhnbh)

	#prop_arr = read_last_few_lines(filename, N=N)
	#Nbh, Nbherr = np.mean(prop_arr[:,2]), np.std(prop_arr[:,2])
	#Nbh_sin, Nbh_sinerr = np.mean(prop_arr[:,3]), np.std(prop_arr[:,3])
	#Nbh_bin, Nbh_binerr = np.mean(prop_arr[:,4]), np.std(prop_arr[:,4])
        #Nbh_bin_bhbh, Nbh_bin_bhbherr = np.mean(prop_arr[:,5]), np.std(prop_arr[:,5])
        #Nbh_bin_bhnbh, Nbh_bin_bhnbherr = np.mean(prop_arr[:,6]), np.std(prop_arr[:,6])

	return Nbh, Nbherr, Nbh_sin, Nbh_sinerr, Nbh_bin, Nbh_binerr, Nbh_bin_bhbh, Nbh_bin_bhbherr, Nbh_bin_bhnbh, Nbh_bin_bhnbherr


def collect_rc_rh_rhoc_vcrms(filestring):
	#relevant header
	##1:t #2:Dt #3:tcount #4:N #5:M #6:VR #7:N_c #8:r_c #9:r_max #10:Etot #11:KE #12:PE #13:Etot_int #14:Etot_bin #15:E_cenma #16:Eesc #17:Ebesc #18:Eintesc #19:Eoops #20:Etot+Eoops #21:r_h #22:rho_0 #23:rc_spitzer #24:v0_rms #25:rc_nb #26.DMse(MSUN) #27.DMrejuv(MSUN) #28.N_c_nb
	filename = filestring+'.dyn.dat'
	lastarr, lastline = read_lastline(filename)
	rc, rh, vc, rhoc, t = float(lastarr[7]), float(lastarr[20]), float(lastarr[23]), float(lastarr[21]), float(lastarr[0])
	return rc, rh, vc, rhoc, t


def collect_rc_rh_rhoc_vcrms_stat(filestring, N=10):
	#relevant header
	##1:t #2:Dt #3:tcount #4:N #5:M #6:VR #7:N_c #8:r_c #9:r_max #10:Etot #11:KE #12:PE #13:Etot_int #14:Etot_bin #15:E_cenma #16:Eesc #17:Ebesc #18:Eintesc #19:Eoops #20:Etot+Eoops #21:r_h #22:rho_0 #23:rc_spitzer #24:v0_rms #25:rc_nb #26.DMse(MSUN) #27.DMrejuv(MSUN) #28.N_c_nb
	filename = filestring+'.dyn.dat'
	prop_arr = read_last_few_lines(filename, N=N)
	rc, rcerr = np.mean(prop_arr[:,7]), np.std(prop_arr[:,7])
	rh, rherr = np.mean(prop_arr[:,20]), np.std(prop_arr[:,20])
	vc, vcerr = np.mean(prop_arr[:,23]), np.std(prop_arr[:,23])
	rhoc, rhocerr = np.mean(prop_arr[:,21]), np.std(prop_arr[:,21])
	t, terr = np.mean(prop_arr[:,0]), np.std(prop_arr[:,0])

	return rc, rcerr, rh, rherr, vc, vcerr, rhoc, rhocerr, t, terr


def collect_final_properties_all_sims(fileloclist, writefilename='final_dynamical_properties.dat'):
	writefile = open(writefilename, 'w')
	writefile.write("#1.t(Myr) 2.Nbh 3.Nbh_sin 4.Nbh_bin 5.Nbh_bin_bhbh 6.Nbh_bin_bhnbh 7.rc(pc) 8.rh(pc) 9.vc(km/s) 10.rhoc(Msun/pc^3) 11.filelocation\n")
	for i in range(len(fileloclist)):
		print fileloclist[i]
		filestring = fileloclist[i]+'/initial'
		print 'units'
		units = scripts.read_units(filestring)
		print 'rc, rh'
		rc, rh, vc, rhoc, t = collect_rc_rh_rhoc_vcrms(filestring)
		print 'Nbh'
		Nbh, Nbh_sin, Nbh_bin, Nbh_bin_bhbh, Nbh_bin_bhnbh = collect_nbh_final(filestring)
		rc = rc * units[0]['l_pc']
		rh = rh * units[0]['l_pc']
		vc = vc * 1e-5*units[0]['l_cgs']/units[0]['nbt_cgs']
		t = t * units[0]['t_myr']
		print 'writing'
		writefile.write("%g %d %d %d %d %d %g %g %g %g %s\n" %(t, Nbh, Nbh_sin, Nbh_bin, Nbh_bin_bhbh, Nbh_bin_bhnbh, rc, rh, vc, rhoc, fileloclist[i]))
	writefile.close()

		



def make_2D_projection_all_data(filestring, snapno, seedy=100, proj=(0,1)):
	np.random.seed(seedy)
	units = scripts.read_units(filestring)
	lpc = units[0]['l_pc']
	t_myr = scripts.find_t_myr(filestring, snapno) 

	writefilename=filestring+'.snap'+snapno+'.2Dproj_alldata.dat'
	writefile=open(writefilename, 'w')
	#writefile.write("#t=%g\n#1.r2D(pc) 2.Ltot(Lsun) 3.binflag 4.startype 5.L(Lsun) 6.startype0 7.startype1 8.L0(Lsun) 9.L1(Lsun) 10.Mtot(Msun) 11.M0(Msun) 12.M1(Msun)\n" %(t_myr))

	#read the snapfile
	snapfile = filestring+'.snap'+snapno+'.dat.gz'
	colnos = (2, 7, 14, 15, 17, 18, 19, 20, 1, 8, 9,)
	#0-r, 1-binflag 2-startype 3-L 4-startype0 5-startype1 6-L0 7-L1 8-Mtot 9-M0 10-M1
	data = np.loadtxt(snapfile, usecols=colnos)
	r2darr, binflagarr, startypearr, Larr, startype0arr, startype1arr, L0arr, L1arr, mtotarr, m0arr, m1arr = [], [], [], [], [], [], [], [], [], [], []
	valid_line = 1
	for i in range(len(data)):
		try:
			for j in range(len(data[i])):
				if str(data[i,j])=='nan' or str(data[i,j])=='inf':
					valid_line = 0
					raise StopIteration()
				else:
					valid_line = 1
		except StopIteration:
			print "there was a nan or inf"
			pass
		
		if valid_line==1:
			rs = scripts.make_3D(data[i,0])
			r2d = 0
			for j in proj:
				r2d += rs[j]**2. 
			r2d = r2d**0.5
			r2darr.append(r2d)
			binflagarr.append(data[i,1])
			startypearr.append(data[i,2])
			Larr.append(data[i,3])
			startype0arr.append(data[i,4])
			startype1arr.append(data[i,5])
			L0arr.append(data[i,6])
			L1arr.append(data[i,7])
			mtotarr.append(data[i,8])
			m0arr.append(data[i,9])	
			m1arr.append(data[i,10])
			
		

	r2d_sorted = np.argsort(r2darr)
	for i in range(len(r2d_sorted)):
		if binflagarr[r2d_sorted[i]]==1.:
			tempL = L0arr[r2d_sorted[i]] + L1arr[r2d_sorted[i]]
		elif binflagarr[r2d_sorted[i]]==0.:
			tempL = Larr[r2d_sorted[i]]
		writefile.write("%g %g %d %d %g %d %d %g %g %g %g %g\n" %( r2darr[r2d_sorted[i]]*units[0]['l_pc'], tempL, binflagarr[r2d_sorted[i]], startypearr[r2d_sorted[i]], Larr[r2d_sorted[i]], startype0arr[r2d_sorted[i]], startype1arr[r2d_sorted[i]], L0arr[r2d_sorted[i]], L1arr[r2d_sorted[i]], mtotarr[r2d_sorted[i]], m0arr[r2d_sorted[i]], m1arr[r2d_sorted[i]]))
	writefile.close()





def binary_bh_prop(filestrings, plotfilename, runstring='initial'):
    import scripts
    #plotfilename = 't_vs_Nbh_comparison.pdf'
    for i in range(len(filestrings)):
	tempstring = filestrings[i]+'/initial'
        units = scripts.read_units(tempstring)
	globstring = tempstring+'.bhinfo*.dat.gz'
	globlist = glob.glob(globstring)
	sorted_globlist = np.sort(globlist)
	filename = sorted_globlist[-1]
	snapno = filename.split('bhinfo')[1].split('.dat.gz')[0]
	t_myr = find_t_myr(tempstring, 'bhinfo', snapno)	
        print 'reading from', filename, snapno, t_myr
        #data = np.loadtxt(filename, usecols=(1,2,3,4,5,6,))
	data = np.loadtxt(filename)
	#[t, m1, m2, a, e, k1, k2, id1, id2]
	bhbh_m1, bhbh_m2 = [], []
	bhnbh_m1, bhnbh_m2 = [], [] 
	for i in range(len(data)):
		if data[i,7]==1: #binaries
			if data[i,17]==14. and data[i,18]==14.: #both BHs
				bhbh_m1.append(data[i,8])
				bhbh_m2.append(data[i,9])
			if (data[i,17]==14. and data[i,18]!=14.) or (data[i,17]!=14. and data[i,18]==14.): #BH-nBH
				bhnbh_m1.append(data[i,8])
				bhnbh_m2.append(data[i,9])
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.hist(bhbh_m1, bins=100, histtype='step', color='black', ls='solid')
	ax.hist(bhbh_m2, bins=100, histtype='step', color='blue', ls='solid')

	plt.show()



def copy_necessary_files(list_of_sims_file, filehandle):
	f = open(list_of_sims_file, 'r')
	f.seek(0)
	filestringlist = []
	for line in f:
		if line.rfind('#')==-1:
			filestringlist.append(line.split()[0])
	string = '/Volumes/SC_pspt2'
	filestringlist1 = []
	for i in range(len(filestringlist)):
		temp = filestringlist[i].split('CMC_results')[1]
		temp1 = string+temp
		filestringlist1.append(temp1)
		create_directory(temp1)
		#dynfile
		#sourcefile = filestringlist[i]+'/initial.'+filehandle
		sourcefile = filestringlist[i]+'/'+filehandle
		destfile = filestringlist1[i]
		copy_files(sourcefile, destfile)
	f.close()


def copy_necessary_lastsnapfiles(list_of_sims_file, filehandle):
	f = open(list_of_sims_file, 'r')
	f.seek(0)
	filestringlist = []
	for line in f:
		if line.rfind('#')==-1:
			filestringlist.append(line.split()[0])
	string = '/Volumes/SC_pspt2'
	filestringlist1 = []
	for i in range(len(filestringlist)):
		temp = filestringlist[i].split('CMC_results')[1]
		temp1 = string+temp
		filestringlist1.append(temp1)
		create_directory(temp1)
		#dynfile
		sourcefilelist = glob.glob(filestringlist[i]+'/initial.'+filehandle)
		sorted_sourcefilelist = np.sort(sourcefilelist)
		sourcefile = sorted_sourcefilelist[-1]
		destfile = filestringlist1[i]
		copy_files(sourcefile, destfile)


def copy_necessary_allsnapfiles(list_of_sims_file, filehandle):
	f = open(list_of_sims_file, 'r')
	f.seek(0)
	filestringlist = []
	for line in f:
		if line.rfind('#')==-1:
			filestringlist.append(line.split()[0])
	string = '/Volumes/SC_pspt2'
	filestringlist1 = []
	for i in range(len(filestringlist)):
		temp = filestringlist[i].split('CMC_results')[1]
		temp1 = string+temp
		filestringlist1.append(temp1)
		create_directory(temp1)
		#dynfile
		sourcefilelist = glob.glob(filestringlist[i]+'/initial.'+filehandle)
		sorted_sourcefilelist = np.sort(sourcefilelist)
		for j in range(len(sorted_sourcefilelist)):
			sourcefile = sorted_sourcefilelist[j]
			destfile = filestringlist1[i]
			copy_files(sourcefile, destfile)




def move_files(sourcefile, destfile):
	print "copying file", sourcefile, "to", destfile
	dataout,dataerr=subprocess.Popen([r"mv",sourcefile,destfile],stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()
	print "dataout", dataout
	print "dataerr", dataerr


def delete_files(sourcefile):
	print "deleting file", sourcefile
	dataout,dataerr=subprocess.Popen([r"rm",sourcefile],stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()
	print "dataout", dataout
	print "dataerr", dataerr	 

def create_directory(directory_name):
	print "creating directory:", directory_name
	dataout,dataerr=subprocess.Popen([r"mkdir","-p",directory_name],stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()
	print "dataout", dataout
	print "dataerr", dataerr


def copy_files(sourcefile, destfile):
	print "copying file", sourcefile, "to", destfile
	dataout,dataerr=subprocess.Popen([r"cp",sourcefile,destfile],stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()
	print "dataout", dataout
	print "dataerr", dataerr
		
			
			
def find_vsigmac(filename2dproj, rcobs):
	data = np.loadtxt(filename2dproj)
	varr = []
	for i in range(len(data)):
		if data[i,0]<rcobs:
			if data[i,2]==0.:
				if data[i,3]<10.:
					temp = (data[i,12]**2. + data[i,13]**2.)**0.5
					varr.append(temp)
			elif data[i,2]==1.:
				if data[i,5]<10. or data[i,6]<10.:
					temp = (data[i,12]**2. + data[i,13]**2.)**0.5
					varr.append(temp)
	varr = np.array(varr)
	vsigmac = np.std(varr)
	return vsigmac


def find_McoverLc(filename2dproj, rcobs):
	data = np.loadtxt(filename2dproj)
	marr = []
	ltot, mtot = 0., 0.
	count = 0
	try:
		while data[count,0]<rcobs:
			ltot += data[count,1]
			mtot += data[count,9]
			count += 1
			#print ltot, mtot, data[count, 0]
	except IndexError:
		ltot, mtot = ltot, mtot
		pass
	return ltot, mtot
			
			
				

def find_Nbssc_obs(filename2dproj, t_myr, rcobs, z=0.001, mcut=1.05):
	data = np.loadtxt(filename2dproj)
	N, Ns, Nb = 0, 0, 0
	mguess = scripts.find_MS_turnoff(t_myr) 
	m = scripts.find_MS_TO(t_myr, z, mguess)
	m_cut = mcut*m
	for i in range(len(data)):
		if data[i,0]<rcobs:
			if data[i,2]==0.:
				if data[i,3]<2. and data[i,9]>=m_cut:
					N += 1
					Ns += 1
			elif data[i,2]==1.:
				if data[i,5]<2. and data[i,10]>=m_cut:
					N += 1
					Nb += 1
				if data[i,6]<2. and data[i,11]>=m_cut:
					N += 1
					Nb += 1
	return N, Ns, Nb		
			
				

def collect_important_obs_data(fileloclist, rcobsfilename, rhobsfilename, writefilename):
	print 'getting rcobs'
	rcobs, rhocobs = [], []
	f = open(rcobsfilename, 'r')
	for i in range(len(fileloclist)):
		f.seek(0)
		temprc, temprhoc = 0., 0.
		for line in f:
			if line.rfind('#')==-1:
				a = line.split()
				#if line.rfind(fileloclist[i])>-1:
				if a[-1]==fileloclist[i]:
					temprc = float(a[2])
					temprhoc = float(a[0])
					rcobs.append(temprc)
					rhocobs.append(temprhoc)
	f.close()
	print 'obtained rcobs'
	#rh
	print 'getting rhobs'
	rhobs = []
	f1 = open(rhobsfilename, 'r')
	for i in range(len(fileloclist)):
		f1.seek(0)
		for line in f1:
			if line.rfind('#')==-1:
				a = line.split()
				if a[-1] == fileloclist[i]>-1:
					temprh = float(a[2])
					rhobs.append(temprh)
	f1.close()
	print 'obtained rhobs'
	#vsigmac, Nbssc
	print 'getting vsigmac'
	vsigmac, Nbssc, Nbssc_sing, Nbssc_bin = [], [], [], []
	Nbh, Nbh_sin, Nbh_bin, Nbh_bin_bhbh, Nbh_bin_bhnbh = [], [], [], [], []
	t_myr = []	
	for i in range(len(fileloclist)):
		print i, rcobs[i], fileloclist[i]
		temp = fileloclist[i]+'/initial.snap*.2Dproj.dat'
		temp1 = glob.glob(temp)
		temp2 = np.sort(temp1)
		filename2dproj = temp2[-1]
		tempvsigmac = find_vsigmac(filename2dproj, rcobs[i])
		vsigmac.append(tempvsigmac)
		temptmyr = read_time(filename2dproj)
		t_myr.append(temptmyr)
		print 't_myr=', temptmyr
		N, Ns, Nb = find_Nbssc_obs(filename2dproj, temptmyr, rcobs[i], z=0.001, mcut=1.05)
		Nbssc.append(N)
		Nbssc_sing.append(Ns)
		Nbssc_bin.append(Nb)
		filestring = fileloclist[i]+'/initial'
		tempNbh, tempNbh_sin, tempNbh_bin, tempNbh_bin_bhbh, tempNbh_bin_bhnbh = collect_nbh_final(filestring) 
		Nbh.append(tempNbh)
		Nbh_sin.append(tempNbh_sin)
		Nbh_bin.append(tempNbh_bin)
		Nbh_bin_bhbh.append(tempNbh_bin_bhbh)
		Nbh_bin_bhnbh.append(tempNbh_bin_bhnbh)
	
	writefile = open(writefilename, 'w')
	writefile.write("#1.rcobs(pc) 2.rhocobs(Lsun/pc^2) 3.rhobs(pc) 4.vsigmac(km/s) 5.Nbssc 6.Nbssc_sing 7.Nbssc_bin 8.Nbh 9.Nbh_sing 10.Nbh_bin 11.Nbh_bin_bhbh 12.Nbh_bin_bhnbh 13.t_myr 14.fileloc\n")
	for i in range(len(fileloclist)):
		writefile.write("%g %g %g %g %d %d %d %d %d %d %d %d %f %s\n" %(rcobs[i], rhocobs[i], rhobs[i], vsigmac[i], Nbssc[i], Nbssc_sing[i], Nbssc_bin[i], Nbh[i], Nbh_sin[i], Nbh_bin[i], Nbh_bin_bhbh[i], Nbh_bin_bhnbh[i], t_myr[i], fileloclist[i]))
	writefile.close()
	
	
	
	return rcobs, rhocobs, rhobs, vsigmac, Nbssc, Nbssc_sing, Nbssc_bin, Nbh, Nbh_sin, Nbh_bin, Nbh_bin_bhbh, Nbh_bin_bhnbh, t_myr


def collect_runfiles(list_of_runs_file):
	f = open(list_of_runs_file, 'r')
	f.seek(0)
	runlist = []
	for line in f:
		if line.rfind('#')==-1:
			a = line.split()
			runlist.append(a[0])
	f.close()	
	return runlist

def collect_runnames(list_of_runname_file):
	f = open(list_of_runname_file)
	f.seek(0)
	runnames = []
	for line in f:
		if line.rfind('#')==-1:
			a = line.split()
			runnames.append(a[0])
	f.close()
	return runnames

def find_last_snap(file_location, filestring):
	globstring = file_location+'/initial.snap*.dat.gz'
	#print globstring
	globfiles = glob.glob(globstring)
	sorted_globfiles = np.sort(globfiles)
	last_snapfile = sorted_globfiles[-1]
	last_snapno = last_snapfile.split('snap')[1].split('.dat.gz')[0]
	return last_snapfile, last_snapno


def get_dyn_props(file_location, filestring):
	f=file_location+'/'+filestring+'.dyn.dat'
	data=np.loadtxt(f)
	N_i, N_f = data[0,3], np.mean(data[-6:-1,3])
	M_i, M_f = data[0,4], np.mean(data[-6:-1,4])
	rc_i, rc_f = data[0,7], np.mean(data[-6:-1,7])
	rh_i, rh_f = data[0,20], np.mean(data[-6:-1,20])
	rhoc_i, rhoc_f = data[0,21], np.mean(data[-6:-1,21])
	t = data[-1,0]
	v0i, v0f = data[0,23], np.mean(data[-6:-1,23])
	return N_i, N_f, M_i, M_f, rc_i, rc_f, rh_i, rh_f, rhoc_i, rhoc_f, t, v0i, v0f


def get_M_sim(file_location, filestring):
	f=file_location+'/'+filestring+'.dyn.dat'
	data=np.loadtxt(f, usecols=(4,))
	convfile=file_location+'/'+filestring
	units=scripts.read_units(convfile)
	M_i, M_f = data[0]*units[0]['m_msun'], np.mean(data[-6:-1])*units[0]['m_msun']
	return M_i, M_f

def get_fb_sim(file_location, filestring):
	f=file_location+'/'+filestring+'.bin.dat'
	try:
		data=np.loadtxt(f, usecols=(10, 11,))
		fb_i, fb_f = data[0,1], np.mean(data[-6:-1,1])
		fbc_i, fbc_f = data[0,0], np.mean(data[-6:-1,0])
	except IOError:
		fb_i, fb_f = 0, 0
		fbc_i, fbc_f = 0, 0
	return fb_i, fb_f, fbc_i, fbc_f


def collect_cluster_props_table(list_of_runs_file, list_of_runname_file, filestring='initial', TWINDOW=5., DYNAMICS_SOURCE=0, WRITEFILENAME='new_plots/table_clusterprop.dat'):
	runnames = collect_runnames(list_of_runname_file)
	temp = collect_runfiles(list_of_runs_file)
	runlist = []
	if DYNAMICS_SOURCE==0:
		for i in range(len(temp)):
			temp1 = temp[i].split('CMC_results')[1]
			temp2 = '/Volumes/SC_pspt2'
			temp3 = temp2+temp1
			runlist.append(temp3)
	else:
		runlist = temp	
		
	writefile = open(WRITEFILENAME, 'w')
	count = 0
	for i in range(len(runlist)):
		count += 1
		#print runlist[i], filestring
		#last_snapfile, last_snapno = find_last_snap(runlist[i], filestring)
		#t_myr = find_t_myr(runlist[i], filestring, 'snap', last_snapno)
		convfile=runlist[i]+'/'+filestring
		units=scripts.read_units(convfile)
		N_i, N_f, M_i, M_f, rc_i, rc_f, rh_i, rh_f, rhoc_i, rhoc_f, t, v0i, v0f = get_dyn_props(runlist[i], filestring)
		M_i, M_f = M_i*units[0]['m_msun'], M_f*units[0]['m_msun']
		rc_i, rc_f, rh_i, rh_f = rc_i*units[0]['l_pc'], rc_f*units[0]['l_pc'], rh_i*units[0]['l_pc'], rh_f*units[0]['l_pc']
		dens = units[0]['m_msun']/(units[0]['l_pc'])**3.
		rhoc_i, rhoc_f = rhoc_i*dens, rhoc_f*dens
		t_myr = t*units[0]['t_myr']
		kms = units[0]['l_cgs']/units[0]['nbt_cgs']/1e5
		v0i, v0f = v0i*kms, v0f*kms
		print runlist[i], N_f/1e4, M_f/1e4, rc_f, rh_f, rhoc_f, t_myr, v0i, v0f
		fb_i, fb_f, fbc_i, fbc_f = get_fb_sim(runlist[i], filestring)


		writefile.write("%d & %s & %.1f & %.0f & %.0f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.1f & %.1f & %.1f & %.2f & %.1f \\\\ \n" %(count, runnames[i], t_myr/1e3, N_f/1e4, M_f/1e4, fb_f, fbc_f, rc_f, 0., rh_f, 0., rhoc_f/1e3, 0., v0f, 0., 0.))
	writefile.close()



def collect_rcobs_rhobs_sigmacobs_table(propertyfilename, runfilestring):
	rcobs, rhobs, sigmacobs = 0., 0., 0.
	f = open(propertyfilename, 'r')
	f.seek(0)
	try:
		for line in f:
			if line.rfind('#')==-1 and line.rfind(runfilestring)>-1:
				a = line.split()
				rcobs = float(a[0])
				rhobs = float(a[2])
				sigmacobs = float(a[1])
				raise StopIteration()
	except StopIteration:
		pass 
				
	return rcobs, rhobs, sigmacobs	

def collect_cluster_props_table1(list_of_runs_file, propertyfilename='new_plots/new_collected_obs_properties1.dat', filestring='initial', TWINDOW=5., DYNAMICS_SOURCE=0, WRITEFILENAME='new_plots/table2_clusterprop.dat'):
	#runnames = collect_runnames(list_of_runname_file)
	temp = collect_runfiles(list_of_runs_file)
	runlist = []
	if DYNAMICS_SOURCE==0:
		for i in range(len(temp)):
			temp1 = temp[i].split('CMC_results')[1]
			temp2 = '/Volumes/SC_pspt2'
			temp3 = temp2+temp1
			runlist.append(temp3)
	else:
		runlist = temp	
		
	writefile = open(WRITEFILENAME, 'w')
	count = 0
	for i in range(len(runlist)):
		count += 1
		#print runlist[i], filestring
		#last_snapfile, last_snapno = find_last_snap(runlist[i], filestring)
		#t_myr = find_t_myr(runlist[i], filestring, 'snap', last_snapno)
		convfile=runlist[i]+'/'+filestring
		units=scripts.read_units(convfile)
		N_i, N_f, M_i, M_f, rc_i, rc_f, rh_i, rh_f, rhoc_i, rhoc_f, t, v0i, v0f = get_dyn_props(runlist[i], filestring)
		M_i, M_f = M_i*units[0]['m_msun'], M_f*units[0]['m_msun']
		rc_i, rc_f, rh_i, rh_f = rc_i*units[0]['l_pc'], rc_f*units[0]['l_pc'], rh_i*units[0]['l_pc'], rh_f*units[0]['l_pc']
		dens = units[0]['m_msun']/(units[0]['l_pc'])**3.
		rhoc_i, rhoc_f = rhoc_i*dens, rhoc_f*dens
		t_myr = t*units[0]['t_myr']
		kms = units[0]['l_cgs']/units[0]['nbt_cgs']/1e5
		v0i, v0f = v0i*kms, v0f*kms
		print runlist[i], N_f/1e4, M_f/1e4, rc_f, rh_f, rhoc_f, t_myr, v0i, v0f
		fb_i, fb_f, fbc_i, fbc_f = get_fb_sim(runlist[i], filestring)

		rcobs, rhobs, sigmacobs = collect_rcobs_rhobs_sigmacobs_table(propertyfilename, runlist[i])

		if rcobs>0.:
			globstring = runlist[i]+'/initial*2Dproj*.dat'
			globfiles = glob.glob(globstring)
			#print globfiles
			filename2dproj = np.sort(globfiles)[-1]
			vsigmac = find_vsigmac(filename2dproj, rcobs)
			ltot, mtot = find_McoverLc(filename2dproj, rcobs)
			moverl = mtot/ltot	
			
		else:
			vsigmac, moverl = 0., 0


		writefile.write("%d & %.1f & %.0f & %.0f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.1f & %.1f & %.1f & %.2f & %.1f \\\\ \n" %(count, t_myr/1e3, N_f/1e4, M_f/1e4, fb_f, fbc_f, rc_f, rcobs, rh_f, rhobs, rhoc_f/1e3, sigmacobs/1e3, v0f, vsigmac, moverl))
		writefile.write("%s\n" %(runlist[i]))
	writefile.close()
		
		 
def collect_similar_runs(filelistname, writefilename):
	writefile = open(writefilename, 'w')
	f = open(filelistname, 'r')
	f.seek(0)
	similar_filelist = []
	rc, rhoc, nbh = [], [], []
	for line in f:
		if line.rfind('#')>-1:
			writefile.write("%s" %(line))
		elif line.rfind('RG3')==-1 and line.rfind('vink')==-1 and line.rfind('wind3')>-1 and line.rfind('#')==-1:
			a = line.split()
			if float(a[-3])>11000.:
				writefile.write("%s" %(line))
	f.close()
	writefile.close()


def collect_escaped_BHB_properties(filelocation, namestring):
	convfilename = filelocation+'/'+namestring
	units = scripts.read_units(convfilename)
	escfilename = filelocation+'/'+namestring+'.esc.dat'
	data = np.genfromtxt(escfilename)
	nbhbh = 0
	bhbh_a, bhbh_e, bhbh_m1, bhbh_m2, tesc, tinsp, tmerge = [], [], [], [], [], [], []
	bhbh_id1, bhbh_id2 = [], []
	for i in range(len(data)):
		if data[i,14]==1.:
			if data[i,22]==14. and data[i,23]==14.:
				nbhbh += 1
				bhbh_a.append(data[i,19]) 
				bhbh_e.append(data[i,20])
				tmpm1, tmpm2, tmpid1, tmpid2 = data[i,15], data[i,16], data[i,17], data[i,18]
				if tmpm1>=tmpm2:
					tmpmp = tmpm1
					tmpmpid = tmpid1
					tmpms = tmpm2
					tmpmsid = tmpid2
				else:
					tmpmp = tmpm2
					tmpmpid = tmpid2
					tmpms = tmpm1
					tmpmsid = tmpid1


				bhbh_m1.append(tmpmp)
				bhbh_m2.append(tmpms)
				bhbh_id1.append(tmpmpid)
				bhbh_id2.append(tmpmsid)
				temp_tesc = data[i,1]*units[0]['t_myr']
				tesc.append(temp_tesc)
				temp_tinsp = inspiral_time_peters(data[i,19], data[i,20], data[i,15], data[i,16])
				tinsp.append(temp_tinsp*1e3)
				temp_tmerge = temp_tesc+temp_tinsp
				tmerge.append(temp_tmerge)
	#snapstring = filelocation+'/'+namestring+'.snap*.dat.gz'
	#snaplist = glob.glob(snapstring)
	#sorted_snaplist = np.sort(snaplist)
	#last_snapfile = sorted_snaplist[-1]
	#last_snapno = last_snapfile.split('snap')[1].split('.dat')[0]
	#tfinal = find_t_myr(filelocation, namestring, 'snap', last_snapno)
					
	#bhbh_a = np.array(bhbh_a)
	#bhbh_e = np.array(bhbh_e)
	#bhbh_m1 = np.array(bhbh_m1)
	#bhbh_m2 = np.array(bhbh_m2)

	return nbhbh, bhbh_a, bhbh_e, bhbh_m1, bhbh_m2, tesc, tinsp, tmerge, bhbh_id1, bhbh_id2


def extract_metallicity(filename, Z=0.001):
	f = open(filename, 'r')
	f.seek(0)
	tmpZ = 0.
	try:
		for line in f:
			if line.rfind('OVERWRITE_Z')>-1:
				#print line
				a = line.split()
				tmpZ = float(a[1])
				if tmpZ==0.:
					tmpZ = Z
				raise StopIteration()
	except StopIteration:
		pass
	#print filename
	#print tmpZ
	return tmpZ
			

def collect_escaped_BHB_properties_all_runs(filelocationlist, jsondumpfilename, namestring='initial', namestringlist=[], Z=0.001):
	import json
	BHBHprops = {}
	for i in range(len(filelocationlist)):
		print filelocationlist[i]
		if len(namestringlist)==0:
			tmpnamestring = namestring
		else:
			tmpnamestring = namestringlist[i]
		convfilename = filelocationlist[i]+'/'+tmpnamestring
		units = scripts.read_units(convfilename)
		
		parsedfilename = filelocationlist[i]+'/'+tmpnamestring+'.cmc.parsed'
		Zvalue = extract_metallicity(parsedfilename, Z=Z)

    		nbhbh, bhbh_a, bhbh_e, bhbh_m1, bhbh_m2, tesc, tinsp, tmerge, bhbh_id1, bhbh_id2 = collect_escaped_BHB_properties(filelocationlist[i], tmpnamestring)
		N_i, N_f, M_i, M_f, rc_i, rc_f, rh_i, rh_f, rhoc_i, rhoc_f, t, v0_i, v0_f = get_dyn_props(filelocationlist[i], tmpnamestring)
		Mi, Mf = units[0]['m_msun']*M_i, units[0]['m_msun']*M_f
		rci, rcf = rc_i*units[0]['l_pc'], rc_f* units[0]['l_pc']
		rhi, rhf = rh_i*units[0]['l_pc'], rh_f* units[0]['l_pc']
		dens = units[0]['m_msun']/units[0]['l_pc']**3.
		rhoci, rhocf = rhoc_i*dens, rhoc_f*dens
		tfinal = t*units[0]['t_myr']
		kms = units[0]['l_cgs']/units[0]['nbt_cgs']/1e5
		v0i, v0f = v0_i*kms, v0_f*kms
		BHBHprops[i] = {'a': bhbh_a,
				'e': bhbh_e,
				'm1': bhbh_m1, 
				'm2': bhbh_m2, 
				'nbhbh': nbhbh,
				'tesc': tesc,
				'tinsp': tinsp, 
				'tmerge': tmerge,
				'loc': filelocationlist[i],
				'namestring': tmpnamestring,
				'id1': bhbh_id1,
				'id2': bhbh_id2,
				'tfinal': tfinal,
				'Ncl': [N_i, N_f],
				'Mcl': [Mi, Mf],
				'rc': [rci, rcf], 	
				'rh': [rhi, rhf], 
				'rhoc': [rhoci, rhocf],
				'v0': [v0i, v0f], 
				'Z': Zvalue, 
				'rv': units[0]['l_pc'] 
				}
	BHBHpropstring = json.dumps(BHBHprops)
	jsondumpfile = open(jsondumpfilename, 'w')
	jsondumpfile.write("%s" %(BHBHpropstring))
	jsondumpfile.close()
	return BHBHprops


def collect_merged_inside_cluster_BHB_properties(filelocation, namestring):
	convfilename = filelocation+'/'+namestring
	units = scripts.read_units(convfilename)
	filename = filelocation+'/'+namestring+'.semergedisrupt.log'
	f = open(filename, 'r')
	f.seek(0)
	t, bhbhm1, bhbhm2, bhbhid1, bhbhid2 = [], [], [], [], []
	nmerged = 0
	for line in f:
		if line.rfind('#')==-1:
			if line.rfind('type1=14')>-1 and line.rfind('type2=14')>-1 and line.rfind(':')>-1:
				tempt = float(line.split()[0].split('=')[1])
				tempt = tempt*units[0]['t_myr']
				t.append(tempt)
				tempm1 = float(line.split('m1=')[1].split(')')[0])
				tempm2 = float(line.split('m2=')[1].split(')')[0])
				tempid1 = int(line.split()[3].split('id1=')[1].split('(')[0])
				tempid2 = int(line.split()[3].split('id2=')[1].split('(')[0])

				if tempm1>=tempm2:
					tempmp = tempm1
					tempmpid = tempid1
					tempms = tempm2
					tempmsid = tempid2
				else:
					tempmp = tempm2
					tempmpid = tempid2
					tempms = tempm1
					tempmsid = tempid1

				bhbhm1.append(tempmp)
				bhbhm2.append(tempms)
				bhbhid1.append(tempmpid)
				bhbhid2.append(tempmsid)

				nmerged += 1
				
	f.close()
	return bhbhm1, bhbhm2, t, nmerged, bhbhid1, bhbhid2
	
				
def collect_mergerd_inside_cluster_BHB_properties_all_runs(filelocationlist, jsondumpfilename, namestring='initial', namestringlist=[], Z=0.001):
	import json
	BHBHprops = {}
	for i in range(len(filelocationlist)):
		print filelocationlist[i]
		if len(namestringlist)==0:
			tmpnamestring = namestring
		else:
			tmpnamestring = namestringlist[i]
		convfilename = filelocationlist[i]+'/'+tmpnamestring
		units = scripts.read_units(convfilename)

		parsedfilename = filelocationlist[i]+'/'+tmpnamestring+'.cmc.parsed'
		Zvalue = extract_metallicity(parsedfilename, Z=Z)

    		bhbhm1, bhbhm2, tmerge, nmerged, bhbhid1, bhbhid2 = collect_merged_inside_cluster_BHB_properties(filelocationlist[i], tmpnamestring)

		N_i, N_f, M_i, M_f, rc_i, rc_f, rh_i, rh_f, rhoc_i, rhoc_f, t, v0_i, v0_f = get_dyn_props(filelocationlist[i], tmpnamestring)
		Mi, Mf = units[0]['m_msun']*M_i, units[0]['m_msun']*M_f
		rci, rcf = rc_i*units[0]['l_pc'], rc_f* units[0]['l_pc']
		rhi, rhf = rh_i*units[0]['l_pc'], rh_f* units[0]['l_pc']
		dens = units[0]['m_msun']/units[0]['l_pc']**3.
		rhoci, rhocf = rhoc_i*dens, rhoc_f*dens
		tfinal = t*units[0]['t_myr']
		kms = units[0]['l_cgs']/units[0]['nbt_cgs']/1e5
		v0i, v0f = v0_i*kms, v0_f*kms

		BHBHprops[i] = {'m1': bhbhm1, 
				'm2': bhbhm2, 
				'nmerged': nmerged,
				'tmerge': tmerge,
				'loc': filelocationlist[i],
				'namestring': tmpnamestring,
				'id1': bhbhid1,
				'id2': bhbhid2,
				'tfinal': tfinal,
				'Ncl': [N_i, N_f],
				'Mcl': [Mi, Mf],
				'rc': [rci, rcf], 	
				'rh': [rhi, rhf], 
				'rhoc': [rhoci, rhocf],
				'v0': [v0i, v0f], 
				'Z': Zvalue, 
				'rv': units[0]['l_pc']
				}
	BHBHpropstring = json.dumps(BHBHprops)
	jsondumpfile = open(jsondumpfilename, 'w')
	jsondumpfile.write("%s" %(BHBHpropstring))
	jsondumpfile.close()
	return BHBHprops
				

	
def inspiral_time_peters(a0,e0,m1,m2):
    	"""
    	Computes the inspiral time, in Gyr, for a binary
	    a0 in Au, and masses in solar masses
    	"""

    	coef = 6.086768e-11 #G^3 / c^5 in au, gigayear, solar mass units
    	beta = (64./5.) * coef * m1 * m2 * (m1+m2)

    	if e0 == 0:
        	return a0**4 / (4*beta)

    	c0 = a0 * (1.-e0**2.) * e0**(-12./19.) * (1.+(121./304.)*e0**2.)**(-870./2299.)

    	time_integrand = lambda e: e**(29./19.)*(1.+(121./304.)*e**2.)**(1181./2299.) / (1.-e**2.)**1.5
    	integral,abserr = scipy.integrate.quad(time_integrand,0,e0)

    	return integral * (12./19.) * c0**4. / beta


def beta(m1,m2):
	coef = 6.086768e-11 #G^3 / c^5 in au, gigayear, solar mass units	
	#coef =  6.086768e-20 #G^3 / c^5 in au, year, solar mass units	
	G = 3.96511851e-14 #G in au, solar mass, second units 
	return (64./5.) * coef * m1 * m2 * (m1+m2)

def c0(a0,e0):
	return a0 * (1.-e0**2.) * e0**(-12./19.) * (1.+(121./304.)*e0**2.)**(-870./2299.)

def time_integrand(e):
	return e**(29./19.)*(1.+(121./304.)*e**2.)**(1181./2299.) / (1.-e**2.)**1.5				
				
				
		
def extract_statistics(xarr, yarr, xbinedges = [0,1,2,3,4,5,6,7,8,9,10,11,12]):
	#if len(xbinedges)==0:
	#	xbinedges = np.linspace(np.min(xarr), np.max(xarr), NBINS)
	tempdict = {}
	count = 0
	for i in range(len(xbinedges)-1):
		tempy, n = [], 0
		for j in range(len(yarr)):
			if xbinedges[i]<=xarr[j]<xbinedges[i+1]:
				tempy.append(yarr[j])
				n+=1
		tempdict[count] = {'x': (xbinedges[i+1]+xbinedges[i])/2.,
				'y': tempy,
				'n': n
				}
		count += 1
	return tempdict
			
			
def plot_t_vs_nbhbh_averaged(filestrings, avewidths, plotfilename='t_vs_nbhbh_averaged.pdf', YSCALE='linear', XSCALE='linear', XLIM=[0.,13.], YLIM=[0.,40.]):
    import scripts
    #plotfilename = 'rc_comparison.pdf'
    fig = plt.figure()
    fig.subplots_adjust(hspace=0.,wspace=0., bottom=0.15, top=0.99)
    #XLIM = [0., 12.]
    #YLIM = [0., 4.]
    colors = ['black', 'red', 'blue', 'green', 'magenta']
    for i in range(len(filestrings)):
        print 'reading from', filestrings[i]
        units = scripts.read_units(filestrings[i])
	
	filename = filestrings[i]+'.bh.dat'
        print 'reading from', filename
	data = np.loadtxt(filename)
        Nbh = data[:,2]
	Nbh_sin = data[:,3]
	Nbh_bin = data[:,4]
        Nbh_bin_bhbh = data[:,5]
        Nbh_bin_bhnbh = data[:,6]
	t = data[:,1]*units[0]['t_myr']/1e3

	tave, terr = running_average(t, averagewidth=avewidths[i])
	bhbh, bhbherr = running_average(Nbh_bin_bhbh, averagewidth=avewidths[i])
	bhnbh, bhnbherr = running_average(Nbh_bin_bhnbh, averagewidth=avewidths[i])
	
        ax = fig.add_subplot(2, 1, 1)
	#ax.yaxis.tick_right()
	#ax.yaxis.set_label_position("right")

        ax.plot(tave, bhbh, ls='solid', lw=1.5, color=colors[i])
	low = np.zeros(len(tave))
	high = np.zeros(len(tave))
	for j in range(len(tave)):
		low[j] = bhbh[j]-bhbherr[j]
		if low[j]<=0.:
			low[j] = 1e-15
		high[j] = bhbh[j]+bhbherr[j]

	ax.fill_between(tave, low, high, facecolor=colors[i], alpha=0.2)
	ax.set_ylabel(r'$N_{\rm{BH-BH}}$', size=18)
	
	
	ax.set_xlim(XLIM)
    	ax.set_ylim(YLIM)
    	ax.set_xscale(XSCALE)
    	ax.set_yscale(YSCALE)
	ax.set_xticklabels([])
	ax.set_xticks([])

	
	ax1 = fig.add_subplot(2, 1, 2)
	#ax.yaxis.tick_right()
	#ax.yaxis.set_label_position("right")

        ax1.plot(tave, bhnbh, ls='solid', lw=1.5, color=colors[i])
	low1 = np.zeros(len(tave))
	high1 = np.zeros(len(tave))
	for j in range(len(tave)):
		low1[j] = bhnbh[j]-bhnbherr[j]
		if low1[j]<=0.:
			low1[j] = 1e-15
		high1[j] = bhnbh[j]+bhnbherr[j]
	ax1.fill_between(tave, low1, high1, facecolor=colors[i], alpha=0.2)
	#ax1.fill_between(tave, np.max([bhnbh-bhnbherr,0.]), np.max([bhnbh+bhnbherr,0.]), facecolor=colors[i], alpha=0.2)
	ax1.set_ylabel(r'$N_{\rm{BH-nBH}}$', size=18)
	#ax.set_xticklabels([])

	ax1.set_xlim(XLIM)
    	ax1.set_ylim(YLIM)
    	ax1.set_xscale(XSCALE)
    	ax1.set_yscale(YSCALE)

    ax1.set_xlabel(r'$t\ ({\rm Gyr})$', size=18)
    #ax.set_ylabel(r'$N_{\rm{BH-BH}}$', size=18)
    	

    plt.savefig(plotfilename)
    plt.show()



def extract_property(filename, runfilename):
	f = open(filename, 'r')
	f.seek(0)
	try:
		for line in f:
			if line.rfind(runfilename) and line.rfind('#')==-1:
				a = line.split()
				nbh = float(a[7])
				nbhbh = float(a[10])
				nbhnbh = float(a[11])
				raise StopIteration()
	except StopIteration:
		print 'found nbh props:', runfilename
		pass
	f.close()
	return nbh, nbhbh, nbhnbh

def extract_property_from_json_escaped(jsonfilename, runfilename, z = [4058.3457426100003, 1477.3043937000007,]):
	f = open(jsonfilename, 'r')
	escline = f.readline()
	bhbhprop_esc = json.loads(escline)
	f.close()
	n1, n2 = 0, 0
	m1, m2 = [], []
	mp1, mp2 = [], []
	ms1, ms2 = [], []
	try:
		for i in bhbhprop_esc.keys():
			if bhbhprop_esc[i]['loc'] == runfilename:
				for j in range(len(bhbhprop_esc[i]['tmerge'])):
					if z[0]<=bhbhprop_esc[i]['tmerge'][j]<=12000.:
						n1 += 1
						tempm = bhbhprop_esc[i]['m1'][j] + bhbhprop_esc[i]['m2'][j]
						m1.append(tempm)
						mp1.append( bhbhprop_esc[i]['m1'][j] )
						ms1.append( bhbhprop_esc[i]['m1'][j] )
					if z[1]<=bhbhprop_esc[i]['tmerge'][j]<=12000.:
						n2 += 1
						tempm = bhbhprop_esc[i]['m1'][j] + bhbhprop_esc[i]['m2'][j]
						m2.append(tempm)
						mp2.append( bhbhprop_esc[i]['m1'][j] )
						ms2.append( bhbhprop_esc[i]['m1'][j] )
				raise StopIteration()
	except StopIteration:
		#print 'found esc props:', runfilename
		pass

	
	return n1, m1, n2, m2, mp1, mp2, ms1, ms2	
				


def extract_property_from_json_bound(jsonfilename, runfilename, z = [4058.3457426100003, 1477.3043937000007,]):
	f = open(jsonfilename, 'r')
	boundline = f.readline()
	bhbhprop_bound = json.loads(boundline)
	f.close()
	n1, n2 = 0, 0
	m1, m2 = [], []
	mp1, mp2 = [], []
	ms1, ms2 = [], []
	try:
		for i in bhbhprop_bound.keys():
			if bhbhprop_bound[i]['loc'] == runfilename:
				for j in range(len(bhbhprop_bound[i]['t'])):
					if z[0]<=bhbhprop_bound[i]['t'][j]<=12000.:
						n1 += 1
						tempm = bhbhprop_bound[i]['m1'][j] + bhbhprop_bound[i]['m2'][j]
						m1.append(tempm)
						mp1.append( bhbhprop_bound[i]['m1'][j] )
						ms1.append( bhbhprop_bound[i]['m1'][j] )
					if z[1]<=bhbhprop_bound[i]['t'][j]<=12000.:
						n2 += 1
						tempm = bhbhprop_bound[i]['m1'][j] + bhbhprop_bound[i]['m2'][j]
						m2.append(tempm)
						mp2.append( bhbhprop_bound[i]['m1'][j] )
						ms2.append( bhbhprop_bound[i]['m1'][j] )
				raise StopIteration()
	except StopIteration:
		#print 'found bound props:', runfilename
		pass

	
	return n1, m1, n2, m2, mp1, mp2, ms1, ms2	



def combine_and_find_stat(x, y):
	comb = []
	for i in range(len(x)):
		comb.append(x[i])
	for i in range(len(y)):
		comb.append(y[i])
	comb = np.array(comb)
	comb_med, comb_2siglow, comb_2sighigh = np.percentile(comb, 50.), np.percentile(comb,2.23), np.percentile(comb, 97.7)
	return comb_med, comb_2siglow, comb_2sighigh


#def make_binary_bh_table_data(list_of_runsfilename, propertyfilename, jsonfileesc, jsonfilebound, writefilename):
def collect_bound_bh_numbers(runlistname):
	filelocations, nbh, nbhbh, nbhnbh = [], [], [], []
	f = open(runlistname, 'r')
	f.seek(0)
	for line in f:
		if line.rfind('#')==-1:
			if line.rfind('CMC_results')>-1:
				a = line.split()[0]
				tmp = a.split('CMC_results')[1]
				tmp1 = '/Volumes/SC_pspt2'
				tmp2 = tmp1+tmp
				filelocations.append(tmp2)
			else:
				filelocations.append(line.split()[0])
	f.close()
	for i in range(len(filelocations)):
		#print filelocations[i]
		#tmp = filelocations[i].split('CMC_results')[1]
		#tmp1 = '/Volumes/SC_pspt2'
		#filename = tmp1+tmp+'/initial.esc.bh.dat'
		filename = filelocations[i]+'/initial.bh.dat'
		lastarr, lastline = read_lastline(filename)
		#data = np.loadtxt(filename)
		nbh.append(int(lastarr[2]))
		nbhbh.append(int(lastarr[5]))
		nbhnbh.append(int(lastarr[6]))

	return filelocations, nbh, nbhbh, nbhnbh


def collect_esc_bh_numbers(runlistname):
	filelocations, nbhesc, nbhbhesc, nbhnbhesc = [], [], [], []
	f = open(runlistname, 'r')
	f.seek(0)
	for line in f:
		if line.rfind('#')==-1:
			if line.rfind('CMC_results')>-1:
				a = line.split()[0]
				tmp = a.split('CMC_results')[1]
				tmp1 = '/Volumes/SC_pspt2'
				tmp2 = tmp1+tmp
				filelocations.append(tmp2)
			else:
				filelocations.append(line.split()[0])
	f.close()
	
	for i in range(len(filelocations)):
		#print filelocations[i]
		#tmp = filelocations[i].split('CMC_results')[1]
		#tmp1 = '/Volumes/SC_pspt2'
		#filename = tmp1+tmp+'/initial.esc.bh.dat'
		filename = filelocations[i]+'/initial.esc.bh.dat'
		lastarr, lastline = read_lastline(filename)
		#data = np.loadtxt(filename)
		nbhesc.append(int(lastarr[2]))
		nbhbhesc.append(int(lastarr[5]))
		nbhnbhesc.append(int(lastarr[6]))

	return filelocations, nbhesc, nbhbhesc, nbhnbhesc
#

def collect_bh_merger_prop(runlistname, jsonfileesc, jsonfilebound):
	filelocations = []
	f = open(runlistname, 'r')
	f.seek(0)
	for line in f:
		if line.rfind('#')==-1:
			a = line.split()[0]
			tmp = a.split('CMC_results')[1]
			tmp1 = '/Volumes/SC_pspt2'
			tmp2 = tmp1+tmp
			#filelocations.append(a[13])
			filelocations.append(tmp2)
	f.close()
	n1escarr, n2escarr = [], []
	n1boundarr, n2boundarr = [], []
	m1_medarr, m2_medarr = [], []
	m1_2siglowarr, m2_2siglowarr = [], []
	m1_2sighigharr, m2_2sighigharr = [], []
	for i in range(len(filelocations)):
		n1esc, m1esc, n2esc, m2esc, mp1esc, mp2esc, ms1esc, ms2esc = extract_property_from_json_escaped(jsonfileesc, filelocations[i], z = [4058.3457426100003, 1477.3043937000007,])
		n1bound, m1bound, n2bound, m2bound, mp1bound, mp2bound, ms1bound, ms2bound = extract_property_from_json_bound(jsonfilebound, filelocations[i], z = [4058.3457426100003, 1477.3043937000007,])
		
		if n1esc+n1bound>0:	
			m1_med, m1_2siglow, m1_2sighigh = combine_and_find_stat(m1esc, m1bound)
		else:	
			m1_med, m1_2siglow, m1_2sighigh = 0, 0, 0
			print 'looking at z=1:', filelocations[i], n1esc, n1bound
		if n2bound+n2esc>0:
			m2_med, m2_2siglow, m2_2sighigh = combine_and_find_stat(m2esc, m2bound)
		else:
			m2_med, m2_2siglow, m2_2sighigh = 0, 0, 0
			print 'looking at z=2:', filelocations[i], n2esc, n2bound

		n1escarr.append(n1esc)
		n2escarr.append(n2esc)
		n1boundarr.append(n1bound)
		n2boundarr.append(n2bound)
		m1_medarr.append(m1_med)
		m2_medarr.append(m2_med)
		m1_2siglowarr.append(m1_2siglow)
		m2_2siglowarr.append(m2_2siglow)
		m1_2sighigharr.append(m1_2sighigh)
		m2_2sighigharr.append(m2_2sighigh)

	return filelocations, n1escarr, n2escarr, n1boundarr, n2boundarr, m1_medarr, m2_medarr, m1_2siglowarr, m2_2siglowarr, m1_2sighigharr, m2_2sighigharr


def read_json_file(jsonfilename):
	f = open(jsonfilename, 'r')
	line = f.readline()
	bhbhprop = json.loads(line)
	f.close()
	return bhbhprop



def collect_stat_percentiles(t, m, tmerge=[0., 0.1, 0.5, 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13.,], PERCENTILE_LIST=[2.3, 97.7]): 
	medarr, lowarr, higharr = [], [], []
	tarr = []
	for i in range(len(tmerge)-1):
	#for i in range(1):
		tmpm, tmpt = [], []
		for j in range(len(t)):
			if tmerge[i]<=t[j]<=tmerge[i+1]:
				print tmerge[i], tmerge[i+1], t[j]
				tmpm.append(m)
		if len(tmpm)>0:
			#print tmerge[i], tmerge[i+1], len(tmpm)
			tmpmed = np.median(tmpm)
			tmplow = np.percentile(tmpm, PERCENTILE_LIST[0])
			tmphigh = np.percentile(tmpm, PERCENTILE_LIST[1])
			medarr.append(tmpmed)
			lowarr.append(tmplow)
			higharr.append(tmphigh)
			tarr.append( (tmerge[i]+tmerge[i+1])/2. )
	tarr = np.array(tarr)
	medarr = np.array(medarr)
	higharr = np.array(higharr)
	lowarr = np.array(lowarr)
	
	return tarr, medarr, higharr, lowarr 
		
		

def tmerger_histogram_one_run(jsonescfilename, jsonboundfilename, STRING='new_plots/tmerge_hist'):
	bhbhpropesc = read_json_file(jsonescfilename)
	bhbhpropbound = read_json_file(jsonboundfilename)
	
	import matplotlib.pyplot as plt

	NBIN = 20
	for i in bhbhpropesc.keys():
		print bhbhpropesc[i]['loc']
		tm_e, tm_b, tm_all = [], [], []
		for j in range(len(bhbhpropesc[i]['tmerge'])):
			tm_all.append(bhbhpropesc[i]['tmerge'][j]/1e3)
			tm_e.append(bhbhpropesc[i]['tmerge'][j]/1e3)
		for j in range(len(bhbhpropbound[i]['t'])):
			tm_all.append(bhbhpropbound[i]['t'][j]/1e3)
			tm_b.append(bhbhpropbound[i]['t'][j]/1e3)

		if len(tm_all)>0:

			plotfilename = STRING+'_'+str(i)+'.pdf'
			fig = plt.figure()
			ax = fig.add_subplot(111)
			
			tm_e = np.array(tm_e)
			tm_all = np.array(tm_all)
			tm_b = np.array(tm_b)
	
			ltm_e = np.log10(tm_e)
			ltm_all = np.log10(tm_all)
			ltm_b = np.log10(tm_b)
	
			ax.hist(ltm_all, bins=NBIN, histtype='step', range=[-3,np.log10(20)], color='black', ls='solid', label='all')
			if len(ltm_e)>0:
				ax.hist(ltm_e, bins=NBIN, histtype='step', range=[-3,np.log10(20)], color='blue', ls='solid', label='out')
			if len(ltm_b)>0:
				ax.hist(ltm_b, bins=NBIN, histtype='step', range=[-3,np.log10(20)], color='red', ls='solid', label='in')
		
			ax.legend(loc='best', numpoints=1, frameon=0)
			ax.set_title('run %s' %(str(i)))
			ax.set_xlabel(r'$\log(t/\rm{Gyr})$', size=14)
			ax.set_ylabel(r'$dN_{\rm{merger}}/d\log(t/\rm{Gyr})$', size=14)

			plt.savefig(plotfilename)
		
		
		
		

		
			
