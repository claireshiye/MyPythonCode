import numpy as np
import matplotlib
#matplotlib.use('PDF')
import matplotlib.pyplot as plt
from glob import glob
import collections
from collections import Counter
import os, sys
import re, gzip
import scipy.stats as ss
import StringIO
import scripts
import ns
import history_cmc_modified_bss as hbss
import dynamics as dyn

savepath='/projects/b1011/syr904/projects/BSS/'


path=np.genfromtxt('/projects/b1011/syr904/projects/BSS/rvgrid_path.dat', dtype='str')
#path=np.genfromtxt('/projects/b1011/syr904/projects/PULSAR/kickgrid_runs/paper/data/path_newmodel.dat', dtype='str')


##Find parameters needed for finding BSS
def find_para(filepath):
    #filepath=str(data[k])
    snaps=np.sort(glob(filepath+'/'+'*.snap*.dat.gz'))
    firstsnap=snaps[0]
    lastsnap=snaps[-1]
    
    #Find prefix
    x=firstsnap.replace('.snap0000.dat.gz','')
    prefix=x.replace(filepath,'')
    
    #Find time
    t_conv=ns.conv('t',filepath+prefix+'.conv.sh')
    time=ns.get_time(lastsnap)*t_conv    ## in Myr
    
    #Find mass of MS
    #datalast=np.genfromtxt(lastsnap)
    #binflag=datalast[:,7]; m=datalast[:,1]; m0=datalast[:,8]; m1=datalast[:,9]
    #kstar=datalast[:,14]; k0=datalast[:,17]; k1=datalast[:,18]
    #mms=[]
    #for j in range(len(binflag)):
    #    if binflag[j]==0:
    #        if kstar[j]==0 or kstar[j]==1: mms.append(m[j])
    #    if binflag[j]==1:
    #        if k0[j]==0 or k0[j]==1: mms.append(m0[j])
    #        if k1[j]==0 or k1[j]==1: mms.append(m1[j])
    
    return time, prefix, lastsnap  #, mms


##Find BSS
def find_BSS(lastsnap, mto):
    bif=[]; m0bss=[]; m1bss=[]; k0bss=[]; k1bss=[]; rbss=[]; id0bss=[]; id1bss=[]; abss=[]; ebss=[]
    
    ##Memory free version
    with gzip.open(lastsnap, 'r') as f:
        for _ in xrange(2):
            next(f)
        for line in f:
            #print line
            datalast=line.split()
            if int(datalast[7])!=1:
                if int(datalast[14])==0 or int(datalast[14])==1:
                    if float(datalast[1])>=1.05*mto:
                        bif.append(0); m0bss.append(float(datalast[1])); m1bss.append(-100)
                        k0bss.append(int(datalast[14])); k1bss.append(-100); rbss.append(float(datalast[2]))
                        id0bss.append(int(datalast[0])); id1bss.append(-100); abss.append(-100); ebss.append(-100)
            if int(datalast[7])==1:
                if int(datalast[17])==0 or int(datalast[17])==1:
                    if float(datalast[8])>=1.05*mto:
                        bif.append(1); m0bss.append(float(datalast[8])); m1bss.append(float(datalast[9]))
                        k0bss.append(int(datalast[17])); k1bss.append(int(datalast[18])); rbss.append(float(datalast[2]))
                        id0bss.append(int(datalast[10])); id1bss.append(int(datalast[11])); abss.append(float(datalast[12])); ebss.append(float(datalast[13]))
                        
                if int(datalast[18])==0 or int(datalast[18])==1:
                    if float(datalast[9])>=1.05*mto:
                        bif.append(1); m0bss.append(float(datalast[9])); m1bss.append(float(datalast[8]))
                        k0bss.append(int(datalast[18])); k1bss.append(int(datalast[17])); rbss.append(float(datalast[2]))
                        id0bss.append(int(datalast[11])); id1bss.append(int(datalast[10])); abss.append(float(datalast[12])); ebss.append(float(datalast[13]))
                            
        
    return bif, m0bss, m1bss, k0bss, k1bss, rbss, id0bss, id1bss, abss, ebss



##Find BH in the last timestep
def find_NBH_NTOT(filestring):     
    #datalast=np.genfromtxt()
    #binflag=datalast[:,7]; kstar=datalast[:,14]; k0=datalast[:,17]; k1=datalast[:,18]
    #nbh=0
    #for j in range(len(binflag)):
    #    if binflag[j]==0:
    #        if kstar[j]==14: nbh+=1
    #    if binflag[j]==1:
    #        if k0[j]==14: nbh+=1
    #        if k1[j]==14: nbh+=1
    filebh=filestring+'.bh.dat'
    filedyn=filestring+'.dyn.dat'
    with open(filebh, 'r') as fbh:
	for line in fbh:pass
	lastbh=line
    databh=lastbh.split()
    nbh=float(databh[2])

    with open(filedyn, 'r') as fdyn:
	for line in fdyn:pass
	lastdyn=line
    datadyn=lastdyn.split()
    ntot=float(datadyn[3])
	   
    return nbh, ntot


##Find Nbss
def find_NBSS(filestring):
    bssfile=filestring+'.BSS.dat'
    classfile=filestring+'.BSSclass.dat'
    n_sin=0; n_bin=0; n_coll=0; n_mtb=0; n_se=0
    n_bin_si=0; n_bin_bi=0; n_coll_si=0; n_coll_bi=0; n_se_si=0; n_se_bi=0
    with open(bssfile, 'r') as fbss:
        next(fbss)
        for line in fbss:
            databss=line.split()
            if int(databss[0])==0: n_sin+=1.; binflag=0
            if int(databss[0])==1: n_bin+=1.; binflag=1
    
    with open(classfile, 'r') as fclass:
        next(fclass)
        for line in fclass:
            dataclass=line.split()
            if int(dataclass[1])==1: n_coll+=1.
            if int(dataclass[9])==1: n_mtb+=1.
            if int(dataclass[5])==1: n_se+=1.


    return n_sin, n_bin, n_coll, n_mtb, n_se



def find_z(filepath):
    meta=-100
    #filepath=str(data[i])
    for fname in os.listdir(filepath):
        if fname.endswith('.cmc'):
            cmcfile=glob(filepath+'/'+'*.cmc')
            #print cmcfile
            thecmcfi=cmcfile[0]
            with open(thecmcfi) as fi:
                for line in fi:
                    if re.findall('OVERWRITE_Z', line):
                        l=re.findall('\d+\.\d+', line)
                        meta=float(l[0])
         
    if meta==-100:
        for fname in os.listdir(filepath):
            if fname.endswith('.sh'):
                shfile=glob(filepath+'*.sh')
                for i in range(len(shfile)):
                    s=shfile[i]
                    s=s.replace(filepath, '')
                    if s[:2]=='ge': 
                        theshfi=filepath+s; print theshfi
                        with open(theshfi) as fish:
                            for line in fish:
                                if re.findall('-Z', line):
                                    l=re.findall('-Z ([\d.]+)', line)
                                    meta=float(l[0])
                                    break                

    return meta



##Printout Nbss-Nbh-Ntot of All Models
def printout_Nbss_Nbh():
    #handle=StringIO.StringIO()
    #sys.stdout=handle
    fhandle=open('/Users/shiye/Documents/ClusterGroup/BSSproject/Num1.dat', 'a', 0)
    for k in range(577, len(data)):
        filepath=str(data[k])
        t, pref, ls=find_MS(filepath)
        mtoguess=scripts.find_MS_turnoff(t)
        z=dataz[k]
        mtotrue=scripts.find_MS_TO(t, z, mtoguess)
            
        Nbss=int(find_BSS(ls, mtotrue))
            
        if os.path.isfile(filepath+pref+'.bh.dat'):
            with open(filepath+pref+'.bh.dat') as fi:
                for line in fi: pass
                databh=line.split()
            Nbh=int(databh[2])
        else:
            Nbh=int(find_BH(ls))
            
        with open(filepath+pref+'.dyn.dat') as fo:
            for line in fo: pass
            datatot=line.split()
        Ntot=int(datatot[3])
            
            
        #fhandle.write(handle.getvalue())
        #print Nbss, Nbh, Ntot
        fhandle.write('%d %d %d\n'%(Nbss, Nbh, Ntot))
        #sys.stdout.close()



##Find hrdiag_L_T
def print_hrdiag_LT(sourcedir):
    pref='initial'
    filepath=sourcedir
    filestr=filepath+'/'+pref
    snaps=np.sort(glob(filestr+'.snap*.dat.gz'))
    lastno=len(snaps)-1
    snapno=str(lastno).zfill(4)
    print snapno
    scripts.hrdiag_L_T(filestr, snapno)



##Print out BSS of all models
def printout_BSS(start, end):
    for k in range(start, end):
        #filepath=path[k]
        filepath='/projects/b1011/syr904/cmc/cmc-mpi-09/rundir/kickgrid/kickgrid_1.0'
        t, pref, ls=find_para(filepath)
        l_conv=ns.conv('l',filepath+'/'+pref+'.conv.sh')
        mtoguess=scripts.find_MS_turnoff(t)
        #z=dataz[k]
        z=0.001
        mtotrue=scripts.find_MS_TO(t, z, mtoguess)
        strnum=str(k).zfill(4)
        bf, m0_bss, m1_bss, k0_bss, k1_bss, r_bss, id0_bss, id1_bss, a_bss, e_bss=find_BSS(ls, mtotrue)
        r_bsspc = [x * l_conv for x in r_bss]
        np.savetxt(filepath+'/'+pref+'.BSS.dat', np.c_[bf, id0_bss, id1_bss, m0_bss, m1_bss, k0_bss, k1_bss, r_bsspc, a_bss, e_bss], fmt ='%d %d %d %f %f %d %d %f %f %f', delimiter= ' ', header = '1.binflag, 2.id0, 3.id1, 4.m0[msun], 5.m1[msun], 6.k0, 7.k1, 8.r[pc], 9.a[AU], 10.e', comments = '#')
        print k
        #print strnum



##Classify BSS
def class_bss(start, end):
    bssfile=np.sort(glob('/Users/shiye/Documents/ClusterGroup/BSSproject/BSS_mcut1.05_lastsnap/BSS*.dat'))
    #COLL_SS=[]; COLL_BS=[]; COLL_BB=[]; SE=[]; SE_MERGER=[]; SE_DISRUPT=[]; MTB=[]
    fhandle=open('/Users/shiye/Documents/ClusterGroup/BSSproject/'+'BSSclass_500more.dat', 'a', 0)
    #fhandle.write('#1.coll, 2.coll_ss, 3.coll_bs, 4.coll_bb, 5.se, 6.se_merger, 7.se_disrupt, 8.se_binint, 9.mtb, 10.mtb_pure, 11.mtb_binint\n')
    for k in range(start, end):
        COLL=0; COLL_SS=0; COLL_BS=0; COLL_BB=0
        SE=0; SE_MERGER=0; SE_DISRUPT=0; SE_BININT=0
        MTB=0; MTB_PURE=0; MTB_BININT=0
        
        filepath=path[k]
        pref=prefixstring[k]
        filestr=filepath+pref
        mcut=1.05*mturnoff[k]; tnow=tlastcode[k]; zmodel=z[k]
        
        binintstring=filestr+'.binint.log'
        binint=glob(binintstring)
        binary=1
        if len(binint)==0: binary=0
            
        with open(bssfile[k], 'r') as fi:
            for _ in xrange(2):
                next(fi)
            for line in fi:
                databss=line.split()
                theid=[long(databss[1])]
                hdict=hbss.history_maker(theid, [1], pref, filepath, binary)
                bss_coll, bss_ss_coll, bss_bs_coll, bss_bb_coll, bss_se, bss_se_merger, bss_se_disruption, bss_had_binint, bss_mtb, bss_mtb_pure, bss_se_binint, bss_mtb_binint, actual_t, actual_position, star_id=hbss.classifying_BSS(hdict, long(databss[1]), binary, mcut, tnow, filestr, zmodel)
                #print bss_coll, bss_ss_coll, bss_bs_coll, bss_bb_coll, bss_se, bss_se_merger, bss_se_disruption, bss_had_binint, bss_mtb, bss_mtb_pure, bss_se_binint, bss_mtb_binint, actual_t, actual_position, star_id
                
                COLL+=bss_coll; COLL_SS+=bss_ss_coll; COLL_BS+=bss_bs_coll; COLL_BB+=bss_bb_coll
                SE+=bss_se; SE_MERGER+=bss_se_merger; SE_DISRUPT+=bss_se_disruption; SE_BININT+=bss_se_binint
                MTB+=bss_mtb; MTB_PURE+=bss_mtb_pure; MTB_BININT+=bss_mtb_binint          

        fhandle.write('%d %d %d %d %d %d %d %d %d %d %d\n'%(COLL, COLL_SS, COLL_BS, COLL_BB, SE, SE_MERGER, SE_DISRUPT, SE_BININT, MTB, MTB_PURE, MTB_BININT))
        
        print k
        


##Extract semimajor axis and eccentricity from history dictionary
def find_binint_ae(hdict, theid, comid, stringnum):
    ain=0; ein=0; aout=0; eout=0
    for j in hdict[theid]['binint']['binint'].keys():
        for i in hdict[theid]['binint']['binint'][j]['interaction']['input'].keys():
            idlenin=len(hdict[theid]['binint']['binint'][j]['interaction']['input'][i]['ids'])
            if idlenin==2:
                #print str(hdict[theid]['binint']['binint'][j]['interaction']['input'][i]['ids'][0]), str(hdict[theid]['binint']['binint'][j]['interaction']['input'][i]['ids'][1])
                if str(hdict[theid]['binint']['binint'][j]['interaction']['input'][i]['ids'][0]).rfind(':')<=-1 and str(hdict[theid]['binint']['binint'][j]['interaction']['input'][i]['ids'][1]).rfind(':')<=-1:
                    if long(hdict[theid]['binint']['binint'][j]['interaction']['input'][i]['ids'][0])==theid and long(hdict[theid]['binint']['binint'][j]['interaction']['input'][i]['ids'][1])==comid:
                        ain=float(hdict[theid]['binint']['binint'][j]['interaction']['input'][i]['a'])
                        ein=float(hdict[theid]['binint']['binint'][j]['interaction']['input'][i]['e'])
                        minp=float(hdict[theid]['binint']['binint'][j]['interaction']['input'][i]['m'][0])
                        #print ain, ein, minp
                
                    if long(hdict[theid]['binint']['binint'][j]['interaction']['input'][i]['ids'][1])==theid and long(hdict[theid]['binint']['binint'][j]['interaction']['input'][i]['ids'][0])==comid:
                        ain=float(hdict[theid]['binint']['binint'][j]['interaction']['input'][i]['a'])
                        ein=float(hdict[theid]['binint']['binint'][j]['interaction']['input'][i]['e'])
                        minp=float(hdict[theid]['binint']['binint'][j]['interaction']['input'][i]['m'][1])
                        #print ain, ein, minp


        if ain!=0:
            for o in hdict[theid]['binint']['binint'][j]['interaction']['output'].keys():
                idlenout=len(hdict[theid]['binint']['binint'][j]['interaction']['output'][o]['ids'])
                if idlenout==2:
                    if str(hdict[theid]['binint']['binint'][j]['interaction']['output'][o]['ids'][0]).rfind(':')<=-1 and str(hdict[theid]['binint']['binint'][j]['interaction']['output'][o]['ids'][1]).rfind(':')<=-1:
                        if long(hdict[theid]['binint']['binint'][j]['interaction']['output'][o]['ids'][0])==theid and long(hdict[theid]['binint']['binint'][j]['interaction']['output'][o]['ids'][1])==comid:
                            aout=float(hdict[theid]['binint']['binint'][j]['interaction']['output'][o]['a'][0])
                            eout=float(hdict[theid]['binint']['binint'][j]['interaction']['output'][o]['e'][0])
                            mout=float(hdict[theid]['binint']['binint'][j]['interaction']['output'][o]['m'][0])
                            #print aout, eout, mout

                        if long(hdict[theid]['binint']['binint'][j]['interaction']['output'][o]['ids'][1])==theid and long(hdict[theid]['binint']['binint'][j]['interaction']['output'][o]['ids'][0])==comid:
                            aout=float(hdict[theid]['binint']['binint'][j]['interaction']['output'][o]['a'][0])
                            eout=float(hdict[theid]['binint']['binint'][j]['interaction']['output'][o]['e'][0])
                            mout=float(hdict[theid]['binint']['binint'][j]['interaction']['output'][o]['m'][1])
                            #print aout, eout, mout

                if idlenout==3:
                    if str(hdict[theid]['binint']['binint'][j]['interaction']['output'][o]['ids'][0]).rfind(':')<=-1 and str(hdict[theid]['binint']['binint'][j]['interaction']['output'][o]['ids'][1]).rfind(':')<=-1:
                        if long(hdict[theid]['binint']['binint'][j]['interaction']['output'][o]['ids'][0])==theid and long(hdict[theid]['binint']['binint'][j]['interaction']['output'][o]['ids'][1])==comid:
                            aout=float(hdict[theid]['binint']['binint'][j]['interaction']['output'][o]['a'][0])
                            eout=float(hdict[theid]['binint']['binint'][j]['interaction']['output'][o]['e'][0])
                            mout=float(hdict[theid]['binint']['binint'][j]['interaction']['output'][o]['m'][0])
                            #print aout, eout, mout

                        if long(hdict[theid]['binint']['binint'][j]['interaction']['output'][o]['ids'][1])==theid and long(hdict[theid]['binint']['binint'][j]['interaction']['output'][o]['ids'][0])==comid:
                            aout=float(hdict[theid]['binint']['binint'][j]['interaction']['output'][o]['a'][0])
                            eout=float(hdict[theid]['binint']['binint'][j]['interaction']['output'][o]['e'][0])
                            mout=float(hdict[theid]['binint']['binint'][j]['interaction']['output'][o]['m'][1])
                            #print aout, eout, mout

                if aout!=0:
                        typeint=hdict[theid]['binint']['binint'][j]['interaction']['type']['type']
                        timeint=float(hdict[theid]['binint']['binint'][j]['interaction']['type']['time'])



        if ain!=0 and aout!=0:
            fbinint=open('/Users/shiye/Documents/ClusterGroup/BSSproject/BSS_class/'+'BSSint'+stringnum+'.dat', 'a+',0)
            fbinint.write('%d %f %s %f %f %f %f %f %f\n'%(theid, timeint, typeint, ain, ein, minp, aout, eout, mout))

                  



##Classify BSS
def printout_class_bss(start, end):
    #bssfile=np.sort(glob('/Users/shiye/Documents/ClusterGroup/BSSproject/BSS_mcut1.05_lastsnap/BSS*.dat'))
    for k in range(start, end):
        filepath=path[k]
        #filepath='/projects/b1011/syr904/cmc/cmc-mpi-09/rundir/kickgrid/kickgrid_1.0'
        #pref=prefixstring[k]
        t, pref, ls=find_para(filepath)
        filestr=filepath+'/'+pref

        mtoguess=scripts.find_MS_turnoff(t)
        #z=dataz[k]
        z=0.001
        mtotrue=scripts.find_MS_TO(t, z, mtoguess)
	
        snaps=np.sort(glob(filestr+'.snap*.dat.gz'))
        lastsnap=snaps[-1]
        tnow=ns.get_time(lastsnap)
        zmodel=z 
        mcut=1.05*mtotrue

        #mcut=1.05*mturnoff[k]; tnow=tlastcode[k]; zmodel=z[k]
        
        l_conv=ns.conv('l',filestr+'.conv.sh')
        
        strnum=str(k).zfill(4)
        #fhandle=open('/Users/shiye/Documents/ClusterGroup/BSSproject/BSS_class/'+'BSSclass'+strnum+'.dat', 'a+', 0)
        fhandle=open(filepath+'/'+'initial.BSSclass.dat', 'w+', 0)
        fhandle.write('#1.star_id, 2.bss_coll, 3.bss_ss_coll, 4.bss_bs_coll, 5.bss_bb_coll, 6.bss_se, 7.bss_se_merger, 8.bss_se_disruption, 9.bss_had_binint, 10.bss_mtb, 11.bss_mtb_pure, 12.bss_se_binint, 13.bss_mtb_binint, 14.actual_t, 15.actual_position, 16.primordial_binary\n')
        binintstring=filestr+'.binint.log'
        binint=glob(binintstring)
        binary=1
        if len(binint)==0: binary=0
        print binary

        bssfile=filestr+'.BSS.dat'            

        with open(bssfile, 'r') as fbss:
            next(fbss)
            for line in fbss:
                databss=line.split()
                theid=long(databss[1]); poscode=float(databss[7])/l_conv
                hd=hbss.history_maker([theid], [poscode], pref, filepath, binary)
                bss_coll, bss_ss_coll, bss_bs_coll, bss_bb_coll, bss_se, bss_se_merger, bss_se_disruption, bss_had_binint, bss_mtb, bss_mtb_pure, bss_se_binint, bss_mtb_binint, actual_t, actual_position, star_id=hbss.classifying_BSS(hd, long(databss[1]), binary, mcut, tnow, filestr, zmodel)
                #print bss_coll, bss_ss_coll, bss_bs_coll, bss_bb_coll, bss_se, bss_se_merger, bss_se_disruption, bss_had_binint, bss_mtb, bss_mtb_pure, bss_se_binint, bss_mtb_binint, actual_t, actual_position, star_id
                
                pribin=0
                if bss_mtb==1:
                    firstsnap=filestr+'.snap0000.dat.gz'
                    with gzip.open(firstsnap, 'r') as fsnap:
                        for _ in xrange(2):
                            next(fsnap)
                        for line in fsnap:
                            datasnap=line.split()
                            if int(datasnap[7])==1:
                                if long(datasnap[10])==theid: pribin=1; compid=long(datasnap[11])
                                if long(datasnap[11])==theid: pribin=1; compid=long(datasnap[10])
                
                #if pribin==1 and bss_had_binint==1:
                #    find_binint_ae(hd, theid, compid, strnum)
                                      
                fhandle.write('%d %d %d %d %d %d %d %d %d %d %d %d %d %f %f %d\n'%(star_id, bss_coll, bss_ss_coll, bss_bs_coll, bss_bb_coll, bss_se, bss_se_merger, bss_se_disruption, bss_had_binint, bss_mtb, bss_mtb_pure, bss_se_binint, bss_mtb_binint, actual_t, actual_position, pribin))
                
        
        print k



##print out Nbh, Ntot, Rc
def printout_numbers(start, end):
    pref='initial'
    NBH=[]; NTOT=[]; RC=[]; RHL=[]; T=[]
    for i in range(start, end):
	filepath=path[i]
	filestr=filepath+'/'+pref
	snapobs=np.sort(glob(filestr+'.snap*.obs_params.dat'))
	lastsnapobs=snapobs[-1]
	Nbh, Ntot=find_NBH_NTOT(filestr)
	Rc, Rhl, T_Gyr, Nbh=dyn.find_rcrh(lastsnapobs)
	NBH.append(Nbh); NTOT.append(Ntot); RC.append(Rc); RHL.append(Rhl); T.append(T_Gyr)
	
	print i

    np.savetxt('/projects/b1011/syr904/projects/BSS/kickgrid_property.dat', np.c_[T, NTOT, NBH, RC, RHL], fmt ='%f %d %d %f %f', delimiter= ' ', header = '1.t_Gyr, 2.Ntot, 3.Nbh, 4.rc[pc], 5.rhl[pc]', comments = '#')


##Extract number of BSS and number of different BSS classes for each model from BSS.dat and BSSclass.dat file
def get_BSSnum(pathlist, start, end):
    databss=np.genfromtxt(savepath+'bss_num.dat')
    Mcut=databss[:,-1]

    sourcedir=np.genfromtxt(pathlist, dtype=str)
    Model=[]; Nbss=[]; Nbh=[]; Ntot=[]; Mtot=[]; Nbss_coll=[]; Nbss_se=[]; Nbss_mtb=[]; Npribin=[]; Rc=[]; Rh=[]; Rhoc=[]; Nbssc=[]
    Mc=[]; Mc_wobh=[]

    for i in range(start, end):
        filestr=sourcedir[i]+'/initial'      

        nbh, ntot=dyn.find_NBH_NTOT_last(filestr)
        mtot, rc, rh, rho0=dyn.find_rcrh_mtotrho0(filestr)

        mcore, mcore_wobh=dyn.find_mc(filestr,rc)

        l_conv=ns.conv('l', filestr+'.conv.sh')
        m_conv=ns.conv('m', filestr+'.conv.sh')
        mtot=mtot*m_conv; rc=rc*l_conv; rh=rh*l_conv; rho0=rho0*m_conv/(l_conv**3)
        Mc.append(m_conv*mcore); Mc_wobh.append(m_conv*mcore_wobh)

        numbss=0; ncoll=0; nse=0; nmtb=0; npb=0
        numbsscore=0; ncollcore=0; nsecore=0; nmtbcore=0; 
        with open(filestr+'.BSS.dat', 'r') as fbss:
            next(fbss)
            for line in fbss: 
                numbss+=1
                databss=line.split()
                if float(databss[7])<=rc:
                    numbsscore+=1


        with open(filestr+'.BSSclass.dat','r') as fclass:
            next(fclass)
            for line in fclass:
                dataclass=line.split()
                if int(dataclass[1])==1: ncoll+=1
                if int(dataclass[5])==1: nse+=1
                if int(dataclass[9])==1: nmtb+=1
                if int(dataclass[15])==1: npb+=1


        l_conv=ns.conv('l', filestr+'.conv.sh')
        m_conv=ns.conv('m', filestr+'.conv.sh')
        mtot=mtot*m_conv; rc=rc*l_conv; rh=rh*l_conv; rho0=rho0*m_conv/(l_conv**3)

        Nbss.append(numbss); Nbh.append(nbh); Ntot.append(ntot); Nbss_coll.append(ncoll); Nbss_se.append(nse); Nbss_mtb.append(nmtb); Npribin.append(npb); Nbssc.append(numbsscore)
        Mtot.append(mtot); Rc.append(rc); Rh.append(rh); Rhoc.append(rho0)
        Model.append(i)

        print i


    np.savetxt(savepath+'bss_num_new.dat', np.c_[Model, Nbss, Nbssc, Nbh, Ntot, Mtot, Nbss_coll, Nbss_se, Nbss_mtb, Npribin, Rc, Rh, Rhoc, Mcut, Mc, Mc_wobh], fmt='%d %d %d %d %d %f %d %d %d %d %f %f %f %f %e %e', delimiter='', header='1.Model 2.Nbss 3.Nbssc 4.Nbh 5.Ntot 6.Mtot 7.Coll 8.SE 9.MTB 10.PriBin 11.rc_3D(pc) 12.rh_3D(pc) 13.rhoc(M_sun/pc^3 14.mcut(Msun) 15.Mc(Msun) 16.Mc_wobh(Msun))', comments='#')


def get_BSSnum_append(pathlist, start, end):
    Mc=[]; Mc_wobh=[]
    sourcedir=np.genfromtxt(pathlist, dtype=str)
    for i in range(start, end):
        filestr=sourcedir[i]+'/initial'
        mcore, mcore_wobh=dyn.find_mc(filestr)

        m_conv=dyn.conv('m', filestr+'.conv.sh')

        Mc.append(m_conv*mcore); Mc_wobh.append(m_conv*mcore_wobh)

        print i
    Mc=np.array([Mc]); Mc_wobh=np.array([Mc_wobh])

    databss=np.genfromtxt(savepath+'bss_num.dat')
    databss=np.array(databss)

    databss=np.concatenate((databss, Mc.T), axis=1)
    databss=np.concatenate((databss. Mc_wobh.T), axis=1)

    np.savetxt(savepath+'bss_num_new.dat', np.c_[databss[:,0], databss[:,1], databss[:,2], databss[:,3], databss[:,4], databss[:,5], databss[:,6], databss[:,7], databss[:,8], databss[:,9], databss[:,10], databss[:,11], databss[:,12], databss[:,13]], fmt='%d %d %d %d %f %d %d %d %d %f %f %f %f %e', delimiter='', header='1.Model 2.Nbss 3.Nbh 4.Ntot 5.Mtot 6.Coll 7.SE 8.MTB 9.PriBin 10.rc_3D(pc) 11.rh_3D(pc) 12.rhoc(M_sun/pc^3), 13.mcut(M_sun), 14.Mc(M_sun)', comments='#')



def get_MSWD(pathlist, start, end): ##Maybe they are failed BSSs because BSE doesn't treat common envelope correctly
     sourcedir=np.genfromtxt(pathlist, dtype=str)
     for i in range(len(sourcedir)):
        id0=[]; id1=[]; k0=[]; k1=[]; a=[]; ecc=[]; radrol0=[]; radrol1=[]; m0=[]; m1=[]
        filestr=sourcedir[i]+'/initial'
        snaps=np.sort(glob(filestr+'.snap*.dat.gz'))
        lastsnap=snaps[-1]

        with gzip.open(lastsnap, 'r') as flast:
            next(flast)
            next(flast)
            for line in flast:
                datalast=line.split()
                if (int(datalast[17])==0 or int(datalast[17])==1) and float(datalast[8])<0.88:
                    if int(datalast[18])>=10 and int(datalast[18])<=12:
                        id0.append(int(datalast[10])); id1.append(int(datalast[11])); k0.append(int(datalast[17])); k1.append(int(datalast[18])); m0.append(float(datalast[8])); m1.append(float(datalast[9])); a.append(float(datalast[12])); ecc.append(float(datalast[13])); radrol0.append(float(datalast[43])); radrol1.append(float(datalast[44]))

                if (int(datalast[18])==0 or int(datalast[18])==1) and float(datalast[9])<0.88:
                    if int(datalast[17])>=10 and int(datalast[17])<=12:
                        id0.append(int(datalast[11])); id1.append(int(datalast[10])); k0.append(int(datalast[18])); k1.append(int(datalast[17])); m0.append(float(datalast[9])); m1.append(float(datalast[8])); a.append(float(datalast[12])); ecc.append(float(datalast[13])); radrol0.append(float(datalast[44])); radrol1.append(float(datalast[43]))


        np.savetxt(savepath+'MSWDbinaries/MSWD'+str(i).zfill(2)+'.dat', np.c_[id0, id1, m0, m1, k0, k1, a, ecc, radrol0, radrol1], fmt='%d %d %f %f %d %d %f %f %f %f', delimiter='', header='1.id0 2.id1 3.m0(Msun) 4.m1(Msun) 5.k0 6.k1 7.a(AU) 8.ecc 9.radrol0 10.radrol1', comments='#')

        print i


def get_MSWD_progenitors(sourcedir, start, end):
    pathlist=np.genfromtxt(sourcedir+'/rvgrid_path.dat', dtype=str)
    mswdfiles=np.sort(glob(sourcedir+'/MSWDbinaries/*.dat'))
    for i in range(start, end):
        filestr=pathlist[i]+'/initial'
        datamswd=np.genfromtxt(mswdfiles[i])
        ID0=datamswd[:,0]; ID1=datamswd[:,1]; K1=datamswd[:,5]

        #id0=[]; id1=[]; k0=[]; k1=[]; a=[]; ecc=[]; radrol0=[]; radrol1=[]; m0=[]; m1=[]; time=[]
        fmswd=open(savepath+'MSWDbinaries/MSWD_projenitor_new'+str(i).zfill(2)+'.dat', 'a+', 0)
        fmswd.write('#1.Time 2.id0 3.id1 4.m0(Msun) 5.m1(Msun) 6.k0 7.k1 8.a(AU) 9.ecc 10.radrol0 11.radrol1 12.massc0(Msun) 13.massc1(Msun)\n')
        snaps=np.sort(glob(filestr+'.snap*.dat.gz'))
        t_conv=ns.conv('t',filestr+'.conv.sh')
        for j in range(len(ID0)):
            theid=int(ID0[j]); theidcom=int(ID1[j])
            print theid
            #checkpri=0
            #with gzip.open(snaps[1], 'r') as fsnap:
            #        for _ in xrange(2): next(fsnap)
            #        for line in fsnap:
            #            datasnap=line.split()
            #            if (int(datasnap[10])==theid and int(datasnap[11])==theidcom) or (int(#datasnap[11])==theid and int(datasnap[10])==theidcom):
            #                print 'primordial'
            #                id0=theid; id1=theidcom; k0=0; k1=0; m0=0; m1=0; a=0; ecc=0; radrol0=0; #radrol1=0
            #                time=ns.get_time(snaps[0])*t_conv
            #                checkpri=1
            #if checkpri==1: 
            #    fmswd.write('%f %d %d %f %f %d %d %f %f %f %f\n'%(time, id0, id1, m0, m1, k0, k1, a, #ecc, radrol0, radrol1))
            #    continue

            for k in range(len(snaps)-1, -1, -1):
                check=0
                t=ns.get_time(snaps[k])
                with gzip.open(snaps[k], 'r') as fsnap:
                    for _ in xrange(2): next(fsnap)
                    for line in fsnap:
                        datasnap=line.split()
                        if int(datasnap[10])==theid and int(datasnap[11])==theidcom and int(datasnap[18])<10:
                            print 'yes'
                            id0=int(datasnap[10]); id1=int(datasnap[11]); k0=int(datasnap[17]); k1=int(datasnap[18]); m0=float(datasnap[8]); m1=float(datasnap[9]); a=float(datasnap[12]); ecc=float(datasnap[13]); radrol0=float(datasnap[43]); radrol1=float(datasnap[44])
                            mc0=float(datasnap[31]); mc1=float(datasnap[32])
                            time=t*t_conv
                            check=1

                        if int(datasnap[11])==theid and int(datasnap[10])==theidcom and int(datasnap[17])<10:
                            print 'yes'
                            id0=int(datasnap[11]); id1=int(datasnap[10]); k0=int(datasnap[18]); k1=int(datasnap[17]); m0=float(datasnap[9]); m1=float(datasnap[8]); a=float(datasnap[12]); ecc=float(datasnap[13]); radrol0=float(datasnap[44]); radrol1=float(datasnap[43])
                            mc0=float(datasnap[32]); mc1=float(datasnap[31])
                            time=t*t_conv
                            check=1

                        if int(datasnap[10])==theid and int(datasnap[11])!=theidcom:
                            print 'no'
                            id0=theid; id1=int(datasnap[11]); k0=-100; k1=-100; m0=-100; m1=-100; a=-100; ecc=-100; radrol0=-100; radrol1=-100; mc0=-100; mc1=-100
                            time=t*t_conv
                            check=1

                        if int(datasnap[11])==theid and int(datasnap[10])!=theidcom:
                            print 'no'
                            id0=theid; id1=int(datasnap[10]); k0=-100; k1=-100; m0=-100; m1=-100; a=-100; ecc=-100; radrol0=-100; radrol1=-100; mc0=-100; mc1=-100
                            time=t*t_conv
                            check=1

                        if int(datasnap[0])==theid:
                            print 'no'
                            id0=theid; id1=-100; k0=-100; k1=-100; m0=-100; m1=-100; a=-100; ecc=-100; radrol0=-100; radrol1=-100; mc0=-100; mc1=-100
                            time=t*t_conv
                            check=1
                print k

                if check==1: 
                    fmswd.write('%f %d %d %f %f %d %d %f %f %f %f %f %f\n'%(time, id0, id1, m0, m1, k0, k1, a, ecc, radrol0, radrol1, mc0, mc1))
                    break

            print j


        #np.savetxt(savepath+'MSWDbinaries/MSWD_projenitor'+str(i).zfill(2)+'.dat', np.c_[time, id0, id1, m0, m1, k0, k1, a, ecc, radrol0, radrol1], fmt='%f %d %d %f %f %d %d %f %f %f %f', delimiter='', header='1.Time 2.id0 3.id1 4.m0(Msun) 5.m1(Msun) 6.k0 7.k1 8.a(AU) 9.ecc 10.radrol0 11.radrol1', comments='#')

        fmswd.close()

        print i


##Find out which MS-WD binaries went through CE evolution
def find_CE_MSWD(sourcedir, start, end):
    progenitorfiles=np.sort(glob(sourcedir+'/MSWDbinaries/MSWD_projenitor_new*.dat'))
    for i in range(len(progenitorfiles)):
        datapro=np.genfromtxt(progenitorfiles[i])
        time=datapro[:,0]; m0=datapro[:,3]; m1=datapro[:,4]; mc0=datapro[:,11]; mc1=datapro[:,12]
        k0=datapro[:,5]; k1=datapro[:,6]

        nce=0

        T=[]; ID0=[]; ID1=[]; M0=[]; M1=[]; K0=[]; K1=[]; A=[]; ECC=[]; ROL0=[]; ROL1=[]; MC0=[]; MC1=[]
        for j in range(len(time)):  ##q_crit from BSE
            if k1[j]==2:
                qc=4.0
                if m1[j]/m0[j]>qc: 
                    nce+=1
                    T.append(datapro[j,0]); ID0.append(datapro[j,1]); ID1.append(datapro[j,2]); M0.append(datapro[j,3]); M1.append(datapro[j,4]); K0.append(datapro[j,5]); K1.append(datapro[j,6]); A.append(datapro[j,7]); ECC.append(datapro[j,8]); ROL0.append(datapro[j,9]); ROL1.append(datapro[j,10]); MC0.append(datapro[j,11]); MC1.append(datapro[j,12])
            elif k1[j]==3 or k1[j]==5 or k1[j]==6:
                qc=0.362+(3.0*(1.0-mc1[j]/m1[j]))**-1
                if m1[j]/m0[j]>qc:
                    nce+=1
                    T.append(datapro[j,0]); ID0.append(datapro[j,1]); ID1.append(datapro[j,2]); M0.append(datapro[j,3]); M1.append(datapro[j,4]); K0.append(datapro[j,5]); K1.append(datapro[j,6]); A.append(datapro[j,7]); ECC.append(datapro[j,8]); ROL0.append(datapro[j,9]); ROL1.append(datapro[j,10]); MC0.append(datapro[j,11]); MC1.append(datapro[j,12])
            elif k1[j]==8 or k1[j]==9:
                qc=0.784
                if m1[j]/m0[j]>qc:
                    nce+=1
                    T.append(datapro[j,0]); ID0.append(datapro[j,1]); ID1.append(datapro[j,2]); M0.append(datapro[j,3]); M1.append(datapro[j,4]); K0.append(datapro[j,5]); K1.append(datapro[j,6]); A.append(datapro[j,7]); ECC.append(datapro[j,8]); ROL0.append(datapro[j,9]); ROL1.append(datapro[j,10]); MC0.append(datapro[j,11]); MC1.append(datapro[j,12])
            else:
                qc=3.0
                if m1[j]/m0[j]>qc:
                    nce+=1
                    T.append(datapro[j,0]); ID0.append(datapro[j,1]); ID1.append(datapro[j,2]); M0.append(datapro[j,3]); M1.append(datapro[j,4]); K0.append(datapro[j,5]); K1.append(datapro[j,6]); A.append(datapro[j,7]); ECC.append(datapro[j,8]); ROL0.append(datapro[j,9]); ROL1.append(datapro[j,10]); MC0.append(datapro[j,11]); MC1.append(datapro[j,12])


        print nce
        np.savetxt(sourcedir+'/MSWDbinaries/'+'MSWD_CE'+str(i).zfill(2)+'.dat', np.c_[T, ID0, ID1, M0, M1, K0, K1, A, ECC, ROL0, ROL1, MC0, MC1], fmt='%f %d %d %f %f %d %d %f %f %f %f %f %f', header='1.Time 2.id0 3.id1 4.m0(Msun) 5.m1(Msun) 6.k0 7.k1 8.a(AU) 9.ecc 10.radrol0 11.radrol1 12.massc0(Msun) 13.massc1(Msun)', delimiter=' ', comments='#')

        print i



            








##Find metallicity, turnoff mass and lastsnap time for Models
#z=[]; mto=[]; tlast=[]; tlastcode=[]; prefix=[]
#for i in range(len(data)):
#    filepath= data[i]
#    
#    ##Metallicity
#    #z.append(find_z(filepath))
#    
#    ##tlast, turnoff mass and prefix
#    t, pref, ls=find_MS(filepath)
#    #tcode=ns.get_time(ls)
#    #mtoguess=scripts.find_MS_turnoff(t)
#    #z=dataz[i]
#    #mtotrue=scripts.find_MS_TO(t, z, mtoguess)
#    #tlast.append(t); mto.append(mtotrue); tlastcode.append(tcode)
#    prefix.append(pref)
#    
#np.savetxt('/Users/shiye/Documents/ClusterGroup/BSSproject/prefix.dat', np.c_[prefix], fmt='%s')
#np.savetxt('/Users/shiye/Documents/ClusterGroup/BSSproject/Modelgt12_prop_appenx.dat', np.c_[tlastcode, tlast, mto], fmt='%f %f %f', header='1.tlastcode, 2.tlast, 3.mto', delimiter=' ', comments='#')
#np.savetxt('/Users/shiye/Documents/ClusterGroup/Modelgt12_prop.dat', np.c_[z], fmt='%f')




#printout_Nbss_Nbh()

#plot_Nbh_Nbss(0, 32)

#hrdiag_LT()

#plot_hrdiag('/Volumes/homes/sourav/CMC_results/BH_variations/kickoutputtest_variations/rundir/wind1/z0.001/normalkick/')

#printout_BSS(100, 150)

#printout_Nbss_sinbin()

#class_bss(579,597)

#plot__Nbss_sinbin()

#plot_Nbss_class()

#printout_BSS(16,32)
#printout_class_bss(16, 32)
