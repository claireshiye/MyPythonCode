
# coding: utf-8

# In[ ]:


#get_ipython().magic(u'matplotlib inline')

import numpy as np
import matplotlib
matplotlib.use('PDF')
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


# In[ ]:


data=np.genfromtxt('/Users/shiye/Documents/ClusterGroup/BSSproject/Modelgt12withpref.dat', dtype='str')
path=data[:,0]; prefixstring=data[:,1]
#dataz=np.genfromtxt('/Users/shiye/Documents/ClusterGroup/BSSproject/Modelgt12_prop.dat')
dataprop=np.genfromtxt('/Users/shiye/Documents/ClusterGroup/BSSproject/Modelgt12properties.dat')
tlastcode=dataprop[:,0]; mturnoff=dataprop[:,2]; z=dataprop[:,3]


# In[ ]:


##Find Main Sequence Star
def find_MS(filepath):
    #filepath=str(data[k])
    snaps=np.sort(glob(filepath+'*.snap*.dat.gz'))
    firstsnap=snaps[0]
    lastsnap=snaps[-1]
    
    #Find prefix
    x=firstsnap.replace('.snap0000.dat.gz','')
    prefix=x.replace(filepath,'')
    
    #Find time
    t_conv=ns.conv('t',filepath+prefix+'.conv.sh')
    time=ns.get_time(lastsnap)*t_conv
    
    
    
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


# In[ ]:


##Find BSS
def find_BSS(lastsnap, mto):
    #datalast=np.genfromtxt(lastsnap)
    #binflag=datalast[:,7]; m=datalast[:,1]; m0=datalast[:,8]; m1=datalast[:,9]
    #kstar=datalast[:,14]; k0=datalast[:,17]; k1=datalast[:,18]
    #sid=datalast[:,0]; id0=datalast[:,10]; id1=datalast[:,11]; rposition=datalast[:,2]
    ##nbss=0
    bif=[]; m0bss=[]; m1bss=[]; k0bss=[]; k1bss=[]; rbss=[]; id0bss=[]; id1bss=[]
    #for j in range(len(binflag)):
    #    if binflag[j]==0:
    #        if kstar[j]==0 or kstar[j]==1:
    #            if m[j]>=1.05*mto: 
    #                #nbss+=1
    #                bif.append(0); m0bss.append(m[j]); m1bss.append(-100)
    #                k0bss.append(kstar[j]); k1bss.append(-100); rbss.append(rposition[j])
    #                id0bss.append(sid[j]); id1bss.append(-100)
    #    if binflag[j]==1:
    #        if k0[j]==0 or k0[j]==1:
    #            if m0[j]>=1.05*mto: 
    #                #nbss+=1
    #                bif.append(1); m0bss.append(m0[j]); m1bss.append(m1[j])
    #                k0bss.append(k0[j]); k1bss.append(k1[j]); rbss.append(rposition[j])
    #                id0bss.append(id0[j]); id1bss.append(id1[j])
    #        if k1[j]==0 or k1[j]==1:
    #            if m1[j]>=1.05*mto: 
    #                #nbss+=1
    #                bif.append(1); m0bss.append(m1[j]); m1bss.append(m0[j])
    #                k0bss.append(k1[j]); k1bss.append(k0[j]); rbss.append(rposition[j])
    #                id0bss.append(id1[j]); id1bss.append(id0[j])
    #                
    #return nbss
    
    ##Memory free version
    with gzip.open(lastsnap, 'r') as f:
        for _ in xrange(2):
            next(f)
        for line in f:
            #print line
            datalast=line.split()
            if int(datalast[7])==0:
                if int(datalast[14])==0 or int(datalast[14])==1:
                    if float(datalast[1])>=1.05*mto:
                        bif.append(0); m0bss.append(float(datalast[1])); m1bss.append(-100)
                        k0bss.append(int(datalast[14])); k1bss.append(-100); rbss.append(float(datalast[2]))
                        id0bss.append(int(datalast[0])); id1bss.append(-100)
            if int(datalast[7])==1:
                if int(datalast[17])==0 or int(datalast[17])==1:
                    if float(datalast[8])>=1.05*mto:
                        bif.append(1); m0bss.append(float(datalast[8])); m1bss.append(float(datalast[9]))
                        k0bss.append(int(datalast[17])); k1bss.append(int(datalast[18])); rbss.append(float(datalast[2]))
                        id0bss.append(int(datalast[10])); id1bss.append(int(datalast[11]))
                        
                if int(datalast[18])==0 or int(datalast[18])==1:
                    if float(datalast[9])>=1.05*mto:
                        bif.append(1); m0bss.append(float(datalast[9])); m1bss.append(float(datalast[8]))
                        k0bss.append(int(datalast[18])); k1bss.append(int(datalast[17])); rbss.append(float(datalast[2]))
                        id0bss.append(int(datalast[11])); id1bss.append(int(datalast[10]))
                            
        
    return bif, m0bss, m1bss, k0bss, k1bss, rbss, id0bss, id1bss


# In[ ]:


##Find BH
def find_BH(lastsnap):
    datalast=np.genfromtxt(lastsnap)
    binflag=datalast[:,7]; kstar=datalast[:,14]; k0=datalast[:,17]; k1=datalast[:,18]
    nbh=0
    for j in range(len(binflag)):
        if binflag[j]==0:
            if kstar[j]==14: nbh+=1
        if binflag[j]==1:
            if k0[j]==14: nbh+=1
            if k1[j]==14: nbh+=1
        
    return nbh


# In[ ]:


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


# In[ ]:


##Plot Mass Distribution
def plot_massdist():
    for k in range(0, 600, 100):
        filepath=str(data[k])
        time, mms=find_MS(filepath)
        mtoguess=scripts.find_MS_turnoff(time)
        print mtoguess
        z=dataz[k]
        mtotrue=scripts.find_MS_TO(time, z, mtoguess)
        print mtotrue
        
        plt.figure()
        plt.yscale('log')
        plt.hist(mms, bins=50,color='orange')
        plt.axvline(x=mtotrue, color='b', linestyle='--')
        plt.xlabel(r'$mass(M_{\odot})$')
        #plt.title(filepath)
        plt.show()


# In[ ]:


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
        


# In[ ]:


##Plot Nbh-Nbss
def plot_Nbh_Nbss():
    datan=np.genfromtxt('/Users/shiye/Documents/ClusterGroup/BSSproject/Num.dat')
    Nbss=np.array(datan[:,0])
    NBSS=Nbss.astype(float)
    Nbh=np.array(datan[:,1])
    NBH=Nbh.astype(float)
    Ntot=np.array(datan[:,2])
    NTOT=Ntot.astype(float)
    #print np.log(NBSS/NTOT), np.log(NBH/NTOT)
    
    rho, p=ss.spearmanr(np.log(NBH/NTOT), np.log(NBSS/NTOT))
    print rho, p    

    plt.figure()
    plt.scatter(np.log(NBH/NTOT), np.log(NBSS/NTOT), marker='.')
    #plt.xlim(-10., -1.)
    #plt.ylim(-10., -1.)
    #plt.xscale('symlog')
    #plt.yscale('symlog')
    plt.xlabel('Nbh/Ntot')
    plt.ylabel('Nbss/Ntot')
    #plt.show()
    plt.savefig('/Users/shiye/Documents/ClusterGroup/BSSproject/nbhnbss_log.png')


# In[ ]:


##Printout Nbss-single-binary
def printout_Nbss_sinbin():
    NBSS_sin=[]; NBSS_bin=[]
    bssfile=np.sort(glob('/Users/shiye/Documents/ClusterGroup/BSSproject/BSS_mcut1.05_lastsnap/BSS*.dat'))
    for k in range(len(bssfile)):
        nbss_sin=0; nbss_bin=0
        with open(bssfile[k], 'r') as fi:
            for _ in xrange(2):
                next(fi)
            for line in fi:
                databss=line.split()
                if int(databss[0])==0:
                    nbss_sin+=1
                if int(databss[0])==1:
                    nbss_bin+=1
        
        NBSS_sin.append(nbss_sin)
        NBSS_bin.append(nbss_bin)
        
    np.savetxt('/Users/shiye/Documents/ClusterGroup/BSSproject/Num_sin_bin.dat', np.c_[NBSS_sin, NBSS_bin], fmt='%d %d', delimiter=' ', header='1.Nbss_single, 2. Nbss_binary', comments='#')


# In[ ]:


##Plot Nbh-Nbss-single-binary
def plot__Nbss_sinbin():
    datan=np.genfromtxt('/Users/shiye/Documents/ClusterGroup/BSSproject/bssnum.dat')
    NBSS=datan[:,0]; NBH=datan[:,1]; NTOT=datan[:,2]; NBSSSI=datan[:,3]; NBSSBI=datan[:,4]
    
    rho, p=ss.spearmanr(np.log(NBH/NTOT), np.log(NBSS/NTOT))
    rhos, ps=ss.spearmanr(np.log(NBH/NTOT), np.log(NBSSSI/NTOT))
    rhob, pb=ss.spearmanr(np.log(NBH/NTOT), np.log(NBSSBI/NTOT))
    
    print rho, p
    print rhos, ps
    print rhob, pb
    
    plt.figure(1)
    plt.scatter(np.log(NBH/NTOT), np.log(NBSSSI/NTOT), color='purple', label='single', s=8)
    plt.scatter(np.log(NBH/NTOT), np.log(NBSSBI/NTOT), color='orange', label='binary', s=5, alpha=0.7)
    plt.xlabel('log(Nbh/Ntot)')
    plt.ylabel('log(Nbss/Ntot)')
    plt.legend(loc='lower left')
    
    #plt.figure(2)
    #plt.scatter(np.log(NBH/NTOT), np.log(NBSSBI/NTOT), color='orange', label='binary', s=8)
    #plt.xlabel('Nbh/Ntot')
    #plt.ylabel('Nbss/Ntot')
    #plt.legend(loc='lower left')
    
    #plt.subplot(133)
    #plt.scatter(np.log(NBH/NTOT), np.log(NBSS/NTOT), marker='.')
    #plt.xlabel('Nbh/Ntot')
    #plt.ylabel('Nbss/Ntot')
    
    plt.savefig('/Users/shiye/Documents/ClusterGroup/BSSproject/nbhnbss_all.png', dpi = 300)
    #plt.show()
    


# In[ ]:


##Find hrdiag_L_T
def hrdiag_LT():
    filestring='/Volumes/homes/sourav/CMC_results/BH_variations/kickoutputtest_variations/rundir/wind1/z0.001/normalkick/initial'
    snapno='0096'
    scripts.hrdiag_L_T(filestring, snapno)


# In[ ]:


##Plot hrdiag
def plot_hrdiag(filepath):
    t, pref, ls=find_MS(filepath)
    mtoguess=scripts.find_MS_turnoff(t)
    z=dataz[48]
    mtotrue=scripts.find_MS_TO(t, z, mtoguess)
    
    datahrd=np.genfromtxt(filepath+pref+'.snap0096.hrdiag.dat')
    binflag=datahrd[:,0]; k0=datahrd[:,1]; k1=datahrd[:,2]; m0=datahrd[:,6]; m1=datahrd[:,7]; Teff=datahrd[:,15]; Leff=datahrd[:,16]
    si_bssL=[]; bi_bssL=[]; starL=[]; si_bssT=[]; bi_bssT=[]; starT=[]
    for k in range(len(binflag)):
        if binflag[k]==1:
            if k0[k]==0 or k0[k]==1:
                if m0[k]>=1.05*mtotrue:
                    bi_bssT.append(Teff[k]); bi_bssL.append(Leff[k])
                else:
                    starT.append(Teff[k]); starL.append(Leff[k])
                    
            elif k1[k]==0 or k1[k]==1:
                if m1[k]>=1.05*mtotrue:
                    bi_bssT.append(Teff[k]); bi_bssL.append(Leff[k])
                else:
                    starT.append(Teff[k]); starL.append(Leff[k])
            else:
                starT.append(Teff[k]); starL.append(Leff[k])
        else:
            if k0[k]==0 or k0[k]==1:
                if m0[k]>=1.05*mtotrue:
                    si_bssT.append(Teff[k]); si_bssL.append(Leff[k])
                else:
                    starT.append(Teff[k]); starL.append(Leff[k])
            else:
                starT.append(Teff[k]); starL.append(Leff[k])
    
    plt.figure()
    plt.xlabel(r'$log(T_{eff}/K)$')
    plt.ylabel(r'$log(L_{eff}/L_{\odot})$')
    plt.xlim(3.3, 4.2)
    plt.ylim(-3.5, 4)
    plt.scatter(starT, starL, marker='.', s=3, edgecolors='none', facecolor='k')
    plt.scatter(si_bssT, si_bssL, marker='o', s=5, facecolors='none', edgecolors='r', label='single')
    plt.scatter(bi_bssT, bi_bssL, marker='^', s=5, facecolors='none', edgecolors='b', label='binary')
    plt.gca().invert_xaxis()
    plt.legend()
    #plt.show()
    plt.savefig('/Users/shiye/Documents/ClusterGroup/BSSproject/hrd.png', dpi = 300)


# In[ ]:


##Print out BSS of all models
def printout_BSS(start, end):
    for k in range(start, end):
        filepath=str(data[k])
        t, pref, ls=find_MS(filepath)
        l_conv=ns.conv('l',filepath+pref+'.conv.sh')
        mtoguess=scripts.find_MS_turnoff(t)
        z=dataz[k]
        mtotrue=scripts.find_MS_TO(t, z, mtoguess)
        strnum=str(k).zfill(4)
        fhandle=open('/Users/shiye/Documents/ClusterGroup/BSSproject/BSS_mcut1.05_lastsnap/'+'BSS'+strnum+'.dat', 'a+', 0)
        bf, m0_bss, m1_bss, k0_bss, k1_bss, r_bss, id0_bss, id1_bss=find_BSS(ls, mtotrue)
        r_bsspc = [x * l_conv for x in r_bss]
        fhandle.write('#%s\n'%(filepath))
        np.savetxt(fhandle, np.c_[bf, id0_bss, id1_bss, m0_bss, m1_bss, k0_bss, k1_bss, r_bsspc], fmt ='%d %d %d %f %f %d %d %f', delimiter= ' ', header = '1.binflag, 2.id0, 3.id1, 4.m0[msun], 5.m1[msun], 6.k0, 7.k1, 8.r[pc]', comments = '#')
        #print k
        print strnum


# In[ ]:


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
        


# In[ ]:


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

                  




# In[ ]:


##Classify BSS
def printout_class_bss(start, end):
    bssfile=np.sort(glob('/Users/shiye/Documents/ClusterGroup/BSSproject/BSS_mcut1.05_lastsnap/BSS*.dat'))
    for k in range(start, end):
        filepath=path[k]
        pref=prefixstring[k]
        filestr=filepath+pref
        mcut=1.05*mturnoff[k]; tnow=tlastcode[k]; zmodel=z[k]
        
        l_conv=ns.conv('l',filestr+'.conv.sh')
        
        strnum=str(k).zfill(4)
        fhandle=open('/Users/shiye/Documents/ClusterGroup/BSSproject/BSS_class/'+'BSSclass'+strnum+'.dat', 'a+', 0)
        fhandle.write('#1.star_id, 2.bss_coll, 3.bss_ss_coll, 4.bss_bs_coll, 5.bss_bb_coll, 6.bss_se, 7.bss_se_merger, 8.bss_se_disruption, 9.bss_had_binint, 10.bss_mtb, 11.bss_mtb_pure, 12.bss_se_binint, 13.bss_mtb_binint, 14.actual_t, 15.actual_position, 16.primordial_binary\n')
        binintstring=filestr+'.binint.log'
        binint=glob(binintstring)
        binary=1
        if len(binint)==0: binary=0
            
        with open(bssfile[k], 'r') as fi:
            for _ in xrange(2):
                next(fi)
            for line in fi:
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
                
                if pribin==1 and bss_had_binint==1:
                    find_binint_ae(hd, theid, compid, strnum)
                                      
                fhandle.write('%d %d %d %d %d %d %d %d %d %d %d %d %d %f %f %d\n'%(star_id, bss_coll, bss_ss_coll, bss_bs_coll, bss_bb_coll, bss_se, bss_se_merger, bss_se_disruption, bss_had_binint, bss_mtb, bss_mtb_pure, bss_se_binint, bss_mtb_binint, actual_t, actual_position, pribin))
                
        
        print k


# In[ ]:


##Plot Nbh_Nbss class
def plot_Nbss_class():
    datan=np.genfromtxt('/Users/shiye/Documents/ClusterGroup/BSSproject/bssnum.dat')
    NBSS=datan[:,0]; NBH=datan[:,1]; NTOT=datan[:,2]; NBSSSI=datan[:,3]; NBSSBI=datan[:,4]
    dataclass=np.genfromtxt('/Users/shiye/Documents/ClusterGroup/BSSproject/BSSclass.dat')
    COLL=dataclass[:,0]; MTB=dataclass[:,8]
    nbh=[]; ntot=[]; coll=[]; mtb=[]
    for i in range(len(COLL)):
        if COLL[i]!=-100:
            nbh.append(NBH[i]); ntot.append(NTOT[i])
            coll.append(COLL[i]); mtb.append(MTB[i])
    
    
    nbh=np.array(nbh); ntot=np.array(ntot); coll=np.array(coll); mtb=np.array(mtb)
    
    rhocoll, pcoll=ss.spearmanr(np.log(nbh/ntot), np.log(coll/ntot))
    rhomtb, pmtb=ss.spearmanr(np.log(nbh/ntot), np.log(mtb/ntot))
    print rhocoll, pcoll
    print rhomtb, pmtb
    
    plt.figure()
    plt.scatter(np.log(nbh/ntot), np.log(coll/ntot), color='blue', s=8, label='coll')
    plt.scatter(np.log(nbh/ntot), np.log(mtb/ntot), color='red', alpha=0.7, s=8, label='mtb')
    plt.legend(loc='lower left')
    plt.xlabel('log(Nbh/Ntot)')
    plt.ylabel('log(Nbss/Ntot)')
    
    plt.savefig('/Users/shiye/Documents/ClusterGroup/BSSproject/nbhnbss_class.png', dpi = 300)
    #plt.show()


# In[ ]:


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


# In[ ]:


#printout_Nbss_Nbh()


# In[ ]:


#plot_Nbh_Nbss()


# In[ ]:


#hrdiag_LT()


# In[ ]:


#plot_hrdiag('/Volumes/homes/sourav/CMC_results/BH_variations/kickoutputtest_variations/rundir/wind1/z0.001/normalkick/')


# In[ ]:


#printout_BSS(100, 150)


# In[ ]:


#printout_Nbss_sinbin()


# In[ ]:


#class_bss(579,597)


# In[ ]:


#plot__Nbss_sinbin()


# In[ ]:


#plot_Nbss_class()


# In[ ]:


printout_class_bss(308, 350)

