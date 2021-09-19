import numpy as np
import pandas as pd
from glob import glob
import collections
from collections import Counter
import os,sys
import subprocess
#import seaborn as sns
import gzip
import math
import re
#import history_maker_full4 as hi4
import history_cmc as hic
import dynamics as dyn
import scripts3
import scripts1
import scripts2
import ns_tidalcapture as ntc

sys.path.insert(1, '/projects/b1095/syr904/MyCodes/cmctoolkit')
import cmctoolkit as cmct


yearsc=31557600.
twopi=6.283185307179586
Gconst=6.674*10**-8 ##cm3*g-1*s-2
clight=3*10**10 ##cm/s
Msun=2*10**33 ##gram
AU=1.496*10**13  ##cm
PC=3.086*10**18  ##cm
Kconst=9.87*10**-48 ##yr/G^2
Lsun=4.02*10**16 ##mJy*kpc^2


def conv_dict(): return {'l':15, 't':19, 'm':7}    #?


def conv(unit,filepath):   # Returns the unit conversion multiplier given a simulation's *.conv.sh file and a unit (either 'l' or 't')
    dict = conv_dict()
    from re import findall
    with open(filepath,'r') as f:
        head = [next(f) for x in range(24)]
        ##print head[dict[unit]]
    return float(head[dict[unit]].strip().split('=')[-1])
    #return float(findall('\d+[\.]?\d*',head[dict[unit]])[0])


def read_keys(thekey):
    return re.findall(r'\d+\.\d+|\d+', thekey)



def print_Nns_snap(modelpath):
    all_keys = np.genfromtxt(modelpath+'snap_keys.txt', dtype = str)

    t_conv = conv('t', modelpath+'initial.conv.sh')

    try:
        fh = open(modelpath+'initial.ns.dat', 'r')
    except:
        #if True:
        fhandle=open(modelpath+'initial.ns.dat', 'a+')
        fhandle.write('#1.Totaltime, 2.Nns,tot, 3.Nns,single, 4.Nns,binary, 5.Nns,mtb, 6.Npulsar, 7.Nmsp, 8.Nns-ns, 9.Nns-bh, 10.Nns-wd, 11.Nns-ms, 12.Nns-postms\n')
        for xx in range(len(all_keys[:,0])):
            keys_str = '/'+all_keys[:,0][xx]+'(t='+all_keys[:,1][xx]+')'

            snap = cmct.Snapshot(fname=modelpath+'initial.snapshots.h5', 
                                snapshot_name=keys_str, conv=modelpath+'initial.conv.sh', 
                                dist=4.52, # distance to cluster in kpc
                                z=0.0038)

            T = float(all_keys[:,1][xx])*t_conv

            binflag = snap.data['binflag']
            ktype = snap.data['startype']
            #print(ktype)
            k0 = snap.data['bin_startype0']; k1 = snap.data['bin_startype1']
            ospin = snap.data['ospin']
            ospin0 = snap.data['ospin0']; ospin1 = snap.data['ospin1']
            B = snap.data['B']
            B0 = snap.data['B0']; B1 = snap.data['B1']
            rad0 = snap.data['radrol0']; rad1 = snap.data['radrol1']

            N_NS=0; N_NS_SIN=0; N_NS_BIN=0; N_NS_MTB=0; N_PULS=0; N_MSP=0; N_NSNS=0; N_NSBH=0; N_NSWD=0; N_NSMS=0; N_NSPOSTMS=0

            for kk in range(len(binflag)): 
                if binflag[kk]!=1:
                    #print('yes1')
                    if ktype[kk]==13: 
                        #print('yes2')
                        N_NS+=1; N_NS_SIN+=1
                        spin=twopi*yearsc/ospin[kk]
                        deathcut=(spin**2)*(0.17*10**12)
                        if deathcut<B[kk]: 
                            N_PULS+=1
                            if spin<=0.03: N_MSP+=1
                else:
                    if k0[kk]==13:
                        #print('yes3')
                        N_NS+=1; N_NS_BIN+=1
                        spin0=twopi*yearsc/ospin0[kk]
                        deathcut0=(spin0**2)*(0.17*10**12)
                        if rad1[kk]>=1: N_NS_MTB+=1
                        if deathcut0<B0[kk]: 
                            N_PULS+=1   
                            if spin0<=0.03: N_MSP+=1

                        if k1[kk]<2: N_NSMS+=1
                        elif k1[kk]>=10 and k1[kk]<=12: N_NSWD+=1
                        elif k1[kk]==13: N_NSNS+=1
                        elif k1[kk]==14: N_NSBH+=1
                        else: N_NSPOSTMS+=1

                    if k1[kk]==13:
                        #print('yes4')
                        N_NS+=1; N_NS_BIN+=1
                        spin1=twopi*yearsc/ospin1[kk]
                        deathcut1=(spin1**2)*(0.17*10**12)
                        if rad0[kk]>=1: N_NS_MTB+=1
                        if deathcut1<B1[kk]: 
                            N_PULS+=1
                            if spin1<=0.03: N_MSP+=1

                        if k0[kk]<2: N_NSMS+=1
                        elif k0[kk]>=10 and k0[kk]<=12: N_NSWD+=1
                        elif k0[kk]==13: print('already counted')
                        elif k0[kk]==14: N_NSBH+=1
                        else: N_NSPOSTMS+=1


            fhandle.write('%f %d %d %d %d %d %d %d %d %d %d %d\n'%(T, N_NS, N_NS_SIN, N_NS_BIN, N_NS_MTB, N_PULS, N_MSP, N_NSNS, N_NSBH, N_NSWD, N_NSMS, N_NSPOSTMS))

            print(xx)

        fhandle.close()
        


def get_allpsr_snapshot(snapshot, mspflag):
    id0_snap=[]; id1_snap=[]; m0_snap=[]; m1_snap=[]; 
    k0_snap=[]; k1_snap=[]; B_snap=[]; P_snap=[]; 
    a_snap=[]; ecc_snap=[]; FC_snap=[]; dmdt0_snap=[]; dmdt1_snap=[]
    radrol0_snap=[]; radrol1_snap=[]; r_snap=[]

    binflag = snapshot.data['binflag']
    ktype = snapshot.data['startype']
    k0 = snapshot.data['bin_startype0']; k1 = snapshot.data['bin_startype1']
    m = snapshot.data['m_MSUN']
    m0 = snapshot.data['m0_MSUN']; m1 = snapshot.data['m1_MSUN']
    B = snapshot.data['B']; B0 = snapshot.data['B0']; B1 = snapshot.data['B1']
    ospin = snapshot.data['ospin']
    ospin0 = snapshot.data['ospin0']; ospin1 = snapshot.data['ospin1']
    ID = snapshot.data['id']; ID0 = snapshot.data['id0']; ID1 = snapshot.data['id1']
    SN = snapshot.data['formation']
    SN0 = snapshot.data['formation0']; SN1 = snapshot.data['formation1']
    dmdt0 = snapshot.data['dmdt0']; dmdt1 = snapshot.data['dmdt1']
    rad0 = snapshot.data['radrol0']; rad1 = snapshot.data['radrol1']
    sma = snapshot.data['a_AU']; ecc = snapshot.data['e']
    r = snapshot.data['r']


    for ii in range(len(binflag)):
        if binflag[ii] != 1:
            if ktype[ii]==13:
                spin=twopi*yearsc/ospin[ii]
                deathcut=(spin**2)*(0.17*10**12)
                if mspflag=='PSR':
                    if spin>0.03 and B[ii]>=deathcut:
                        id0_snap.append(ID[ii]); id1_snap.append(-100)
                        m0_snap.append(m[ii]); m1_snap.append(-100)
                        k0_snap.append(ktype[ii]); k1_snap.append(-100)
                        FC_snap.append(SN[ii]); B_snap.append(B[ii]); P_snap.append(spin)
                        a_snap.append(-100); ecc_snap.append(-100)
                        dmdt0_snap.append(-100); dmdt1_snap.append(-100)
                        radrol0_snap.append(-100); radrol1_snap.append(-100)
                        r_snap.append(r[ii]) ##code unit

                else:
                    if spin<=0.03 and B[ii]>=deathcut:
                        id0_snap.append(ID[ii]); id1_snap.append(-100)
                        m0_snap.append(m[ii]); m1_snap.append(-100)
                        k0_snap.append(ktype[ii]); k1_snap.append(-100)
                        FC_snap.append(SN[ii]); B_snap.append(B[ii]); P_snap.append(spin)
                        a_snap.append(-100); ecc_snap.append(-100)
                        dmdt0_snap.append(-100); dmdt1_snap.append(-100)
                        radrol0_snap.append(-100); radrol1_snap.append(-100)
                        r_snap.append(r[ii]) ##code unit
                        

        else:
            if k0[ii]==13:
                spin0=twopi*yearsc/ospin0[ii]
                deathcut0=(spin0**2)*(0.17*10**12)
                if mspflag=='PSR':
                    if spin0>0.03 and B0[ii]>=deathcut0:
                        id0_snap.append(ID0[ii]); id1_snap.append(ID1[ii])
                        m0_snap.append(m0[ii]); m1_snap.append(m1[ii])
                        k0_snap.append(k0[ii]); k1_snap.append(k1[ii])
                        FC_snap.append(SN0[ii]); B_snap.append(B0[ii]); P_snap.append(spin0)
                        a_snap.append(sma[ii]); ecc_snap.append(ecc[ii])
                        dmdt0_snap.append(dmdt0[ii]); dmdt1_snap.append(dmdt1[ii])
                        radrol0_snap.append(rad0[ii]); radrol1_snap.append(rad1[ii])
                        r_snap.append(r[ii]) ##code unit

                else:
                    if spin0<=0.03 and B0[ii]>=deathcut0:
                        id0_snap.append(ID0[ii]); id1_snap.append(ID1[ii])
                        m0_snap.append(m0[ii]); m1_snap.append(m1[ii])
                        k0_snap.append(k0[ii]); k1_snap.append(k1[ii])
                        FC_snap.append(SN0[ii]); B_snap.append(B0[ii]); P_snap.append(spin0)
                        a_snap.append(sma[ii]); ecc_snap.append(ecc[ii])
                        dmdt0_snap.append(dmdt0[ii]); dmdt1_snap.append(dmdt1[ii])
                        radrol0_snap.append(rad0[ii]); radrol1_snap.append(rad1[ii])
                        r_snap.append(r[ii]) ##code unit


            if k1[ii]==13:
                spin1=twopi*yearsc/ospin1[ii]
                deathcut1=(spin1**2)*(0.17*10**12)
                if mspflag=='PSR':
                    if spin1>0.03 and B1[ii]>=deathcut1:
                        id0_snap.append(ID1[ii]); id1_snap.append(ID0[ii])
                        m0_snap.append(m1[ii]); m1_snap.append(m0[ii])
                        k0_snap.append(k1[ii]); k1_snap.append(k0[ii])
                        FC_snap.append(SN1[ii]); B_snap.append(B1[ii]); P_snap.append(spin1)
                        a_snap.append(sma[ii]); ecc_snap.append(ecc[ii])
                        dmdt0_snap.append(dmdt1[ii]); dmdt1_snap.append(dmdt0[ii])
                        radrol0_snap.append(rad1[ii]); radrol1_snap.append(rad0[ii])
                        r_snap.append(r[ii]) ##code unit

                else:
                    if spin1<=0.03 and B1[ii]>=deathcut1:
                        id0_snap.append(ID1[ii]); id1_snap.append(ID0[ii])
                        m0_snap.append(m1[ii]); m1_snap.append(m0[ii])
                        k0_snap.append(k1[ii]); k1_snap.append(k0[ii])
                        FC_snap.append(SN1[ii]); B_snap.append(B1[ii]); P_snap.append(spin1)
                        a_snap.append(sma[ii]); ecc_snap.append(ecc[ii])
                        dmdt0_snap.append(dmdt1[ii]); dmdt1_snap.append(dmdt0[ii])
                        radrol0_snap.append(rad1[ii]); radrol1_snap.append(rad0[ii])
                        r_snap.append(r[ii]) ##code unit
    
    return B_snap, P_snap, id0_snap, id1_snap, m0_snap, m1_snap, k0_snap, k1_snap, a_snap, ecc_snap, FC_snap, dmdt0_snap, dmdt1_snap, radrol0_snap, radrol1_snap, r_snap


##Find all the pulsars at a snapshot
def get_allpsr_atsnap(modelpath, mspfg, snapno):
    T=[]; BF=[]; S=[]; F=[]; ID_0=[]; ID_1=[]; M_0=[]; M_1=[]; K_0=[]; K_1=[]; Aaxis=[]; E=[]; DMDT0=[]; DMDT1=[]; RAD0=[]; RAD1=[]; RADIUS=[]; TCFLAG=[]

    all_keys = np.genfromtxt(modelpath+'snap_keys.txt', dtype = 'str')
    all_snapno = all_keys[:,0]; all_snaptime = all_keys[:,1]

    thekey = '/'+str(snapno)+'(t='+all_snaptime[snapno]+')'
        
    pref='initial'
    filestr=modelpath+pref
    t_conv=conv('t', filestr+'.conv.sh')
    l_conv=conv('l', filestr+'.conv.sh')

    snap = cmct.Snapshot(fname=modelpath+'initial.snapshots.h5', snapshot_name=thekey, conv=modelpath+'initial.conv.sh', 
                        dist=4.52, # distance to cluster in kpc
                        z=0.0038)

    print('read snap')

    t_age = snap.age

    Bf, Spin, Id0, Id1, M0, M1, K0, K1, A, Ecc, Fc, Dmdt0, Dmdt1, Radrol0, Radrol1, R=get_allpsr_snapshot(snap, mspfg)

    print('get psr data')

    prop_init, prop_finl, prop_des=ntc.find_tc_properties_final(modelpath)
    tc_id0 = prop_init['id0']; tc_id1 = prop_init['id1']; tc_type = prop_init['type']

    tcflag = []
    for i in range(len(Id0)):
        tcflag.append(4)
        for j in range(len(tc_id0)):
            if Id1[i]!=-100:
                if (Id0[i]==tc_id0[j] and Id1[i]==tc_id1[j]) or (Id1[i]==tc_id0[j] and Id0[i]==tc_id1[j]):
                    if tc_type[j]=='SS_COLL_Giant':
                        tcflag[i]=81
                        break
                    if tc_type[j]=='SS_COLL_TC_P':
                        tcflag[i]=91
                        break

                elif Id0[i]==tc_id0[j] or Id0[i]==tc_id1[j]:
                    if tc_type[j]=='SS_COLL_Giant':
                        tcflag[i]=82
                        break
                    if tc_type[j]=='SS_COLL_TC_P':
                        tcflag[i]=92
                        break

            else:
                if Id0[i]==tc_id0[j] or Id0[i]==tc_id1[j]:
                    if tc_type[j]=='SS_COLL_Giant':
                        tcflag[i]=83
                        break
                    if tc_type[j]=='SS_COLL_TC_P':
                        tcflag[i]=93
                        break


    Time=list(np.full_like(Id0, t_age, dtype=np.double))
    #print(Time, Model, Mst)

    T=T+Time; BF=BF+Bf; S=S+Spin; F=F+Fc; ID_0=ID_0+Id0; ID_1=ID_1+Id1; M_0=M_0+M0; M_1=M_1+M1; K_0=K_0+K0; K_1=K_1+K1; Aaxis=Aaxis+A; E=E+Ecc
    DMDT0=DMDT0+Dmdt0; DMDT1=DMDT1+Dmdt1; RAD0=RAD0+Radrol0; RAD1=RAD1+Radrol1; RADIUS=RADIUS+R
    TCFLAG=TCFLAG+tcflag

    RADIUS=np.array(RADIUS)*l_conv
    print('done')

    print(len(T), len(RADIUS), len(BF), len(S), len(DMDT0), len(DMDT1), len(RAD0), len(RAD1), len(M_0), len(M_1), len(ID_0), len(ID_1), len(K_0), len(K_1), len(Aaxis), len(E), len(F), len(TCFLAG))
    np.savetxt(modelpath+mspfg+str(snapno)+'.dat', np.c_[T, RADIUS, BF, S, DMDT0, DMDT1, RAD0, RAD1, M_0, M_1, ID_0, ID_1, K_0, K_1, Aaxis, E, F, TCFLAG], fmt ='%f %f %e %f %f %f %f %f %f %f %d %d %d %d %f %f %d %d', delimiter= ' ', header = '1.Time(Myr) 2.r(pc) 3.B(G) 4.P(sec) 5.dmdt0(Msun/yr) 6.dmdt1(Msun/yr) 7.rolrad0 8.rolrad1 9.m0(Msun) 10.m1(Msun) 11.ID0 12.ID1 13.k0 14.k1 15.a(AU) 16.ecc 17.Formation 18.TCflag', comments = '#')
