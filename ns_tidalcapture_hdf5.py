import numpy as np 
import pandas as pd
import os,sys
from collections import Counter
import re
import gzip
import scripts
import scripts1
import scripts2
import scripts3
import dynamics as dyn
import unit_convert as uc
import gw_ecc_calc as gwcalc
import LISA_calculations as lisa
import ns_history as nh
import useful_function as uf

sys.path.insert(1, '/fs/lustre/cita/claireshiye/cmctoolkit')
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


##Types of tidal captures in the tidal capture file
def find_tc_properties(filepath):
    filestr=filepath+'initial'
    tcfile=filestr+'.tidalcapture.log'
    t=[]; types=[]; id0i=[]; id1i=[]; m0i=[]; m1i=[]; k0i=[]; k1i=[]; r0i=[]; r1i=[]; rperi=[]
    a=[]; e=[]; id0f=[]; id1f=[]; m0f=[]; m1f=[]; k0f=[]; k1f=[]; r0f=[]; r1f=[]
    v_inf = []; r_cm = []
    mc0=[]; mc1=[]; rc0=[]; rc1=[]

    t_des=[]; type_des=[]; rperi_des=[]; idf_des=[]; mf_des=[]; kf_des=[]
    id0_des=[]; id1_des=[]; m0_des=[]; m1_des=[]; k0_des=[]; k1_des=[]; r0_des=[]; r1_des=[]
    v_inf_des=[]; r_cm_des = []
    mc0_des=[]; mc1_des=[]; rc0_des=[]; rc1_des=[]

    n_tc_merged = 0; n_giant_coll = 0
    with open(tcfile, 'r') as ftc:
        next(ftc)
        for line in ftc:
            data=line.split()
            #print(data)
            #numstr=re.findall(r"\d*\.\d+|\d+", data[2])


            if data[0] == 'coll_CE_debug': continue
            if data[1]=='SS_COLL_TC_P_FAILED' or data[1]=='SS_COLL_GW':
                continue
            if data[2][-1] == 'd' and data[1] == 'SS_COLL_TC_P': 
                n_tc_merged+=1
                numstr=re.split(',|\(|\)|->|\[|\]', data[2])
                numstr=list(filter(None, numstr))
                numstr=list(filter(lambda x: x!='+', numstr))
                #print(numstr)

                numstr9 = numstr[9].split('+')
                numstr9=list(filter(None, numstr9))
                #print(numstr9)
                #numstr[9] = numstr[9][1:]
                numstr90 = numstr9[0]; numstr91 = numstr9[1]

                t_des.append(float(data[0])); type_des.append(data[1])
                id0_des.append(int(numstr[0])); id1_des.append(int(numstr[3]))
                m0_des.append(float(numstr[1])); m1_des.append(float(numstr[4]))
                k0_des.append(int(numstr[2])); k1_des.append(int(numstr[5]))
                r0_des.append(float(numstr[6])); r1_des.append(float(numstr[7]))
                rperi_des.append(float(numstr[8]));
                v_inf_des.append(float(numstr90)); r_cm_des.append(float(numstr91))
                mc0_des.append(-100); mc1_des.append(-100)
                rc0_des.append(-100); rc1_des.append(-100)

                continue

            if data[2][-1] == 'd' and data[1] == 'SS_COLL_Giant': 
                n_giant_coll+=1
                numstr=re.split(',|\(|\)|->|\[|\]', data[2])
                numstr=list(filter(None, numstr))
                numstr=list(filter(lambda x: x!='+', numstr))
                #print(numstr)

                numstr9 = numstr[9].split('+')
                numstr9=list(filter(None, numstr9))
                #print(numstr9)
                #numstr[9] = numstr[9][1:]
                numstr90 = numstr9[0]; numstr91 = numstr9[1]
                #print(numstr[9])

                t_des.append(float(data[0])); type_des.append(data[1])
                id0_des.append(int(numstr[0])); id1_des.append(int(numstr[3]))
                m0_des.append(float(numstr[1])); m1_des.append(float(numstr[4]))
                k0_des.append(int(numstr[2])); k1_des.append(int(numstr[5]))
                r0_des.append(float(numstr[6])); r1_des.append(float(numstr[7]))
                rperi_des.append(float(numstr[8])) 
                v_inf_des.append(float(numstr90)); r_cm_des.append(float(numstr91))
                mc0_des.append(float(numstr[10])); mc1_des.append(float(numstr[11]))
                rc0_des.append(float(numstr[12])); rc1_des.append(float(numstr[13]))

                continue

            if data[2][-1] != 'd' and data[1] == 'SS_COLL_TC_P':
                numstr=re.split(',|\(|\)|->|\[|\]', data[2])
                numstr=list(filter(None, numstr))
                numstr=list(filter(lambda x: x!='+', numstr))
                #print(numstr)

                numstr9 = numstr[9].split('+')
                numstr9=list(filter(None, numstr9))
                #print(numstr9)
                #numstr[9] = numstr[9][1:]
                numstr90 = numstr9[0]; numstr91 = numstr9[1]
                numstr[13] = numstr[13][1:]; numstr[14] = numstr[14][:-1]
                #print(numstr)


                #print(numstr)
                ##Initial properties
                t.append(float(data[0]))
                types.append(data[1])
                id0i.append(int(numstr[0])); id1i.append(int(numstr[3]))
                m0i.append(float(numstr[1])); m1i.append(float(numstr[4]))
                k0i.append(int(numstr[2])); k1i.append(int(numstr[5]))
                r0i.append(float(numstr[6])); r1i.append(float(numstr[7]))
                rperi.append(float(numstr[8]))
                v_inf.append(float(numstr90)); r_cm.append(float(numstr91))
                mc0.append(-100); mc1.append(-100)
                rc0.append(-100); rc1.append(-100)

                ##Final properties
                id0f.append(int(numstr[10])); id1f.append(int(numstr[15]))
                m0f.append(float(numstr[11])); m1f.append(float(numstr[16]))
                k0f.append(int(numstr[12])); k1f.append(int(numstr[17]))
                a.append(float(numstr[13])); e.append(float(numstr[14]))
                r0f.append(float(numstr[18])); r1f.append(float(numstr[19]))


            if data[2][-1] != 'd' and data[1] == 'SS_COLL_Giant':
                numstr=re.split(',|\(|\)|->|\[|\]', data[2])
                numstr=list(filter(None, numstr))
                numstr=list(filter(lambda x: x!='+', numstr))
                #print(numstr)

                numstr9 = numstr[9].split('+')
                numstr9=list(filter(None, numstr9))
                #print(numstr9)
                #numstr[9] = numstr[9][1:]
                numstr90 = numstr9[0]; numstr91 = numstr9[1]
                numstr[17] = numstr[17][1:]; numstr[18] = numstr[18][:-1]
                #print(numstr)

                #print(numstr)
                ##Initial properties
                t.append(float(data[0]))
                types.append(data[1])
                id0i.append(int(numstr[0])); id1i.append(int(numstr[3]))
                m0i.append(float(numstr[1])); m1i.append(float(numstr[4]))
                k0i.append(int(numstr[2])); k1i.append(int(numstr[5]))
                r0i.append(float(numstr[6])); r1i.append(float(numstr[7]))
                rperi.append(float(numstr[8]))
                v_inf.append(float(numstr90)); r_cm.append(float(numstr91))
                mc0.append(float(numstr[10])); mc1.append(float(numstr[11]))
                rc0.append(float(numstr[12])); rc1.append(float(numstr[13]))

                ##Final properties
                id0f.append(int(numstr[14])); id1f.append(int(numstr[19]))
                m0f.append(float(numstr[15])); m1f.append(float(numstr[20]))
                k0f.append(int(numstr[16])); k1f.append(int(numstr[21]))
                a.append(float(numstr[17])); e.append(float(numstr[18]))
                r0f.append(float(numstr[22])); r1f.append(float(numstr[23]))

    Prop_init = {'time':t, 'type':types, 'id0': id0i, 'id1': id1i, 'm0': m0i, 'm1': m1i, 'k0': k0i, 'k1': k1i, 'r0': r0i, 'r1': r1i, 'rperi': rperi, 'vinf': v_inf, 'mc0': mc0, 'mc1': mc1, 'rc0': rc0, 'rc1': rc1}
    Prop_finl = {'id0': id0f, 'id1': id1f, 'm0': m0f, 'm1': m1f, 'k0': k0f, 'k1': k1f, 'r0': r0f, 'r1': r1f, 'sma': a, 'ecc': e}
    Prop_des = {'time':t_des, 'type':type_des, 'id0': id0_des, 'id1': id1_des, 'm0': m0_des, 'm1': m1_des, 'k0': k0_des, 'k1': k1_des, 'r0': r0_des, 'r1': r1_des, 'rperi': rperi_des, 'vinf': v_inf_des, 'mc0': mc0_des, 'mc1': mc1_des, 'rc0': rc0_des, 'rc1': rc1_des}

    #print('n_giant_coll:', n_giant_coll, 'n_tc_merged:', n_tc_merged)

    return Prop_init, Prop_finl, Prop_des



##Find NS-XX star binaries at a snapshot and check if they are formed in tidal capture/giant collision
def find_NS_XX_atsnap(filepath, lowlim, highlim, snapno, savename, thedist, themetal):
    property_init, property_finl, property_des = find_tc_properties(filepath)
    ID0 = property_init['id0']; ID1 = property_init['id1']; T=property_init['time']
    Types = property_init['type']

    all_keys = np.genfromtxt(filepath+'snap_keys.txt', dtype = 'str')
    all_snapno = all_keys[:,0]; all_snaptime = all_keys[:,1]

    thekey = '/'+str(snapno)+'(t='+all_snaptime[snapno]+')'

    filestr=filepath+'initial'
    snap = cmct.Snapshot(fname=filepath+'initial.snapshots.h5', 
                        snapshot_name=thekey, 
                        conv=filepath+'initial.conv.sh', 
                        dist=thedist, # distance to cluster in kpc
                        z=themetal)

    t_conv=dyn.conv('t', filestr+'.conv.sh')
    l_conv=dyn.conv('l', filestr+'.conv.sh')
    
    #os.system('rm '+savepath)
    fmsb=open(filepath+savename+str(snapno)+'.dat', 'w+')
    fmsb.write('#1.ID0 2.ID1 3.M0 4.M1 5.K0 6.K1 7.a(AU) 8.ecc 9.radrol0 10.radrol1 11.B(G) 12.P(sec) 13.tcflag 14.SN 15.TC_time(Myr) 16.R(pc)\n')


    binflag = snap.data['binflag']
    k0 = snap.data['bin_startype0']; k1 = snap.data['bin_startype1']
    id0 = snap.data['id0']; id1 = snap.data['id1']
    m0 = snap.data['m0_MSUN']; m1 = snap.data['m1_MSUN']
    sma = snap.data['a_AU']; ecc = snap.data['e']
    rad0 = snap.data['radrol0']; rad1 = snap.data['radrol1']
    B0 = snap.data['B0']; B1 = snap.data['B1']
    ospin0 = snap.data['ospin0']; ospin1 = snap.data['ospin1']
    SN0 = snap.data['formation0']; SN1 = snap.data['formation1']
    Rpc = snap.data['r']*l_conv

    for xx in range(len(binflag)):
        if binflag[xx]==1:
            if k0[xx]==13 and lowlim<=k1[xx]<=highlim:
                ID0ms=id0[xx]; ID1ms=id1[xx]
                tcflag=4
                tctime = -100
                for ii in range(len(ID0)):    
                    if (ID0ms==ID0[ii] and ID1ms==ID1[ii]) or (ID1ms==ID0[ii] and ID0ms==ID1[ii]):
                        if Types[ii]=='SS_COLL_Giant':
                            tcflag=81
                            tctime = T[ii]*t_conv
                            break
                        if Types[ii]=='SS_COLL_TC_P':
                            tcflag=91
                            tctime = T[ii]*t_conv
                            break
                    elif ID0ms==ID0[ii] or ID0ms==ID1[ii]:
                        if Types[ii]=='SS_COLL_Giant':
                            tcflag=82
                            tctime = T[ii]*t_conv
                            break
                        if Types[ii]=='SS_COLL_TC_P':
                            tcflag=92
                            tctime = T[ii]*t_conv
                            break

                fmsb.write('%d %d %f %f %d %d %f %f %f %f %e %f %d %d %f %f\n'%(id0[xx], id1[xx], m0[xx], m1[xx], k0[xx], k1[xx], sma[xx], ecc[xx], rad0[xx], rad1[xx], B0[xx], twopi*yearsc/ospin0[xx], tcflag, SN0[xx], tctime, Rpc[xx]))


            if lowlim<=k0[xx]<=highlim and k1[xx]==13:
                ID0ms=id1[xx]; ID1ms=id0[xx]
                tcflag=4
                tctime = -100
                for ii in range(len(ID0)):    
                    if (ID0ms==ID0[ii] and ID1ms==ID1[ii]) or (ID1ms==ID0[ii] and ID0ms==ID1[ii]):
                        if Types[ii]=='SS_COLL_Giant':
                            tcflag=81
                            tctime = T[ii]*t_conv
                            break
                        if Types[ii]=='SS_COLL_TC_P':
                            tcflag=91
                            tctime = T[ii]*t_conv
                            break
                    elif ID0ms==ID0[ii] or ID0ms==ID1[ii]:
                        if Types[ii]=='SS_COLL_Giant':
                            tcflag=82
                            tctime = T[ii]*t_conv
                            break
                        if Types[ii]=='SS_COLL_TC_P':
                            tcflag=92
                            tctime = T[ii]*t_conv    
                            break
                                
                
                fmsb.write('%d %d %f %f %d %d %f %f %f %f %e %f %d %d %f %f\n'%(id1[xx], id0[xx], m1[xx], m0[xx], k1[xx], k0[xx], sma[xx], ecc[xx], rad1[xx], rad0[xx], B1[xx], twopi*yearsc/ospin1[xx], tcflag, SN1[xx], tctime, Rpc[xx]))
    fmsb.close()
