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
import ecc_calc as gwcalc
import LISA_calculations as lisa
import ns_history as nh
import useful_function as uf
import ns_tidalcapture_hdf5 as nhf
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


##Find WD-XX star binaries at a snapshot and check if they are formed in tidal capture/giant collision
def find_WD_XX_atsnap_hdf5(filepath, lowlim, highlim, snapno, savename):
    property_init, property_finl, property_des = nhf.find_tc_properties(filepath)
    ID0 = property_init['id0']; ID1 = property_init['id1']; T=property_init['time']
    Types = property_init['type']

    all_keys = np.genfromtxt(filepath+'snap_keys.txt', dtype = 'str')
    all_snapno = all_keys[:,0]; all_snaptime = all_keys[:,1]

    thekey = '/'+str(snapno)+'(t='+all_snaptime[snapno]+')'

    filestr=filepath+'initial'
    snap = cmct.Snapshot(fname=filepath+'initial.snapshots.h5', 
                        snapshot_name=thekey, 
                        conv=filepath+'initial.conv.sh', 
                        dist=4.52, # distance to cluster in kpc
                        z=0.0038)

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
            if 10<=k0[xx]<=12 and lowlim<=k1[xx]<=highlim:
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


            if lowlim<=k0[xx]<=highlim and 10<=k1[xx]<=12:
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



##Find NS-XX star binaries at a snapshot and check if they are formed in tidal capture/giant collision
def find_WD_XX_atsnap(filepath, lowlim, highlim, snapno, savename):
    property_init, property_finl, property_des = ntc.find_tc_properties_final(filepath)
    ID0 = property_init['id0']; ID1 = property_init['id1']; T=property_init['time']
    Types = property_init['type']

    filestr=filepath+'initial'
    #snaps=dyn.get_snapshots(filestr)
    #lastsnap=snaps[-1]
    snap = filepath+'initial.snap0'+str(snapno)+'.dat.gz'
    t_conv=dyn.conv('t', filestr+'.conv.sh')
    l_conv = dyn.conv('l', filestr+'.conv.sh')
    time=dyn.get_time(snap)*t_conv
    
    #os.system('rm '+savepath)
    fmsb=open(filepath+savename+str(snapno)+'.dat', 'w+')
    fmsb.write('#1.ID0 2.ID1 3.M0 4.M1 5.K0 6.K1 7.a(AU) 8.ecc 9.radrol0 10.radrol1 11.B(G) 12.P(sec) 13.tcflag 14.SN 15.TC_time(Myr) 16.r(pc)\n')

    with gzip.open(snap, 'r') as flast:
        next(flast); next(flast)
        for line in flast:
            datalast=line.split()
            if int(datalast[7])==1:
                if 10<=int(datalast[17])<=12 and lowlim<=int(datalast[18])<=highlim:
                    ID0ms=int(datalast[10]); ID1ms=int(datalast[11])
                    tcflag=4
                    tctime=-100
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

                    fmsb.write('%d %d %f %f %d %d %f %f %f %f %e %f %d %d %f %f\n'%(int(datalast[10]), int(datalast[11]), float(datalast[8]), float(datalast[9]), int(datalast[17]), int(datalast[18]), float(datalast[12]), float(datalast[13]), float(datalast[43]), float(datalast[44]), float(datalast[47]), float(twopi*yearsc/float(datalast[45])), tcflag, int(datalast[49]), tctime, float(datalast[2])*l_conv))


                if lowlim<=int(datalast[17])<=highlim and 10<=int(datalast[18])<=12:
                    ID0ms=int(datalast[11]); ID1ms=int(datalast[10])
                    tcflag=4
                    tctime=-100
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
                                
                    fmsb.write('%d %d %f %f %d %d %f %f %f %f %e %f %d %d %f %f\n'%(int(datalast[11]), int(datalast[10]), float(datalast[9]), float(datalast[8]), int(datalast[18]), int(datalast[17]), float(datalast[12]), float(datalast[13]), float(datalast[44]), float(datalast[43]), float(datalast[48]), float(twopi*yearsc/float(datalast[46])), tcflag, int(datalast[50]),tctime, float(datalast[2])*l_conv))

    fmsb.close()




##Find BH-XX star binaries at a snapshot and check if they are formed in tidal capture/giant collision
def find_BH_XX_atsnap_hdf5(filepath, lowlim, highlim, snapno, savename):
    property_init, property_finl, property_des = nhf.find_tc_properties(filepath)
    ID0 = property_init['id0']; ID1 = property_init['id1']; T=property_init['time']
    Types = property_init['type']

    all_keys = np.genfromtxt(filepath+'snap_keys.txt', dtype = 'str')
    all_snapno = all_keys[:,0]; all_snaptime = all_keys[:,1]

    thekey = '/'+str(snapno)+'(t='+all_snaptime[snapno]+')'

    filestr=filepath+'initial'
    snap = cmct.Snapshot(fname=filepath+'initial.snapshots.h5', 
                        snapshot_name=thekey, 
                        conv=filepath+'initial.conv.sh', 
                        dist=4.52, # distance to cluster in kpc
                        z=0.0038)

    t_conv=dyn.conv('t', filestr+'.conv.sh')
    l_conv=dyn.conv('l', filestr+'.conv.sh')
    
    #os.system('rm '+savepath)
    fmsb=open(filepath+savename+str(snapno)+'.dat', 'a+')
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
            if k0[xx]==14 and lowlim<=k1[xx]<=highlim:
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


            if lowlim<=k0[xx]<=highlim and k1[xx]==14:
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



##Find BH-XX star binaries at a snapshot and check if they are formed in tidal capture/giant collision
def find_BH_XX_atsnap(filepath, lowlim, highlim, snapno, savename):
    property_init, property_finl, property_des = ntc.find_tc_properties_final(filepath)
    ID0 = property_init['id0']; ID1 = property_init['id1']; T=property_init['time']
    Types = property_init['type']

    filestr=filepath+'initial'
    #snaps=dyn.get_snapshots(filestr)
    #lastsnap=snaps[-1]
    snap = filepath+'initial.snap0'+str(snapno)+'.dat.gz'
    t_conv=dyn.conv('t', filestr+'.conv.sh')
    l_conv = dyn.conv('l', filestr+'.conv.sh')
    time=dyn.get_time(snap)*t_conv
    
    #os.system('rm '+savepath)
    fmsb=open(filepath+savename+str(snapno)+'.dat', 'w+')
    fmsb.write('#1.ID0 2.ID1 3.M0 4.M1 5.K0 6.K1 7.a(AU) 8.ecc 9.radrol0 10.radrol1 11.B(G) 12.P(sec) 13.tcflag 14.SN 15.TC_time(Myr) 16.r(pc)\n')

    with gzip.open(snap, 'r') as flast:
        next(flast); next(flast)
        for line in flast:
            datalast=line.split()
            if int(datalast[7])==1:
                if int(datalast[17])==14 and lowlim<=int(datalast[18])<=highlim:
                    ID0ms=int(datalast[10]); ID1ms=int(datalast[11])
                    tcflag=4
                    tctime=-100
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

                    fmsb.write('%d %d %f %f %d %d %f %f %f %f %e %f %d %d %f %f\n'%(int(datalast[10]), int(datalast[11]), float(datalast[8]), float(datalast[9]), int(datalast[17]), int(datalast[18]), float(datalast[12]), float(datalast[13]), float(datalast[43]), float(datalast[44]), float(datalast[47]), float(twopi*yearsc/float(datalast[45])), tcflag, int(datalast[49]), tctime, float(datalast[2])*l_conv))


                if lowlim<=int(datalast[17])<=highlim and int(datalast[18])==14:
                    ID0ms=int(datalast[11]); ID1ms=int(datalast[10])
                    tcflag=4
                    tctime=-100
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
                                
                    fmsb.write('%d %d %f %f %d %d %f %f %f %f %e %f %d %d %f %f\n'%(int(datalast[11]), int(datalast[10]), float(datalast[9]), float(datalast[8]), int(datalast[18]), int(datalast[17]), float(datalast[12]), float(datalast[13]), float(datalast[44]), float(datalast[43]), float(datalast[48]), float(twopi*yearsc/float(datalast[46])), tcflag, int(datalast[50]),tctime, float(datalast[2])*l_conv))

    fmsb.close()