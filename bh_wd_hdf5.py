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
def find_WD_XX_atsnap_hdf5(filepath, lowlim, highlim, snapno, savename, thedist, themetal):
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


##Check where some massive WDs come from
def check_massive_WD(pathlist):
    ntot_col = 0; ntot_mer = 0
    ncol_no = 0; nmer_no = 0
    for xx in range(len(pathlist)):
        ns_wdwd_col = np.genfromtxt(pathlist[xx]+'ns_wdwd_coll.dat')
        k0_col = ns_wdwd_col[:,10]; k1_col = ns_wdwd_col[:,11]
        m0_col = ns_wdwd_col[:,5]; m1_col = ns_wdwd_col[:,6]
        id_col = ns_wdwd_col[:,2]
        idcol_check1 = id_col[m0_col>=1.3]
        idcol_check2 = id_col[m1_col>=1.3]
        #idcol_check = np.concatenate((idcol_check1,idcol_check2))
        idcol_check = [idcol_check1, idcol_check2]
        k0col_check = k0_col[m0_col>=1.3]
        k1col_check = k1_col[m1_col>=1.3]
        #kcol_check = np.concatenate((k0col_check,k1col_check))
        kcol_check = [k0col_check, k1col_check]
        print('coll', len(k0col_check)+len(k1col_check))
        ntot_col+=len(k0col_check)+len(k1col_check)

        ns_wdwd_mer = np.genfromtxt(pathlist[xx]+'ns_wdwd_merger.dat')
        k0_mer = ns_wdwd_mer[:,8]; k1_mer = ns_wdwd_mer[:,9]
        m0_mer = ns_wdwd_mer[:,5]; m1_mer = ns_wdwd_mer[:,6]
        id_mer = ns_wdwd_mer[:,2]
        idmer_check1 = id_mer[m0_mer>=1.3]
        idmer_check2 = id_mer[m1_mer>=1.3]
        #idmer_check = np.concatenate((idmer_check1, idmer_check2))
        idmer_check = [idmer_check1,idmer_check2]
        k0mer_check = k0_mer[m0_mer>=1.3]
        k1mer_check = k1_mer[m1_mer>=1.3]
        #kmer_check = np.concatenate((k0mer_check, k1mer_check))
        kmer_check = [k0mer_check, k1mer_check]
        print('merge', len(k0mer_check)+len(k1mer_check))
        ntot_mer+=len(k0mer_check)+len(k1mer_check)


        coldata = scripts1.collision(pathlist[xx]+'initial.collision.log')
        collkeys = list(coldata.keys())
        merdata = scripts2.readmergefile(pathlist[xx]+'initial.semergedisrupt.log')

        progen_idcol = []; progen_mcol = []; progen_kcol = []
        progen_idmer = []; progen_mmer = []; progen_kmer = []
        for ii in range(len(idcol_check)):
            for kk in range(len(idcol_check[ii])):
                par = coldata[int(idcol_check[ii][kk])]['parents']
                if len(np.array(par['IDs'])[np.array(par['types'])==kcol_check[ii][kk]])>1:
                    theidcol = np.array(par['IDs'])[np.array(par['types'])==kcol_check[ii][kk]][ii]
                else:
                    theidcol = np.array(par['IDs'])[np.array(par['types'])==kcol_check[ii][kk]][0]
                progen_idcol.append(theidcol)

                checkcol = 0
    
                if theidcol in collkeys:
                    gpar = coldata[theidcol]['parents']
                    progen_mcol.append(gpar['masses'])
                    progen_kcol.append(gpar['types'])
                    checkcol=1
                else:
                    for yy in range(len(merdata)):
                        merline = merdata[yy].split()
                        if len(merline)==12 and (theidcol == int(merline[2])):
                            progen_mcol.append([float(merline[5]), float(merline[7])])
                            progen_kcol.append([int(merline[-2]), int(merline[-1])])
                            checkcol=1
                if checkcol==0:
                    progen_mcol.append(-100); progen_kcol.append(-100)
                    ncol_no+=1


        for jj in range(len(idmer_check)):
            for ll in range(len(idmer_check[jj])):
                for mm in range(len(merdata)):
                    merline = merdata[mm].split()
                    if int(merline[2]) ==  int(idmer_check[jj][ll]):
                        if int(merline[jj-2])==kmer_check[jj][ll]:
                            theidmer = int(merline[4+2*jj])
    
                progen_idmer.append(theidmer)
    
                checkmer = 0
    
                if theidmer in collkeys:
                    gpar = coldata[theidmer]['parents']
                    progen_mmer.append(gpar['masses'])
                    progen_kmer.append(gpar['types'])
                    checkmer=1
                else:
                    for zz in range(len(merdata)):
                        merline = merdata[zz].split()
                        if len(merline)==12 and (theidmer == int(merline[2])):
                            progen_mmer.append([float(merline[5]), float(merline[7])])
                            progen_kmer.append([int(merline[-2]), int(merline[-1])])
                            checkmer=1
                if checkmer==0:
                    progen_mmer.append(-100); progen_kmer.append(-100)
                    nmer_no+=1
        

        print(progen_idcol, progen_mcol, progen_kcol)
        print(progen_idmer, progen_mmer, progen_kmer)


    print(ntot_col, ntot_mer)
    print(ncol_no, nmer_no)


##Check if some WD binaries are from primordial binaries
def check_primordial(pathlist):
    nmer_tot = np.zeros(len(pathlist))
    nmer_tot_early = np.zeros(len(pathlist))
    nmer_tot_late = np.zeros(len(pathlist))
    n_prim = np.zeros(len(pathlist))
    n_prim_early = np.zeros(len(pathlist))
    n_prim_late = np.zeros(len(pathlist))
    for xx in range(len(pathlist)):
        data_wdmer = np.genfromtxt(pathlist[xx]+'ns_wdwd_merger.dat')
        t_mer = data_wdmer[:,1]  ##in Myr
        id_mer = data_wdmer[:,2]
        nmer_tot[xx]+=len(t_mer)
        nmer_tot_early[xx]+=len(t_mer[t_mer<=6000.])
        nmer_tot_late[xx]+=len(t_mer[t_mer>6000.])

        snap = cmct.Snapshot(fname=pathlist[xx]+'initial.window.snapshots.h5', snapshot_name='0(t=0Gyr)', 
                         conv=pathlist[xx]+'initial.conv.sh', 
                         dist=4.125, # distance to cluster in kpc
                         z=0.0002)
        binflag_snap = np.array(snap.data['binflag'])
        id0_snap = np.array(snap.data['id0'])[binflag_snap==1]
        id1_snap = np.array(snap.data['id1'])[binflag_snap==1]

        merdata = scripts2.readmergefile(pathlist[xx]+'initial.semergedisrupt.log')
        for yy in range(len(id_mer)):
            for kk in range(len(merdata)):
                merline = merdata[kk].split()
                if int(merline[2])==int(id_mer[yy]):
                    if int(merline[-3])==13 and 10<=int(merline[-2])<=12 and 10<=int(merline[-1])<=12:
                        prog_ids = [int(merline[4]), int(merline[6])]
                        break
                
            for ii in range(len(id0_snap)):
                if prog_ids[0] == id0_snap[ii] and prog_ids[1] == id1_snap[ii]:
                    n_prim[xx] +=1
                    if t_mer[yy]<=6000:
                        n_prim_early[xx]+=1
                    if t_mer[yy]>6000:
                        n_prim_late[xx]+=1
                if prog_ids[1] == id0_snap[ii] and prog_ids[0] == id1_snap[ii]:
                    n_prim[xx] +=1
                    if t_mer[yy]<=6000:
                        n_prim_early[xx]+=1
                    if t_mer[yy]>6000:
                        n_prim_late[xx]+=1

    print(nmer_tot)
    print(n_prim)
    print(nmer_tot_early)
    print(n_prim_early)
    print(nmer_tot_late)
    print(n_prim_late)
    print(sum(nmer_tot), sum(n_prim), sum(nmer_tot_early), sum(n_prim_early), sum(nmer_tot_late), sum(n_prim_late))



##Find BH-XX star binaries at a snapshot and check if they are formed in tidal capture/giant collision
def find_BH_XX_atsnap_hdf5(filepath, lowlim, highlim, snapno, savename, thedist, themetal):
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
                        dist=thedist, # distance to cluster in kpc
                        z=themetal)

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