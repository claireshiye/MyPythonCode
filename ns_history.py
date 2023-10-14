import numpy as np
from glob import glob
import collections
from collections import Counter
import os,sys
import subprocess
import gzip
import math
import re
import history_cmc as hic
import dynamics as dyn
import scripts3
import scripts1
import scripts2
import ns
import ns_tidalcapture_hdf5 as ntc_hdf5
#from scipy import stats

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


##Find MSP tidal capture/giant collision history
def get_tcgc_msp(modelpath, mspfile):
    msp_data = np.genfromtxt(modelpath+mspfile+'.dat')
    id0 = msp_data[:,10]; id1 = msp_data[:,11]
    tcflag = msp_data[:,17]
    m1 = msp_data[:,9]; k1 = msp_data[:,13]

    prop_init, prop_finl, prop_des=ntc.find_tc_properties_final(modelpath)
    tc_id0 = np.array(prop_init['id0']); tc_id1 = np.array(prop_init['id1']); tc_type = prop_init['type']
    tc_k0 = np.array(prop_init['k0']); tc_k1 = np.array(prop_init['k1'])
    tc_m0 = np.array(prop_init['m0']); tc_m1 = np.array(prop_init['m1'])
    tc_r0 = np.array(prop_init['r0']); tc_r1 = np.array(prop_init['r1'])


    k0_prog = []; k1_prog = []; m0_prog = []; m1_prog = []; r0_prog = []; r1_prog = []
    id0_prog = []; id1_prog = []
    ucxb_flag = []
    tcflag_prog = []

    n_81 = 0; n_82 = 0; n_91 = 0; n_92 = 0
    for ii in range(len(tcflag)):
        ucxb = 0 
        if tcflag[ii] == 81:
            n_81+=1
            if msp_data[:,9][ii]<=0.01 and msp_data[:,13][ii]>10:
                n_81-=1
                ucxb = 1
        
        if tcflag[ii] == 82 or tcflag[ii] == 83:
            n_82+=1
            if msp_data[:,9][ii]<=0.01 and msp_data[:,13][ii]>10:
                n_82-=1
                ucxb = 1
        
        if tcflag[ii] == 91:
            n_91+=1
            if msp_data[:,9][ii]<=0.01 and msp_data[:,13][ii]>10:
                n_91-=1
                ucxb = 1
        
        if tcflag[ii] == 92 or tcflag[ii] == 93:
            n_92+=1
            if msp_data[:,9][ii]<=0.01 and msp_data[:,13][ii]>10:
                n_92-=1

        if tcflag[ii]!=4:
            for xx in range(len(tc_id0)):
                if int(tc_id0[xx])==int(id0[ii]):
                    k0_prog.append(tc_k0[tc_id0 == id0[ii]][0]); k1_prog.append(tc_k1[tc_id0 == id0[ii]][0])
                    m0_prog.append(tc_m0[tc_id0 == id0[ii]][0]); m1_prog.append(tc_m1[tc_id0 == id0[ii]][0])
                    r0_prog.append(tc_r0[tc_id0 == id0[ii]][0]); r1_prog.append(tc_r1[tc_id0 == id0[ii]][0])
                    id0_prog.append(id0[ii]); id1_prog.append(id1[ii])
                    tcflag_prog.append(tcflag[ii])
                    ucxb_flag.append(ucxb)
                elif int(tc_id1[xx]) == int(id0[ii]):
                    k0_prog.append(tc_k1[tc_id1 == id0[ii]][0]); k1_prog.append(tc_k0[tc_id1 == id0[ii]][0])
                    m0_prog.append(tc_m1[tc_id1 == id0[ii]][0]); m1_prog.append(tc_m0[tc_id1 == id0[ii]][0])
                    r0_prog.append(tc_r1[tc_id1 == id0[ii]][0]); r1_prog.append(tc_r0[tc_id1 == id0[ii]][0])
                    id0_prog.append(id0[ii]); id1_prog.append(id1[ii])
                    ucxb_flag.append(ucxb)
                    tcflag_prog.append(tcflag[ii])

        
    np.savetxt(modelpath+mspfile+'_progenitor.dat', np.c_[id0_prog, id1_prog, k0_prog, k1_prog, m0_prog, m1_prog, r0_prog, r1_prog, ucxb_flag, tcflag_prog], fmt = '%d %d %d %d %f %f %f %f %d %d', header = '1.ID0 2.ID1 3.K0 4.K1 5.M0 6.M1 7.R0 8.R1 9.UCXB 10.TCFLAG', comments = '#', delimiter = ' ')

    print(n_81, n_82, n_91, n_92)


##Find the history of an ID in the psrfile
def get_history_inpsrfile(ids, sourcedir):
    pref='initial'
    filestr=sourcedir+pref
    t_conv = dyn.conv('t',filestr+'.conv.sh')
    fname=filestr+'.morepulsars.dat'
    Age=[]; B0=[]; B1=[]; P0=[]; P1=[]; M0=[]; M1=[]; K0=[]; K1=[]; a=[]; ecc=[]
    id0=[]; id1=[]; radrol0=[]; radrol1=[]
    with open(fname, 'r') as fpsr:
        next(fpsr)
        for line in fpsr:
            datapsr=line.split()
            if int(datapsr[2])!=1:
                if int(datapsr[3])==ids:
                    Age.append(float(datapsr[1])); id0.append(int(datapsr[3])); id1.append(-100)
                    M0.append(float(datapsr[5])); M1.append(-100)
                    B0.append(float(datapsr[7])); B1.append(-100)
                    P0.append(float(datapsr[9])); P1.append(-100)
                    K0.append(int(datapsr[11])); K1.append(-100)
                    a.append(-100); ecc.append(-100)
                    radrol0.append(-100); radrol1.append(-100)
            if int(datapsr[2])==1:
                if int(datapsr[3])==ids:
                    Age.append(float(datapsr[1])); id0.append(int(datapsr[3])); id1.append(int(datapsr[4]))
                    M0.append(float(datapsr[5])); M1.append(float(datapsr[6]))
                    B0.append(float(datapsr[7])); B1.append(float(datapsr[8]))
                    P0.append(float(datapsr[9])); P1.append(float(datapsr[10]))
                    K0.append(int(datapsr[11])); K1.append(int(datapsr[12]))
                    a.append(float(datapsr[13])); ecc.append(float(datapsr[14]))
                    radrol0.append(float(datapsr[15])); radrol1.append(float(datapsr[16]))
                if int(datapsr[4])==ids:
                    Age.append(float(datapsr[1])); id0.append(int(datapsr[4])); id1.append(int(datapsr[3]))
                    M0.append(float(datapsr[6])); M1.append(float(datapsr[5]))
                    B0.append(float(datapsr[8])); B1.append(float(datapsr[7]))
                    P0.append(float(datapsr[10])); P1.append(float(datapsr[9]))
                    K0.append(int(datapsr[12])); K1.append(int(datapsr[11]))
                    a.append(float(datapsr[13])); ecc.append(float(datapsr[14]))
                    radrol0.append(float(datapsr[16])); radrol1.append(float(datapsr[15]))

    return Age, B0, B1, P0, P1, M0, M1, K0, K1, a, ecc, id0, id1, radrol0, radrol1


##Find the history of a NS.
def find_star_history(theid, modelpath):
    hist_dict=hic.history_maker([theid],[1],'initial', modelpath, 1.0)
    np.save('/projects/b1095/syr904/projects/PULSAR2/tc_comparison/interact_dict/'+str(theid)+'_dict.npy', hist_dict)

    t, b0, b1, p0, p1, m0, m1, k0, k1, sma, e, i0, i1, rad0, rad1=ns.get_history_inpsrfile(theid, modelpath)
    np.savetxt('/projects/b1095/syr904/projects/PULSAR2/tc_comparison/psr_hist/'+str(theid)+'_hist.dat', np.c_[t, i0, i1, m0, m1, k0, k1, b0, b1, p0, p1, sma, e, rad0, rad1], fmt='%f %d %d %f %f %d %d %e %e %f %f %f %f %f %f', header='1.Time 2.ID0 3.ID1 4.M0 5.M1 6.K0 7.K1 8.B0 9.B1 10.P0 11.P1 12.a(AU) 13.ecc 14.RADROL0 15.RADROL1', comments='#', delimiter='')


##Find the history of multiple NSs
def find_histories():
    tc_bin=np.genfromtxt('/projects/b1095/syr904/projects/PULSAR2/tc_comparison/tc_binarypsr.dat')
    tc_sin=np.genfromtxt('/projects/b1095/syr904/projects/PULSAR2/tc_comparison/tc_singlepsr.dat')

    for i in range(len(tc_bin[:,0])):
        if int(tc_bin[:,5][i])==13:
            find_star_history(int(tc_bin[:,1][i]), '/projects/b1095/syr904/cmc/cmc-mpi-tidalcapture/rundir/8e5rv0.5rg8z0.002/')
        if int(tc_bin[:,6][i])==13:
            find_star_history(int(tc_bin[:,2][i]), '/projects/b1095/syr904/cmc/cmc-mpi-tidalcapture/rundir/8e5rv0.5rg8z0.002/')

        print(i)

    for j in range(len(tc_sin[:,0])):
        if int(tc_sin[:,5][j])==13:
            find_star_history(int(tc_sin[:,1][j]), '/projects/b1095/syr904/cmc/cmc-mpi-tidalcapture/rundir/8e5rv0.5rg8z0.002/')
        if int(tc_sin[:,6][j])==13:
            find_star_history(int(tc_sin[:,2][j]), '/projects/b1095/syr904/cmc/cmc-mpi-tidalcapture/rundir/8e5rv0.5rg8z0.002/')

        print(j)


##Find how accretion-induced collapse NSs are formed
def find_aic_accretion(pathlist, psrfile):
    paths = np.genfromtxt(pathlist, dtype=str)
    data_msp = np.genfromtxt(psrfile)
    models = data_msp[:,0]; ID0 = data_msp[:,12]; ID1 = data_msp[:,13]
    K0 = data_msp[:,14]; K1 = data_msp[:,15]; sn = data_msp[:,18]; tcflag = data_msp[:,19]

    n_ce = 0; n_dwd = 0; n_star = 0
    n_tot = 0

    f = open('/projects/b1095/syr904/cmc/cmc-mpi-tidalcapture/rvgrid/standard_models_tcon/aic_ns_accretion.dat', 'a+')
    for kk in range(len(paths)):
        binint = scripts3.read_binint(paths[kk]+'initial.binint.log')
        #print('binint read')
        #fbinint=open(paths[kk]+'initial.binint.log','r')
        #print 'inputfile %s' %(filename)
        #positions=scripts3.find_positions(fbinint)
        #for mm in range(len(positions)-1):
            #print positions[i]
            #binint=scripts3.read_segment(fbinint,positions[mm])

        prop_init, prop_finl, prop_des = ntc.find_tc_properties_final(paths[kk])
        tcid0 = prop_init['id0']; tcid1 = prop_init['id1']
        tck0 = prop_init['k0']; tck1 = prop_init['k1']
        tck0_fnl = prop_finl['k0']; tck1_fnl = prop_finl['k1']
        #print('tcfile read')


        for ii in range(len(models)):
            if int(models[ii]) == kk:
                check = 0
                id_ns = int(ID0[ii]); id_comp = int(ID1[ii])
                if sn[ii] == 4:
                    n_tot+=1
                    if tcflag[ii] == 81 or tcflag[ii] == 91:            
                        for xx in range(len(tcid0)):
                            if int(tcid0[xx]) == id_ns and int(tcid1[xx]) == id_comp:
                                if tck0[xx] == 13.:
                                    break
                                elif tck0[xx] == 12.:
                                    if 2.<=tck1[xx]<=9. and 10.<=tck1_fnl[xx]<=12.:
                                        n_ce+=1; n_dwd+=1
                                        print(ii, models[ii], id_ns, id_comp, 'DWD, Need CE accretion', file = f)
                                        check = 1
                                        break
                                    elif 2.<=tck1[xx]<=9. and tck1_fnl[xx]==7.:
                                        n_star +=1
                                        print(ii, models[ii], id_ns, id_comp, 'Maybe mass transfer', file = f)
                                        check = 1
                                        break
                                else:
                                    break

                            elif int(tcid1[xx]) == id_ns and int(tcid0[xx]) == id_comp:
                                #print('yes')
                                if tck1[xx] == 13.:
                                    break
                                elif tck1[xx] == 12.:
                                    if 2.<=tck0[xx]<=9. and 10.<=tck0_fnl[xx]<=12.:
                                        n_ce+=1; n_dwd+=1
                                        print(ii, models[ii], id_ns, id_comp, 'DWD, Need CE accretion', file = f)
                                        check = 1
                                        break
                                elif 2.<=tck0[xx]<=9. and tck0_fnl[xx]==7.:
                                        n_star +=1
                                        print(ii, models[ii], id_ns, id_comp, 'Maybe mass transfer', file = f)
                                        check = 1
                                        break
                                else:
                                    break

                    else:
                        for yy in range(len(binint)):
                            outputs = binint[yy]['output']
                            for m in range(len(outputs)):
                                if int(outputs[m]['no'])==2 and not outputs[m]['merge']['ids']:
                                    if int(outputs[m]['ids'][0])==id_ns and int(outputs[m]['ids'][1]) == id_comp:
                                        if int(outputs[m]['startype'][0])==12 and 10<=int(outputs[m]['startype'][1])<=12:
                                            n_dwd +=1
                                            print(ii, models[ii], id_ns, id_comp, 'DWD', file = f) 
                                            check = 1

                                    if int(outputs[m]['ids'][1])==id_ns and int(outputs[m]['ids'][0]) == id_comp:
                                        if int(outputs[m]['startype'][1])==12 and 10<=int(outputs[m]['startype'][0])<=12:
                                            n_dwd +=1
                                            print(ii, models[ii],id_ns, id_comp, 'DWD', file = f) 
                                            check = 1

                                if int(outputs[m]['no'])==3 and not outputs[m]['merge']['ids']:
                                    if int(outputs[m]['ids'][0])==id_ns and int(outputs[m]['ids'][1]) == id_comp:
                                        if int(outputs[m]['startype'][0])==12 and 10<=int(outputs[m]['startype'][1])<=12:
                                            n_dwd +=1
                                            print(ii, models[ii], id_ns, id_comp, 'DWD', file = f) 
                                            check = 1

                                    if int(outputs[m]['ids'][1])==id_ns and int(outputs[m]['ids'][0]) == id_comp:
                                        if int(outputs[m]['startype'][1])==12 and 10<=int(outputs[m]['startype'][0])<=12:
                                            n_dwd +=1
                                            print(ii, models[ii], id_ns, id_comp, 'DWD', file = f) 
                                            check = 1

                            if check == 1:
                                break


                    if check == 0:
                        print(ii, models[ii], ID0[ii], ID1[ii], "Can't find", file = f)

    print(n_tot, n_ce, n_dwd, n_star)
    f.close()



##How many pulsar-wd binaries that could have become psr-ms binaries but not?
def find_psr_wd_evol(filepath, modelpath):
    data_msp=np.genfromtxt(filepath+'MSP_last.dat')
    data_psr=np.genfromtxt(filepath+'PSR_last.dat')
    id_psrs=np.concatenate((data_msp[:,12],data_psr[:,12]))
    id_coms=np.concatenate((data_msp[:,13],data_psr[:,13]))

    psrfile=modelpath+'initial.morepulsars.dat'
    for i in range(len(id_psrs)):
        linepsr=[]
        with open(psrfile, 'r') as fpsr:
            next(fpsr)
            for line in fpsr:
                datapsr=line.split()
                if int(datapsr[2])==1:
                    if int(datapsr[3])==id_psrs[i] and int(datapsr[4])==id_coms[i] and 0<=int(datapsr[12])<=1:
                        linepsr.append(datapsr)
                    if int(datapsr[4])==id_psrs[i] and int(datapsr[3])==id_coms[i] and 0<=int(datapsr[11])<=1:
                        linepsr.append(datapsr)
        
        if linepsr:
            print(id_psrs[i])
            print(linepsr[-1])



##Find how MSPs accrete mass
def MSP_accretion(filepath, modelpath):
    data_msp=np.genfromtxt(filepath+'MSP_last.dat')
    ids=data_msp[:,12]
    
    fswitch=open(filepath+'MSP_switch.dat', 'a+')
    fswitch.write('#1.Time 2.Time(Myr) 3.ID0 4.ID1 5.K0 6.K1 7.B(G) 8.P(sec) 9.A(AU) 10.Ecc 11.radrol_com 12.M0 13.M1\n')

    psrfile=modelpath+'initial.morepulsars.dat'
    #size = os.path.getsize(psrfile)-1
    #print(size)
    t_conv=dyn.conv('t', modelpath+'initial.conv.sh')
    for i in range(len(ids)):
        lineid=[]; Pid=[]; pos=[]

        with open(psrfile, 'r') as fpsr:
            next(fpsr)
            for line in fpsr:
                #size -= len(line)
                #if not size:
                #    break
                datapsr=line.split()
                #print(datapsr)
                #next(fpsr)
                #datanext=line.split()
                #print(datanext)
                t_myr=float(datapsr[1])*t_conv
                if int(datapsr[3])==ids[i]:
                    lineid.append(datapsr)
                    Pid.append(float(datapsr[9])); pos.append(0)
                if int(datapsr[4])==ids[i]:
                    lineid.append(datapsr)
                    Pid.append(float(datapsr[10])); pos.append(1)

                if len(Pid)>1 and Pid[-1]<=0.03:
                    if pos[-2]==0:
                        fswitch.write('%f %f %d %d %d %d %e %f %f %f %f %f %f\n'%(float(lineid[-2][1]), t_myr, int(lineid[-2][3]), int(lineid[-2][4]), int(lineid[-2][11]), int(lineid[-2][12]), float(lineid[-2][7]), float(lineid[-2][9]), float(lineid[-2][13]), float(lineid[-2][14]), float(lineid[-2][16]), float(lineid[-2][5]), float(lineid[-2][6])))
                    if pos[-2]==1:
                        fswitch.write('%f %f %d %d %d %d %e %f %f %f %f %f %f\n'%(float(lineid[-2][1]), t_myr, int(lineid[-2][4]), int(lineid[-2][3]), int(lineid[-2][12]), int(lineid[-2][11]), float(lineid[-2][8]), float(lineid[-2][10]), float(lineid[-2][13]), float(lineid[-2][14]), float(lineid[-2][15]), float(lineid[-2][6]), float(lineid[-2][5])))

                    if pos[-1]==0:
                        fswitch.write('%f %f %d %d %d %d %e %f %f %f %f %f %f\n'%(float(lineid[-1][1]), t_myr, int(lineid[-1][3]), int(lineid[-1][4]), int(lineid[-1][11]), int(lineid[-1][12]), float(lineid[-1][7]), float(lineid[-1][9]), float(lineid[-1][13]), float(lineid[-1][14]), float(lineid[-1][16]), float(lineid[-1][5]), float(lineid[-1][6])))
                    if pos[-1]==1:
                        fswitch.write('%f %f %d %d %d %d %e %f %f %f %f %f %f\n'%(float(lineid[-1][1]), t_myr, int(lineid[-1][4]), int(lineid[-1][3]), int(lineid[-1][12]), int(lineid[-1][11]), float(lineid[-1][8]), float(lineid[-1][10]), float(lineid[-1][13]), float(lineid[-1][14]), float(lineid[-1][15]), float(lineid[-1][6]), float(lineid[-1][5])))
                    break


        print(i)


##Find how a binary is disrupted or if it is disurupted
def find_binary_disruipt(bin_id0, bin_id1, modelpath):
    filestr=modelpath+'initial'
    collfile=filestr+'.collision.log'
    sefile=filestr+'.semergedisrupt.log'
    escfile=filestr+'.esc.dat'

    colldata=scripts1.readcollfile(collfile)
    sedata=scripts2.readmergefile(sefile)

    snaps = dyn.get_snapshots(filestr)
    lastsnap = snaps[-1]

    check_disrupt = -100
    time_disrupt = -100
    type_disrupt = -100

    with gzip.open(lastsnap, 'r') as flast:
        next(flast)
        next(flast)
        for line in flast:
            datalast = line.split()
            if int(datalast[7]) == 1:
                if (int(datalast[10])==bin_id0 and int(datalast[11])==bin_id1) or (int(datalast[11])==bin_id0 and int(datalast[10])==bin_id1):
                    check_disrupt = 'binintact'
                    time_disrupt = dyn.get_time(lastsnap)
                    type_disrupt = '-'

    if check_disrupt==-100:
        for i in range(len(sedata)):
            line=sedata[i].split()
            #print(line)
            if int(line[1])<3:
                if (int(line[4])==int(bin_id0) and int(line[6])==int(bin_id1)) or (int(line[6])==int(bin_id0) and int(line[4])==int(bin_id1)):
                    check_disrupt='binmerge'
                    time_disrupt = float(line[0])
                    type_disrupt = '-'
                    break
            else:
                if (int(line[4])==int(bin_id0) and int(line[6])==int(bin_id1)) or (int(line[6])==int(bin_id0) and int(line[4])==int(bin_id1)):
                    check_disrupt='bindisrupt'
                    time_disrupt = float(line[0])
                    type_disrupt = '-'
                    break

    if check_disrupt==-100:
        for j in range(len(colldata)):
            line=colldata[j].split()
            #print(line)
            if int(line[2])==2:
                collids=[int(line[5]), int(line[7])]
            elif len(line)==13:
                collids=[int(line[5]), int(line[7])]
            elif len(line)==16:
                collids=[int(line[5]), int(line[7]), int(line[9])]
            elif len(line)==19:
                collids=[int(line[5]), int(line[7]), int(line[9]),int(line[11])]

            if int(bin_id0) in collids and int(bin_id1) in collids:
                check_disrupt='dyncollision'
                time_disrupt = float(line[0])
                type_disrupt = line[1]
                break

    if check_disrupt==-100:
        with open(escfile, 'r') as fesc:
            next(fesc)
            for line in fesc:
                data=line.split()
                if int(data[14])==1:
                    if (int(data[17])==bin_id0 and int(data[18])==bin_id1) or (int(data[18])==bin_id0 and int(data[17])==bin_id1):
                        check_disrupt='ejected'
                        time_disrupt = float(data[1])
                        type_disrupt = '-'


    if check_disrupt==-100:
        check_disrupt='dynexchange'
        time_disrupt = '-'
        type_disrupt = '-'

    print(check_disrupt)
    return check_disrupt, time_disrupt, type_disrupt


##Find how a binary is formed
def find_binary_form(bin_id0, bin_id1, modelpath):
    filestr=modelpath+'initial'

    check_form=-100

    snaps=dyn.get_snapshots(filestr)
    firstsnap=snaps[0]
    with gzip.open(firstsnap, 'r') as fsnap:
        next(fsnap); next(fsnap)
        for line in fsnap:
            data=line.split()
            if int(data[7])==1:
                priids=[int(data[10]), int(data[11])]
                #print(priids)
                if int(bin_id0) in priids and int(bin_id1) in priids:
                    check_form='primordial'
                    break

    if check_form==-100:
        tc_id0, tc_id1, tc_m0, tc_m1, tc_k0, tc_k1, tc_a, tc_e, tc_t=tc.find_tc_properties(modelpath)
        #print(tc_id0, tc_id1)
        if (int(bin_id0) in tc_id0 and int(bin_id1) in tc_id1) or (int(bin_id1) in tc_id0 and int(bin_id0) in tc_id1):
            check_form='tidalcapture'

    if check_form==-100:
        check_form='dynexchange'

    print(check_form)
    return check_form


##Find binary encounters between two timesteps
##Timeinterval is a list with two values
def find_binary_encounter(bin_id0, bin_id1, modelpath, timeinterval):
    filestr=modelpath+'initial'
    binintfile=filestr+'.binint.log'
    
    t_conv=dyn.conv('t', filestr+'.conv.sh')

    binint=scripts3.read_binint(binintfile)

    #noenc=0
    anyenc='no'
    for i in range(len(binint)):
        if timeinterval[0]<=t_conv*float(binint[i]['type']['time'])<=timeinterval[1]:
            bininput=binint[i]['input']
            for j in range(len(bininput)):
                if int(bininput[j]['no'])==2:
                    if (int(bininput[j]['ids'][0])==bin_id0 and int(bininput[j]['ids'][1])==bin_id1) or (int(bininput[j]['ids'][0])==bin_id1 and int(bininput[j]['ids'][1])==bin_id0):
                        #noenc+=1
                        anyenc='yes'

            if anyenc==1: break

    print(anyenc)
    return anyenc


##Find binary encounters between two timesteps for a MSP
##Timeinterval is a list with two values
def find_binary_encounter_id(theid, binintpath, timeinterval, tconv):
    anyenc=0

    f=open(binintpath,'r')
    #print 'inputfile %s' %(filename)
    positions=scripts3.find_positions(f)
    for i in range(len(positions)-1):
        #print positions[i]
        binint=scripts3.read_segment(f,positions[i])
        if timeinterval[0]<=tconv*float(binint['type']['time'])<=timeinterval[1]:
            bininput = binint['input']
            binoutput = binint['output']
            #print('yes')
            for j in range(len(bininput)):
                if int(bininput[j]['no'])==1 and int(bininput[j]['ids'][0])==theid:
                    anyenc = 1
                    print('yes1')
                    #for k in range(len(binoutput)):
                    #    if int(binoutput[k]['no'])==2 and (int(binoutput[k]['ids'][0])==theid or int(binoutput[k]['ids'][1])==theid):
                    #        anyenc = 3

                if int(bininput[j]['no'])==2 and (int(bininput[j]['ids'][0])==theid or int(bininput[j]['ids'][1])==theid):
                    anyenc = 2
                    print('yes2')
                    #for k in range(len(binoutput)):
                    #    if int(binoutput[k]['no'])==2 and (int(binoutput[k]['ids'][0])==theid or int(binoutput[k]['ids'][1])==theid):
                    #        anyenc = 4

            if anyenc!=0: break

        elif tconv*float(binint['type']['time'])>timeinterval[1]:
            break
    print('finish binint search')

    f.close()

    return anyenc


##Find binary encounters throughout the cluster evolution
##theids is a list with two values
##Unfinished
#def find_binary_encounter_id(theids, binintpath, tconv):
#    anyenc=0
#    timeall = []; time = []; k0 = []; k1 = []; m0 = []; m1 = []; sma = []; ecc = []; id0 = []; id1 = []
#
#    f=open(binintpath,'r')
#    #print 'inputfile %s' %(filename)
#    positions=scripts3.find_positions(f)
#    for ii in range(len(positions)-1):
#        #print positions[i]
#        binint=scripts3.read_segment(f,positions[ii])
#        #if timeinterval[0]<=tconv*float(binint['type']['time'])<=timeinterval[1]:
#        bininput = binint['input']
#        binoutput = binint['output']
#            #print('yes')
#        for j in range(len(bininput)):
#            if int(bininput[j]['no'])==1 and int(bininput[j]['ids'][0])==theid:
#                anyenc = 1
#                timeall.append(tconv*float(binint['type']['time']))
#                for k in range(len(binoutput)):
#                    if int(binoutput[k]['no'])==1 and (int(binoutput[k]['ids'][0])==theid:
#                        k0.append(int(binoutput[k]['startype'][0])); k1.append(-100)
#                        m0.append(float(binoutput[k]['m'][0])); m1.append(-100)
#                        id0.append(int(binoutput[k]['ids'][0])); ids.append(-100)
#                        sma.append(-100); ecc.append(-100)
#                        timeall.append(tconv*float(binint['type']['time']))
#
#
#                    if int(binoutput[k]['no'])==2 and (int(binoutput[k]['ids'][0])==theid:
#
#
#                    if int(binoutput[k]['no'])==2 and (int(binoutput[k]['ids'][1])==theid:
#
#
#
#            if int(bininput[j]['no'])==2 and (int(bininput[j]['ids'][0])==theid or int(bininput[j]['ids'][1])==theid):
#                anyenc = 2
#                    for k in range(len(binoutput)):
#                        if int(binoutput[k]['no'])==2 and (int(binoutput[k]['ids'][0])==theid or int(binoutput[k]['ids'][1])==theid):


#    print('finish binint search')

#    f.close()

#    return anyenc


##Find MSPs at the time of their formation
def find_msp_atbirth(pathlist):
    sourcedir = np.genfromtxt(pathlist, dtype=str)
    paths = sourcedir[:,0]; status = sourcedir[:,1]

    #f = open('/projects/b1095/syr904/projects/GCE/catalog/msp_catalog_atbirth.dat', 'a+')
    #f.write('#1.Model 2.ID0 3.ID1 4.M0 5.M1 6.K0 7.K1 8.B(G) 9.Ps(sec) 10.SMA(AU) 11.ECC\n')

    ID0 = []; ID1 = []; SPIN = []; MFIELD = []; MODEL = []; TMYR = []; ST = []
    M0 = []; M1 = []; SMA = []; ECC = []; K0 = []; K1 = []
    for ii in range(len(paths)):
        s=paths[ii].split('/')
        n_star=float(s[-2])
        z=float(s[-3][1:])
        rg=int(s[-4][2:])
        rv=float(s[-5][2:])

        filestr = paths[ii]+'initial'

        id_lim = int(n_star*1.05)

        t_conv = ns.conv('t', filestr+'.conv.sh')

        id0 = []; id1 = []; spin = []; mfield = []; model = []; tmyr = []; st = []
        m0 = []; m1 = []; sma = []; ecc = []; k0 = []; k1 = []
        with open(filestr+'.morepulsars.dat', 'r') as fpsr:
            next(fpsr)
            for line in fpsr:
                data = line.split()
                if int(data[2])!=1:
                    Pspin=float(data[9])  ##in sec
                    B=float(data[7])
                    deathcut=(Pspin**2)*(0.17*10**12)
                    if deathcut<B and Pspin<=0.03 and int(data[3]) not in id0:
                        id0.append(int(data[3])); id1.append(-100)
                        m0.append(float(data[5])); m1.append(-100)
                        k0.append(int(data[11])); k1.append(-100)
                        spin.append(float(data[9])); mfield.append(float(data[7]))
                        sma.append(-100); ecc.append(-100)
                        model.append(ii); tmyr.append(float(data[1])*t_conv)
                        st.append(int(status[ii]))

                else:
                    if int(data[11]) == 13:
                        Pspin=float(data[9])  ##in sec
                        B=float(data[7])
                        deathcut=(Pspin**2)*(0.17*10**12)
                        if deathcut<B and Pspin<=0.03 and int(data[3]) not in id0:
                            id0.append(int(data[3])); id1.append(int(data[4]))
                            m0.append(float(data[5])); m1.append(float(data[6]))
                            k0.append(int(data[11])); k1.append(int(data[12]))
                            spin.append(float(data[9])); mfield.append(float(data[7]))
                            sma.append(float(data[13])); ecc.append(float(data[14]))
                            model.append(ii); tmyr.append(float(data[1])*t_conv)
                            st.append(int(status[ii]))

                    if int(data[12]) == 13:
                        Pspin=float(data[10])  ##in sec
                        B=float(data[8])
                        deathcut=(Pspin**2)*(0.17*10**12)
                        if deathcut<B and Pspin<=0.03 and int(data[4]) not in id0:
                            id0.append(int(data[4])); id1.append(int(data[3]))
                            m0.append(float(data[6])); m1.append(float(data[5]))
                            k0.append(int(data[12])); k1.append(int(data[11]))
                            spin.append(float(data[10])); mfield.append(float(data[8]))
                            sma.append(float(data[13])); ecc.append(float(data[14]))
                            model.append(ii); tmyr.append(float(data[1])*t_conv)
                            st.append(int(status[ii]))


        ID0 = ID0+id0; ID1 = ID1+id1; SPIN = SPIN+spin; MFIELD = MFIELD+mfield
        MODEL = MODEL+model; TMYR = TMYR+tmyr; ST = ST+st
        M0 = M0+m0; M1 = M1+m1; SMA = SMA+sma; ECC = ECC+ecc; K0 = K0+k0; K1 = K1+k1
        print(ii)

    np.savetxt('/projects/b1095/syr904/projects/GCE/catalog/msp_catalog_atbirth.dat', np.c_[MODEL, TMYR, ID0, ID1, M0, M1, K0, K1, MFIELD, SPIN, SMA, ECC, ST], fmt = '%d %f %d %d %f %f %d %d %e %f %f %f %d', header = '#1.Model 2.Time(Myr) 3.ID0 4.ID1 5.M0 6.M1 7.K0 8.K1 9.B(G) 10.Ps(sec) 11.SMA(AU) 12.ECC 13.Model_status')


##Find MSPs at the time of their formation for one simulation and check if they are primordial 
##or if they are from tidal capture
def find_msp_atbirth_single(modelpath):
    filestr = modelpath+'initial'
    t_conv = ns.conv('t', filestr+'.conv.sh')

    property_init, property_finl, property_des = ntc_hdf5.find_tc_properties(modelpath)
    ID0_tc = property_init['id0']; ID1_tc = property_init['id1']; T_tc=property_init['time']
    Types_tc = property_init['type']

    snap0 = cmct.Snapshot(fname=modelpath+'initial.snapshots.h5', 
                                snapshot_name='/0(t=0)', conv=modelpath+'initial.conv.sh', 
                                dist=4.52, # distance to cluster in kpc
                                z=0.0038)
    binflag = np.array(snap0.data['binflag'])
    id0_bin = np.array(snap0.data['id0'])[binflag==1]
    id1_bin = np.array(snap0.data['id1'])[binflag==1]


    id0 = []; id1 = []; spin = []; mfield = []; model = []; tmyr = []
    m0 = []; m1 = []; sma = []; ecc = []; k0 = []; k1 = []; primordial = []
    tcflag = []
    with open(filestr+'.morepulsars.dat', 'r') as fpsr:
        next(fpsr)
        for line in fpsr:
            data = line.split()
            if int(data[2])!=1:
                Pspin=float(data[9])  ##in sec
                B=float(data[7])
                deathcut=(Pspin**2)*(0.17*10**12)
                if deathcut<B and Pspin<=0.03 and int(data[3]) not in id0:
                    id0.append(int(data[3])); id1.append(-100)
                    m0.append(float(data[5])); m1.append(-100)
                    k0.append(int(data[11])); k1.append(-100)
                    spin.append(float(data[9])); mfield.append(float(data[7]))
                    sma.append(-100); ecc.append(-100)
                    tmyr.append(float(data[1])*t_conv)
                    primordial.append(-100)
                    tcflag.append(-100)

            else:
                if int(data[11]) == 13:
                    Pspin=float(data[9])  ##in sec
                    B=float(data[7])
                    deathcut=(Pspin**2)*(0.17*10**12)
                    if deathcut<B and Pspin<=0.03 and int(data[3]) not in id0:
                        id0.append(int(data[3])); id1.append(int(data[4]))
                        m0.append(float(data[5])); m1.append(float(data[6]))
                        k0.append(int(data[11])); k1.append(int(data[12]))
                        spin.append(float(data[9])); mfield.append(float(data[7]))
                        sma.append(float(data[13])); ecc.append(float(data[14]))
                        tmyr.append(float(data[1])*t_conv)
                        primordial.append(0)
                        tcflag.append(0)

                        for xx in range(len(id0_bin)):
                            if int(data[3])==int(id0_bin[xx]) and int(data[4])==int(id1_bin[xx]):
                                primordial[-1] = 1
                            elif int(data[3])==int(id1_bin[xx]) and int(data[4])==int(id0_bin[xx]):
                                primordial[-1] = 1

                        for yy in range(len(ID0_tc)):
                            if int(data[3])==int(ID0_tc[yy]) and int(data[4])==int(ID1_tc[yy]):
                                if Types_tc[yy]=='SS_COLL_Giant':
                                    tcflag[-1] = 1
                                if Types_tc[yy]=='SS_COLL_TC_P':
                                    tcflag[-1] = 2
                            elif int(data[3])==int(ID1_tc[yy]) and int(data[4])==int(ID0_tc[yy]):
                                if Types_tc[yy]=='SS_COLL_Giant':
                                    tcflag[-1] = 1
                                if Types_tc[yy]=='SS_COLL_TC_P':
                                    tcflag[-1] = 2


                if int(data[12]) == 13:
                    Pspin=float(data[10])  ##in sec
                    B=float(data[8])
                    deathcut=(Pspin**2)*(0.17*10**12)
                    if deathcut<B and Pspin<=0.03 and int(data[4]) not in id0:
                        id0.append(int(data[4])); id1.append(int(data[3]))
                        m0.append(float(data[6])); m1.append(float(data[5]))
                        k0.append(int(data[12])); k1.append(int(data[11]))
                        spin.append(float(data[10])); mfield.append(float(data[8]))
                        sma.append(float(data[13])); ecc.append(float(data[14]))
                        tmyr.append(float(data[1])*t_conv)
                        primordial.append(0)
                        tcflag.append(0)

                        for xx in range(len(id0_bin)):
                            if int(data[3])==int(id0_bin[xx]) and int(data[4])==int(id1_bin[xx]):
                                primordial[-1] = 1
                            elif int(data[3])==int(id1_bin[xx]) and int(data[4])==int(id0_bin[xx]):
                                primordial[-1] = 1

                        for yy in range(len(ID0_tc)):
                            if int(data[3])==int(ID0_tc[yy]) and int(data[4])==int(ID1_tc[yy]):
                                if Types_tc[yy]=='SS_COLL_Giant':
                                    tcflag[-1] = 1
                                if Types_tc[yy]=='SS_COLL_TC_P':
                                    tcflag[-1] = 2
                            elif int(data[3])==int(ID1_tc[yy]) and int(data[4])==int(ID0_tc[yy]):
                                if Types_tc[yy]=='SS_COLL_Giant':
                                    tcflag[-1] = 1
                                if Types_tc[yy]=='SS_COLL_TC_P':
                                    tcflag[-1] = 2


    np.savetxt(modelpath+'allmsp_atbirth.dat', np.c_[tmyr, id0, id1, m0, m1, k0, k1, mfield, spin, sma, ecc, primordial, tcflag], fmt = '%f %d %d %f %f %d %d %e %f %f %f %d %d', header = '#1.Time(Myr) 2.ID0 3.ID1 4.M0 5.M1 6.K0 7.K1 8.B(G) 9.Ps(sec) 10.SMA(AU) 11.ECC 12.Primordial_flag 13.TC_flag')


##Find if MSPs are from primordial binaries or from dynamical interactions
def find_msp_pridyn(pathlist, msp_birth_file):
    sourcedir = np.genfromtxt(pathlist, dtype = str)
    paths = sourcedir[:,0]; status = sourcedir[:,1]

    msp_birth = np.genfromtxt(msp_birth_file)
    models = msp_birth[:,0]; time = msp_birth[:,1]; id0 = msp_birth[:,2]; id1 = msp_birth[:,3]

    allkey = list(Counter(models).keys())
    #print(allkey)

    priflag = []
    for ii in range(len(allkey)):
        mno = int(allkey[ii])
        print(mno)
        filestr = paths[mno]+'initial'

        firstsnap = filestr+'.snap0000.dat.gz'
        psrfile = filestr+'.morepulsars.dat'
        #semerge = filestr+'.semergedisrupt.log'
        binint = filestr+'.binint.log'
        #datasemer = scripts2.readmergefile(semerge)
        #databinint = scripts3.read_binint(binint)
        #print('read binint')

        t_conv = ns.conv('t', filestr+'.conv.sh')

        bid0 = []; bid1 = []
        with gzip.open(firstsnap, 'r') as fsnap:
            next(fsnap); next(fsnap)
            for line in fsnap:
                datasnap = line.split()
                if int(datasnap[7])==1:
                    bid0.append(int(datasnap[10])); bid1.append(int(datasnap[11]))

        for xx in range(len(models)):
            if int(models[xx])==mno:
                check = 0
                priflag.append(0)
                if id1[xx]!=-100:
                    for h in range(len(bid0)):
                        if (id0[xx]==bid0[h] and id1[xx]==bid1[h]) or (id0[xx]==bid1[h] and id1[xx]==bid0[h]):
                            check = 1
                            encs = find_binary_encounter_id(id0[xx], binint, [0., time[xx]], t_conv)
                            if encs == 0:
                                priflag[xx]=1
                            else:
                                priflag[xx]=2
                            break

                        elif id0[xx]==bid0[h] or id0[xx]==bid1[h]:
                            check = 2
                            encs = find_binary_encounter_id(id0[xx], binint, [0., time[xx]], t_conv)
                            if encs == 0:
                                priflag[xx]=1
                                break
                            else:
                                with open(psrfile, 'r') as fpsr:
                                    next(fpsr)
                                    for line in fpsr:
                                        datapsr = line.split()
                                        if int(datapsr[3])==id0[xx] and int(datapsr[4])==-100 and float(datapsr[9])>0.03:
                                            priflag[xx]=3
                                            break


                    if check != 1 and check != 2: 
                        priflag[xx]=3


                else:
                    for h in range(len(bid0)):
                        if id0[xx]==bid0[h] or id0[xx]==bid1[h]:
                            check = 2
                            encs = find_binary_encounter_id(id0[xx], binint, [0., time[xx]], t_conv)
                            if encs == 0:
                                priflag[xx]=1
                                break
                            else:
                                with open(psrfile, 'r') as fpsr:
                                    next(fpsr)
                                    for line in fpsr:
                                        datapsr = line.split()
                                        if int(datapsr[3])==id0[xx] and int(datapsr[4])==-100 and float(datapsr[9])>0.03:
                                            priflag[xx]=3
                                            break

                    if check != 2: 
                        priflag[xx]=3


                print(xx, priflag[xx])

    np.savetxt('/projects/b1095/syr904/projects/GCE/catalog/msp_pridyn.dat', np.c_[models, time, id0, id1, priflag], fmt = '%d %f %d %d %d', header = '#1.Model 2.Time(Myr) 3.ID0 4.ID1 5.Priflag (0:undeci, 1:pri 2:pri+dyn 3:dyn)', comments = '#')




##Separate primordial and dynamical msps at all times for multiple models
def separate_msp_primordial(pathlist, msp_pridyn_file):
    sourcedir = np.genfromtxt(pathlist, dtype = str)
    paths = sourcedir[:,0]; status = sourcedir[:,1]

    data = np.genfromtxt(msp_pridyn_file)
    model = data[:,0]; id0 = data[:,2]; id1 = data[:,3]; priflag = data[:,4]

    for ii in range(len(paths)):
        filestr = paths[ii]+'initial'
        firstsnap = filestr+'.snap0000.dat.gz'
        psrfile = filestr+'.morepulsars.dat'
        #semerge = filestr+'.semergedisrupt.log'
        #datasemer = scripts2.readmergefile(semerge)

        t_conv = dyn.conv('t', filestr+'.conv.sh')

        #bid0 = []; bid1 = []
        #with gzip.open(firstsnap, 'r') as fsnap:
        #    next(fsnap); next(fsnap)
        #    for line in fsnap:
        #        datasnap = line.split()
        #        if int(datasnap[7])==1:
        #            bid0.append(int(datasnap[10])); bid1.append(int(datasnap[11]))
    
        ntot = [0]
        nsin = [0]; nspri = [0]; nsdyn = [0]
        nbin = [0]; nbpri = [0]; nbdyn = [0]
        xx = 0
        time = []
        with open(filestr+'.morepulsars.dat', 'r') as fpsr:
            next(fpsr)
            line0 = fpsr.readline().strip()
            if not line0:
                print(ii, 'empty')
                continue
            data0 = line0.split()
            t0 = float(data0[1])
            #print(t0)
            time.append(t0*t_conv)
            for line in fpsr:
                data = line.split()
                if float(data[1])==t0:
                    if int(data[2])!=1:
                        Pspin=float(data[9])  ##in sec
                        B=float(data[7])
                        deathcut=(Pspin**2)*(0.17*10**12)
                        if deathcut<B and Pspin<=0.03:
                            ntot[xx]+=1; nsin[xx]+=1
                            check=0
                            for jj in range(len(id0)):
                                if int(model[jj])==ii and int(data[3])==id0[jj]:
                                    check=1
                                    if priflag[jj]==1.:
                                        nspri[xx]+=1
                                    else:
                                        nsdyn[xx]+=1
                                    break

                            if check==0:
                                print(ii, int(data[3]))


                    else:
                        if int(data[11]) == 13:
                            Pspin=float(data[9])  ##in sec
                            B=float(data[7])
                            deathcut=(Pspin**2)*(0.17*10**12)
                            if deathcut<B and Pspin<=0.03:
                                ntot[xx]+=1; nbin[xx]+=1
                                check = 0
                                for jj in range(len(id0)):
                                    if int(model[jj])==ii and int(data[3])==id0[jj]:
                                        check=1
                                        if priflag[jj]==1.:
                                            nbpri[xx]+=1
                                        else:
                                            nbdyn[xx]+=1
                                        break

                                if check==0:
                                    print(ii, int(data[3]))
                                

                        if int(data[12]) == 13:
                            Pspin=float(data[10])  ##in sec
                            B=float(data[8])
                            deathcut=(Pspin**2)*(0.17*10**12)
                            if deathcut<B and Pspin<=0.03:
                                ntot[xx]+=1; nbin[xx]+=1
                                check = 0
                                for jj in range(len(id0)):
                                    if int(model[jj])==ii and int(data[4])==id0[jj]:
                                        check=1
                                        if priflag[jj]==1.:
                                            nbpri[xx]+=1
                                        else:
                                            nbdyn[xx]+=1
                                        break

                                if check==0:
                                    print(ii, int(data[4]))


                else:
                    ntot.append(0)
                    nsin.append(0); nspri.append(0); nsdyn.append(0) 
                    nbin.append(0); nbpri.append(0); nbdyn.append(0)
                    xx+=1
                    t0 = float(data[1])
                    time.append(t0*t_conv)
                    if int(data[2])!=1:
                        Pspin=float(data[9])  ##in sec
                        B=float(data[7])
                        deathcut=(Pspin**2)*(0.17*10**12)
                        if deathcut<B and Pspin<=0.03:
                            ntot[xx]+=1; nsin[xx]+=1
                            check = 0
                            for jj in range(len(id0)):
                                if int(model[jj])==ii and int(data[3])==id0[jj]:
                                    check=1
                                    if priflag[jj]==1.:
                                        nspri[xx]+=1
                                    else:
                                        nsdyn[xx]+=1
                                    break

                            if check==0:
                                print(ii, int(data[3]))


                    else:
                        if int(data[11]) == 13:
                            Pspin=float(data[9])  ##in sec
                            B=float(data[7])
                            deathcut=(Pspin**2)*(0.17*10**12)
                            if deathcut<B and Pspin<=0.03:
                                ntot[xx]+=1; nbin[xx]+=1
                                check = 0
                                for jj in range(len(id0)):
                                    if int(model[jj])==ii and int(data[3])==id0[jj]:
                                        check=1
                                        if priflag[jj]==1.:
                                            nbpri[xx]+=1
                                        else:
                                            nbdyn[xx]+=1
                                        break

                                if check==0:
                                    print(ii, int(data[3]))

                        if int(data[12]) == 13:
                            Pspin=float(data[10])  ##in sec
                            B=float(data[8])
                            deathcut=(Pspin**2)*(0.17*10**12)
                            if deathcut<B and Pspin<=0.03:
                                ntot[xx]+=1; nbin[xx]+=1
                                check = 0
                                for jj in range(len(id0)):
                                    if int(model[jj])==ii and int(data[4])==id0[jj]:
                                        check=1
                                        if priflag[jj]==1.:
                                            nspri[xx]+=1
                                        else:
                                            nsdyn[xx]+=1
                                        break

                                if check==0:
                                    print(ii, int(data[4]))


        np.savetxt('/projects/b1095/syr904/projects/GCE/catalog/msp_primordial_dynamical_catalog/ns_pridyn_'+str(ii)+'.dat', np.c_[time, ntot, nsin, nspri, nsdyn, nbin, nbpri, nbdyn], fmt = '%f %d %d %d %d %d %d %d', header = '1.Time(Myr) 2.Nmsp_tot 3.Nmsp_sin 4.Nmsp_sin_pri 5.Nmsp_sin_dyn 6.Nmsp_bin 7.Nmsp_bin_pri 8.Nmsp_bin_dyn', comments='#')
        print(ii)









