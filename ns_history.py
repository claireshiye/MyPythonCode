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
import ns_tidalcapture as ntc
#from scipy import stats

yearsc=31557600.
twopi=6.283185307179586
Gconst=6.674*10**-8 ##cm3*g-1*s-2
clight=3*10**10 ##cm/s
Msun=2*10**33 ##gram
AU=1.496*10**13  ##cm
PC=3.086*10**18  ##cm
Kconst=9.87*10**-48 ##yr/G^2
Lsun=4.02*10**16 ##mJy*kpc^2


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













