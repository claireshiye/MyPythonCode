import numpy as np 
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
def find_tc_properties_final(filepath):
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



##Types of tidal captures in the tidal capture file
def find_tc_properties(filepath):
    filestr=filepath+'initial'
    tcfile=filestr+'.tidalcapture.log'
    t=[]; types=[]; id0i=[]; id1i=[]; m0i=[]; m1i=[]; k0i=[]; k1i=[]; r0i=[]; r1i=[]; rperi=[]
    a=[]; e=[]; id0f=[]; id1f=[]; m0f=[]; m1f=[]; k0f=[]; k1f=[]; r0f=[]; r1f=[]
    v_inf = []; rcm = []
    mc0=[]; mc1=[]; rc0=[]; rc1=[]

    t_des=[]; type_des=[]; rperi_des=[]; idf_des=[]; mf_des=[]; kf_des=[]
    id0_des=[]; id1_des=[]; m0_des=[]; m1_des=[]; k0_des=[]; k1_des=[]; r0_des=[]; r1_des=[]
    v_inf_des=[]; rcm_des = []
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

                #if len(numstr)>19: 
                #    expo=float(numstr.pop(14))
                #    numstr[13]=str(float(numstr[13])**(-expo))
                numstr[9] = numstr[9][1:]

                t_des.append(float(data[0])); type_des.append(data[1])
                id0_des.append(int(numstr[0])); id1_des.append(int(numstr[3]))
                m0_des.append(float(numstr[1])); m1_des.append(float(numstr[4]))
                k0_des.append(int(numstr[2])); k1_des.append(int(numstr[5]))
                r0_des.append(float(numstr[6])); r1_des.append(float(numstr[7]))
                rperi_des.append(float(numstr[8])); v_inf_des.append(float(numstr[9]))
                #mc0_des = list(np.full_like(id0_des, -100)); mc1_des = mc0_des
                #rc0_des = mc0_des; rc1_des = mc0_des
                continue

            if data[2][-1] == 'd' and data[1] == 'SS_COLL_Giant': 
                n_giant_coll+=1
                numstr=re.split(',|\(|\)|->|\[|\]', data[2])
                numstr=list(filter(None, numstr))
                numstr=list(filter(lambda x: x!='+', numstr))

                #if len(numstr)>19: 
                #    expo=float(numstr.pop(14))
                #    numstr[13]=str(float(numstr[13])**(-expo))
                numstr[9] = numstr[9][1:]
                #print(numstr[9])

                t_des.append(float(data[0])); type_des.append(data[1])
                id0_des.append(int(numstr[0])); id1_des.append(int(numstr[3]))
                m0_des.append(float(numstr[1])); m1_des.append(float(numstr[4]))
                k0_des.append(int(numstr[2])); k1_des.append(int(numstr[5]))
                r0_des.append(float(numstr[6])); r1_des.append(float(numstr[7]))
                rperi_des.append(float(numstr[8])); v_inf_des.append(float(numstr[9]))
                mc0_des.append(float(numstr[10])); mc1_des.append(float(numstr[11]))
                rc0_des.append(float(numstr[12])); rc1_des.append(float(numstr[13]))

                continue


            numstr=re.split(',|\(|\)|->|\[|\]', data[2])
            numstr=list(filter(None, numstr))
            numstr=list(filter(lambda x: x!='+', numstr))

            #if len(numstr)>19: 
            #    expo=float(numstr.pop(14))
            #    numstr[13]=str(float(numstr[13])**(-expo))
            numstr[9] = numstr[9][1:]
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
            v_inf.append(float(numstr[9]))

            ##Final properties
            id0f.append(int(numstr[10])); id1f.append(int(numstr[15]))
            m0f.append(float(numstr[11])); m1f.append(float(numstr[16]))
            k0f.append(int(numstr[12])); k1f.append(int(numstr[17]))
            a.append(float(numstr[13])); e.append(float(numstr[14]))
            r0f.append(float(numstr[18])); r1f.append(float(numstr[19]))

    Prop_init = {'time':t, 'type':types, 'id0': id0i, 'id1': id1i, 'm0': m0i, 'm1': m1i, 'k0': k0i, 'k1': k1i, 'r0': r0i, 'r1': r1i, 'rperi': rperi, 'vinf': v_inf}
    Prop_finl = {'id0': id0f, 'id1': id1f, 'm0': m0f, 'm1': m1f, 'k0': k0f, 'k1': k1f, 'r0': r0f, 'r1': r1f, 'sma': a, 'ecc': e}
    Prop_des = {'time':t_des, 'type':type_des, 'id0': id0_des, 'id1': id1_des, 'm0': m0_des, 'm1': m1_des, 'k0': k0_des, 'k1': k1_des, 'r0': r0_des, 'r1': r1_des, 'rperi': rperi_des, 'vinf': v_inf_des, 'mc0': mc0_des, 'mc1': mc1_des, 'rc0': rc0_des, 'rc1': rc1_des}

    print('n_giant_coll:', n_giant_coll, 'n_tc_merged:', n_tc_merged)

    return Prop_init, Prop_finl, Prop_des


##Find how many tidal capture binaries have become NS binaries
def find_tc_ns(filepath, savepath):
    property_init, property_finl, property_des = find_tc_properties(filepath)
    ID0 = property_init['id0']; ID1 = property_init['id1']; T=property_init['time']
    filestr=filepath+'initial'
    psrfile=filestr+'.morepulsars.dat'

    f1=open(savepath+'tc_binarypsr.dat', 'a+')##The tc binary is intact and one of the stars becomes a NS
    f2=open(savepath+'tc_singlepsr.dat', 'a+')##The tc binary is not intact but one of the stars becomes a NS
    f1.write('#1.Time 2.ID0 3.ID1 4.M0 5.M1 6.K0 7.K1 8.A(AU) 9.ECC 10.B0 11.B1 12.P0 13.P1\n')
    f2.write('#1.Time 2.ID0 3.ID1 4.M0 5.M1 6.K0 7.K1 8.A(AU) 9.ECC 10.B0 11.B1 12.P0 13.P1\n')
    for i in range(len(ID0)):
        check0=0; check1=0
        with open(psrfile, 'r') as fpsr:
            next(fpsr)
            for line in fpsr:
                datapsr=line.split()
                if float(datapsr[1])>T[i]:
                    if (int(datapsr[3])==ID0[i] and int(datapsr[4])==ID1[i]) or (int(datapsr[3])==ID1[i] and int(datapsr[4])==ID0[i]):
                        f1.write('%f %d %d %f %f %d %d %f %f %e %e %f %f\n'%(float(datapsr[1]), int(datapsr[3]), int(datapsr[4]), float(datapsr[5]), float(datapsr[6]), int(datapsr[11]), int(datapsr[12]), float(datapsr[13]), float(datapsr[14]), float(datapsr[7]), float(datapsr[8]), float(datapsr[9]), float(datapsr[10])))
                        break

                    elif check0==0 and (int(datapsr[3])==ID0[i] or int(datapsr[4])==ID0[i]):
                        f2.write('%f %d %d %f %f %d %d %f %f %e %e %f %f\n'%(float(datapsr[1]), int(datapsr[3]), int(datapsr[4]), float(datapsr[5]), float(datapsr[6]), int(datapsr[11]), int(datapsr[12]), float(datapsr[13]), float(datapsr[14]), float(datapsr[7]), float(datapsr[8]), float(datapsr[9]), float(datapsr[10])))
                        check0=1


                    elif check1==0 and (int(datapsr[3])==ID1[i] and int(datapsr[4])==ID1[i]):
                        f2.write('%f %d %d %f %f %d %d %f %f %e %e %f %f\n'%(float(datapsr[1]), int(datapsr[3]), int(datapsr[4]), float(datapsr[5]), float(datapsr[6]), int(datapsr[11]), int(datapsr[12]), float(datapsr[13]), float(datapsr[14]), float(datapsr[7]), float(datapsr[8]), float(datapsr[9]), float(datapsr[10])))
                        check1=0

                if check0==1 and check1==1: break

        print(i)

    f1.close()
    f2.close()


##Find tc binaries or singles that are still pulsars at the end
def find_tc_psr():
    data_tcbin=np.genfromtxt('/projects/b1095/syr904/projects/PULSAR2/tc_comparison/tc_binarypsr.dat')
    data_tcsin=np.genfromtxt('/projects/b1095/syr904/projects/PULSAR2/tc_comparison/tc_singlepsr.dat')
    data_msp=np.genfromtxt('/projects/b1095/syr904/projects/PULSAR2/tc_comparison/MSP_last.dat')
    data_psr=np.genfromtxt('/projects/b1095/syr904/projects/PULSAR2/tc_comparison/PSR_last.dat')

    id_tcbin=[]; id_tcsin=[]; id_msp=[]; id_psr=[]

    id0_tc=np.concatenate((data_tcbin[:,1],data_tcsin[:,1]))
    id1_tc=np.concatenate((data_tcbin[:,2],data_tcsin[:,2]))
    id_psrs=np.concatenate((data_msp[:,12],data_psr[:,12]))

    print(np.intersect1d(id0_tc, id_psrs))
    print(np.intersect1d(id1_tc, id_psrs))



##Find NS-XX star binaries at the last snapshot and check if they are formed in tidal capture/giant collision
def find_NS_XX_last(filepath, savepath, lowlim, highlim, loopnum):
    property_init, property_finl, property_des = find_tc_properties_final(filepath)
    ID0 = property_init['id0']; ID1 = property_init['id1']; T=property_init['time']
    Types = property_init['type']

    filestr=filepath+'initial'
    snaps=dyn.get_snapshots(filestr)
    lastsnap=snaps[-1]
    t_conv=dyn.conv('t', filestr+'.conv.sh')
    time=dyn.get_time(lastsnap)*t_conv
    model=loopnum

    #os.system('rm '+savepath)
    fmsb=open(savepath, 'a+')
    fmsb.write('#1.Model 2.Time 3.ID0 4.ID1 5.M0 6.M1 7.K0 8.K1 9.a(AU) 10.ecc 11.radrol0 12.radrol1 13.B(G) 14.P(sec) 15.tcflag\n')

    with gzip.open(lastsnap, 'r') as flast:
        next(flast); next(flast)
        for line in flast:
            datalast=line.split()
            if int(datalast[7])==1:
                if int(datalast[17])==13 and lowlim<=int(datalast[18])<=highlim:
                    ID0ms=int(datalast[10]); ID1ms=int(datalast[11])
                    tcflag=4
                    for ii in range(len(ID0)):    
                        if (ID0ms==ID0[ii] and ID1ms==ID1[ii]) or (ID1ms==ID0[ii] and ID0ms==ID1[ii]):
                            if Types[ii]=='SS_COLL_Giant':
                                tcflag=81
                                break
                            if Types[ii]=='SS_COLL_TC_P':
                                tcflag=91
                                break
                        elif ID0ms==ID0[ii] or ID0ms==ID1[ii]:
                            if Types[ii]=='SS_COLL_Giant':
                                tcflag=82
                                break
                            if Types[ii]=='SS_COLL_TC_P':
                                tcflag=92
                                break

                    fmsb.write('%d %f %d %d %f %f %d %d %f %f %f %f %e %f %d\n'%(model, time, int(datalast[10]), int(datalast[11]), float(datalast[8]), float(datalast[9]), int(datalast[17]), int(datalast[18]), float(datalast[12]), float(datalast[13]), float(datalast[43]), float(datalast[44]), float(datalast[47]), float(twopi*yearsc/float(datalast[45])), tcflag))


                if lowlim<=int(datalast[17])<=highlim and int(datalast[18])==13:
                    ID0ms=int(datalast[11]); ID1ms=int(datalast[10])
                    tcflag=4
                    for ii in range(len(ID0)):    
                        if (ID0ms==ID0[ii] and ID1ms==ID1[ii]) or (ID1ms==ID0[ii] and ID0ms==ID1[ii]):
                            if Types[ii]=='SS_COLL_Giant':
                                tcflag=81
                                break
                            if Types[ii]=='SS_COLL_TC_P':
                                tcflag=91
                                break
                        elif ID0ms==ID0[ii] or ID0ms==ID1[ii]:
                            if Types[ii]=='SS_COLL_Giant':
                                tcflag=82
                                break
                            if Types[ii]=='SS_COLL_TC_P':
                                tcflag=92
                                break
                                
                    fmsb.write('%d %f %d %d %f %f %d %d %f %f %f %f %e %f %d\n'%(model, time, int(datalast[11]), int(datalast[10]), float(datalast[9]), float(datalast[8]), int(datalast[18]), int(datalast[17]), float(datalast[12]), float(datalast[13]), float(datalast[44]), float(datalast[43]), float(datalast[48]), float(twopi*yearsc/float(datalast[46])), tcflag))

    fmsb.close()


##Print out all 13+X binaries from tidal capture and giant collisions
def find_NS_XX_tc_giantcoll(pathlist, savepath):
    #paths = np.genfromtxt(pathlist, dtype=str)
    paths = ['/projects/b1095/syr904/cmc/47Tuc/rundir/47Tuc/best_fits/MOCHA47Tuc_elson_rv4_3e6_tcon/']
    ftc = open(savepath + 'all_tc_NSXX.dat', 'a+')
    fcoll = open(savepath + 'all_coll_NSXX.dat', 'a+')
    ftc.write('#1.Model 2.Time(Myr) 3.ID0 4.ID1 5.K0_init 6.K1_init 7.M0_init 8.M1_init 9.R0_init 10.R1_init 11.R_peri 12.V_INF 13.K0_fnl 14.K1_fnl 15.M0_fnl 16.M1_fnl 17.SMA 18.ECC\n')
    fcoll.write('#1.Model 2.Time(Myr) 3.ID0 4.ID1 5.K0_init 6.K1_init 7.M0_init 8.M1_init 9.R0_init 10.R1_init 11.R_peri 12.V_INF 13.K0_fnl 14.K1_fnl 15.M0_fnl 16.M1_fnl 17.SMA 18.ECC 19.Period(days)\n')

    for ii in range(len(paths)):
        t_conv = dyn.conv('t', paths[ii]+'initial.conv.sh')

        prop_init, prop_finl, prop_des = find_tc_properties_final(paths[ii])

        for xx in range(len(prop_finl['id0'])):
            if prop_finl['k0'][xx]==13 and prop_init['type'][xx] == 'SS_COLL_TC_P':
                period = uc.au_to_period(prop_finl['sma'][xx], prop_finl['m0'][xx], prop_finl['m0'][xx])
                ftc.write('%d %f %d %d %d %d %f %f %f %f %f %f %d %d %f %f %f %f %f\n'%(ii, prop_init['time'][xx]*t_conv, prop_init['id0'][xx], prop_init['id1'][xx], 
                      prop_init['k0'][xx], prop_init['k1'][xx], 
                      prop_init['m0'][xx], prop_init['m1'][xx],
                      prop_init['r0'][xx], prop_init['r1'][xx],
                      prop_init['rperi'][xx], prop_init['vinf'][xx],
                      prop_finl['k0'][xx], prop_finl['k1'][xx], 
                      prop_finl['m0'][xx], prop_finl['m1'][xx],
                      prop_finl['sma'][xx], prop_finl['ecc'][xx],
                      period))
            if prop_finl['k1'][xx]==13 and prop_init['type'][xx] == 'SS_COLL_TC_P':
                period = uc.au_to_period(prop_finl['sma'][xx], prop_finl['m0'][xx], prop_finl['m0'][xx])
                ftc.write('%d %f %d %d %d %d %f %f %f %f %f %f %d %d %f %f %f %f %f\n'%(ii, prop_init['time'][xx]*t_conv, prop_init['id1'][xx], prop_init['id0'][xx], 
                      prop_init['k1'][xx], prop_init['k0'][xx], 
                      prop_init['m1'][xx], prop_init['m0'][xx],
                      prop_init['r1'][xx], prop_init['r0'][xx],
                      prop_init['rperi'][xx], prop_init['vinf'][xx],
                      prop_finl['k1'][xx], prop_finl['k0'][xx], 
                      prop_finl['m1'][xx], prop_finl['m0'][xx],
                      prop_finl['sma'][xx], prop_finl['ecc'][xx],
                      period))

            if prop_finl['k0'][xx]==13 and prop_init['type'][xx] == 'SS_COLL_Giant':
                period = uc.au_to_period(prop_finl['sma'][xx], prop_finl['m0'][xx], prop_finl['m0'][xx])
                fcoll.write('%d %f %d %d %d %d %f %f %f %f %f %f %d %d %f %f %f %f %f\n'%(ii, prop_init['time'][xx]*t_conv, prop_init['id0'][xx], prop_init['id1'][xx], 
                      prop_init['k0'][xx], prop_init['k1'][xx], 
                      prop_init['m0'][xx], prop_init['m1'][xx],
                      prop_init['r0'][xx], prop_init['r1'][xx],
                      prop_init['rperi'][xx], prop_init['vinf'][xx],
                      prop_finl['k0'][xx], prop_finl['k1'][xx], 
                      prop_finl['m0'][xx], prop_finl['m1'][xx],
                      prop_finl['sma'][xx], prop_finl['ecc'][xx],
                      period))
            if prop_finl['k1'][xx]==13 and prop_init['type'][xx] == 'SS_COLL_Giant':
                period = uc.au_to_period(prop_finl['sma'][xx], prop_finl['m0'][xx], prop_finl['m0'][xx])
                fcoll.write('%d %f %d %d %d %d %f %f %f %f %f %f %d %d %f %f %f %f %f\n'%(ii, prop_init['time'][xx]*t_conv, prop_init['id1'][xx], prop_init['id0'][xx], 
                      prop_init['k1'][xx], prop_init['k0'][xx], 
                      prop_init['m1'][xx], prop_init['m0'][xx],
                      prop_init['r1'][xx], prop_init['r0'][xx],
                      prop_init['rperi'][xx], prop_init['vinf'][xx],
                      prop_finl['k1'][xx], prop_finl['k0'][xx], 
                      prop_finl['m1'][xx], prop_finl['m0'][xx],
                      prop_finl['sma'][xx], prop_finl['ecc'][xx],
                      period))

        print(ii)

    ftc.close()
    fcoll.close()

##find what happens to binaries, especially tc and giant collision binaries
##Here inputfile is /projects/b1095/syr904/cmc/cmc-mpi-tidalcapture/rvgrid/standard_models_tcon/all_tc_NSXX.dat etc.
def find_binary_interaction(inputfile, pathlist):
    paths = np.genfromtxt(pathlist, dtype=str)
    data = np.genfromtxt(inputfile)
    models = data[:,0]; id0 = data[:,2]; id1 = data[:,3]

    f = open('/projects/b1095/syr904/cmc/cmc-mpi-tidalcapture/rvgrid/standard_models_tcon/coll_interactions.dat', 'a+')
    print('#1.Model 2.ID0 3.ID1 4.Time(Myr) 5.Status', file = f)

    for ii in range(len(paths)):
        t_conv = dyn.conv('t', paths[ii]+'initial.conv.sh')
        data_binint = scripts3.read_binint(paths[ii]+'initial.binint.log')
        #print('binint read', file = f)

        for xx in range(len(models)):
            if models[xx] == ii:
                id_ns = int(id0[xx]); id_comp = int(id1[xx])
                check = 0

                for yy in range(len(data_binint)-1, -1, -1):
                    bininput = data_binint[yy]['input']
                    binoutput = data_binint[yy]['output']
                    time = data_binint[yy]['type']['time']*t_conv

                    for m in range(len(bininput)):
                        #print(bininput[m]['no'], bininput[m]['ids'])
                        if int(bininput[m]['no'])==2 and ((int(bininput[m]['ids'][0])==id_ns and int(bininput[m]['ids'][1]) == id_comp) or (int(bininput[m]['ids'][0])==id_comp and int(bininput[m]['ids'][1]) == id_ns)):
                            check = 1
                            check_out = 0
                            for n in range(len(binoutput)):
                                if binoutput[n]['merge']['ids']:    
                                    #print(binoutput[n]['merge']['merge'].index(1))
                                    merge_index = int(binoutput[n]['merge']['merge'].index(1))

                                if binoutput[n]['merge']['ids'] and (str(id_ns) in binoutput[n]['merge']['ids'][merge_index] and str(id_comp) in binoutput[n]['merge']['ids'][merge_index]):
                                    print(ii, id_ns, id_comp, time, 'binary_merged', file = f)
                                    check_out = 1
                                    break

                                elif binoutput[n]['merge']['ids'] and (str(id_ns) in binoutput[n]['merge']['ids'][merge_index] or str(id_comp) in binoutput[n]['merge']['ids'][merge_index]):
                                    print(ii, id_ns, id_comp, time, 'one_merged', file = f)
                                    check_out = 1
                                    break

                                elif int(binoutput[n]['no'])==2 and str(id_ns) in binoutput[n]['ids'] and str(id_comp) in binoutput[n]['ids']:
                                    print(ii, id_ns, id_comp, time, binoutput[n], file = f)
                                    check_out = 1
                                    break

                                elif int(binoutput[n]['no'])==3 and str(id_ns) in binoutput[n]['ids'][:2] and str(id_comp) in binoutput[n]['ids'][:2]:
                                    print(ii, id_ns, id_comp, time, binoutput[n], file = f)
                                    check_out = 1
                                    break

                            if check_out == 0:    
                                print(ii, id_ns, id_comp, time, 'disrupted', file = f)
                                break


                    if check: break

                if not check:
                    print(ii, id_ns, id_comp, time, 'no_interactions', file = f)

        print(ii)

    f.close()



##Check if the pulsars at the last snapshot is formed from tidal capture
def check_psr_tc_last(filepath):
    datamsp = np.genfromtxt(filepath+'msp_last.dat')
    #datapsr = np.genfromtxt(filepath+'normalpsr_last.dat')

    datansms = np.genfromtxt(filepath+'NS_MS_last.dat')
    datanswd = np.genfromtxt(filepath+'NS_WD_last.dat')

    msp_id0 = datamsp[:,12]; msp_id1 = datamsp[:,13]
    msp_k1 = datamsp[:,15]

    nsms_id0 = datansms[:,2]; nsms_id1 = datansms[:,3]; nsmsflag = datansms[:,14]
    nswd_id0 = datanswd[:,2]; nswd_id1 = datanswd[:,3]; nswdflag = datanswd[:,14]

    prop_init, prop_finl, prop_des=find_tc_properties_final(filepath)
    tc_id0 = prop_init['id0']; tc_id1 = prop_init['id1']

    msp_tcflag = []
    for i in range(len(msp_id0)):
        if msp_id1[i]!=-100 and msp_k1[i]<=1:
            for j in range(len(nsms_id0)):
                if msp_id0[i]==nsms_id0[j] and msp_id1[i]==nsms_id1[j]:
                    msp_tcflag.append(int(nsmsflag[j]))
                    break
        
        if msp_id1[i]!=-100 and 10<=msp_k1[i]<=12:
            for j in range(len(nswd_id0)):
                if msp_id0[i]==nswd_id0[j] and msp_id1[i]==nswd_id1[j]:
                    msp_tcflag.append(int(nswdflag[j]))
                    break

        if msp_id1[i]!=-100 and msp_k1[i]>=13:
            msp_tcflag.append(3)
            for j in range(len(tc_id0)):
                if (msp_id0[i]==tc_id0[j] and msp_id1[i]==tc_id1[j]) or (msp_id1[i]==tc_id0[j] and msp_id0[i]==tc_id1[j]):
                    msp_tcflag[i]=1
                    break
                elif msp_id0[i]==tc_id0[j] or msp_id0[i]==tc_id1[j]:
                    msp_tcflag[i]=2
                    break

        if msp_id1[i]==-100:
            msp_tcflag.append(3)
            for j in range(len(tc_id0)):
                if msp_id0[i]==tc_id0[j] or msp_id0[i]==tc_id1[j]:
                    msp_tcflag[i]=2
                    break
                
    return msp_tcflag


##Check if the pulsars at the last snapshot is formed from tidal capture in multiple models
def check_psr_tc_models(pathlist):
    paths = np.genfromtxt(pathlist, dtype=str)

    for kk in range(len(paths)):
        print(kk)
        tc_flags = check_psr_tc_last(paths[kk])
        print()
        uf.add_column(paths[kk]+'msp_last.dat', paths[kk]+'msp_last.dat', tc_flags)



##Find NS-MS star binaries at all times and check if they are formed in tidal capture
def find_NSMS_atalltime(filepath, savepath):
    property_init, property_finl, property_des = find_tc_properties(filepath)
    ID0 = property_init['id0']; ID1 = property_init['id1']; T=property_init['time']

    filestr=filepath+'initial'
    t_conv=dyn.conv('t', filestr+'.conv.sh')
    l_conv=dyn.conv('l', filestr+'.conv.sh')
    
    fnsms=open(savepath+'/NS_MS_alltimes.dat', 'w+')
    fnsms.write('#1:TotalTime(Myr) 2:id0 3:id1 4:m0[MSUN] 5:m1[MSUN] 6:B[G] 7:P[sec] 8:startype0 9:startype1 10:Porb(days) 11:ecc 12:radrol0 13:radrol1 14:r 15:tcflag 16:rbflag\n')

    fnsms_selected=open(savepath+'/rb_progenitor_alltimes.dat', 'w+')
    fnsms_selected.write('#1:TotalTime(Myr) 2:id0 3:id1 4:m0[MSUN] 5:m1[MSUN] 6:B[G] 7:P[sec] 8:startype0 9:startype1 10:Porb(days) 11:ecc 12:radrol0 13:radrol1 14:r 15:tcflag 16:rbflag\n')
    with open(filestr+'.morepulsars.dat', 'r') as fpsr:
        next(fpsr)
        for line in fpsr:
            data=line.split()
            if int(data[2])==1:
                if int(data[11])==13 and int(data[12])<=1:
                    if (int(data[3]) in ID0 and int(data[4]) in ID1) or (int(data[4]) in ID0 and int(data[3]) in ID1):
                        tcflag=1
                    elif int(data[3]) in ID0 or int(data[3]) in ID1:
                        tcflag=2
                    else:
                        tcflag=3

                    porb=uc.au_to_period(float(data[13]), float(data[5]), float(data[6]))
                    if float(data[6])<=2.0 and porb<=2.0:
                        rbflag=1
                    else:
                        rbflag=0

                    fnsms.write('%f %d %d %f %f %e %f %d %d %f %f %f %f %f %d %d\n'%(float(data[1])*t_conv, int(data[3]), int(data[4]), float(data[5]), float(data[6]), float(data[7]), float(data[9]), int(data[11]), int(data[12]), porb, float(data[14]), float(data[15]), float(data[16]), float(data[19])*l_conv, tcflag, rbflag))

                    if rbflag==1:
                        fnsms_selected.write('%f %d %d %f %f %e %f %d %d %f %f %f %f %f %d %d\n'%(float(data[1])*t_conv, int(data[3]), int(data[4]), float(data[5]), float(data[6]), float(data[7]), float(data[9]), int(data[11]), int(data[12]), porb, float(data[14]), float(data[15]), float(data[16]), float(data[19])*l_conv, tcflag, rbflag))


                if int(data[12])==13 and int(data[11])<=1:
                    if (int(data[3]) in ID0 and int(data[4]) in ID1) or (int(data[4]) in ID0 and int(data[3]) in ID1):
                        tcflag=1
                    elif int(data[4]) in ID0 or int(data[4]) in ID1:
                        tcflag=2
                    else:
                        tcflag=3

                    porb=uc.au_to_period(float(data[13]), float(data[5]), float(data[6]))
                    if float(data[5])<=2.0 and porb<=2.0:
                        rbflag=1
                    else:
                        rbflag=0

                    fnsms.write('%f %d %d %f %f %e %f %d %d %f %f %f %f %f %d %d\n'%(float(data[1])*t_conv, int(data[4]), int(data[3]), float(data[6]), float(data[5]), float(data[8]), float(data[10]), int(data[12]), int(data[11]), porb, float(data[14]), float(data[16]), float(data[15]), float(data[19])*l_conv, tcflag, rbflag))

                    if rbflag==1:
                       fnsms_selected.write('%f %d %d %f %f %e %f %d %d %f %f %f %f %f %d %d\n'%(float(data[1])*t_conv, int(data[4]), int(data[3]), float(data[6]), float(data[5]), float(data[8]), float(data[10]), int(data[12]), int(data[11]), porb, float(data[14]), float(data[16]), float(data[15]), float(data[19])*l_conv, tcflag, rbflag)) 

    fnsms.close()
    fnsms_selected.close()



##Find the unique redback progenitors at their first and last timesteps from the rb_progenitor_alltimes.dat file
def find_rbprogen_Unique(savepath, modelpath):
    data=np.genfromtxt(savepath+'rb_progenitor_alltimes.dat')
    times=data[:,0]; id0=data[:,1]; id1=data[:,2]
    alltimes=list(Counter(times).keys())
    #print(len(alltimes))

    id0_unique_end=[]; id1_unique_end=[]; time_unique_end=[]
    idstr_hold1_end=[str(0)]
    idstr_hold2_end=[str(0)]
    for i in range(len(alltimes)-1, 0, -1):
        timeno=alltimes[i]

        for j in range(len(times)-1, 0, -1):
            if times[j]==timeno:
                idstr=str(int(id0[j]))+str(int(id1[j]))
                check=1
                for k in range(len(idstr_hold1_end)):
                    if idstr==idstr_hold1_end[k] or idstr==idstr_hold2_end[k]:
                        check=0
    
                if check==1:
                    time_unique_end.append(times[j]); id0_unique_end.append(id0[j]); id1_unique_end.append(id1[j])

                    idstr_hold1_end.append(idstr)
                    idstr_hold2_end.append(str(int(id1[j]))+str(int(id0[j])))

        #print(i)
    print(id0_unique_end)

    id0_unique_begin=[]; id1_unique_begin=[]; time_unique_begin=[]
    idstr_hold1_begin=[str(0)]
    idstr_hold2_begin=[str(0)]
    for i in range(len(alltimes)):
        timeno=alltimes[i]

        for j in range(len(times)):
            if times[j]==timeno:
                idstr=str(int(id0[j]))+str(int(id1[j]))
                check=1
                for k in range(len(idstr_hold1_begin)):
                    if idstr==idstr_hold1_begin[k] or idstr==idstr_hold2_begin[k]:
                        check=0
    
                if check==1:
                    time_unique_begin.append(times[j]); id0_unique_begin.append(id0[j]); id1_unique_begin.append(id1[j])

                    idstr_hold1_begin.append(idstr)
                    idstr_hold2_begin.append(str(int(id1[j]))+str(int(id0[j])))
    print(id0_unique_begin)

    funique=open(savepath+'rb_prog_unique.dat', 'w+')
    funique.write('#1:TotalTime(Myr) 2:id0 3:id1 4:m0[MSUN] 5:m1[MSUN] 6:B[G] 7:P[sec] 8:startype0 9:startype1 10:Porb(days) 11:ecc 12:radrol0 13:radrol1 14:r 15:tcflag 16:rbflag 17.Formation 18.Disruption 29.Nenc\n')
    frball=open(savepath+'rb_progenitor_alltimes.dat', 'r')
    datarball=frball.readlines()
    #print(datarball)

    for m in range(len(id0_unique_begin)):
        for x in range(1, len(datarball)):
            data=datarball[x].split()
            if float(data[0])==time_unique_begin[m] and int(data[1])==id0_unique_begin[m] and int(data[2])==id1_unique_begin[m]:
                theline=datarball[x].strip('\n')+' '+str(-100)+' '+str(-100)+' '+str(-100)+'\n'
                funique.write(theline)
                print(time_unique_begin[m])

        for n in range(len(id0_unique_end)):
            if id0_unique_begin[m]==id0_unique_end[n]:
                print(id0_unique_begin[m], id0_unique_end[n])
                disruption=nh.find_binary_disruipt(id0_unique_begin[m], id1_unique_begin[m], modelpath)
                formation=nh.find_binary_form(id0_unique_begin[m], id1_unique_begin[m], modelpath)
                any_enc=nh.find_binary_encounter(id0_unique_begin[m], id1_unique_begin[m], modelpath, [time_unique_begin[m], time_unique_end[n]])
                for x in range(1, len(datarball)):
                    data=datarball[x].split()
                    if float(data[0])==time_unique_end[n] and int(data[1])==id0_unique_end[n] and int(data[2])==id1_unique_end[n]:
                        theline=datarball[x].strip('\n')+' '+formation+' '+disruption+' '+any_enc+'\n'
                        funique.write(theline)
    
    frball.close()
    funique.close()



##Find NS-MS binaries that contain a MSP in many models
def find_msp_NSMS_9to14Gyr(savepath):
    data=np.genfromtxt(savepath+'msp_9to14Gyr.dat')
    model=data[:,0]; id0=data[:,12]; id1=data[:,13]; k0=data[:,14]; k1=data[:,15]
    model_unique=list(Counter(model).keys())

    f_all=open(savepath+'msp_9to14Gyr.dat')
    data_all=f_all.readlines()


    #f=open(savepath+'msp_nsms_9to14Gyr.dat', 'w+')
    #f.write('#1.Model 2.Time(Myr) 3.Status 4.r(pc) 5.B(G) 6.P(sec) 7.dmdt0(Msun/yr) 8.dmdt1(Msun/yr) 9.rolrad0 10.rolrad1 11.m0(Msun) 12.m1(Msun) 13.ID0 14.ID1 15.k0 16.k1 17.a(AU) 18.ecc 19.Formation\n')
    #for i in range(1, len(data_all)):
    #    line=data_all[i].split()
    #    if int(line[14])==13 and 0<=int(line[15])<=1:
    #        f.write(data_all[i])

    #f.close()


##Find NSs that exchange into a binary with low-mass main-sequence star companions.
def find_NS_MS_exchange(pathlist, savepath):
    paths = np.genfromtxt(pathlist, dtype=str)

    f = open(savepath, 'a+')
    f.write('#1.Model 2.Time(Myr) 3.ID0 4.ID1 5.M0 6.M1 7.K0 8.K1 9.SMA[AU] 10.ECC 11.Period[days]\n')

    for ii in range(len(paths)):
        filestr=paths[ii]+'initial'
        t_conv=dyn.conv('t', filestr+'.conv.sh')

        binintfile=filestr+'.binint.log'

        binint=scripts3.read_binint(binintfile)
        print('binint read')
        for xx in range(len(binint)):
            inputs=binint[xx]['input']
            outputs=binint[xx]['output']
            for m in range(len(outputs)):
                if int(outputs[m]['no'])==2: 
                    if int(outputs[m]['startype'][1])==13 and 0<=int(outputs[m]['startype'][0])<=1 and float(outputs[m]['m'][0])<=2.:
                        nsid = int(outputs[m]['ids'][1]); compid = int(outputs[m]['ids'][0])
                        mns = float(outputs[m]['m'][1]); mcomp = float(outputs[m]['m'][0])
                        time = binint[xx]['type']['time']*t_conv

                        for n in range(len(inputs)):
                            if int(inputs[n]['no'])==1 and int(inputs[n]['ids'][0])== nsid and float(outputs[m]['a'][0])<=1.:
                                period = uc.au_to_period(float(outputs[m]['a'][0]), mns, mcomp)
                                f.write('%d %f %d %d %f %f %d %d %f %f %f\n'%(ii, time, nsid, compid, mns, mcomp, int(outputs[m]['startype'][1]), int(outputs[m]['startype'][0]), float(outputs[m]['a'][0]), float(outputs[m]['e'][0]), period))
                                break
                            elif int(inputs[n]['no'])==2:
                                if nsid in inputs[n]['ids'] and compid in inputs[n]['ids']:
                                    break
                                elif nsid in inputs[n]['ids'] and float(outputs[m]['a'][0])<=1.:
                                    period = uc.au_to_period(float(outputs[m]['a'][0]), mns, mcomp)
                                    f.write('%d %f %d %d %f %f %d %d %f %f %f\n'%(ii, time, nsid, compid, mns, mcomp, int(outputs[m]['startype'][1]), int(outputs[m]['startype'][0]), float(outputs[m]['a'][0]), float(outputs[m]['e'][0]), period))
                                    break
         

                    if int(outputs[m]['startype'][0])==13 and 0<=int(outputs[m]['startype'][1])<=1 and float(outputs[m]['m'][1])<=2.:
                        nsid = int(outputs[m]['ids'][0]); compid = int(outputs[m]['ids'][1])
                        mns = float(outputs[m]['m'][0]); mcomp = float(outputs[m]['m'][1])
                        time = binint[xx]['type']['time']*t_conv

                        for n in range(len(inputs)):
                            if int(inputs[n]['no'])==1 and int(inputs[n]['ids'][0])== nsid and float(outputs[m]['a'][0])<=1.:
                                period = uc.au_to_period(float(outputs[m]['a'][0]), mns, mcomp)
                                f.write('%d %f %d %d %f %f %d %d %f %f %f\n'%(ii, time, nsid, compid, mns, mcomp, int(outputs[m]['startype'][0]), int(outputs[m]['startype'][1]), float(outputs[m]['a'][0]), float(outputs[m]['e'][0]), period))
                                break
                            elif int(inputs[n]['no'])==2:
                                if nsid in inputs[n]['ids'] and compid in inputs[n]['ids']:
                                    break
                                elif nsid in inputs[n]['ids'] and float(outputs[m]['a'][0])<=1.:
                                    period = uc.au_to_period(float(outputs[m]['a'][0]), mns, mcomp)
                                    f.write('%d %f %d %d %f %f %d %d %f %f %f\n'%(ii, time, nsid, compid, mns, mcomp, int(outputs[m]['startype'][0]), int(outputs[m]['startype'][1]), float(outputs[m]['a'][0]), float(outputs[m]['e'][0]), period))
                                    break


                if int(outputs[m]['no'])==3:
                    if int(outputs[m]['startype'][1])==13 and 0<=int(outputs[m]['startype'][0])<=1 and float(outputs[m]['m'][0])<=2.:
                        nsid = int(outputs[m]['ids'][1]); compid = int(outputs[m]['ids'][0])
                        mns = float(outputs[m]['m'][1]); mcomp = float(outputs[m]['m'][0])
                        time = binint[xx]['type']['time']*t_conv

                        for n in range(len(inputs)):
                            if int(inputs[n]['no'])==1 and int(inputs[n]['ids'][0])== nsid:
                                period = uc.au_to_period(float(outputs[m]['a'][0]), mns, mcomp)
                                f.write('%d %f %d %d %f %f %d %d %f %f %f\n'%(ii, time, nsid, compid, mns, mcomp, int(outputs[m]['startype'][1]), int(outputs[m]['startype'][0]), float(outputs[m]['a'][0]), float(outputs[m]['e'][0]), period))
                                break
                            elif int(inputs[n]['no'])==2:
                                if nsid in inputs[n]['ids'] and compid in inputs[n]['ids']:
                                    break
                                elif nsid in inputs[n]['ids']:
                                    period = uc.au_to_period(float(outputs[m]['a'][0]), mns, mcomp)
                                    f.write('%d %f %d %d %f %f %d %d %f %f %f\n'%(ii, time, nsid, compid, mns, mcomp, int(outputs[m]['startype'][1]), int(outputs[m]['startype'][0]), float(outputs[m]['a'][0]), float(outputs[m]['e'][0]), period))
                                    break
         

                    if int(outputs[m]['startype'][0])==13 and 0<=int(outputs[m]['startype'][1])<=1 and float(outputs[m]['m'][1])<=2.:
                        nsid = int(outputs[m]['ids'][0]); compid = int(outputs[m]['ids'][1])
                        mns = float(outputs[m]['m'][0]); mcomp = float(outputs[m]['m'][1])
                        time = binint[xx]['type']['time']*t_conv

                        for n in range(len(inputs)):
                            if int(inputs[n]['no'])==1 and int(inputs[n]['ids'][0])== nsid:
                                period = uc.au_to_period(float(outputs[m]['a'][0]), mns, mcomp)
                                f.write('%d %f %d %d %f %f %d %d %f %f %f\n'%(ii, time, nsid, compid, mns, mcomp, int(outputs[m]['startype'][0]), int(outputs[m]['startype'][1]), float(outputs[m]['a'][0]), float(outputs[m]['e'][0]), period))
                                break
                            elif int(inputs[n]['no'])==2:
                                if nsid in inputs[n]['ids'] and compid in inputs[n]['ids']:
                                    break
                                elif nsid in inputs[n]['ids']:
                                    period = uc.au_to_period(float(outputs[m]['a'][0]), mns, mcomp)
                                    f.write('%d %f %d %d %f %f %d %d %f %f %f\n'%(ii, time, nsid, compid, mns, mcomp, int(outputs[m]['startype'][0]), int(outputs[m]['startype'][1]), float(outputs[m]['a'][0]), float(outputs[m]['e'][0]), period))
                                    break                  

        print(ii)

    f.close()


##Find the frequency of all collisions from single-single and from fewbody in the collision file
##Make a table like the BSE collision table
def find_num_allcollision():
    pathlist = np.genfromtxt('/projects/b1095/syr904/projects/tidal_capture/path_allfinished_newruns_maingrid.dat', dtype=str)
    filepaths = pathlist[:,0]; status = pathlist[:,1]
    
    index_delete = []
    for i in range(len(status)):
        if status[i] != '1':
            index_delete.append(i)

    path_finished = list(np.delete(filepaths, index_delete))
    #print(len(path_finished))

    ss = np.zeros([15, 15]); fewbody = np.zeros([15, 15])
    fewbody_34 = np.zeros([15, 15])
    n34 = 0; n34_gg = 0; n34_g = 0
    n_ss_gg = 0; n_fb_gg = 0
    n_ss_g = 0; n_fb_g = 0
    for j in range(len(path_finished)):
        #print(path_finished[j])
        collfile = path_finished[j]+'initial.collision.log'
        colldata = scripts1.collision(collfile)
        print(len(colldata))
        
        IDs = list(colldata.keys())

        for k in range(len(IDs)):
            if colldata[IDs[k]]['interaction'] == 'single-single':
                x_tmp = int(colldata[IDs[k]]['parents']['types'][0])
                y_tmp = int(colldata[IDs[k]]['parents']['types'][1])

                if x_tmp<y_tmp: xvalue = x_tmp; yvalue = y_tmp
                else: xvalue = y_tmp; yvalue = x_tmp

                ss[xvalue, yvalue]+=1

                if (2<=xvalue<=9 and xvalue!=7) and (2<=yvalue<=9 and yvalue!=7):
                    n_ss_gg += 1
                elif (2<=xvalue<=9 and xvalue!=7) or (2<=yvalue<=9 and yvalue!=7):
                    n_ss_g += 1

            elif colldata[IDs[k]]['interaction'] == 'binary-single' or colldata[IDs[k]]['interaction'] == 'binary-binary':
                nopars = len(colldata[IDs[k]]['parents']['types'])
                partypes = colldata[IDs[k]]['parents']['types']
                if nopars < 3:
                    x_tmp = int(colldata[IDs[k]]['parents']['types'][0])
                    y_tmp = int(colldata[IDs[k]]['parents']['types'][1])

                    if x_tmp<y_tmp: xvalue = x_tmp; yvalue = y_tmp
                    else: xvalue = y_tmp; yvalue = x_tmp

                    fewbody[xvalue, yvalue]+=1
                    fewbody_34[xvalue, yvalue]+=1

                    if (2<=xvalue<=9 and xvalue!=7) and (2<=yvalue<=9 and yvalue!=7):
                        n_fb_gg += 1
                    elif (2<=xvalue<=9 and xvalue!=7) or (2<=yvalue<=9 and yvalue!=7):
                        n_fb_g += 1
                else:
                    n34 += 1
                    n_giant = 0
                    for x in range(nopars):
                        if 2<=int(partypes[x])<=9 and int(partypes[x])!=7:
                            n_giant += 1

                    if n_giant>=2: 
                        n34_gg += 1
                    elif n_giant>=1:
                        n34_g += 1

                    pt_sort = np.sort(partypes)
                    if nopars == 3:
                        xvalue = int(pt_sort[0])
                        yvalue = int(pt_sort[1])
                        zvalue = int(pt_sort[2])

                        fewbody_34[xvalue, yvalue]+=1
                        fewbody_34[xvalue, zvalue]+=1
                        fewbody_34[yvalue, zvalue]+=1

                    if nopars == 4:
                        xvalue = int(pt_sort[0])
                        yvalue = int(pt_sort[1])
                        zvalue = int(pt_sort[2])
                        wvalue = int(pt_sort[3])

                        fewbody_34[xvalue, yvalue]+=1
                        fewbody_34[xvalue, zvalue]+=1
                        fewbody_34[xvalue, wvalue]+=1
                        fewbody_34[yvalue, zvalue]+=1
                        fewbody_34[yvalue, wvalue]+=1
                        fewbody_34[zvalue, wvalue]+=1

                    #print(colldata[IDs[k]])


    sscoll_sum = np.sum(ss)
    fewbodycoll_sum = np.sum(fewbody)
    tot_sum = sscoll_sum+fewbodycoll_sum

    ss_frac = ss/tot_sum
    fewbody_frac = fewbody/tot_sum

    ss_permodel = ss/len(path_finished)
    fewbody_permodel = fewbody/len(path_finished)
    fewbody34_permodel = fewbody_34/len(path_finished)


    #np.savetxt('/projects/b1095/syr904/projects/tidal_capture/sscoll_frac_maingrid.txt', ss_frac, fmt = '%f', delimiter=" ")
    #np.savetxt('/projects/b1095/syr904/projects/tidal_capture/fewbodycoll_frac_maingrid.txt', fewbody_frac, fmt = '%f', delimiter=" ")

    np.savetxt('/projects/b1095/syr904/projects/tidal_capture/sscoll_permodel_maingrid.txt', ss_permodel, fmt = '%f', delimiter=" ")
    np.savetxt('/projects/b1095/syr904/projects/tidal_capture/fewbodycoll_permodel_maingrid.txt', fewbody_permodel, fmt = '%f', delimiter=" ")
    np.savetxt('/projects/b1095/syr904/projects/tidal_capture/fewbody34coll_permodel_maingrid.txt', fewbody34_permodel, fmt = '%f', delimiter=" ")


    #print(ss)
    #print(fewbody)
    print(sscoll_sum, fewbodycoll_sum, tot_sum, n34)
    print(n_ss_gg, n_fb_gg, n34_gg)
    print(n_ss_g, n_fb_g, n34_g)


##Find how the tidal capture binaries/giant collision binaries are disrupted or not
#def find_howbin_disrupt():



##Find the number of NS/BH collisions in the models
##Startype=13 or 14
def get_NS_collision(pathlist, start, end, startype):
    sourcedir=np.genfromtxt(pathlist, dtype=str)
    filepaths=sourcedir[:,0]; status=sourcedir[:,1]

    #model=[]; model_status=[]; mm=[]; mcom=[]; ktypem=[]; kcom=[]; timem=[]; idm=[]; rm=[]; colltype=[]
    fcoll=open('/projects/b1095/syr904/projects/PULSAR2/newruns/tidal_capture/BH_coll_all.dat', 'a+')
    fcoll.write('#1.Model 2.Time(Myr) 3.IDcoll 4.Radius(pc) 5.Mcoll 6.M0 7.M1 8.M2 9.M3 10.kcoll 11.k0 12.k1 13.k2 14.k3 15.model_status 16.COLLTYPE\n')
    for i in range(len(filepaths)):
        filestr=filepaths[i]+'initial'

        t_conv=dyn.conv('t', filestr+'.conv.sh')
        l_conv=dyn.conv('l', filestr+'.conv.sh')

        collfile=filestr+'.collision.log'
        collfile2=filestr+'2.collision.log'
        colldata=scripts1.readcollfile(collfile)
        if os.path.isfile(collfile2) and os.path.getsize(collfile2) > 0:
            colldata2=scripts1.readcollfile(collfile2)
            colldata=colldata+colldata2

        for j in range(len(colldata)):
            line=colldata[j].split()
            if line[1]=='single-single':  ##Single-single star collision
                colltype='SS'
                if int(line[11])==startype or int(line[12])==startype:
                    model=i; model_status=int(status[i]); timem=t_conv*float(line[0])
                    mm=float(line[4]); m0=float(line[6]); m1=float(line[8]); m2=-100; m3=-100
                    ktypem=int(line[10]); ktype0=int(line[11]); ktype1=int(line[12]); ktype2=-100; ktype3=-100
                    idm=int(line[3]); rm=float(line[9])*l_conv

                    fcoll.write('%d %f %d %f %f %f %f %f %f %d %d %d %d %d %d %s\n'%(model, timem, idm, rm, mm, m0, m1, m2, m3, ktypem, ktype0, ktype1, ktype2, ktype3, model_status, colltype))


            if line[1]=='binary-single':   ##Binary-single collision
                colltype='BS'
                if int(line[2])==2:
                    if int(line[11])==startype or int(line[12])==startype:
                        model=i; model_status=int(status[i]); timem=t_conv*float(line[0])
                        mm=float(line[4]); m0=float(line[6]); m1=float(line[8]); m2=-100; m3=-100
                        ktypem=int(line[10]); ktype0=int(line[11]); ktype1=int(line[12]); ktype2=-100; ktype3=-100
                        idm=int(line[3]); rm=float(line[9])*l_conv

                        fcoll.write('%d %f %d %f %f %f %f %f %f %d %d %d %d %d %d %s\n'%(model, timem, idm, rm, mm, m0, m1, m2, m3, ktypem, ktype0, ktype1, ktype2, ktype3, model_status, colltype))

                if int(line[2])==3:
                    if int(line[13])==startype or int(line[14])==startype or int(line[15])==startype:
                        model=i; model_status=int(status[i]); timem=t_conv*float(line[0])
                        mm=float(line[4]); m0=float(line[6]); m1=float(line[8]); m2=float(line[10]); m3=-100
                        ktypem=int(line[12]); ktype0=int(line[13]); ktype1=int(line[14]); ktype2=int(line[15]); ktype3=-100
                        idm=int(line[3]); rm=float(line[11])*l_conv

                        fcoll.write('%d %f %d %f %f %f %f %f %f %d %d %d %d %d %d %s\n'%(model, timem, idm, rm, mm, m0, m1, m2, m3, ktypem, ktype0, ktype1, ktype2, ktype3, model_status, colltype))


            if line[1]=='binary-binary':   ##Binary-binary collision
                colltype='BB'
                if int(line[2])==2:
                    if int(line[11])==startype or int(line[12])==startype:
                        model=i; model_status=int(status[i]); timem=t_conv*float(line[0])
                        mm=float(line[4]); m0=float(line[6]); m1=float(line[8]); m2=-100; m3=-100
                        ktypem=int(line[10]); ktype0=int(line[11]); ktype1=int(line[12]); ktype2=-100; ktype3=-100
                        idm=int(line[3]); rm=float(line[9])*l_conv

                        fcoll.write('%d %f %d %f %f %f %f %f %f %d %d %d %d %d %d %s\n'%(model, timem, idm, rm, mm, m0, m1, m2, m3, ktypem, ktype0, ktype1, ktype2, ktype3, model_status, colltype))

                if int(line[2])==3:
                    if int(line[13])==startype or int(line[14])==startype or int(line[15])==startype:
                        model=i; model_status=int(status[i]); timem=t_conv*float(line[0])
                        mm=float(line[4]); m0=float(line[6]); m1=float(line[8]); m2=float(line[10]); m3=-100
                        ktypem=int(line[12]); ktype0=int(line[13]); ktype1=int(line[14]); ktype2=int(line[15]); ktype3=-100
                        idm=int(line[3]); rm=float(line[11])*l_conv

                        fcoll.write('%d %f %d %f %f %f %f %f %f %d %d %d %d %d %d %s\n'%(model, timem, idm, rm, mm, m0, m1, m2, m3, ktypem, ktype0, ktype1, ktype2, ktype3, model_status, colltype))

                if int(line[2])==4:
                    if int(line[15])==startype or int(line[16])==startype or int(line[17])==startype or int(line[18])==startype:
                        model=i; model_status=int(status[i]); timem=t_conv*float(line[0])
                        mm=float(line[4]); m0=float(line[6]); m1=float(line[8]); m2=float(line[10]); m3=float(line[12])
                        ktypem=int(line[14]); ktype0=int(line[15]); ktype1=int(line[16]); ktype2=int(line[17]); ktype3=int(line[18])
                        idm=int(line[3]); rm=float(line[13])*l_conv

                        fcoll.write('%d %f %d %f %f %f %f %f %f %d %d %d %d %d %d %s\n'%(model, timem, idm, rm, mm, m0, m1, m2, m3, ktypem, ktype0, ktype1, ktype2, ktype3, model_status, colltype))


        print(i)

    fcoll.close()



##Total numbers of different collision systems
def get_num_collision(nscollfile, typeinput):
    data=np.genfromtxt(nscollfile, dtype=None)
    
    count_ms=0; count_giant=0; count_wd=0; count_nsbh=0 #; count_bh=0
    ntot=0
    for i in range(len(data)):
        if data[i][14]==1: ntot+=1
        if data[i][14]==1 and data[i][-1].decode("utf-8")==typeinput:
            klist=[data[i][10], data[i][11], data[i][12], data[i][13]]
            knum=Counter(klist)
            for j in range(15):
                if j<=1:
                    if knum[j]>=1:
                        count_ms+=1
                        break

                elif 2<=j<=9:
                   if knum[j]>=1:
                        count_giant+=1
                        break

                elif 10<=j<=12:
                    if knum[j]>=1:
                        count_wd+=1
                        break

                else:
                    if knum[j]>=1:
                        count_nsbh+=1
                        break

    print(count_ms, count_giant, count_wd, count_nsbh)
    print(ntot)


##Find initial properties of a NS binary that goes into binary interaction and leads to DNS and DNS merger. Return the initial properties of sma and m0 and m1.
def get_NS_binint(pathlist, therv):
    sourcedir=np.genfromtxt(pathlist, dtype=str)
    #filepaths=sourcedir[:,0]; status=sourcedir[:,1]
    filepaths=['/projects/b1095/syr904/cmc/extreme_model/N8e5fbh100rv0.5_NSkick20_BHkick_300_IMF20/']; status=[1]

    id0=[]; id1=[]; m0=[]; m1=[]; k0=[]; k1=[]; a=[]; e=[]; time=[]; model=[]
    count_allbinint=[[],[],[],[]]; count_alldns=[[],[],[],[]]; count_allmerge=[[],[],[],[]]
    sma=[[], [], [], []]; m0_sma=[[],[],[],[]]; m1_sma=[[],[],[],[]]


    limit_low=[0, 2, 10, 13]; limit_high=[1, 9, 12, 14]
    for i in range(len(filepaths)):
        count_ns=[0, 0, 0, 0]; count_dns=[0, 0, 0, 0]; count_dns_merge=[0, 0, 0, 0]
        #n_ini, metal, r_g, r_v = uc.find_init_conditions(filepaths[i])
        r_v=0.5

        filestr=filepaths[i]+'initial'
        t_conv=dyn.conv('t', filestr+'.conv.sh')

        binintfile=filestr+'.binint.log'
        binintfile2=filestr+'2.binint.log'

        if int(status[i])==1 and r_v==therv:
            if os.path.isfile(binintfile2) and os.path.getsize(binintfile2) > 0:
                binint=scripts3.read_binint(binintfile2)
                for j in range(len(binint)):
                    bininput=binint[j]['input']
                    binoutput=binint[j]['output']
                    for k in range(len(bininput)):
                        if int(bininput[k]['no'])==2: 
                            for x in range(4):
                                if (int(bininput[k]['startype'][1])==13 and limit_low[x]<=int(bininput[k]['startype'][0])<=limit_high[x]) or (int(bininput[k]['startype'][0])==13 and limit_low[x]<=int(bininput[k]['startype'][1])<=limit_high[x]):
                                    count_ns[x]+=1
                                    #print bininput[k]['a']
                                    sma[x].append(float(bininput[k]['a']))
                                    m0_sma[x].append(float(bininput[k]['m'][0])); m1_sma[x].append(float(bininput[k]['m'][1]))
                                    for l in range(len(binoutput)):
                                        if int(binoutput[l]['no'])==2:
                                            if binoutput[l]['ids'][0].find(':')==-1 and binoutput[l]['ids'][1].find(':')==-1:
                                                if int(binoutput[l]['startype'][0])==13 and int(binoutput[l]['startype'][1])==13:
                                                    count_dns[x]+=1
                                                    time.append(float(binint[j]['type']['time'])*t_conv)
                                                    id0.append(int(binoutput[l]['ids'][0])); id1.append(int(binoutput[l]['ids'][1]))
                                                    m0.append(float(binoutput[l]['m'][0])); m1.append(float(binoutput[l]['m'][1]))
                                                    k0.append(int(binoutput[l]['startype'][0])); k1.append(int(binoutput[l]['startype'][1]))
                                                    a.append(float(binoutput[l]['a'])); e.append(float(binoutput[l]['e']))
                                                    model.append(i)
                                                    t_inspiral=gwcalc.t_inspiral_2(float(binoutput[l]['a']), float(binoutput[l]['e']), float(binoutput[l]['m'][0]), float(binoutput[l]['m'][1]), 0, 0, 0, 1100)/10**6 ##in Myr
                                                    if t_inspiral+float(binint[j]['type']['time'])*t_conv<14000.0:
                                                        count_dns_merge[x]+=1


            binint=scripts3.read_binint(binintfile)
            for j in range(len(binint)):
                bininput=binint[j]['input']
                binoutput=binint[j]['output']
                for k in range(len(bininput)):
                    if int(bininput[k]['no'])==2:
                        for x in range(4):
                            if (int(bininput[k]['startype'][1])==13 and limit_low[x]<=int(bininput[k]['startype'][0])<=limit_high[x]) or (int(bininput[k]['startype'][0])==13 and limit_low[x]<=int(bininput[k]['startype'][1])<=limit_high[x]):
                                count_ns[x]+=1
                                #print bininput[k]['a']
                                sma[x].append(float(bininput[k]['a']))
                                m0_sma[x].append(float(bininput[k]['m'][0])); m1_sma[x].append(float(bininput[k]['m'][1]))
                                for l in range(len(binoutput)):
                                    if int(binoutput[l]['no'])==2:
                                        if binoutput[l]['ids'][0].find(':')==-1 and binoutput[l]['ids'][1].find(':')==-1:
                                            if int(binoutput[l]['startype'][0])==13 and int(binoutput[l]['startype'][1])==13:
                                                count_dns[x]+=1
                                                time.append(float(binint[j]['type']['time'])*t_conv)
                                                id0.append(int(binoutput[l]['ids'][0])); id1.append(int(binoutput[l]['ids'][1]))
                                                m0.append(float(binoutput[l]['m'][0])); m1.append(float(binoutput[l]['m'][1]) )
                                                k0.append(int(binoutput[l]['startype'][0])); k1.append(int(binoutput[l]['startype'][1]))
                                                a.append(float(binoutput[l]['a'][0])); e.append(float(binoutput[l]['e'][0]))
                                                model.append(i)
                                                t_inspiral=lisa.inspiral_time_peters(float(binoutput[l]['a'][0]), float(binoutput[l]['e'][0]), float(binoutput[l]['m'][0]), float(binoutput[l]['m'][1]))*10**3 ##in Myr
                                                if t_inspiral+float(binint[j]['type']['time'])*t_conv<14000.0:
                                                    count_dns_merge[x]+=1
                                                    


            for y in range(4):
                count_allbinint[y].append(count_ns[y])
                count_alldns[y].append(count_dns[y])
                count_allmerge[y].append(count_dns_merge[y])


        print(i)

    print(np.sum(count_allbinint[0]), np.sum(count_alldns[0]), np.sum(count_allmerge[0]))
    print(np.sum(count_allbinint[1]), np.sum(count_alldns[1]), np.sum(count_allmerge[1]))
    print(np.sum(count_allbinint[2]), np.sum(count_alldns[2]), np.sum(count_allmerge[2]))

    return sma, m0_sma, m1_sma



            


