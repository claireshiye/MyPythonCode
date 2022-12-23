import numpy as np
import pandas as pd
import collections
from collections import Counter
import dynamics as dyn
from glob import glob


def extract_psr_binary_evol(mspflag):
    sourcecdir = np.genfromtxt('/projects/b1091/CMC_Grid_March2019/rundir/path_final.dat', dtype = 'str')
    allpath = sourcecdir[:,0]; allstatus = sourcecdir[:,1]

    paths = allpath[allstatus == '1']

    for kk in range(len(paths)):
        psrfile = np.genfromtxt(paths[kk]+'initial.morepulsars.dat')
        t_conv = dyn.conv('t', paths[kk]+'initial.conv.sh')
        l_conv = dyn.conv('l', paths[kk]+'initial.conv.sh')

        ##Extract all pulsar IDs
        allmsp_id = []
        with open(psrfile, 'r') as fpsr:
            next(fpsr)
            for line in fpsr:
                data = line.split()
                if int(data[2]) == 1:
                    deathcut1=(float(data[9])**2)*(0.17*10**12)
                    if deathcut1<=float(data[7]): 
                        if mspflag == 'msp' and float(data[9])<=0.03: 
                           allmsp_id.append(data[3])

                    deathcut2=(float(data[10])**2)*(0.17*10**12)
                    if deathcut1<=float(data[8]): 
                        if mspflag == 'msp' and float(data[10])<=0.03: 
                           allmsp_id.append(data[4])

                else:
                    deathcut=(float(data[9])**2)*(0.17*10**12)
                    if deathcut1<=float(data[7]): 
                        if mspflag == 'msp' and float(data[9])<=0.03: 
                           allmsp_id.append(data[3])

        ##Select unique pulsar IDs
        unique_IDs = Counter(allmsp_id).keys()

        psr_info = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
        #time, id0, id1, m0, m1, k0, k1, rad0, rad1, dmt0, dmt1, bacc0, bacc1, a, e, B, P, formation, r
        ##Extract all pulsar info
        with open(psrfile, 'r') as fpsr:
            next(fpsr)
            for line in fpsr:
                data = line.split()
                for xx in range(len(unique_IDs)):
                    if int(data[2])==1:
                        if data[3] == unique_IDs[xx]:
                            psr_info[1].append(data[3]); psr_info[2].append(data[4])
                            psr_info[3].append(float(data[5])); psr_info[4].append(float(data[6]))
                            psr_info[5].append(float(data[11])); psr_info[6].append(float(data[12]))
                            psr_info[7].append(float(data[15])); psr_info[8].append(float(data[16]))
                            psr_info[9].append(float(data[17])); psr_info[10].append(float(data[18]))
                            psr_info[11].append(float(data[22])); psr_info[12].append(float(data[23]))
                            psr_info[15].append(float(data[7]))
                            psr_info[16].append(float(data[9]))
                            psr_info[17].append(float(data[26]))

                            psr_info[0].append(float(data[1])*t_conv)
                            psr_info[13].append(float(data[13])); psr_info[14].append(float(data[14]))
                            psr_info[18].append(float(data[19])*l_conv)


                        if data[4] == unique_IDs[xx]:
                            psr_info[1].append(data[4]); psr_info[2].append(data[3])
                            psr_info[3].append(float(data[6])); psr_info[4].append(float(data[5]))
                            psr_info[5].append(float(data[12])); psr_info[6].append(float(data[11]))
                            psr_info[7].append(float(data[16])); psr_info[8].append(float(data[15]))
                            psr_info[9].append(float(data[18])); psr_info[10].append(float(data[17]))
                            psr_info[11].append(float(data[23])); psr_info[12].append(float(data[22]))
                            psr_info[15].append(float(data[8]))
                            psr_info[16].append(float(data[10]))
                            psr_info[17].append(float(data[27]))

                            psr_info[0].append(float(data[1])*t_conv)
                            psr_info[13].append(float(data[13])); psr_info[14].append(float(data[14]))
                            psr_info[18].append(float(data[19])*l_conv)


        s=paths[kk].split('/')
        n_star=s[-2]
        z=s[-3][1:]
        rg=s[-4][2:]
        rv=s[-5][2:]
        filename  = 'N'+n_star+'_rv'+rv+'_rg'+rg+'_z'+z+'.hdf5'
        columns = ['Time_myr', 'ID0', 'ID1', 'M0_msun', 'M1_msun', 'K0', 'K1', 'Radrol0', 'Radrol1', 'Dmdt0', 'Dmdt1', 'bacc0', 'bacc1', 'a_AU', 'e', 'B_G', 'P_sec', 'formation', 'r_pc']

        ##Extract separate pulsar system info
        for jj in range(len(unique_IDs)):
            theid = unique_IDs[jj]
            
            allIDs = np.array(psr_info[1])

            psr_info_time = np.array(psr_info[0])[allIDs == theid]
            print(len(psr_info_time))
            psr_info_time = psr_info_time.reshape(len(psr_info_time),1)
            print(len(psr_info_time))

            psr_info_id0 = np.array(psr_info[1])[allIDs == theid]
            psr_info_id0 = psr_info_id0.reshape(len(psr_info_id0),1)
            
            unique_psr_info = np.hstack((psr_info_time, psr_info_id0))
            print(unique_psr_info)

            for yy in range(2, len(psr_info)):
                psr_info_yy = np.array(psr_info[yy])[allIDs == theid]
                psr_info_yy = psr_info_yy.reshape(len(psr_info_yy),1)

                unique_psr_info = np.hstack((unique_psr_info, psr_info_yy))
                
        
        #df_unique = pd.DataFrame(data=unique_psr_info,columns=columns)
        #store = pd.HDFStore(filename)
        #store.put(, df_unique)
        #store.close()


def extract_NS_binary_presentday(typeflag, klowlim, khighlim, lowtimelim, hightimelim, savename):
    sourcecdir = np.genfromtxt('/projects/b1095/syr904/projects/PULSAR_Catalog/newruns/path_allfinished_newruns_maingrid.dat', dtype = 'str')
    allpath = sourcecdir[:,0]; allstatus = sourcecdir[:,1]
    allpathno = np.linspace(0, len(allpath)-1, len(allpath))

    paths = allpath[allstatus == '1']
    pathnos = allpathno[allstatus == '1']


    f = open('/projects/b1095/syr904/projects/PULSAR_Catalog/newruns/'+savename+'.dat', 'w+')
    f.write('#1.Model, 2.Time[Myr], 3.R[PC], 4.ID0, 5.ID1, 6.M0[Msun], 7.M1[Msun], 8.K0, 9.K1, 10.Rad0, 11.Rad1, 12.dmt0, 13.dmt1, 14.bacc0, 15.bacc1, 16.a[AU], 17.e, 18.B0[G], 19.P0[sec], 20.B1[G], 21.P1[sec]\n')
    #model, time, r, id0, id1, m0, m1, k0, k1, rad0, rad1, dmt0, dmt1, bacc0, bacc1, a, e, B0, P0, B1, P1, formation0, formation1
    
    for kk in range(len(paths)):
        print(paths[kk])
        psrfile = paths[kk]+'initial.morepulsars.dat'
        #print(psrfile)
        t_conv = dyn.conv('t', paths[kk]+'initial.conv.sh')
        l_conv = dyn.conv('l', paths[kk]+'initial.conv.sh')

        ##Extract all NS info
        with open(psrfile, 'r') as fpsr:
            next(fpsr)
            for line in fpsr:
                data = line.split()
                if typeflag == 'NSNS':
                    if int(data[2])==1 and lowtimelim <= float(data[1])*t_conv <= hightimelim:
                        if int(data[11])==13 and int(data[12])==13:
                            f.write('%d %f %f %s %s %f %f %d %d %f %f %f %f %f %f %f %f %e %f %e %f\n'%(pathnos[kk], float(data[1])*t_conv, float(data[19])*l_conv, data[3], data[4], float(data[5]), float(data[6]), int(data[11]), int(data[12]), float(data[15]), float(data[16]), float(data[17]), float(data[18]), float(data[22]), float(data[23]), float(data[13]), float(data[14]),float(data[7]), float(data[9]), float(data[8]), float(data[10])))

                if typeflag == 'MSP':
                    if lowtimelim <= float(data[1])*t_conv <= hightimelim:
                        if int(data[2]) == 1:
                            deathcut1=(float(data[9])**2)*(0.17*10**12)
                            deathcut2=(float(data[10])**2)*(0.17*10**12)
                            if deathcut1<=float(data[7]): 
                                if float(data[9])<=0.03: 
                                   f.write('%d %f %f %s %s %f %f %d %d %f %f %f %f %f %f %f %f %e %f %e %f\n'%(pathnos[kk], float(data[1])*t_conv, float(data[19])*l_conv, data[3], data[4], float(data[5]), float(data[6]), int(data[11]), int(data[12]), float(data[15]), float(data[16]), float(data[17]), float(data[18]), float(data[22]), float(data[23]), float(data[13]), float(data[14]),float(data[7]), float(data[9]), float(data[8]), float(data[10])))

                            elif deathcut2<=float(data[8]): 
                                if float(data[10])<=0.03:
                                    f.write('%d %f %f %s %s %f %f %d %d %f %f %f %f %f %f %f %f %e %f %e %f\n'%(pathnos[kk], float(data[1])*t_conv, float(data[19])*l_conv, data[4], data[3], float(data[6]), float(data[5]), int(data[12]), int(data[11]), float(data[16]), float(data[15]), float(data[18]), float(data[17]), float(data[23]), float(data[22]), float(data[13]), float(data[14]),float(data[8]), float(data[10]), float(data[7]), float(data[9])))

                            
                        else:
                            deathcut=(float(data[9])**2)*(0.17*10**12)
                            if deathcut<=float(data[7]): 
                                if float(data[9])<=0.03: 
                                    f.write('%d %f %f %s %s %f %f %d %d %f %f %f %f %f %f %f %f %e %f %e %f\n'%(pathnos[kk], float(data[1])*t_conv, float(data[19])*l_conv, data[3], data[4], float(data[5]), float(data[6]), int(data[11]), int(data[12]), float(data[15]), float(data[16]), float(data[17]), float(data[18]), float(data[22]), float(data[23]), float(data[13]), float(data[14]),float(data[7]), float(data[9]), float(data[8]), float(data[10])))


    f.close()


def extract_NS_binary_attime(typeflag, klowlim, khighlim, timelim, savename):
    sourcecdir = np.genfromtxt('/projects/b1095/syr904/projects/PULSAR_Catalog/newruns/path_allfinished_newruns_maingrid.dat', dtype = 'str')
    allpath = sourcecdir[:,0]; allstatus = sourcecdir[:,1]
    allpathno = np.linspace(0, len(allpath)-1, len(allpath))

    paths = allpath[allstatus == '1']
    pathnos = allpathno[allstatus == '1']


    f = open('/projects/b1095/syr904/projects/isolated_MSP/'+savename+'.dat', 'w+')
    f.write('#1.Model, 2.Time[Myr], 3.R[PC], 4.ID0, 5.ID1, 6.M0[Msun], 7.M1[Msun], 8.K0, 9.K1, 10.Rad0, 11.Rad1, 12.dmt0, 13.dmt1, 14.bacc0, 15.bacc1, 16.a[AU], 17.e, 18.B0[G], 19.P0[sec], 20.B1[G], 21.P1[sec]\n')
    #model, time, r, id0, id1, m0, m1, k0, k1, rad0, rad1, dmt0, dmt1, bacc0, bacc1, a, e, B0, P0, B1, P1, formation0, formation1
    
    Mtot = []; Rc = []; Rh = []; MSP_sin = []; MSP_bin = []; Model = []; Time = []
    Rc_obs = []; Rh_obs = []
    for kk in range(len(paths)):
        print(paths[kk])

        #print(psrfile)
        t_conv = dyn.conv('t', paths[kk]+'initial.conv.sh')
        l_conv = dyn.conv('l', paths[kk]+'initial.conv.sh')
        m_conv = dyn.conv('m', paths[kk]+'initial.conv.sh')

        psrfile = paths[kk]+'initial.morepulsars.dat'
        paramfile = np.sort(glob(paths[kk]+'initial.snap*.cluster_params.dat'))
        #print(paramfile)

        for xx in range(len(paramfile)):
            data_param = np.genfromtxt(paramfile[xx])
            if len(data_param)==0:
                continue
            if data_param[0,0]>=timelim:
                Rc_obs.append(data_param[0,9])
                Rh_obs.append(data_param[0,10])
                break


        Model.append(pathnos[kk])
        MSP_sin.append(0); MSP_bin.append(0)

        thetime_old = 0
        thetime_code = 0
        ##Extract all NS info
        with open(psrfile, 'r') as fpsr:
            next(fpsr)
            for line in fpsr:
                data = line.split()
                thetime = float(data[1])*t_conv
                if thetime > thetime_old > timelim:
                    psr_time = thetime_code
                    break
                if typeflag == 'NSNS':
                    if int(data[2])==1 and timelim <= thetime:
                        if int(data[11])==13 and int(data[12])==13:
                            f.write('%d %f %f %s %s %f %f %d %d %f %f %f %f %f %f %f %f %e %f %e %f\n'%(pathnos[kk], float(data[1])*t_conv, float(data[19])*l_conv, data[3], data[4], float(data[5]), float(data[6]), int(data[11]), int(data[12]), float(data[15]), float(data[16]), float(data[17]), float(data[18]), float(data[22]), float(data[23]), float(data[13]), float(data[14]),float(data[7]), float(data[9]), float(data[8]), float(data[10])))

                if typeflag == 'MSP':
                    if timelim <= thetime:
                        if int(data[2]) == 1:
                            deathcut1=(float(data[9])**2)*(0.17*10**12)
                            deathcut2=(float(data[10])**2)*(0.17*10**12)
                            if deathcut1<=float(data[7]): 
                                if float(data[9])<=1.: 
                                   f.write('%d %f %f %s %s %f %f %d %d %f %f %f %f %f %f %f %f %e %f %e %f\n'%(pathnos[kk], float(data[1])*t_conv, float(data[19])*l_conv, data[3], data[4], float(data[5]), float(data[6]), int(data[11]), int(data[12]), float(data[15]), float(data[16]), float(data[17]), float(data[18]), float(data[22]), float(data[23]), float(data[13]), float(data[14]),float(data[7]), float(data[9]), float(data[8]), float(data[10])))
                                   MSP_bin[-1]+=1

                            elif deathcut2<=float(data[8]): 
                                if float(data[10])<=1.:
                                    f.write('%d %f %f %s %s %f %f %d %d %f %f %f %f %f %f %f %f %e %f %e %f\n'%(pathnos[kk], float(data[1])*t_conv, float(data[19])*l_conv, data[4], data[3], float(data[6]), float(data[5]), int(data[12]), int(data[11]), float(data[16]), float(data[15]), float(data[18]), float(data[17]), float(data[23]), float(data[22]), float(data[13]), float(data[14]),float(data[8]), float(data[10]), float(data[7]), float(data[9])))
                                    MSP_bin[-1]+=1

                            
                        else:
                            deathcut=(float(data[9])**2)*(0.17*10**12)
                            if deathcut<=float(data[7]): 
                                if float(data[9])<=1.: 
                                    f.write('%d %f %f %s %s %f %f %d %d %f %f %f %f %f %f %f %f %e %f %e %f\n'%(pathnos[kk], float(data[1])*t_conv, float(data[19])*l_conv, data[3], data[4], float(data[5]), float(data[6]), int(data[11]), int(data[12]), float(data[15]), float(data[16]), float(data[17]), float(data[18]), float(data[22]), float(data[23]), float(data[13]), float(data[14]),float(data[7]), float(data[9]), float(data[8]), float(data[10])))
                                    MSP_sin[-1]+=1

                thetime_old = thetime
                thetime_code = float(data[1])

        if len(np.sort(glob(paths[kk]+'*.dyn.dat')))>1:
            psrfile2 = paths[kk]+'initial2.morepulsars.dat'

            with open(psrfile2, 'r') as fpsr2:
                next(fpsr2)
                for line in fpsr2:
                    data = line.split()
                    thetime = float(data[1])*t_conv
                    if thetime > thetime_old > timelim:
                        psr_time = thetime_code
                        break
                    if typeflag == 'NSNS':
                        if int(data[2])==1 and timelim <= thetime:
                            if int(data[11])==13 and int(data[12])==13:
                                f.write('%d %f %f %s %s %f %f %d %d %f %f %f %f %f %f %f %f %e %f %e %f\n'%(pathnos[kk], float(data[1])*t_conv, float(data[19])*l_conv, data[3], data[4], float(data[5]), float(data[6]), int(data[11]), int(data[12]), float(data[15]), float(data[16]), float(data[17]), float(data[18]), float(data[22]), float(data[23]), float(data[13]), float(data[14]),float(data[7]), float(data[9]), float(data[8]), float(data[10])))
    
                    if typeflag == 'MSP':
                        if timelim <= thetime:
                            if int(data[2]) == 1:
                                deathcut1=(float(data[9])**2)*(0.17*10**12)
                                deathcut2=(float(data[10])**2)*(0.17*10**12)
                                if deathcut1<=float(data[7]): 
                                    if float(data[9])<=1.: 
                                       f.write('%d %f %f %s %s %f %f %d %d %f %f %f %f %f %f %f %f %e %f %e %f\n'%(pathnos[kk], float(data[1])*t_conv, float(data[19])*l_conv, data[3], data[4], float(data[5]), float(data[6]), int(data[11]), int(data[12]), float(data[15]), float(data[16]), float(data[17]), float(data[18]), float(data[22]), float(data[23]), float(data[13]), float(data[14]),float(data[7]), float(data[9]), float(data[8]), float(data[10])))
                                       MSP_bin[-1]+=1
    
                                elif deathcut2<=float(data[8]): 
                                    if float(data[10])<=1.:
                                        f.write('%d %f %f %s %s %f %f %d %d %f %f %f %f %f %f %f %f %e %f %e %f\n'%(pathnos[kk], float(data[1])*t_conv, float(data[19])*l_conv, data[4], data[3], float(data[6]), float(data[5]), int(data[12]), int(data[11]), float(data[16]), float(data[15]), float(data[18]), float(data[17]), float(data[23]), float(data[22]), float(data[13]), float(data[14]),float(data[8]), float(data[10]), float(data[7]), float(data[9])))
                                        MSP_bin[-1]+=1
    
                                
                            else:
                                deathcut=(float(data[9])**2)*(0.17*10**12)
                                if deathcut<=float(data[7]): 
                                    if float(data[9])<=1.: 
                                        f.write('%d %f %f %s %s %f %f %d %d %f %f %f %f %f %f %f %f %e %f %e %f\n'%(pathnos[kk], float(data[1])*t_conv, float(data[19])*l_conv, data[3], data[4], float(data[5]), float(data[6]), int(data[11]), int(data[12]), float(data[15]), float(data[16]), float(data[17]), float(data[18]), float(data[22]), float(data[23]), float(data[13]), float(data[14]),float(data[7]), float(data[9]), float(data[8]), float(data[10])))
                                        MSP_sin[-1]+=1
    
                    thetime_old = thetime
                    thetime_code = float(data[1])



        ## Extract rc and rh as well
        filedyn=paths[kk]+'initial.dyn.dat'
        check = 0

        with open(filedyn, 'r') as fdyn:
            next(fdyn)
            next(fdyn)
            for line in fdyn:
                datadyn=line.split()
                if float(datadyn[0])>=psr_time: 
                    Mtot.append(float(datadyn[4])*m_conv)
                    Rc.append(float(datadyn[7])*l_conv); Rh.append(float(datadyn[20])*l_conv)
                    check = 1
                    break

        if len(np.sort(glob(paths[kk]+'*.dyn.dat')))>1 and check != 1:
            filedyn2=paths[kk]+'initial2.dyn.dat'

            with open(filedyn2, 'r') as fdyn2:
                next(fdyn2)
                next(fdyn2)
                for line in fdyn2:
                    datadyn=line.split()
                    if float(datadyn[0])>=psr_time: 
                        Mtot.append(float(datadyn[4])*m_conv)
                        Rc.append(float(datadyn[7])*l_conv); Rh.append(float(datadyn[20])*l_conv)
                        break


        Time.append(psr_time)

    print(len(Model), len(Time), len(Mtot), len(Rc), len(Rh), len(MSP_sin), len(MSP_bin))
    
    f.close()

    np.savetxt('/projects/b1095/syr904/projects/isolated_MSP/nmsp1000ms_rc_rh_'+str(timelim)+'myr.dat', np.c_[Model, Time, Mtot, Rc, Rh, Rc_obs, Rh_obs, MSP_sin, MSP_bin], fmt = '%d %f %f %f %f %f %f %d %d', header = '1.Model, 2.Time[code unit], 3.Mtot[Msun], 4.Rc[pc], 5.Rh[pc], 6.Rc_obs[pc], 7.Rhl_obs[pc], 8.MSP_sin, 9.MSP_bin', delimiter = ' ', comments = '#')
        











