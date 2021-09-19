import numpy as np
from glob import glob
import collections
from collections import Counter
import os,sys
import subprocess
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
#import seaborn as sns
import gzip
import math
import re
import random
import history_cmc as hic
import dynamics as dyn
import scripts1
import scripts2
import scripts3
import LISA_calculations as lisa
import ecc_calc as gwcalc
import unit_convert
import single_psr_evolv as psrev

savepath='/projects/b1095/syr904/projects/SGRB'

yearsc=31557600.
twopi=6.283185307179586
Gconst=6.674*10**-8 ##cm3*g-1*s-2
clight=3*10**10 ##cm/s
Msun=2*10**33 ##gram
AU=1.496*10**13  ##cm
PC=3.086*10**18  ##cm
Kconst=9.87*10**-48 ##yr/G^2
Lsun=4.02*10**16 ##mJy*kpc^2


##Colllecting useful models from the set of new models
def find_models():
    allmodels=np.genfromtxt(savepath+'/newruns/finaldata/path_final.dat', dtype=str)
    status=allmodels[:,1]; paths=allmodels[:,0]

    sourcedir=[]; st=[]
    path_nondissolved=[]
    for i in range(len(paths)):
        if status[i]=='1' or status[i]=='2' or status[i]=='3': sourcedir.append(paths[i]); st.append(status[i])
        if status[i]=='1': path_nondissolved.append(paths[i])

    np.savetxt(savepath+'/newruns/finaldata/path_nondissolved_newruns.dat', np.c_[path_nondissolved], fmt='%s')
    fpath=open(savepath+'/newruns/finaldata/path_allfinished_newruns.dat', 'w+', 0)
    fpath.write('#1.Path 2.Status(1-done; 2&3-dissolved)\n')
    for j in range(len(sourcedir)):
        fpath.write('%s %s\n'%(sourcedir[j], st[j]))

    fpath.close()


def find_escmerger_nsns(pathlist, start, end):
    sourcedir=np.genfromtxt(pathlist, dtype=str)

    model_esc=[]; timeesc=[]; timeesc_myr=[]; tins=[]; m0=[]; m1=[]; id0=[]; id1=[]; a=[]; ecc=[]
    for i in range(start, end):
        filestr=sourcedir[i]+'/initial'
        escfile=filestr+'.esc.dat'

        t_conv=dyn.conv('t', filestr+'.conv.sh')

        nesc=0; nescmerger=0

        with open(escfile, 'r') as fesc:
            next(fesc)
            for line in fesc:
                dataesc=line.split()
                if int(dataesc[14])==1:
                    if int(dataesc[22])==13 and int(dataesc[23])==13:
                        nesc+=1
                        t_inspiral=lisa.inspiral_time_peters(float(dataesc[19]), float(dataesc[20]), float(dataesc[15]), float(dataesc[16]))
                        tesc=float(dataesc[1])*t_conv/1000.   ##In Gyr
                        if t_inspiral+tesc<=14.:
                            nescmerger+=1
                            #print int(dataesc[17]), int(dataesc[18]), tesc, t_inspiral
                            model_esc.append(i)
                            timeesc_myr.append(float(dataesc[1])*t_conv); tins.append(t_inspiral*1000.)
                            timeesc.append(float(dataesc[1])); m0.append(float(dataesc[15])); m1.append(float(dataesc[16])); id0.append(int(dataesc[17])); id1.append(int(dataesc[18])); a.append(float(dataesc[19])); ecc.append(float(dataesc[20]))


        print(nesc, nescmerger)     
        #print i

    np.savetxt(savepath+'/extrememodels/Escmerger.dat', np.c_[model_esc, timeesc, timeesc_myr, tins, m0, m1, id0, id1, a, ecc], fmt='%d %f %f %f %f %f %d %d %f %f', header='1.Model 2.Time_esc(code) 3.Time_esc(Myr) 4.T_insp(Myr) 5.M0 6.M1 7.ID0 8.ID1 9.A 10.ECC', delimiter='', comments='#')



def find_formationtime_DNS_new(ids, filestring, flag):
    binintfile=filestring+'.binint.log'
    binintfile2=filestring+'2.binint.log'

    t_conv=dyn.conv('t', filestring+'.conv.sh')

    #snaps=dyn.get_snapshots(filestring)
    #firstsnap=snaps[0]
    sys.path.insert(1, '/projects/b1095/syr904/MyCodes/cmctoolkit')
    import cmctoolkit as cmct

    path = filestring.rsplit('/', 1)[0]
    firstsnap = cmct.Snapshot(fname=path+'/initial.snapshots.h5', 
        snapshot_name='/0(t=0)', conv=path+'/initial.conv.sh', 
                     dist=4.52, # distance to cluster in kpc
                     z=0.0038)
    id0_snap0 = firstsnap.data['id0']; id1_snap0 = firstsnap.data['id1']


    tform=-100; snapno=-100
    tinteract_last=0
    primordial=0; dynamics_eject=0; anydyn=0
    #with gzip.open(firstsnap, 'r') as ffirst:
    #    next(ffirst)
    #    next(ffirst)
    #    for line in ffirst:
    #        datafirst=line.split()
    #        if int(datafirst[7])==1:
    #            if (int(datafirst[10])==ids[0] and int(datafirst[11])==ids[1]) or (int(datafirst[10])==ids[1] and int(datafirst[11])==ids[0]):
    #                primordial=1
    #                break
    for ii in range(len(id0_snap0)):
        if (ids[0] == int(id0_snap0[ii]) and ids[1] == int(id1_snap0[ii])) or (ids[0] == int(id1_snap0[ii]) and ids[1] == int(id0_snap0[ii])):
            primordial = 1
            break


    if os.path.isfile(binintfile2) and os.path.getsize(binintfile2) > 0:
        binint=scripts3.read_binint(binintfile2)
        for i in range(len(binint)-1, -1, -1):
            bininput=binint[i]['input']
            if primordial==1:   ##Check if the binary is dynamically ejected or not

                if tinteract_last>0: break

                for j in range(len(bininput)):
                    if int(bininput[j]['no'])==2:
                        if (int(bininput[j]['ids'][0])==ids[0] and int(bininput[j]['ids'][1])==ids[1]) or (int(bininput[j]['ids'][0])==ids[1] and int(bininput[j]['ids'][1])==ids[0]):
                            anydyn=1
                            list_star=[int(bininput[j]['startype'][0]), int(bininput[j]['startype'][1])]
                            num_star=Counter(list_star)
                            if (flag=='DNS' and num_star[13]==2) or (flag=='NSBH' and num_star[13]==1 and num_star[14]==1) or (flag=='BBH' and num_star[14]==2):
                                dynamics_eject=1
                                tinteract_last=float(binint[i]['type']['time'])
                                break


            binoutput=binint[i]['output']
            if primordial==0:
                anydyn=1

                if tinteract_last>0: break

                for j in range(len(binoutput)):
                    if int(binoutput[j]['no'])==2:
                        if binoutput[j]['ids'][0].find(':')==-1 and binoutput[j]['ids'][1].find(':')==-1:
                            if (int(binoutput[j]['ids'][0])==ids[0] and int(binoutput[j]['ids'][1])==ids[1]) or (int(binoutput[j]['ids'][0])==ids[1] and int(binoutput[j]['ids'][1])==ids[0]):
                                list_star=[int(binoutput[j]['startype'][0]), int(binoutput[j]['startype'][1])]
                                num_star=Counter(list_star)
                                if (flag=='DNS' and num_star[13]==2) or (flag=='NSBH' and num_star[13]==1 and num_star[14]==1) or (flag=='BBH' and num_star[14]==2):
                                    tinteract_last=float(binint[i]['type']['time'])
                                    dynamics_eject=1
                                    break

    if tinteract_last==0:
        binint=scripts3.read_binint(binintfile)
        #fbinint=open(binintfile, 'r')
        #positions=scripts3.find_positions(fbinint)
        #print positions
        #for i in range(len(positions)-1, 0, -1):
        for i in range(len(binint)-1, -1, -1):
            #binint_line=scripts3.read_segment(fbinint,positions[i])
            #print i
            #print binint[i]
            bininput=binint[i]['input']
            if primordial==1:   ##Check if the binary is dynamically ejected or not

                if tinteract_last>0: break

                for j in range(len(bininput)):
                    if int(bininput[j]['no'])==2:
                        if (int(bininput[j]['ids'][0])==ids[0] and int(bininput[j]['ids'][1])==ids[1]) or (int(bininput[j]['ids'][0])==ids[1] and int(bininput[j]['ids'][1])==ids[0]):
                            anydyn=1
                            list_star=[int(bininput[j]['startype'][0]), int(bininput[j]['startype'][1])]
                            num_star=Counter(list_star)
                            if (flag=='DNS' and num_star[13]==2) or (flag=='NSBH' and num_star[13]==1 and num_star[14]==1) or (flag=='BBH' and num_star[14]==2):
                            #if int(bininput[j]['startype'][0])==13 and int(bininput[j]['startype'][1])==13:
                                dynamics_eject=1
                                tinteract_last=float(binint[i]['type']['time'])
                                break


            binoutput=binint[i]['output']
            if primordial==0:
                anydyn=1

                if tinteract_last>0: break

                for j in range(len(binoutput)):
                    if int(binoutput[j]['no'])==2:
                        if binoutput[j]['ids'][0].find(':')==-1 and binoutput[j]['ids'][1].find(':')==-1:
                            if (int(binoutput[j]['ids'][0])==ids[0] and int(binoutput[j]['ids'][1])==ids[1]) or (int(binoutput[j]['ids'][0])==ids[1] and int(binoutput[j]['ids'][1])==ids[0]):
                                list_star=[int(binoutput[j]['startype'][0]), int(binoutput[j]['startype'][1])]
                                num_star=Counter(list_star)
                                if (flag=='DNS' and num_star[13]==2) or (flag=='NSBH' and num_star[13]==1 and num_star[14]==1) or (flag=='BBH' and num_star[14]==2):
                                #if int(binoutput[j]['startype'][0])==13 and int(binoutput[j]['startype'][1])==13:  
                                    tinteract_last=float(binint[i]['type']['time'])
                                    dynamics_eject=1
                                    break


    if tinteract_last>0: tform=tinteract_last

    return tform, snapno, primordial, dynamics_eject, anydyn


def find_formationtime_randombin(ids, filestring):
    binintfile=filestring+'.binint.log'
    binintfile2=filestring+'2.binint.log'

    snaps=dyn.get_snapshots(filestring)
    t_conv=dyn.conv('t', filestring+'.conv.sh')
    firstsnap=snaps[0]

    tform=-100; snapno=-100
    tinteract_last=0
    primordial=0; dynamics_eject=0; anydyn=0
    with gzip.open(firstsnap, 'r') as ffirst:
        for _ in range(2): next(ffirst)
        for line in ffirst:
            datafirst=line.split()
            if int(datafirst[7])==1:
                if (int(datafirst[10])==ids[0] and int(datafirst[11])==ids[1]) or (int(datafirst[10])==ids[1] and int(datafirst[11])==ids[0]):
                    primordial=1
                    break


    if os.path.isfile(binintfile2) and os.path.getsize(binintfile2) > 0:
        binint=scripts3.read_binint(binintfile2)
        for i in range(len(binint)-1, -1, -1):
            bininput=binint[i]['input']
            if primordial==1:   ##Check if the binary is dynamically ejected or not

                if tinteract_last>0: break

                for j in range(len(bininput)):
                    if int(bininput[j]['no'])==2:
                        if (int(bininput[j]['ids'][0])==ids[0] and int(bininput[j]['ids'][1])==ids[1]) or (int(bininput[j]['ids'][0])==ids[1] and int(bininput[j]['ids'][1])==ids[0]):
                            anydyn=1
                            #list_star=[int(bininput[j]['startype'][0]), int(bininput[j]['startype'][1])]
                            #num_star=Counter(list_star)
                            #if (flag=='DNS' and num_star[13]==2) or (flag=='NSBH' and num_star[13]==1 and num_star[14]==1) or (flag=='BBH' and num_star[14]==2):
                            dynamics_eject=1   ##Not trusted here
                            tinteract_last=float(binint[i]['type']['time'])
                            break


            binoutput=binint[i]['output']
            if primordial==0:
                anydyn=1

                if tinteract_last>0: break

                for j in range(len(binoutput)):
                    if int(binoutput[j]['no'])==2:
                        if binoutput[j]['ids'][0].find(':')==-1 and binoutput[j]['ids'][1].find(':')==-1:
                            if (int(binoutput[j]['ids'][0])==ids[0] and int(binoutput[j]['ids'][1])==ids[1]) or (int(binoutput[j]['ids'][0])==ids[1] and int(binoutput[j]['ids'][1])==ids[0]):
                                #list_star=[int(binoutput[j]['startype'][0]), int(binoutput[j]['startype'][1])]
                                #num_star=Counter(list_star)
                                #if (flag=='DNS' and num_star[13]==2) or (flag=='NSBH' and num_star[13]==1 and num_star[14]==1) or (flag=='BBH' and num_star[14]==2):
                                tinteract_last=float(binint[i]['type']['time'])
                                dynamics_eject=1
                                break

    if tinteract_last==0:
        binint=scripts3.read_binint(binintfile)
        for i in range(len(binint)-1, -1, -1):
            bininput=binint[i]['input']
            if primordial==1:   ##Check if the binary is dynamically ejected or not

                if tinteract_last>0: break

                for j in range(len(bininput)):
                    if int(bininput[j]['no'])==2:
                        if (int(bininput[j]['ids'][0])==ids[0] and int(bininput[j]['ids'][1])==ids[1]) or (int(bininput[j]['ids'][0])==ids[1] and int(bininput[j]['ids'][1])==ids[0]):
                            anydyn=1
                            #list_star=[int(bininput[j]['startype'][0]), int(bininput[j]['startype'][1])]
                            #num_star=Counter(list_star)
                            #if (flag=='DNS' and num_star[13]==2) or (flag=='NSBH' and num_star[13]==1 and num_star[14]==1) or (flag=='BBH' and num_star[14]==2):
                            #if int(bininput[j]['startype'][0])==13 and int(bininput[j]['startype'][1])==13:
                            dynamics_eject=1
                            tinteract_last=float(binint[i]['type']['time'])
                            break


            binoutput=binint[i]['output']
            if primordial==0:
                anydyn=1

                if tinteract_last>0: break

                for j in range(len(binoutput)):
                    if int(binoutput[j]['no'])==2:
                        if binoutput[j]['ids'][0].find(':')==-1 and binoutput[j]['ids'][1].find(':')==-1:
                            if (int(binoutput[j]['ids'][0])==ids[0] and int(binoutput[j]['ids'][1])==ids[1]) or (int(binoutput[j]['ids'][0])==ids[1] and int(binoutput[j]['ids'][1])==ids[0]):
                                #list_star=[int(binoutput[j]['startype'][0]), int(binoutput[j]['startype'][1])]
                                #num_star=Counter(list_star)
                                #if (flag=='DNS' and num_star[13]==2) or (flag=='NSBH' and num_star[13]==1 and num_star[14]==1) or (flag=='BBH' and num_star[14]==2):
                                #if int(binoutput[j]['startype'][0])==13 and int(binoutput[j]['startype'][1])==13:  
                                tinteract_last=float(binint[i]['type']['time'])
                                dynamics_eject=1
                                break


    if tinteract_last>0: tform=tinteract_last

    return tform, snapno, primordial, dynamics_eject, anydyn



def find_formationtime_NSBH_new(ids, filestring): ##Find formation time of merging NSBHs
    binintfile=filestring+'.binint.log'
    snaps=np.sort(glob(filestring+'.snap*.dat.gz'))
    t_conv=dyn.conv('t', filestring+'.conv.sh')

    binint=scripts3.read_binint(binintfile)
    firstsnap=snaps[0]

    tform=-100; snapno=-100
    tinteract_last=0
    primordial=0; dynamics_eject=0
    with gzip.open(firstsnap, 'r') as ffirst:
        for _ in xrange(2): next(ffirst)
        for line in ffirst:
            datafirst=line.split()
            if int(datafirst[7])==1:
                if (int(datafirst[10])==ids[0] and int(datafirst[11])==ids[1]) or (int(datafirst[10])==ids[1] and int(datafirst[11])==ids[0]):
                    primordial=1
                    break


    for i in range(len(binint)-1, 0, -1):

        #print i
        #print binint[i]
        bininput=binint[i]['input']
        if primordial==1:   ##Check if the binary is dynamically ejected or not

            if tinteract_last>0: break

            for j in range(len(bininput)):
                if int(bininput[j]['no'])==2:
                    if (int(bininput[j]['ids'][0])==ids[0] and int(bininput[j]['ids'][1])==ids[1]) or (int(bininput[j]['ids'][0])==ids[1] and int(bininput[j]['ids'][1])==ids[0]):
                        if (int(bininput[j]['startype'][0])==13 and int(bininput[j]['startype'][1])==14) or (int(bininput[j]['startype'][0])==14 and int(bininput[j]['startype'][1])==13):
                            dynamics_eject=1
                            #print binint[i]
                            tinteract_last=float(binint[i]['type']['time'])
                            break


        binoutput=binint[i]['output']
        if primordial==0:

            if tinteract_last>0: break

            for j in range(len(binoutput)):
                if int(binoutput[j]['no'])==2:
                    if binoutput[j]['ids'][0].find(':')==-1 and binoutput[j]['ids'][1].find(':')==-1:
                        if (int(binoutput[j]['ids'][0])==ids[0] and int(binoutput[j]['ids'][1])==ids[1]) or (int(binoutput[j]['ids'][0])==ids[1] and int(binoutput[j]['ids'][1])==ids[0]):
                            if (int(binoutput[j]['startype'][0])==13 and int(binoutput[j]['startype'][1])==14) or (int(binoutput[j]['startype'][0])==14 and int(binoutput[j]['startype'][1])==13):  
                                tinteract_last=float(binint[i]['type']['time'])
                                dynamics_eject=1
                                break


    if tinteract_last>0: tform=tinteract_last

    return tform, snapno, primordial, dynamics_eject



def find_formationtime_BBH(ids, filestring):
    binintfile=filestring+'.binint.log'
    collfile=filestring+'.collision.log'
    snaps=dyn.get_snapshots(filestring)
    t_conv=dyn.conv('t', filestring+'.conv.sh')

    binint=scripts3.read_binint(binintfile)
    #collsion=scripts1.readcollfile(collfile)
    firstsnap=snaps[0]

    tform=-100; snapno=-100
    tinteract_last=0
    primordial=0; dynamics_eject=0
    with gzip.open(firstsnap, 'r') as ffirst:
        for _ in xrange(2): next(ffirst)
        for line in ffirst:
            datafirst=line.split()
            if int(datafirst[7])==1:
                if (int(datafirst[10])==ids[0] and int(datafirst[11])==ids[1]) or (int(datafirst[10])==ids[1] and int(datafirst[11])==ids[0]):
                    primordial=1
                    break


    for i in range(len(binint)-1, 0, -1):
        #print i
        #print binint[i]
        bininput=binint[i]['input']
        if primordial==1:   ##Check if the binary is dynamically ejected or not
            if tinteract_last>0: break
            for j in range(len(bininput)):
                if int(bininput[j]['no'])==2:
                    if (int(bininput[j]['ids'][0])==ids[0] and int(bininput[j]['ids'][1])==ids[1]) or (int(bininput[j]['ids'][0])==ids[1] and int(bininput[j]['ids'][1])==ids[0]):
                        if int(bininput[j]['startype'][0])==14 and int(bininput[j]['startype'][1])==14:
                            dynamics_eject=1
                            tinteract_last=float(binint[i]['type']['time'])


        binoutput=binint[i]['output']
        if primordial==0:
            if tinteract_last>0: break
            for j in range(len(binoutput)):
                if int(binoutput[j]['no'])==2:
                    #if binoutput[j]['ids'][0].find(':')!=-1 and binoutput[j]['ids'][1].find(':')!=-1:
                    #   print 'two-collision'
                    #   continue

                    #elif binoutput[j]['ids'][1].find(':')!=-1:
                    #   if (int(binoutput[j]['ids'][0])==ids[0] or int(binoutput[j]['ids'][0])==ids[1]) and int(binoutput[j]['startype'][0])==14:
                    #       #dynamics_eject=1
                    #       print 'collision'

                    #elif binoutput[j]['ids'][0].find(':')!=-1:
                    #   if (int(binoutput[j]['ids'][1])==ids[0] or int(binoutput[j]['ids'][1])==ids[1]) and int(binoutput[j]['startype'][1])==14:
                    #       #dynamics_eject=1
                    #       print 'collision'

                    if binoutput[j]['ids'][0].find(':')==-1 and binoutput[j]['ids'][1].find(':')==-1:
                        if (int(binoutput[j]['ids'][0])==ids[0] and int(binoutput[j]['ids'][1])==ids[1]) or (int(binoutput[j]['ids'][0])==ids[1] and int(binoutput[j]['ids'][1])==ids[0]):
                            if int(binoutput[j]['startype'][0])==14 and int(binoutput[j]['startype'][1])==14:   
                                tinteract_last=float(binint[i]['type']['time'])
                                dynamics_eject=1
                                break


    #if primordial==1:  ##Check the formation time of the primordial binary
    #   for k in range(1, len(snaps)):
    #       tsnap_code=dyn.get_time(snaps[k])
    #       if tsnap_code>time_esc: break

    #       if tform!=-100: break
    #       with gzip.open(snaps[k], 'r') as fsnap:
    #           for _ in xrange(2): next(fsnap)
    #           for line in fsnap:
    #               datasnap=line.split()
    #               if int(datasnap[7])==1:
    #                   if (int(datasnap[10])==ids[0] and int(datasnap[11])==ids[1]) or (int(datasnap[10])==ids[1] and int(datasnap[11])==ids[0]):
    #                       if int(datasnap[17])==14 and int(datasnap[18])==14:
    #                           tform=tsnap_code
    #                           snapno=k

    #                           if dynamics_eject!=0:
    #                               if tinteract_first<tform:
    #                                   tform=tinteract_first
    #                                   snapno=-100
    #                           break

    if tinteract_last>0: tform=tinteract_last

    return tform, snapno, primordial, dynamics_eject




def find_all_mergers_DNS(pathlist, start, end, typeflag):   ##DNSs both in the cluster and escaped
    #sourcedir=np.genfromtxt(pathlist, dtype=str)
    #status = np.ones(len(sourcedir))
    #status=map(int, sourcedir[:,1]); 
    #sourcedir=sourcedir[:,0]
    sourcedir=['/projects/b1095/syr904/cmc/CMC-COSMIC/master_tc_test/ver_0601/MOCHA47Tuc_elson_rv4_3e6_tcon/']
    status=[1]

    ##Mergers in the clusters
    model_coll=[]; model_merge=[]; status_coll=[]; status_merge=[]

    idcoll=[]; id0coll=[]; id1coll=[]; id2coll=[]; id3coll=[]
    idsemerge=[]; id0merge=[]; id1merge=[]
    typecoll=[]; k0=[]; k1=[]; k2=[]; k3=[]
    timecoll=[]; timesemerge=[]
    timecoll_myr=[]; timesemerge_myr=[]
    mf_merge=[]; m0_merge=[]; m1_merge=[]
    mf_coll=[]; m0_coll=[]; m1_coll=[]; m2_coll=[]; m3_coll=[]
    r_merge=[]; r_coll=[]
    tform_semerge=[]; snapno_semerge=[]
    primordial_bin=[]; dynamics_ejection=[]; any_interact=[]

    ##Mergers outside of the clusters
    #model_esc=[]; timeesc=[]; timeesc_myr=[]; tins=[]; m0=[]; m1=[]; id0=[]; id1=[]; a=[]; ecc=[]
    #tform_esc=[]; snapno_esc=[]


    ##Numbers
    Ncoll2=[]; Ncoll3=[]; Ncoll4=[]; Ncoll=[]
    Nsemerge=[]
    Nesc=[]; Nescmerge=[]
    Models=[]


    ####For mergers that happen in the clusters####
    for i in range(start, end):
        filestr=sourcedir[i]+'initial'
        collfile=filestr+'.collision.log'
        binintfile=filestr+'.binint.log'
        semergefile=filestr+'.semergedisrupt.log'
        collfile2=filestr+'2.collision.log'
        binintfile2=filestr+'2.binint.log'
        semergefile2=filestr+'2.semergedisrupt.log'


        t_conv=dyn.conv('t', filestr+'.conv.sh')

        ncoll2=0; ncoll3=0; ncoll4=0
        nsemerge=0
        nesc=0; nescmerger=0

        ##Check in-cluster merger in the collision file
        colldata=scripts1.readcollfile(collfile)
        if os.path.isfile(collfile2) and os.path.getsize(collfile2) > 0:
            colldata2=scripts1.readcollfile(collfile2)
            colldata=colldata+colldata2
        
        for j in range(len(colldata)):
            line=colldata[j].split()
            if int(line[2])==2:  ##Single-single star collision
                starlist=[int(line[11]), int(line[12])]
                starnum=Counter(starlist)
                if (typeflag=='DNS' and starnum[13]>=2) or (typeflag=='NSBH' and starnum[13]>=1 and starnum[14]>=1) or (typeflag=='BBH' and starnum[14]>=2):
                #if int(line[11])==13 and int(line[12])==13:    
                    model_coll.append(i); status_coll.append(status[i])
                    ncoll2+=1
                    idcoll.append(int(line[3])); id0coll.append(int(line[5])); id1coll.append(int(line[7])); id2coll.append(-100); id3coll.append(-100)
                    typecoll.append(int(line[2]))
                    timecoll.append(float(line[0])); timecoll_myr.append(t_conv*float(line[0]))
                    mf_coll.append(float(line[4])); m0_coll.append(float(line[6])); m1_coll.append(float(line[8])); m2_coll.append(-100); m3_coll.append(-100)
                    r_coll.append(float(line[9]))
                    k0.append(int(line[11])); k1.append(int(line[12])); k2.append(-100); k3.append(-100)
                    

            if int(line[2])==3:  ##Binary-single star collision
                if len(line)==16:  ##Three stars collision
                    starlist=[int(line[13]), int(line[14]), int(line[15])]
                    starnum=Counter(starlist)
                    if (typeflag=='DNS' and starnum[13]>=2) or (typeflag=='NSBH' and starnum[13]>=1 and starnum[14]>=1) or (typeflag=='BBH' and starnum[14]>=2):
                    #if (int(line[13])==13 and int(line[14])==13) or (int(line[14])==13 and int(line[15])==13) or (int(line[13])==13 and int(line[15])==13):
                        model_coll.append(i); status_coll.append(status[i])
                        ncoll3+=1
                        idcoll.append(int(line[3])); id0coll.append(int(line[5])); id1coll.append(int(line[7])); id2coll.append(int(line[9])); id3coll.append(-100)
                        typecoll.append(int(line[2]))
                        timecoll.append(float(line[0])); timecoll_myr.append(t_conv*float(line[0]))
                        mf_coll.append(float(line[4])); m0_coll.append(float(line[6])); m1_coll.append(float(line[8])); m2_coll.append(float(line[10])); m3_coll.append(-100)
                        r_coll.append(float(line[11]))
                        k0.append(int(line[13])); k1.append(int(line[14])); k2.append(int(line[15])); k3.append(-100)
                        
                        
                if len(line)==13: ##Two stars collision
                    starlist=[int(line[11]), int(line[12])]
                    starnum=Counter(starlist)
                    if (typeflag=='DNS' and starnum[13]>=2) or (typeflag=='NSBH' and starnum[13]>=1 and starnum[14]>=1) or (typeflag=='BBH' and starnum[14]>=2):
                    #if int(line[11])==13 and int(line[12])==13:
                        model_coll.append(i); status_coll.append(status[i])
                        ncoll3+=1
                        idcoll.append(int(line[3])); id0coll.append(int(line[5])); id1coll.append(int(line[7])); id2coll.append(-100); id3coll.append(-100)
                        typecoll.append(int(line[2]))
                        timecoll.append(float(line[0])); timecoll_myr.append(t_conv*float(line[0]))
                        mf_coll.append(float(line[4])); m0_coll.append(float(line[6])); m1_coll.append(float(line[8])); m2_coll.append(-100); m3_coll.append(-100)
                        r_coll.append(float(line[9]))
                        k0.append(int(line[11])); k1.append(int(line[12])); k2.append(-100); k3.append(-100)


            if int(line[2])==4:  ##Binary-binary star collision
                if len(line)==19: ##Four stars collision
                    starlist=[int(line[15]), int(line[16]), int(line[17]), int(line[18])]
                    starnum=Counter(starlist)
                    if (typeflag=='DNS' and starnum[13]>=2) or (typeflag=='NSBH' and starnum[13]>=1 and starnum[14]>=1) or (typeflag=='BBH' and starnum[14]>=2):
                    #if (int(line[15])==13 and int(line[16])==13) or (int(line[16])==13 and int(line[17])==13) or (int(line[17])==13 and int(line[18])==13) or (int(line[15])==13 and int(line[17])==13) or (int(line[15])==13 and int(line[18])==13) or (int(line[16])==13 and int(line[18])==13):
                        model_coll.append(i); status_coll.append(status[i])
                        ncoll4+=1
                        idcoll.append(int(line[3])); id0coll.append(int(line[5])); id1coll.append(int(line[7])); id2coll.append(int(line[9])); id3coll.append(int(line[11]))
                        typecoll.append(int(line[2]))
                        timecoll.append(float(line[0])); timecoll_myr.append(t_conv*float(line[0]))
                        mf_coll.append(float(line[4])); m0_coll.append(float(line[6])); m1_coll.append(float(line[8])); m2_coll.append(float(line[10])); m3_coll.append(float(line[12]))
                        r_coll.append(float(line[13]))
                        k0.append(int(line[15])); k1.append(int(line[16])); k2.append(int(line[17])); k3.append(int(line[18]))


                if len(line)==16:  ##Three stars collision
                    starlist=[int(line[13]), int(line[14]), int(line[15])]
                    starnum=Counter(starlist)
                    if (typeflag=='DNS' and starnum[13]>=2) or (typeflag=='NSBH' and starnum[13]>=1 and starnum[14]>=1) or (typeflag=='BBH' and starnum[14]>=2):
                    #if (int(line[13])==13 and int(line[14])==13) or (int(line[14])==13 and int(line[15])==13) or (int(line[13])==13 and int(line[15])==13):
                        model_coll.append(i); status_coll.append(status[i])
                        ncoll4+=1
                        idcoll.append(int(line[3])); id0coll.append(int(line[5])); id1coll.append(int(line[7])); id2coll.append(int(line[9])); id3coll.append(-100)
                        typecoll.append(int(line[2]))
                        timecoll.append(float(line[0])); timecoll_myr.append(t_conv*float(line[0]))
                        mf_coll.append(float(line[4])); m0_coll.append(float(line[6])); m1_coll.append(float(line[8])); m2_coll.append(float(line[10])); m3_coll.append(-100)
                        r_coll.append(float(line[11]))
                        k0.append(int(line[13])); k1.append(int(line[14])); k2.append(int(line[15])); k3.append(-100)
                        
                        
                if len(line)==13: ##Two stars collision
                    starlist=[int(line[11]), int(line[12])]
                    starnum=Counter(starlist)
                    if (typeflag=='DNS' and starnum[13]>=2) or (typeflag=='NSBH' and starnum[13]>=1 and starnum[14]>=1) or (typeflag=='BBH' and starnum[14]>=2):
                    #if int(line[11])==13 and int(line[12])==13:
                        model_coll.append(i); status_coll.append(status[i])
                        ncoll4+=1
                        idcoll.append(int(line[3])); id0coll.append(int(line[5])); id1coll.append(int(line[7])); id2coll.append(-100); id3coll.append(-100)
                        typecoll.append(int(line[2]))
                        timecoll.append(float(line[0])); timecoll_myr.append(t_conv*float(line[0]))
                        mf_coll.append(float(line[4])); m0_coll.append(float(line[6])); m1_coll.append(float(line[8])); m2_coll.append(-100); m3_coll.append(-100)
                        r_coll.append(float(line[9]))
                        k0.append(int(line[11])); k1.append(int(line[12])); k2.append(-100); k3.append(-100)


        print('collfile:', ncoll2, ncoll3, ncoll4)#, typecoll, timecoll
        Ncoll2.append(ncoll2); Ncoll3.append(ncoll3); Ncoll4.append(ncoll4)
        Ncoll.append(ncoll2+ncoll3+ncoll4)

    
        ##Check in-cluster merger in the semerge file
        semergedata=scripts2.readmergefile(semergefile)
        if os.path.isfile(semergefile2) and os.path.getsize(semergefile2)>0:
            semergedata2=scripts2.readmergefile(semergefile2)
            semergedata=semergedata+semergedata2

        for k in range(len(semergedata)):
            line=semergedata[k].split()
            if int(line[1])<3:
                starlist=[int(line[-1]), int(line[-2])]
                starnum=Counter(starlist)
                if (typeflag=='DNS' and starnum[13]==2) or (typeflag=='NSBH' and starnum[13]==1 and starnum[14]==1) or (typeflag=='BBH' and starnum[14]==2):
                #if int(line[-1])==13 and int(line[-2])==13:
                    model_merge.append(i); status_merge.append(status[i])
                    nsemerge+=1
                    timesemerge.append(float(line[0])); timesemerge_myr.append(t_conv*float(line[0]))
                    idsemerge.append(int(line[2])); id0merge.append(int(line[4])); id1merge.append(int(line[6]))
                    mf_merge.append(float(line[3])); m0_merge.append(float(line[5])); m1_merge.append(float(line[7]))
                    r_merge.append(float(line[8]))
                    tf_semerge, sno_semerge, pribin, dyn_eject, dyn_any=find_formationtime_DNS_new([int(line[4]), int(line[6])], filestr, typeflag)
                    if tf_semerge!=-100:
                       tform_semerge.append(tf_semerge*t_conv); snapno_semerge.append(sno_semerge)
                    else:
                       tform_semerge.append(t_conv*float(line[0])); snapno_semerge.append(-100)
                    primordial_bin.append(pribin); dynamics_ejection.append(dyn_eject); any_interact.append(dyn_any)

        print('semerge:', nsemerge)#, idsemerge, timesemerge
        Nsemerge.append(nsemerge)


        fesc_merger=open('/projects/b1095/syr904/cmc/CMC-COSMIC/master_tc_test/ver_0601/MOCHA47Tuc_elson_rv4_3e6_tcon/Esc_'+typeflag+'.dat', 'a+')
        #fescbns.write('#1.Model 2.Time_esc(code) 3.Time_esc(Myr) 4.T_insp(Myr) 5.M0 6.M1 7.ID0 8.ID1 9.A 10.ECC 11.Tform(Myr) 12.Primordial? 13.Dynamically_Ejected? 14.Status\n')
        ####For mergers that happen outside of the clusters####
        escfile=filestr+'.esc.dat'
        with open(escfile, 'r') as fesc:
            next(fesc)
            for line in fesc:
                dataesc=line.split()
                if int(dataesc[14])==1:
                    starlist=[int(dataesc[22]), int(dataesc[23])]
                    starnum=Counter(starlist)
                    if (typeflag=='DNS' and starnum[13]==2) or (typeflag=='NSBH' and starnum[13]==1 and starnum[14]==1) or (typeflag=='BBH' and starnum[14]==2):
                    #if int(dataesc[22])==13 and int(dataesc[23])==13:
                        nesc+=1
                        #print float(dataesc[19]), float(dataesc[20]), float(dataesc[15]), float(dataesc[16])
                        if typeflag=='BBH':
                            t_inspiral=lisa.inspiral_time_peters(float(dataesc[19]), float(dataesc[20]), float(dataesc[15]), float(dataesc[16]))*10**3 ##in Myr
                        else:
                            t_inspiral=gwcalc.t_inspiral_2(float(dataesc[19]), float(dataesc[20]), float(dataesc[15]), float(dataesc[16]), 0, 0, 0, 1100)/10**6 ##in Myr
                        tesc=float(dataesc[1])*t_conv  ##In Myr
                        model_esc=i
                        timeesc_myr=float(dataesc[1])*t_conv; tins=t_inspiral
                        timeesc=float(dataesc[1]); m0=float(dataesc[15]); m1=float(dataesc[16]); id0=int(dataesc[17]); id1=int(dataesc[18]); a=float(dataesc[19]); ecc=float(dataesc[20])
                    
                        tf_esc, sno_esc, pribin, dyn_eject, dyn_any=find_formationtime_DNS_new([int(dataesc[17]), int(dataesc[18])], filestr, typeflag)

                        if tf_esc!=-100:
                            tform_esc=tf_esc*t_conv; snapno_esc=sno_esc
                        else:
                            tform_esc=float(dataesc[1])*t_conv; snapno_esc=-100

                        fesc_merger.write('%d %f %f %f %f %f %d %d %f %f %f %d %d %d\n'%(model_esc, timeesc, timeesc_myr, tins, m0, m1, id0, id1, a, ecc, tform_esc, pribin, dyn_eject, status[i]))
                        ##Turn off dyn_any for BBH


                        if t_inspiral+tesc<=14000.:
                            nescmerger+=1
                            #print int(dataesc[17]), int(dataesc[18]), tesc, t_inspiral
                            
                            #model_esc=i
                            #timeesc_myr=float(dataesc[1])*t_conv; tins=t_inspiral
                            #timeesc=float(dataesc[1]); m0=float(dataesc[15]); m1=float(dataesc[16]); id0=int(dataesc[17]); id1=int(dataesc[18]); a=float(dataesc[19]); ecc=float(dataesc[20])
                            #tf_esc, sno_esc=find_formationtime_DNS([int(dataesc[17]), int(dataesc[18])], filestr)
                            #if tf_esc!=-100:
                            #   tform_esc=tf_esc*t_conv; snapno_esc=sno_esc
                            #else:
                            #   tform_esc=float(dataesc[1])*t_conv; snapno_esc=-100

                            #fescbns.write('%d %f %f %f %f %f %d %d %f %f %f %d\n'%(model_esc, timeesc, timeesc_myr, tins, m0, m1, id0, id1, a, ecc, tform_esc, snapno_esc))

        escfile2=filestr+'2.esc.dat'
        if os.path.isfile(escfile2) and os.path.getsize(escfile2) > 0:
            with open(escfile2, 'r') as fesc:
                for line in fesc:
                    dataesc=line.split()
                    if int(dataesc[14])==1:
                        starlist=[int(dataesc[22]), int(dataesc[23])]
                        starnum=Counter(starlist)
                        if (typeflag=='DNS' and starnum[13]==2) or (typeflag=='NSBH' and starnum[13]==1 and starnum[14]==1) or (typeflag=='BBH' and starnum[14]==2):
                        #if int(dataesc[22])==13 and int(dataesc[23])==13:
                            nesc+=1
                            if typeflag=='BBH':
                                t_inspiral=lisa.inspiral_time_peters(float(dataesc[19]), float(dataesc[20]), float(dataesc[15]), float(dataesc[16]))*10**3 ##in Myr
                            else:
                                t_inspiral=gwcalc.t_inspiral_2(float(dataesc[19]), float(dataesc[20]), float(dataesc[15]), float(dataesc[16]), 0, 0, 0, 1100)/10**6 ##in Myr
                            tesc=float(dataesc[1])*t_conv  ##In Myr
                            model_esc=i
                            timeesc_myr=float(dataesc[1])*t_conv; tins=t_inspiral
                            timeesc=float(dataesc[1]); m0=float(dataesc[15]); m1=float(dataesc[16]); id0=int(dataesc[17]); id1=int(dataesc[18]); a=float(dataesc[19]); ecc=float(dataesc[20])

                            tf_esc, sno_esc, pribin, dyn_eject, dyn_any=find_formationtime_DNS_new([int(dataesc[17]), int(dataesc[18])], filestr, typeflag)

                            if tf_esc!=-100:
                                tform_esc=tf_esc*t_conv; snapno_esc=sno_esc
                            else:
                                tform_esc=float(dataesc[1])*t_conv; snapno_esc=-100

                            fesc_merger.write('%d %f %f %f %f %f %d %d %f %f %f %d %d %d\n'%(model_esc, timeesc, timeesc_myr, tins, m0, m1, id0, id1, a, ecc, tform_esc, pribin, dyn_eject, status[i]))
                            ##Turn off dyn_any for BBH


                            if t_inspiral+tesc<=14000.: nescmerger+=1

        fesc_merger.close()
        print('escaped:', nesc, nescmerger)
        Nesc.append(nesc); Nescmerge.append(nescmerger)


        Models.append(i)


    ##Output files
    np.savetxt('/projects/b1095/syr904/cmc/CMC-COSMIC/master_tc_test/ver_0601/MOCHA47Tuc_elson_rv4_3e6_tcon/GWcap_'+typeflag+'.dat', np.c_[model_coll, timecoll, timecoll_myr, typecoll, idcoll, id0coll, id1coll, id2coll, id3coll, mf_coll, m0_coll, m1_coll, m2_coll, m3_coll, r_coll, k0, k1, k2, k3, status_coll], fmt='%d %f %f %d %d %d %d %d %d %f %f %f %f %f %f %d %d %d %d %d', delimiter='', header='1.Model 2.Time(code) 3.Time(Myr) 4.Type 5.IDM 6.ID0 7.ID1 8.ID2 9.ID3 10.MM 11.M0 12.M1 13.M2 14.M3 15.R 16.K0 17.K1 18.K2 19.K3 20.Status(1-done; 2&3-dissolved)', comments='#')

    np.savetxt('/projects/b1095/syr904/cmc/CMC-COSMIC/master_tc_test/ver_0601/MOCHA47Tuc_elson_rv4_3e6_tcon/Incluster_'+typeflag+'.dat', np.c_[model_merge, timesemerge, timesemerge_myr, idsemerge, id0merge, id1merge, mf_merge, m0_merge, m1_merge, r_merge, tform_semerge, primordial_bin, dynamics_ejection, status_merge], fmt='%d %f %f %d %d %d %f %f %f %f %f %d %d %d', delimiter='', header='1.Model 2.Time(code) 3.Time(Myr) 4.IDM 5.ID0 6.ID1 7.MM 8.M0 9.M1 10.R 11.Tform(Myr) 12.Primordial 13.Dynamics? 14.Status(1-done; 2&3-dissolved)', comments='#')
    ##Turn off dyn_any for BBH

    np.savetxt('/projects/b1095/syr904/cmc/CMC-COSMIC/master_tc_test/ver_0601/MOCHA47Tuc_elson_rv4_3e6_tcon/num_merger_'+typeflag+'.dat', np.c_[Models, Ncoll, Ncoll2, Ncoll3, Ncoll4, Nsemerge, Nesc, Nescmerge, status], fmt='%d %d %d %d %d %d %d %d %d', header='1.Model 2.Ncoll 3.Ncoll2 4.Ncoll3 5.Ncoll4 6.Nsemerge 7.Nesc 8.Nescmerge 9.Status(1-done; 2&3-dissolved)', delimiter='', comments='#')



def find_all_mergers_pnmodels_DNS(pathlist, start, end):   
##Both in the cluster and escaped. Code adjusted for Carl's PN models since their are extension files from reruns.
    sourcedir=np.genfromtxt(pathlist, dtype=str)
    pref=sourcedir[:,1]
    sourcedir=sourcedir[:,0]

    ##Mergers in the clusters
    model_coll=[]; model_merge=[]

    idcoll=[]; id0coll=[]; id1coll=[]; id2coll=[]; id3coll=[]
    idsemerge=[]; id0merge=[]; id1merge=[]
    typecoll=[]; k0=[]; k1=[]; k2=[]; k3=[]
    timecoll=[]; timesemerge=[]
    timecoll_myr=[]; timesemerge_myr=[]
    mf_merge=[]; m0_merge=[]; m1_merge=[]
    mf_coll=[]; m0_coll=[]; m1_coll=[]; m2_coll=[]; m3_coll=[]
    r_merge=[]; r_coll=[]
    tform_semerge=[]; snapno_semerge=[]

    ##Mergers outside of the clusters
    model_esc=[]; timeesc=[]; timeesc_myr=[]; tins=[]; m0=[]; m1=[]; id0=[]; id1=[]; a=[]; ecc=[]
    tform_esc=[]; snapno_esc=[]


    ##Numbers
    Ncoll2=[]; Ncoll3=[]; Ncoll4=[]; Ncoll=[]
    Nsemerge=[]
    Nesc=[]; Nescmerge=[]
    Models=[]


    ####For mergers that happen in the clusters####
    for i in range(start, end):
        filestr=sourcedir[i]+pref[i]
        collfile=filestr+'.collision.log'
        collfileext=filestr+'-extension.collision.log'
        binintfile=filestr+'.binint.log'
        binintfileext=filestr+'-extension.binint.log'
        semergefile=filestr+'.semergedisrupt.log'
        semergefileext=filestr+'-extension.semergedisrupt.log'

        t_conv=dyn.conv('t', filestr+'.conv.sh')

        ncoll2=0; ncoll3=0; ncoll4=0
        nsemerge=0
        nesc=0; nescmerger=0

        ##Check in-cluster merger in the collision file
        colldata=scripts1.readcollfile(collfile)
        if os.path.isfile(collfileext) and os.path.getsize(collfileext) > 0:
            colldataext=scripts1.readcollfile(collfileext)
            colldata=colldata+colldataext

        for j in range(len(colldata)):
            line=colldata[j].split()
            if int(line[2])==2:  ##Single-single star collision
                if int(line[11])==13 and int(line[12])==13:
                    model_coll.append(i)
                    ncoll2+=1
                    idcoll.append(int(line[3])); id0coll.append(int(line[5])); id1coll.append(int(line[7])); id2coll.append(-100); id3coll.append(-100)
                    typecoll.append(int(line[2]))
                    timecoll.append(float(line[0])); timecoll_myr.append(t_conv*float(line[0]))
                    mf_coll.append(float(line[4])); m0_coll.append(float(line[6])); m1_coll.append(float(line[8])); m2_coll.append(-100); m3_coll.append(-100)
                    r_coll.append(float(line[9]))
                    k0.append(int(line[11])); k1.append(int(line[12])); k2.append(-100); k3.append(-100)
                    

            if int(line[2])==3:  ##Binary-single star collision
                if len(line)==16:  ##Three stars collision
                    if (int(line[13])==13 and int(line[14])==13) or (int(line[14])==13 and int(line[15])==13) or (int(line[13])==13 and int(line[15])==13):
                        model_coll.append(i)
                        ncoll3+=1
                        idcoll.append(int(line[3])); id0coll.append(int(line[5])); id1coll.append(int(line[7])); id2coll.append(int(line[9])); id3coll.append(-100)
                        typecoll.append(int(line[2]))
                        timecoll.append(float(line[0])); timecoll_myr.append(t_conv*float(line[0]))
                        mf_coll.append(float(line[4])); m0_coll.append(float(line[6])); m1_coll.append(float(line[8])); m2_coll.append(float(line[10])); m3_coll.append(-100)
                        r_coll.append(float(line[11]))
                        k0.append(int(line[13])); k1.append(int(line[14])); k2.append(int(line[15])); k3.append(-100)
                        
                        
                if len(line)==13: ##Two stars collision
                    if int(line[11])==13 and int(line[12])==13:
                        model_coll.append(i)
                        ncoll3+=1
                        idcoll.append(int(line[3])); id0coll.append(int(line[5])); id1coll.append(int(line[7])); id2coll.append(-100); id3coll.append(-100)
                        typecoll.append(int(line[2]))
                        timecoll.append(float(line[0])); timecoll_myr.append(t_conv*float(line[0]))
                        mf_coll.append(float(line[4])); m0_coll.append(float(line[6])); m1_coll.append(float(line[8])); m2_coll.append(-100); m3_coll.append(-100)
                        r_coll.append(float(line[9]))
                        k0.append(int(line[11])); k1.append(int(line[12])); k2.append(-100); k3.append(-100)


            if int(line[2])==4:  ##Binary-binary star collision
                if len(line)==19: ##Four stars collision
                    if (int(line[15])==13 and int(line[16])==13) or (int(line[16])==13 and int(line[17])==13) or (int(line[17])==13 and int(line[18])==13) or (int(line[15])==13 and int(line[17])==13) or (int(line[15])==13 and int(line[18])==13) or (int(line[16])==13 and int(line[18])==13):
                        model_coll.append(i)
                        ncoll4+=1
                        idcoll.append(int(line[3])); id0coll.append(int(line[5])); id1coll.append(int(line[7])); id2coll.append(int(line[9])); id3coll.append(int(line[11]))
                        typecoll.append(int(line[2]))
                        timecoll.append(float(line[0])); timecoll_myr.append(t_conv*float(line[0]))
                        mf_coll.append(float(line[4])); m0_coll.append(float(line[6])); m1_coll.append(float(line[8])); m2_coll.append(float(line[10])); m3_coll.append(float(line[12]))
                        r_coll.append(float(line[13]))
                        k0.append(int(line[15])); k1.append(int(line[16])); k2.append(int(line[17])); k3.append(int(line[18]))


                if len(line)==16:  ##Three stars collision
                    if (int(line[13])==13 and int(line[14])==13) or (int(line[14])==13 and int(line[15])==13) or (int(line[13])==13 and int(line[15])==13):
                        model_coll.append(i)
                        ncoll4+=1
                        idcoll.append(int(line[3])); id0coll.append(int(line[5])); id1coll.append(int(line[7])); id2coll.append(int(line[9])); id3coll.append(-100)
                        typecoll.append(int(line[2]))
                        timecoll.append(float(line[0])); timecoll_myr.append(t_conv*float(line[0]))
                        mf_coll.append(float(line[4])); m0_coll.append(float(line[6])); m1_coll.append(float(line[8])); m2_coll.append(float(line[10])); m3_coll.append(-100)
                        r_coll.append(float(line[11]))
                        k0.append(int(line[13])); k1.append(int(line[14])); k2.append(int(line[15])); k3.append(-100)
                        
                        
                if len(line)==13: ##Two stars collision
                    if int(line[11])==13 and int(line[12])==13:
                        model_coll.append(i)
                        ncoll4+=1
                        idcoll.append(int(line[3])); id0coll.append(int(line[5])); id1coll.append(int(line[7])); id2coll.append(-100); id3coll.append(-100)
                        typecoll.append(int(line[2]))
                        timecoll.append(float(line[0])); timecoll_myr.append(t_conv*float(line[0]))
                        mf_coll.append(float(line[4])); m0_coll.append(float(line[6])); m1_coll.append(float(line[8])); m2_coll.append(-100); m3_coll.append(-100)
                        r_coll.append(float(line[9]))
                        k0.append(int(line[11])); k1.append(int(line[12])); k2.append(-100); k3.append(-100)


        print('collfile:', ncoll2, ncoll3, ncoll4)#, typecoll, timecoll
        Ncoll2.append(ncoll2); Ncoll3.append(ncoll3); Ncoll4.append(ncoll4)
        Ncoll.append(ncoll2+ncoll3+ncoll4)

    
        ##Check in-cluster merger in the semerge file
        semergedata=scripts2.readmergefile(semergefile)
        if os.path.isfile(semergefileext) and os.path.getsize(semergefileext)>0:
            semergedataext=scripts2.readmergefile(semergefileext)
            semergedata=semergedata+semergedataext

        for k in range(len(semergedata)):
            line=semergedata[k].split()
            if int(line[1])<3:
                if int(line[-1])==13 and int(line[-2])==13:
                    model_merge.append(i)
                    nsemerge+=1
                    timesemerge.append(float(line[0])); timesemerge_myr.append(t_conv*float(line[0]))
                    idsemerge.append(int(line[2])); id0merge.append(int(line[4])); id1merge.append(int(line[6]))
                    mf_merge.append(float(line[3])); m0_merge.append(float(line[5])); m1_merge.append(float(line[7]))
                    r_merge.append(float(line[8]))
                    tf_semerge, sno_semerge=find_formationtime([int(line[4]), int(line[6])], sourcedir[i])
                    if tf_semerge!=-100:
                        tform_semerge.append(tf_semerge*t_conv); snapno_semerge.append(sno_semerge)
                    else:
                        tform_semerge.append(t_conv*float(line[0])); snapno_semerge.append(-100)


        print('semerge:', nsemerge)#, idsemerge, timesemerge
        Nsemerge.append(nsemerge)


        ####For mergers that happen outside of the clusters####
        escfile=filestr+'.esc.dat'
        with open(escfile, 'r') as fesc:
            next(fesc)
            for line in fesc:
                dataesc=line.split()
                if int(dataesc[14])==1:
                    if int(dataesc[22])==13 and int(dataesc[23])==13:
                        nesc+=1
                        #print float(dataesc[19]), float(dataesc[20]), float(dataesc[15]), float(dataesc[16])
                        t_inspiral=gwcalc.t_inspiral_2(float(dataesc[19]), float(dataesc[20]), float(dataesc[15]), float(dataesc[16]), 0, 0, 0, 1100)/10**6 ##in Myr
                        tesc=float(dataesc[1])*t_conv  ##In Myr
                        model_esc.append(i)
                        timeesc_myr.append(float(dataesc[1])*t_conv); tins.append(t_inspiral)
                        timeesc.append(float(dataesc[1])); m0.append(float(dataesc[15])); m1.append(float(dataesc[16])); id0.append(int(dataesc[17])); id1.append(int(dataesc[18])); a.append(float(dataesc[19])); ecc.append(float(dataesc[20]))
                        tf_esc, sno_esc=find_formationtime([int(dataesc[17]), int(dataesc[18])], filestr)
                        if tf_esc!=-100:
                            tform_esc.append(tf_esc*t_conv); snapno_esc.append(sno_esc)
                        else:
                            tform_esc.append(float(dataesc[1])*t_conv); snapno_esc.append(-100)

                        if t_inspiral+tesc<=14000.:
                            nescmerger+=1
                            #print int(dataesc[17]), int(dataesc[18]), tesc, t_inspiral
                            #model_esc.append(i)
                            #timeesc_myr.append(float(dataesc[1])*t_conv); tins.append(t_inspiral)
                            #timeesc.append(float(dataesc[1])); m0.append(float(dataesc[15])); m1.append(float(dataesc[16])); id0.append(int(dataesc[17])); id1.append(int(dataesc[18])); a.append(float(dataesc[19])); ecc.append(float(dataesc[20]))
                            #tf_esc, sno_esc=find_formationtime([int(dataesc[17]), int(dataesc[18])], filestr)
                            #if tf_esc!=-100:
                            #   tform_esc.append(tf_esc*t_conv); snapno_esc.append(sno_esc)
                            #else:
                            #   tform_esc.append(float(dataesc[1])*t_conv); snapno_esc.append(-100)

        escfileext=filestr+'-extension.esc.dat'
        if os.path.isfile(escfileext) and os.path.getsize(escfileext) > 0:
            print('model=', i)
            with open(escfileext, 'r') as fesc:
                next(fesc)
                for line in fesc:
                    dataesc=line.split()
                    if int(dataesc[14])==1:
                        if int(dataesc[22])==13 and int(dataesc[23])==13:
                            nesc+=1
                            #print float(dataesc[19]), float(dataesc[20]), float(dataesc[15]), float(dataesc[16])
                            t_inspiral=gwcalc.t_inspiral_2(float(dataesc[19]), float(dataesc[20]), float(dataesc[15]), float(dataesc[16]), 0, 0, 0, 1100)/10**6 ##in Myr
                            tesc=float(dataesc[1])*t_conv  ##In Myr
                            model_esc.append(i)
                            timeesc_myr.append(float(dataesc[1])*t_conv); tins.append(t_inspiral)
                            timeesc.append(float(dataesc[1])); m0.append(float(dataesc[15])); m1.append(float(dataesc[16])); id0.append(int(dataesc[17])); id1.append(int(dataesc[18])); a.append(float(dataesc[19])); ecc.append(float(dataesc[20]))
                            tf_esc, sno_esc=find_formationtime([int(dataesc[17]), int(dataesc[18])], filestr)
                            if tf_esc!=-100:
                                tform_esc.append(tf_esc*t_conv); snapno_esc.append(sno_esc)
                            else:
                                tform_esc.append(float(dataesc[1])*t_conv); snapno_esc.append(-100)

                            if t_inspiral+tesc<=14000.:
                                nescmerger+=1
                                #print int(dataesc[17]), int(dataesc[18]), tesc, t_inspiral
                                #model_esc.append(i)
                                #timeesc_myr.append(float(dataesc[1])*t_conv); tins.append(t_inspiral)
                                #timeesc.append(float(dataesc[1])); m0.append(float(dataesc[15])); m1.append(float(dataesc[16])); id0.append(int(dataesc[17])); id1.append(int(dataesc[18])); a.append(float(dataesc[19])); ecc.append(float(dataesc[20]))
                                #tf_esc, sno_esc=find_formationtime([int(dataesc[17]), int(dataesc[18])], sourcedir[i])
                                #if tf_esc!=-100:
                                #   tform_esc.append(tf_esc*t_conv); snapno_esc.append(sno_esc)
                                #else:
                                #   tform_esc.append(float(dataesc[1])*t_conv); snapno_esc.append(-100) 


        print('escaped:', nesc, nescmerger)
        Nesc.append(nesc); Nescmerge.append(nescmerger)


        Models.append(i)


    ##Output files
    np.savetxt(savepath+'/pnmodels/GWcap.dat', np.c_[model_coll, timecoll, timecoll_myr, typecoll, idcoll, id0coll, id1coll, id2coll, id3coll, mf_coll, m0_coll, m1_coll, m2_coll, m3_coll, r_coll, k0, k1, k2, k3], fmt='%d %f %f %d %d %d %d %d %d %f %f %f %f %f %f %d %d %d %d', delimiter='', header='1.Model 2.Time(code) 3.Time(Myr) 4.Type 5.IDM 6.ID0 7.ID1 8.ID2 9.ID3 10.MM 11.M0 12.M1 13.M2 14.M3 15.R 16.K0 17.K1 18.K2 19.K3', comments='#')

    np.savetxt(savepath+'/pnmodels/Incluster.dat', np.c_[model_merge, timesemerge, timesemerge_myr, idsemerge, id0merge, id1merge, mf_merge, m0_merge, m1_merge, r_merge, tform_semerge, snapno_semerge], fmt='%d %f %f %d %d %d %f %f %f %f %f %d', delimiter='', header='1.Model 2.Time(code) 3.Time(Myr) 4.IDM 5.ID0 6.ID1 7.MM 8.M0 9.M1 10.R 11.Tform(Myr) 12.Snapno', comments='#')

    #np.savetxt(savepath+'/pnmodels/Escmerger.dat', np.c_[model_esc, timeesc, timeesc_myr, tins, m0, m1, id0, id1, a, ecc, tform_esc, snapno_esc], fmt='%d %f %f %f %f %f %d %d %f %f %f %d', header='1.Model 2.Time_esc(code) 3.Time_esc(Myr) 4.T_insp(Myr) 5.M0 6.M1 7.ID0 8.ID1 9.A 10.ECC 11.Tform(Myr) 12.Snapno', delimiter='', comments='#')

    np.savetxt(savepath+'/pnmodels/num_bnsmerger.dat', np.c_[Models, Ncoll, Ncoll2, Ncoll3, Ncoll4, Nsemerge, Nesc, Nescmerge], fmt='%d %d %d %d %d %d %d %d', header='1.Model 2.Ncoll 3.Ncoll2 4.Ncoll3 5.Ncoll4 6.Nsemerge 7.Nesc 8.Nescmerge', delimiter='', comments='#')

    np.savetxt(savepath+'/pnmodels/Escbns.dat', np.c_[model_esc, timeesc, timeesc_myr, tins, m0, m1, id0, id1, a, ecc, tform_esc, snapno_esc], fmt='%d %f %f %f %f %f %d %d %f %f %f %d', header='1.Model 2.Time_esc(code) 3.Time_esc(Myr) 4.T_insp(Myr) 5.M0 6.M1 7.ID0 8.ID1 9.A 10.ECC 11.Tform(Myr) 12.Snapno', delimiter='', comments='#')



def find_formationchannel_escbns(pathlist, escbnsfile, start, end):
    databns=np.genfromtxt(escbnsfile)
    id0=databns[:,6]; id1=databns[:,7]; model=databns[:,0]
    sourcedir=np.genfromtxt(pathlist, dtype=str)
    ffm=open('/projects/b1011/syr904/projects/SGRB/extrememodels/fmchannel.dat', 'a+', 0)
    ffm.write('#1.Model 2.ID0 3.ID1 4.Primordial 5.Dynamics\n')
    ##6.ExcNS 7.ExcNonNS 8.CollNS 9.EvolNS 10.DynEject 11.SNEject
    for i in range(start, end):
        filestr=sourcedir[i]+'/initial'
        modelno=i
        for j in range(len(model)):
            if model[j]==modelno:
                firstsnap=filestr+'.snap0000.dat.gz'
                dyn=1; pri=0
                with gzip.open(firstsnap, 'r') as ffirst:
                    next(ffirst)
                    next(ffirst)
                    for line in ffirst:
                        datafirst=line.split()
                        if (int(datafirst[10])==id0[j] and int(datafirst[11])==id1[j]) or (int(datafirst[11])==id0[j] and int(datafirst[10])==id1[j]):
                            pri=1; dyn=0
                            break
                
                #if dyn==1: 
                        #hdict=hic.history_maker([long(id0[j])], [1], 'initial', sourcedir[i], 1)
                        ##Can add extra checking of histories here


                ffm.write('%d %d %d %d %d\n'%(modelno, id0[j], id1[j], pri, dyn))

        print(i)

    ffm.close()



##Calculate the maximum event rate from the "best" model
def find_eventrate_maximum(pathlist, totaldraw, deltat):  ##deltat in Myr
    sourcedir=np.genfromtxt(pathlist, dtype=str)
    gwfile=np.genfromtxt('/projects/b1011/syr904/projects/SGRB/extrememodels/GWcap.dat')
    inclusterfile=np.genfromtxt('/projects/b1011/syr904/projects/SGRB/extrememodels/Incluster.dat')
    escbnsfile=np.genfromtxt('/projects/b1011/syr904/projects/SGRB/extrememodels/Escbns_new.dat')
    modelgw=gwfile[:,0]; tmggw=gwfile[:,2] ##All times in Myr
    modelin=inclusterfile[:,0]; tmgin=inclusterfile[:,2]
    modelesc=escbnsfile[:,0]; tesc=escbnsfile[:,2]; tinsp=escbnsfile[:,3]
    tmgesc=np.array(tesc)+np.array(tinsp)

    agedist=np.genfromtxt('/projects/b1011/syr904/projects/SGRB/AgeDistribution/metallicity_0.001.dat')
    agelist=agedist[:,0]

    ##Maximum event rate
    model_maxbns=sourcedir[2]
    filestr_max=model_maxbns+'/initial'
    collmax=filestr_max+'.collision.log'
    semergemax=filestr_max+'.semergedisrupt.log'
    escmax=filestr_max+'.esc.dat'

    draw_coll=[]; draw_semerge=[]; draw_esc=[]
    ndraw=0
    i=0
    while i < totaldraw:
        ncoll=0; nsemerge=0; nesc=0
        ageno=random.randint(0, len(agelist)-1)
        agemin=agelist[ageno]*1000. ##In Myr
        agemax=agemin+deltat  ##in Myr
        #print i
        if agemax>12000.: continue

        i+=1
        ndraw+=1
        #print agemin, agemax
        for j in range(len(modelgw)):
            if modelgw[j]==2:
                if tmggw[j]>=agemin and tmggw[j]<=agemax:
                    ncoll+=1
        for k in range(len(modelin)):
            if modelin[k]==2:
                if tmgin[k]>=agemin and tmgin[k]<=agemax:
                    nsemerge+=1
        for h in range(len(modelesc)):
            if modelesc[h]==2:
                if tmgesc[h]>=agemin and tmgesc[h]<=agemax:
                    nesc+=1

        draw_coll.append(ncoll); draw_semerge.append(nsemerge); draw_esc.append(nesc)

    print(np.mean(draw_coll), np.mean(draw_semerge), np.mean(draw_esc))
    print(ndraw)

    mean_coll=np.mean(draw_coll); mean_semerge=np.mean(draw_semerge); mean_esc=np.mean(draw_esc)

    ntot=mean_coll+mean_semerge+mean_esc

    est1=ntot*0.3*10**9/(deltat*10**6)  ##Gpc^-3 yr^-1
    est2=ntot*0.7*10**9/(deltat*10**6)
    est3=ntot*2.1*10**9/(deltat*10**6)

    print(est1, est2, est3)



##Calculate the average event rate from a small set of models with no weighting on the models
def find_eventrate_average(pathlist, totaldraw, deltat, modeldraw):   ##deltat in Myr
    sourcedir=np.genfromtxt(pathlist, dtype=str)
    pref=sourcedir[:,1]; sourcedir=sourcedir[:,0]
    gwfile=np.genfromtxt('/projects/b1011/syr904/projects/SGRB/extrememodels/GWcap.dat')
    inclusterfile=np.genfromtxt('/projects/b1011/syr904/projects/SGRB/extrememodels/Incluster.dat')
    escbnsfile=np.genfromtxt('/projects/b1011/syr904/projects/SGRB/extrememodels/Escbns.dat')
    modelgw=gwfile[:,0]; tmggw=gwfile[:,2] ##All times in Myr
    modelin=inclusterfile[:,0]; tmgin=inclusterfile[:,2]
    modelesc=escbnsfile[:,0]; tesc=escbnsfile[:,2]; tinsp=escbnsfile[:,3]
    tmgesc=np.array(tesc)+np.array(tinsp)

    agedist=np.genfromtxt('/projects/b1011/syr904/projects/SGRB/AgeDistribution/metallicity_0.001.dat')
    agelist=agedist[:,0]

    mean_coll=[]; mean_semerge=[]; mean_esc=[]

    n=0; nn=0
    while n < modeldraw:
        modelno=random.randint(0, len(sourcedir)-1)
        model=sourcedir[modelno]
        filestr=model+pref[modelno]

        draw_coll=[]; draw_semerge=[]; draw_esc=[]
        ni=0

        n+=1; nn+=1
        i=0
        while i < totaldraw:
            ncoll=0; nsemerge=0; nesc=0
            ageno=random.randint(0, len(agelist)-1)
            agemin=agelist[ageno]*1000. ##In Myr
            agemax=agemin+deltat  ##in Myr
            #print i
            if agemax>12000.: continue

            i+=1
            ni+=1
            #print agemin, agemax
            if len(modelgw)>0:
                for j in range(len(modelgw)):
                    if modelgw[j]==modelno:
                        if tmggw[j]>=agemin and tmggw[j]<=agemax:
                            ncoll+=1.

            if len(modelin)>0:          
                for k in range(len(modelin)):
                    if modelin[k]==modelno:
                        if tmgin[k]>=agemin and tmgin[k]<=agemax:
                            nsemerge+=1.

            if len(modelesc)>0: 
                for h in range(len(modelesc)):
                    if modelesc[h]==modelno:
                        if tmgesc[h]>=agemin and tmgesc[h]<=agemax:
                            nesc+=1.

            draw_coll.append(ncoll); draw_semerge.append(nsemerge); draw_esc.append(nesc)

        #print np.mean(draw_coll), np.mean(draw_semerge), np.mean(draw_esc)
        #print ni

        mean_coll.append(np.mean(draw_coll)); mean_semerge.append(np.mean(draw_semerge)); mean_esc.append(np.mean(draw_esc))

        print(n)

    #print nn


    totmean_coll=np.mean(mean_coll); totmean_semerge=np.mean(mean_semerge); totmean_esc=np.mean(mean_esc)
    print(totmean_coll, totmean_semerge, totmean_esc)
    ntot=totmean_coll+totmean_semerge+totmean_esc

    est1=ntot*0.3*10**9/(deltat*10**6)  ##Gpc^-3 yr^-1
    est2=ntot*0.7*10**9/(deltat*10**6)  ##Gpc^-3 yr^-1
    est3=ntot*2.1*10**9/(deltat*10**6)  ##Gpc^-3 yr^-1

    print(est1, est2, est3)
    



##Calculate the average event rate from new runs with no weighting on the models
def find_eventrate_average_differentZ(pathlist, totaldraw, deltat, modeldraw):   ##deltat in Myr
    sourcedir=np.genfromtxt(pathlist, dtype=str)
    filepath=sourcedir[:,0]; status=sourcedir[:,1]
    pref='initial'

    ##Read the necessary files
    gwfile=np.genfromtxt('/projects/b1095/syr904/projects/SGRB/newruns/GWcap_DNS.dat')
    #print gwfile.shape
    #gwfile=np.atleast_2d(gwfile)   ##Make it a normal list since this file has only one row.
    #print len(gwfile[0,:])
    inclusterfile=np.genfromtxt('/projects/b1095/syr904/projects/SGRB/newruns/Incluster_DNS.dat')
    #inclusterfile=np.atleast_2d(inclusterfile)
    #print len(inclusterfile[0,:])
    escfile=np.genfromtxt('/projects/b1095/syr904/projects/SGRB/newruns/Esc_DNS.dat')


    ##Read the age distribution files
    agedist_z1=np.genfromtxt('/projects/b1095/syr904/projects/SGRB/AgeDistribution/bin2/Metallicity_0.0002.dat')
    agedist_z2=np.genfromtxt('/projects/b1095/syr904/projects/SGRB/AgeDistribution/bin2/Metallicity_0.002.dat')
    agedist_z3=np.genfromtxt('/projects/b1095/syr904/projects/SGRB/AgeDistribution/bin2/Metallicity_0.02.dat')


    ##Read metallicity
    modelinfo=np.genfromtxt('/projects/b1095/syr904/projects/SGRB/newruns/table_allnums.txt')
    zmetal=modelinfo[:,3]

    if len(gwfile[0,:])!=0:
        modelgw=gwfile[:,0]; tmggw=gwfile[:,2] ##All times in Myr
    else:
        modelgw=[]; tmggw=[]

    if len(inclusterfile[0,:])!=0:
        modelin=inclusterfile[:,0]; tmgin=inclusterfile[:,2]
    else:
        modelin=[]; tmgin=[]

    if len(escfile[0,:])!=0:
        modelesc=escfile[:,0]; tesc=escfile[:,2]; tinsp=escfile[:,3]
        tmgesc=np.array(tesc)+np.array(tinsp)
    else:
        modelesc=[]; tmgesc=[]


    age_z1=agedist_z1[:,0]; age_z2=agedist_z2[:,0]; age_z3=agedist_z3[:,0]

    mean_coll=[]; mean_semerge=[]; mean_esc=[]

    n=0; nn=0
    while n < modeldraw:
        modelno=random.randint(0, len(sourcedir)-1)
        model=sourcedir[modelno]
        filestr=model+pref


        draw_coll=[]; draw_semerge=[]; draw_esc=[]
        ni=0

        n+=1; nn+=1
        i=0
        while i < totaldraw:
            ncoll=0; nsemerge=0; nesc=0
            if zmetal[modelno]<=0.00065:
                ageno=random.randint(0, len(age_z1)-1)
                agemin=age_z1[ageno]*1000. ##In Myr
                agemax=agemin+deltat  ##in Myr
            elif zmetal[modelno]<=0.0065:
                ageno=random.randint(0, len(age_z2)-1)
                agemin=age_z2[ageno]*1000. ##In Myr
                agemax=agemin+deltat  ##in Myr
            else:
                ageno=random.randint(0, len(age_z3)-1)
                agemin=age_z3[ageno]*1000. ##In Myr
                agemax=agemin+deltat  ##in Myr

            #print i
            if agemax>12000.: continue

            i+=1
            ni+=1
            #print agemin, agemax
            if len(modelgw)>0:
                for j in range(len(modelgw)):
                    if modelgw[j]==modelno:
                        if tmggw[j]>=agemin and tmggw[j]<=agemax:
                            ncoll+=1.

            if len(modelin)>0:          
                for k in range(len(modelin)):
                    if modelin[k]==modelno:
                        if tmgin[k]>=agemin and tmgin[k]<=agemax:
                            nsemerge+=1.

            if len(modelesc)>0: 
                for h in range(len(modelesc)):
                    if modelesc[h]==modelno:
                        if tmgesc[h]>=agemin and tmgesc[h]<=agemax:
                            nesc+=1.

            draw_coll.append(ncoll); draw_semerge.append(nsemerge); draw_esc.append(nesc)

        #print np.mean(draw_coll), np.mean(draw_semerge), np.mean(draw_esc)
        #print ni

        mean_coll.append(np.mean(draw_coll)); mean_semerge.append(np.mean(draw_semerge)); mean_esc.append(np.mean(draw_esc))

        print(n)

    #print nn


    totmean_coll=np.mean(mean_coll); totmean_semerge=np.mean(mean_semerge); totmean_esc=np.mean(mean_esc)
    print(totmean_coll, totmean_semerge, totmean_esc)
    ntot=totmean_coll+totmean_semerge+totmean_esc

    est1=ntot*0.33*10**9/(deltat*10**6)  ##Gpc^-3 yr^-1
    est2=ntot*0.77*10**9/(deltat*10**6)  ##Gpc^-3 yr^-1
    est3=ntot*2.31*10**9/(deltat*10**6)  ##Gpc^-3 yr^-1

    print(est1, est2, est3)


def find_number_mergers():
    esc_dns=np.genfromtxt('/projects/b1095/syr904/projects/SGRB/newruns/Esc_DNS.dat')
    esc_nsbh=np.genfromtxt('/projects/b1095/syr904/projects/SGRB/newruns/Esc_NSBH.dat')
    in_dns=np.genfromtxt('/projects/b1095/syr904/projects/SGRB/newruns/Incluster_DNS.dat')
    in_nsbh=np.genfromtxt('/projects/b1095/syr904/projects/SGRB/newruns/Incluster_NSBH.dat')
    gw_dns=np.genfromtxt('/projects/b1095/syr904/projects/SGRB/newruns/GWcap_DNS.dat')
    gw_nsbh=np.genfromtxt('/projects/b1095/syr904/projects/SGRB/newruns/GWcap_NSBH.dat')
    
    paths=np.genfromtxt('/projects/b1095/syr904/projects/SGRB/newruns/path_allfinished_newruns.dat')
    st=paths[:,1]

    model_esc_dns=list(esc_dns[:,0]); nesc_dns=len(model_esc_dns)
    model_esc_nsbh=list(esc_nsbh[:,0]); nesc_nsbh=len(model_esc_nsbh)

    model_in_dns=list(in_dns[:,0]); nin_dns=len(model_in_dns)
    model_in_nsbh=list(in_nsbh[:,0]); nin_nsbh=len(model_in_nsbh)

    model_gw_dns=list(gw_dns[:,0]); ngw_dns=len(model_gw_dns)
    model_gw_nsbh=list(gw_nsbh[:,0]); ngw_nsbh=len(model_gw_nsbh)

    ndns=nesc_dns+nin_dns+ngw_dns
    nnsbh=nesc_nsbh+nin_nsbh+ngw_nsbh
    model_dns=model_esc_dns+model_in_dns+model_gw_dns
    model_nsbh=model_esc_nsbh+model_in_nsbh+model_gw_nsbh

    print(ndns, nnsbh)

    ndone=0; ndissolved=0
    for i in range(len(st)):
        if st[i]==2: ndone+=1
        else:
            ndissolved+=1
            print(i)

            for j in range(len(model_dns)):
                if int(model_dns[j])==i:
                    ndns-=1

                    print('dns', i)

            for k in range(len(model_nsbh)):
                if int(model_nsbh[k])==i:
                    nnsbh-=1

                    print('nsbh', i)

    rdns=0.7*float(ndns)/float(ndone)/12.
    rnsbh=0.7*float(nnsbh)/float(ndone)/12.

    print('done', ndone, 'dns', ndns, 'nsbh', nnsbh)
    print(rdns, rnsbh)



def find_DNS_NSBH(pathlist, start, end, timestart, timeend):  ##times in Myr
    sourcedir=np.genfromtxt(pathlist, dtype=str)
    Status=sourcedir[:,1]; sourcedir=sourcedir[:,0]
    pref='initial'
    Model=[]; Time=[]; ID0=[]; ID1=[]; P0=[]; P1=[]; B0=[]; B1=[]; M0=[]; M1=[]; A=[]; E=[]; Type=[]; st=[]
    for i in range(start, end):
        filestr=sourcedir[i]+pref
        t_conv=dyn.conv('t', filestr+'.conv.sh')
        snaps=dyn.get_snapshots(filestr)

        t_last=dyn.get_time(snaps[-1])*t_conv

        if int(Status[i])!=1:  
            ##For dissolved models
            with gzip.open(snaps[-1], 'r') as fsnap:
                next(fsnap)
                next(fsnap)
                for line in fsnap:
                    datasnap=line.split()
                    if int(datasnap[17])==13 and int(datasnap[18])==13:
                        Model.append(i); Time.append(t_last); ID0.append(int(datasnap[10])); ID1.append(int(datasnap[11])); P0.append(twopi*yearsc/float(datasnap[45])); P1.append(twopi*yearsc/float(datasnap[46])); B0.append(float(datasnap[47])); B1.append(float(datasnap[48])); M0.append(float(datasnap[8])); M1.append(float(datasnap[9])); A.append(float(datasnap[12])); E.append(float(datasnap[13])); Type.append(1313); st.append(int(Status[i]))
                    if int(datasnap[17])==13 and int(datasnap[18])==14:
                        Model.append(i); Time.append(t_last); ID0.append(int(datasnap[10])); ID1.append(int(datasnap[11])); P0.append(twopi*yearsc/float(datasnap[45])); P1.append(-100); B0.append(float(datasnap[47])); B1.append(-100); M0.append(float(datasnap[8])); M1.append(float(datasnap[9])); A.append(float(datasnap[12])); E.append(float(datasnap[13])); Type.append(1314); st.append(int(Status[i]))
                    if int(datasnap[17])==14 and int(datasnap[18])==13:
                        Model.append(i); Time.append(t_last); ID0.append(int(datasnap[11])); ID1.append(int(datasnap[10])); P0.append(twopi*yearsc/float(datasnap[46])); P1.append(-100); B0.append(float(datasnap[48])); B1.append(-100); M0.append(float(datasnap[9])); M1.append(float(datasnap[8])); A.append(float(datasnap[12])); E.append(float(datasnap[13])); Type.append(1314); st.append(int(Status[i]))

        else:
            ##Non-dissolved models
            with open(filestr+'.ns.dat', 'r') as fns:
                next(fns)
                for line in fns:
                    datans=line.split()
                    if float(datans[0])*t_conv>=timestart and float(datans[0])*t_conv<=timeend:
                        if int(datans[7])>0:
                            model=i; t=float(datans[0])
                            print(model, t, 'DNS')
                            for j in range(len(snaps)-1, -1, -1):
                                time=dyn.get_time(snaps[j])
                                if round(time, 6)==round(t, 6):
                                    print(model, t, 'DNS')
                                    with gzip.open(snaps[j], 'r') as fsnap:
                                        next(fsnap)
                                        next(fsnap)
                                        for line in fsnap:
                                            datasnap=line.split()
                                            if int(datasnap[17])==13 and int(datasnap[18])==13:
                                                Model.append(i); Time.append(time*t_conv); ID0.append(int(datasnap[10])); ID1.append(int(datasnap[11])); P0.append(twopi*yearsc/float(datasnap[45])); P1.append(twopi*yearsc/float(datasnap[46])); B0.append(float(datasnap[47])); B1.append(float(datasnap[48])); M0.append(float(datasnap[8])); M1.append(float(datasnap[9])); A.append(float(datasnap[12])); E.append(float(datasnap[13])); Type.append(1313); st.append(int(Status[i])) 



                        if int(datans[8])>0: 
                            model=i; t=float(datans[0])
                            print(model, t, 'NSBH')
                            for j in range(len(snaps)-1, -1, -1):
                                time=dyn.get_time(snaps[j])
                                if round(time, 6)==round(t, 6):
                                    print(model, t, 'NSBH')
                                    with gzip.open(snaps[j], 'r') as fsnap:
                                        next(fsnap)
                                        next(fsnap)
                                        for line in fsnap:
                                            datasnap=line.split()
                                            if int(datasnap[17])==13 and int(datasnap[18])==14:
                                                Model.append(i); Time.append(time*t_conv); ID0.append(int(datasnap[10])); ID1.append(int(datasnap[11])); P0.append(twopi*yearsc/float(datasnap[45])); P1.append(-100); B0.append(float(datasnap[47])); B1.append(-100); M0.append(float(datasnap[8])); M1.append(float(datasnap[9])); A.append(float(datasnap[12])); E.append(float(datasnap[13])); Type.append(1314); st.append(int(Status[i]))
                                            if int(datasnap[17])==14 and int(datasnap[18])==13:
                                                Model.append(i); Time.append(time*t_conv); ID0.append(int(datasnap[11])); ID1.append(int(datasnap[10])); P0.append(twopi*yearsc/float(datasnap[46])); P1.append(-100); B0.append(float(datasnap[48])); B1.append(-100); M0.append(float(datasnap[9])); M1.append(float(datasnap[8])); A.append(float(datasnap[12])); E.append(float(datasnap[13])); Type.append(1314); st.append(int(Status[i]))

        #print Model, ID0, ID1, P0, P1, B0, B1, M0, M1, A, E, Type
        print(i)

    np.savetxt('/projects/b1095/syr904/projects/SGRB/newruns/finaldata/DNS_NSBH_06Gyr_maingrid_31.dat', np.c_[Model, Time, ID0, ID1, P0, P1, B0, B1, M0, M1, A, E, Type, st], fmt='%d %f %d %d %f %f %e %e %f %f %f %f %d %d', header='1.Model, 2.Time(Myr), 3.ID0, 4.ID1, 5.P0(sec), 6.P1(sec), 7.B0(G), 8.B1(G), 9.M0(msun), 10.M1(msun), 11.a(AU), 12.ecc, 13.Type 14.Model Status', delimiter='', comments='#')



def find_DNS_NSBH_lastsnap(pathlist, start, end):   ##Unfinished
    sourcedir=np.genfromtxt(pathlist, dtype=str)
    Status=sourcedir[:,1]; sourcedir=sourcedir[:,0]
    pref='initial'
    Model=[]; Time=[]; ID0=[]; ID1=[]; P0=[]; P1=[]; B0=[]; B1=[]; M0=[]; M1=[]; A=[]; E=[]; Type=[]; st=[]; T_insp=[]

    ntot=0   ##Calculate how many dns and nsbh are at the last snapshot in total
    for i in range(start, end):
        filestr=sourcedir[i]+pref
        t_conv=dyn.conv('t', filestr+'.conv.sh')
        snaps=dyn.get_snapshots(filestr)

        t_last=dyn.get_time(snaps[-1])*t_conv
   
        with gzip.open(snaps[-1], 'r') as fsnap:
            next(fsnap)
            next(fsnap)
            for line in fsnap:
                datasnap=line.split()
                if int(datasnap[17])==13 and int(datasnap[18])==13:
                    Model.append(i); Time.append(t_last); ID0.append(int(datasnap[10])); ID1.append(int(datasnap[11])); P0.append(twopi*yearsc/float(datasnap[45])); P1.append(twopi*yearsc/float(datasnap[46])); B0.append(float(datasnap[47])); B1.append(float(datasnap[48])); M0.append(float(datasnap[8])); M1.append(float(datasnap[9])); A.append(float(datasnap[12])); E.append(float(datasnap[13])); Type.append(1313); st.append(int(Status[i]))
                    t_inspiral=gwcalc.t_inspiral_2(float(datasnap[12]), float(datasnap[13]), float(datasnap[8]), float(datasnap[9]), 0, 0, 0, 1100)/10**6 ##in Myr
                    T_insp.append(t_inspiral)
                    if Status[i]==1: ntot+=1

                if int(datasnap[17])==13 and int(datasnap[18])==14:
                    Model.append(i); Time.append(t_last); ID0.append(int(datasnap[10])); ID1.append(int(datasnap[11])); P0.append(twopi*yearsc/float(datasnap[45])); P1.append(-100); B0.append(float(datasnap[47])); B1.append(-100); M0.append(float(datasnap[8])); M1.append(float(datasnap[9])); A.append(float(datasnap[12])); E.append(float(datasnap[13])); Type.append(1314); st.append(int(Status[i]))
                    t_inspiral=gwcalc.t_inspiral_2(float(datasnap[12]), float(datasnap[13]), float(datasnap[8]), float(datasnap[9]), 0, 0, 0, 1100)/10**6 ##in Myr
                    T_insp.append(t_inspiral)
                    if Status[i]==1: ntot+=1

                if int(datasnap[17])==14 and int(datasnap[18])==13:
                    Model.append(i); Time.append(t_last); ID0.append(int(datasnap[11])); ID1.append(int(datasnap[10])); P0.append(twopi*yearsc/float(datasnap[46])); P1.append(-100); B0.append(float(datasnap[48])); B1.append(-100); M0.append(float(datasnap[9])); M1.append(float(datasnap[8])); A.append(float(datasnap[12])); E.append(float(datasnap[13])); Type.append(1314); st.append(int(Status[i]))
                    t_inspiral=gwcalc.t_inspiral_2(float(datasnap[12]), float(datasnap[13]), float(datasnap[8]), float(datasnap[9]), 0, 0, 0, 1100)/10**6 ##in Myr
                    T_insp.append(t_inspiral)
                    if Status[i]==1: ntot+=1


        print(i)

    np.savetxt('/projects/b1095/syr904/projects/SGRB/newruns/finaldata/DNS_NSBH_lastsnap_maingrid.dat', np.c_[Model, Time, ID0, ID1, P0, P1, B0, B1, M0, M1, A, E, Type, T_insp, st], fmt='%d %f %d %d %f %f %e %e %f %f %f %f %d %f %d', header='1.Model, 2.Time(Myr), 3.ID0, 4.ID1, 5.P0(sec), 6.P1(sec), 7.B0(G), 8.B1(G), 9.M0(msun), 10.M1(msun), 11.a(AU), 12.ecc, 13.Type 14.Tinsp(Myr) 15.Model Status', delimiter='', comments='#')




def find_DNS_NSBH_Unique(datfile):
    data=np.genfromtxt(datfile)
    model=data[:,0]; id0=data[:,2]; id1=data[:,3]; types=data[:,12]
    #print np.array(data[0,:])
    allmodel=list(Counter(model).keys())
    print(allmodel)

    #lines=[]; lines_needed=[]
    id0_unique=[]; id1_unique=[]; model_unique=[]; systype=[]
    ndns=0; nnsbh=0
    for i in range(len(allmodel)):
        modelno=int(allmodel[i])

        idstr_hold1=[str(0)]
        idstr_hold2=[str(0)]
        for j in range(len(model)):
            if int(model[j])==modelno:
                idstr=str(int(id0[j]))+str(int(id1[j]))
                check=1
                for k in range(len(idstr_hold1)):
                    if idstr==idstr_hold1[k] or idstr==idstr_hold2[k]:
                        check=0
    
                if check==1:
                    if int(types[j])==1313.:
                        model_unique.append(model[j]); id0_unique.append(id0[j]); id1_unique.append(id1[j]); systype.append(1313)
                        ndns+=1

                        idstr_hold1.append(idstr)
                        idstr_hold2.append(str(int(id1[j]))+str(int(id0[j])))

                    if int(types[j])==1314.:
                        model_unique.append(model[j]); id0_unique.append(id0[j]); id1_unique.append(id1[j]); systype.append(1314)
                        nnsbh+=1

                        idstr_hold1.append(idstr)
                        idstr_hold2.append(str(int(id1[j]))+str(int(id0[j])))

        print(i)
    print(ndns, nnsbh)
    print(id0_unique)

    funique=open('/projects/b1095/syr904/projects/SGRB/newruns/finaldata/DNS_NSBH_Unique_06Gyr_maingrid.dat', 'a+')
    funique.write('#1.Model, 2.Time, 3.ID0, 4.ID1, 5.P0, 6.P1, 7.B0, 8.B1, 9.M0, 10.M1, 11.a, 12.ecc, 13.Type 14.Status\n')
    for m in range(len(id0_unique)):
        with open(datfile, 'r') as fdat:
            next(fdat)
            for line in fdat:
                data=line.split()
                if int(data[0])==model_unique[m] and int(data[2])==id0_unique[m] and int(data[3])==id1_unique[m]:
                    theline=line

        funique.write(theline)

    funique.close()


    #np.savetxt('/projects/b1095/syr904/projects/SGRB/newruns/DNS_NSBH_Unique_nondissolved.dat', np.c_[model_unique, id0_unique, id1_unique, systype], fmt='%d %d %d %d', delimiter='', header='1.Model 2.ID0 3.ID1 4.Types', comments='#')


        #idkey=Counter(idmodel).keys()
        #print idkey 


        #for k in range(len(idkey)):
        #   for l in range(len(model)):
        #       if model[l]==modelno:
        #           ids=str(id0[l])+str(id1[l])
        #           if ids==idkey[k]:
        #               lines.append(l)
        #   lines_needed.append(lines[-1])

        #print modelno

    #print lines_needed

    #with open('/projects/b1011/syr904/projects/PULSAR/kickgrid_runs/paper/data/DNS_NSBH_Unique.dat', 'a+') as out_file:
    #   for m in range(len(lines_needed)):
    #       listm=' '.join(str(e) for e in data[lines_needed[m],:])
    #       out_file.write(listm+'\n')




def DNS_NSBH_mergertime(datfile, pathlist):
    data=np.genfromtxt(datfile)
    modelno=data[:,0]; t=data[:,1]; id0=data[:,2]; id1=data[:,3]; m0=data[:,8]; m1=data[:,9]; a=data[:,10]; ecc=data[:,11]; tp=data[:,12]
    sourcedir=np.genfromtxt(pathlist, dtype=str)
    nnsbhm=0; nnsnsm=0
    
    idstr_hold1=str(0)
    idstr_hold2=str(0)
    for j in range(len(modelno)):
        filestr=sourcedir[int(modelno[j])]+'/initial'
        idstr=str(int(id0[j]))+str(int(id1[j]))
        t_conv=dyn.conv('t', filestr+'.conv.sh')
        time=t[j]*t_conv
        t_inspiral=gwcalc.t_inspiral_2(a[j], ecc[j], m0[j], m1[j], 0, 0, 0, 1100)/10**6 #in Myr

        if t_inspiral+time<=12000. and (tp[j]==1314. or tp[j]==1413.) and idstr!=idstr_hold1 and idstr!=idstr_hold2: 
            print(time, t_inspiral)
            print(idstr)
            nnsbhm+=1
            idstr_hold1=idstr
            idstr_hold2=str(int(id1[j]))+str(int(id0[j]))
        if t_inspiral+time<=12000. and tp[j]==1313. and idstr!=idstr_hold1 and idstr!=idstr_hold2:
            print(time, t_inspiral)
            print(idstr)
            nnsnsm+=1
            idstr_hold1=idstr
            idstr_hold2=str(int(id1[j]))+str(int(id0[j]))

    print(nnsbhm, nnsnsm)




def BNS_NSBH_Appearperiod(datfile,pathlist):
    data=np.genfromtxt(datfile, usecols=(0, 1, 2, 3, 12))
    model=data[:,0]; time=data[:,1]; id0=data[:,2]; id1=data[:,3]; tp=data[:,4]
    allmodel=Counter(model).keys()
    sourcedir=np.genfromtxt(pathlist, dtype=str)

    ids=[]; periods=[]; types=[]
    for i in range(len(allmodel)):
        modelno=int(allmodel[i])
        t_conv=conv('t', sourcedir[modelno]+'/initial.conv.sh')

        idmodel=[]; pmodel=[]; tmodel=[]
        for j in range(len(model)):
            if model[j]==modelno:
                idmodel.append(str(id0[j])+str(id1[j]))
                pmodel.append(time[j])
                tmodel.append(tp[j])

        idkey=Counter(idmodel).keys()
        ids=ids+idkey

        for k in range(len(idkey)):
            p_temp=[]; t_temp=[]
            for l in range(len(idmodel)):
                if idmodel[l]==idkey[k]:
                    p_temp.append(pmodel[l])
                    t_temp.append(tmodel[l])

            #print p_temp

            if len(p_temp)==1: periods.append(30); types.append(t_temp[-1])
            if len(p_temp)>1: periods.append((p_temp[-1]-p_temp[0])*t_conv); types.append(t_temp[-1])


    return ids, periods, types


def find_psr_inmergers(pathlist):
    ##Read the path file
    sourcedir=np.genfromtxt(pathlist, dtype=str)
    status=sourcedir[:,1]; sourcedir=sourcedir[:,0]

    ##Read the merger files
    esc_dns=np.genfromtxt('/projects/b1095/syr904/projects/SGRB/newruns/finaldata/Esc_DNS_maingrid_v2.dat')
    esc_nsbh=np.genfromtxt('/projects/b1095/syr904/projects/SGRB/newruns/finaldata/Esc_NSBH_maingrid_v2.dat')

    gw_dns=np.genfromtxt('/projects/b1095/syr904/projects/SGRB/newruns/finaldata/GWcap_DNS_maingrid.dat')
    gw_nsbh=np.genfromtxt('/projects/b1095/syr904/projects/SGRB/newruns/finaldata/GWcap_NSBH_maingrid.dat')

    in_dns=np.genfromtxt('/projects/b1095/syr904/projects/SGRB/newruns/finaldata/Incluster_DNS_maingrid_v2.dat')
    in_nsbh=np.genfromtxt('/projects/b1095/syr904/projects/SGRB/newruns/finaldata/Incluster_NSBH_maingrid_v2.dat')


    esc_all=[esc_dns, esc_nsbh]
    in_all=[in_dns, in_nsbh]
    gw_all=[gw_dns, gw_nsbh]


    ##Make new file names
    filenames_esc=['Esc_DNS_maingrid_v3.dat', 'Esc_NSBH_maingrid_v3.dat']
    filenames_in=['Incluster_DNS_maingrid_v3.dat', 'Incluster_NSBH_maingrid_v3.dat']
    filenames_gw=['GWcap_DNS_maingrid_v3.dat', 'GWcap_NSBH_maingrid_v3.dat']


    for i in range(2):    ##DNS->NSBH

        models_esc=esc_all[i][:,0]
        id0_esc=esc_all[i][:,6]; id1_esc=esc_all[i][:,7]

        models_in=in_all[i][:,0] 
        idm_in=in_all[i][:,3]; id0_in=in_all[i][:,4]; id1_in=in_all[i][:,5]

        models_gw=gw_all[i][:,0]
        types=gw_all[i][:,3]; idm_gw=gw_all[i][:,4]
        id0_gw=gw_all[i][:,5]; id1_gw=gw_all[i][:,6]; id2_gw=gw_all[i][:,7]; id3_gw=gw_all[i][:,8]
        idgw=[id0_gw, id1_gw, id2_gw, id3_gw]
        k0_gw=gw_all[i][:,15]; k1_gw=gw_all[i][:,16]; k2_gw=gw_all[i][:,17]; k3_gw=gw_all[i][:,18]
        kgw=[k0_gw, k1_gw, k2_gw, k3_gw]

        for j in range(2):     ##ESC->In
            if j==0:
                models=models_esc; ids0=id0_esc; ids1=id1_esc
                filenames=filenames_esc
            if j==1:
                models=models_in; ids0=id0_in; ids1=id1_in
                filenames=filenames_in
            #if j==2:
            #   models=models_gw
            #   filenames=filenames_gw
            #   ids0=[]; ids1=[]
            #   for x in range(len(types)):
            #       kgw_temp=[k0_gw[x], k1_gw[x], k2_gw[x], k3_gw[x]]
            #       index_gw=[i for i, e in enumerate(kgw_temp) if e == 13]
            #       if i==0:
            #           ids0.append(idgw[int(index_gw[0])][x]); ids1.append(idgw[int(index_gw[1])][x])
            #       else:
            #           ids0.append(idgw[int(index_gw[0])][x]); ids1.append(-100)
            #       


            B0=-100; B1=-100; P0=-100; P1=-100; dmdt0=-100; dmdt1=-100

            with open('/projects/b1095/syr904/projects/SGRB/newruns/finaldata/'+filenames[i], 'a+') as f:
                print(filenames[i])
                #f.write('#1.Model 2.Time_esc(code) 3.Time_esc(Myr) 4.T_insp(Myr) 5.M0 6.M1 7.ID0 8.ID1 9.A 10.ECC 11.Tform(Myr) 12.Primordial? 13.Dynamically_Ejected? 14.B0 15.B1 16.P0 17.P1 18.dmdt0 19.dmdt1\n')
                for k in range(len(models)):
                    modelno=int(models[k])
                    filestr=sourcedir[modelno]+'initial'
                    if j>0:   ##for in-cluster files
                        datapsr=np.genfromtxt(filestr+'.morepulsars.dat')
                        id0_psr=datapsr[:,3]; id1_psr=datapsr[:,4]

                        for m in range(len(id0_psr)-1, -1, -1):
                            if id0_psr[m]==ids0[k] and id1_psr[m]==ids1[k]:
                                B0=datapsr[:,7][m]; B1=datapsr[:,8][m]
                                P0=datapsr[:,9][m]; P1=datapsr[:,10][m]
                                dmdt0=datapsr[:,17][m]; dmdt1=datapsr[:,18][m]
                                break

                            if id1_psr[m]==ids0[k] and id0_psr[m]==ids1[k]:
                                B1=datapsr[:,7][m]; B0=datapsr[:,8][m]
                                P1=datapsr[:,9][m]; P0=datapsr[:,10][m]
                                dmdt1=datapsr[:,17][m]; dmdt0=datapsr[:,18][m]
                                break
                                

                        appendstring=list(in_all[i][k,:])+[B0, B1, P0, P1, dmdt0, dmdt1]

                        newline=" ".join(str(x) for x in appendstring)+'\n'
                        f.write(newline)


                    else:   ##for escaped files
                        dataesc=np.genfromtxt(filestr+'.esc.dat')
                        id0_esc=dataesc[:,17]; id1_esc=dataesc[:,18]

                        for l in range(len(id0_esc)):
                            if id0_esc[l]==ids0[k] and id1_esc[l]==ids1[k]:
                                B0=dataesc[:,45][l]; B1=dataesc[:,46][l]
                                P0=twopi*yearsc/dataesc[:,43][l]; P1=twopi*yearsc/dataesc[:,44][l]
                                dmdt0=dataesc[:,39][l]; dmdt1=dataesc[:,40][l]
                                break
                            if id1_esc[l]==ids0[k] and id0_esc[l]==ids1[k]:
                                B1=dataesc[:,45][l]; B0=dataesc[:,46][l]
                                P1=twopi*yearsc/dataesc[:,43][l]; P0=twopi*yearsc/dataesc[:,44][l]
                                dmdt1=dataesc[:,39][l]; dmdt0=dataesc[:,40][l]
                                break

                        appendstring=list(esc_all[i][k,:])+[B0, B1, P0, P1, dmdt0, dmdt1]

                        newline=" ".join(str(x) for x in appendstring)+'\n'
                        f.write(newline)

                    print(k)

            print(j)
    

def selectmergers_lessthanHubbletime():
    ##Read the data files
    esc_dns=np.genfromtxt('/projects/b1095/syr904/projects/SGRB/newruns/Esc_DNS_nondissolved.dat')
    esc_nsbh=np.genfromtxt('/projects/b1095/syr904/projects/SGRB/newruns/Esc_NSBH_nondissolved.dat')
    esc_bbh=np.genfromtxt('/projects/b1095/syr904/projects/SGRB/newruns/Esc_BBH_nondissolved.dat')



    esc_all=[esc_dns, esc_nsbh, esc_bbh]

    filenames_esc=['Esc_DNS_nondissolved_hubble.dat', 'Esc_NSBH_nondissolved_hubble.dat','Esc_BBH_nondissolved_hubble.dat']


    for i in range(3):
        mergertime_esc=np.add(esc_all[i][:,2], esc_all[i][:,3])
        with open(savepath+'/newruns/'+filenames_esc[i], 'a+') as fesc:
            fesc.write('#1.Model 2.Time_esc(code) 3.Time_esc(Myr) 4.T_insp(Myr) 5.M0 6.M1 7.ID0 8.ID1 9.A 10.ECC 11.Tform(Myr) 12.Primordial? 13.Dynamically_Ejected? 14.B0 15.B1 16.P0 17.P1 18.dmdt0 19.dmdt1\n')
            for j in range(len(mergertime_esc)):
                if mergertime_esc[j]<=14000.:
                    newline_esc=" ".join(str(x) for x in esc_all[i][j,:])+'\n'
                    fesc.write(newline_esc)


def num_analysis(filepath):
    ##Read in the files
    #dns_esc=np.genfromtxt(filepath+'Esc_DNS_maingrid.dat')
    #dns_in=np.genfromtxt(filepath+'Incluster_DNS_maingrid.dat')
    #dns_cap=np.genfromtxt(filepath+'GWcap_DNS_maingrid.dat')

    #nsbh_esc=np.genfromtxt(filepath+'Esc_NSBH_maingrid.dat')
    #nsbh_in=np.genfromtxt(filepath+'Incluster_NSBH_maingrid.dat')
    #nsbh_cap=np.genfromtxt(filepath+'GWcap_NSBH_maingrid.dat')

    dns_esc=np.genfromtxt('/projects/b1095/syr904/projects/SGRB/newruns/finaldata/extreme_model/Esc_DNS_extreme.dat')
    dns_in=np.genfromtxt('/projects/b1095/syr904/projects/SGRB/newruns/finaldata/extreme_model/Incluster_DNS_extreme.dat')
    dns_cap=np.genfromtxt('/projects/b1095/syr904/projects/SGRB/newruns/finaldata/extreme_model/GWcap_DNS_extreme.dat')

    nsbh_esc=np.genfromtxt('/projects/b1095/syr904/projects/SGRB/newruns/finaldata/extreme_model/Esc_NSBH_extreme.dat')
    nsbh_in=np.genfromtxt('/projects/b1095/syr904/projects/SGRB/newruns/finaldata/extreme_model/Incluster_NSBH_extreme.dat')
    nsbh_cap=np.genfromtxt('/projects/b1095/syr904/projects/SGRB/newruns/finaldata/extreme_model/GWcap_NSBH_extreme.dat')
    
    bbh_esc=np.genfromtxt(filepath+'Esc_BBH_maingrid.dat')
    bbh_in=np.genfromtxt(filepath+'Incluster_BBH_maingrid.dat')
    bbh_cap=np.genfromtxt(filepath+'GWcap_BBH_maingrid.dat')    

    data_names=['dns_esc', 'dns_in', 'dns_cap', 'nsbh_esc', 'nsbh_in', 'nsbh_cap', 'bbh_esc', 'bbh_in', 'bbh_cap']
    data_objects=[dns_esc, dns_in, dns_cap, nsbh_esc, nsbh_in, nsbh_cap, bbh_esc, bbh_in, bbh_cap]


    ##Find the number of ejected systems that merge within a hubble time
    time_esc_dns=dns_esc[:,2]; time_inspiral_dns=dns_esc[:,3]
    time_esc_nsbh=nsbh_esc[:,2]; time_inspiral_nsbh=nsbh_esc[:,3]
    time_esc_bbh=bbh_esc[:,2]; time_inspiral_bbh=bbh_esc[:,3]

    modtype=[[], [], [], [], [], [], [], [], []]
    for x in range(len(data_objects)):
        modtype[x]=data_objects[x][:,-1]

    nesc_dns=0; nesc_nsbh=0; nesc_bbh=0
    pribinary=[[], [], [], [], [], [], [], [], []]; dynamic_eject=[[], [], [], [], [], [], [], [], []]
    for j in range(len(time_esc_dns)):
        if time_esc_dns[j]+time_inspiral_dns[j]<=14000.0 and modtype[0][j]==1: 
            nesc_dns+=1
            pribinary[0].append(data_objects[0][:,11][j])
            dynamic_eject[0].append(data_objects[0][:,12][j])
    for k in range(len(time_esc_nsbh)):
        if time_esc_nsbh[k]+time_inspiral_nsbh[k]<=14000.0 and modtype[3][k]==1: 
            nesc_nsbh+=1
            pribinary[3].append(data_objects[3][:,11][k])
            dynamic_eject[3].append(data_objects[3][:,12][k])
    for l in range(len(time_esc_bbh)):
        if time_esc_bbh[l]+time_inspiral_bbh[l]<=14000.0 and modtype[6][l]==1: 
            nesc_bbh+=1
            pribinary[6].append(data_objects[6][:,11][l])
            dynamic_eject[6].append(data_objects[6][:,12][l])

    #print pribinary[3], dynamic_eject[3]
    #print pribinary[3], dynamic[3]

    #print nesc_dns, nesc_nsbh


    ##Find the number of all systems
    nums=[0, 0, 0, 0, 0, 0, 0, 0, 0]
    for i in range(len(data_objects)):
        for y in range(len(modtype[i])):
            if modtype[i][y]==1: 
                nums[i]+=1

                if i==2 or i==5 or i==8:
                    pribinary[i].append(0)
                    dynamic_eject[i].append(1)
                if i==1 or i==4 or i==7:
                    pribinary[i].append(data_objects[i][:,11][y])
                    dynamic_eject[i].append(data_objects[i][:,12][y])


    ##Find the number of primordially- and dynamically-formed systems
    num_pri=[]; num_dynej=[]
    for i in range(len(data_objects)):
        #print len(pribinary[i]), len(dynamic[i])
        num_pri.append(np.sum(pribinary[i]))
        #num_dynej.append(np.sum(dynamic_eject[i]))


    #mod_status=[]
    #num_dissolve=0
    #fdiss=open(filepath+'Merger_in_DissolvedCluster.dat', 'a+')
    ###Find the number of mergers in the dissolved clusters
    #for x in range(len(data_objects)):
    #   mod_status.append(data_objects[x][:,-1])
    #   for y in range(len(mod_status[x])):
    #       if mod_status[x][y]!=1:
    #           #print str(data_objects[x][y,:])
    #           newline=" ".join(str(x) for x in list(data_objects[x][y,:]))+'\n'
    #           fdiss.write(newline)
    #           num_dissolve+=1

    #print num_dissolve


    ##If calculating escaped number of mergers within a hubble time
    nums[0]=nesc_dns; nums[3]=nesc_nsbh; nums[6]=nesc_bbh


    num_dns=np.array([nums[0], nums[1], nums[2]])
    num_nsbh=np.array([nums[3], nums[4], nums[5]])
    num_bbh=np.array([nums[6], nums[7], nums[8]])

    ntot_dns=float(sum(num_dns)); ntot_nsbh=float(sum(num_nsbh)); ntot_bbh=float(sum(num_bbh))
    frac_dns=num_dns/ntot_dns; frac_nsbh=num_nsbh/ntot_nsbh; frac_bbh=num_bbh/ntot_bbh

    print('nums:', num_dns, num_nsbh, num_bbh)
    print('fracs:', frac_dns, frac_nsbh, frac_bbh)
    print('ntots:', ntot_dns, ntot_nsbh, ntot_bbh)
    print('num_pri', num_pri)
    #print 'num_dyneject', num_dynej



def psr_analysis(filepath):
    ##Read in the files
    dns_esc=np.genfromtxt(filepath+'Esc_DNS_maingrid_v3.dat')
    dns_in=np.genfromtxt(filepath+'Incluster_DNS_maingrid_v3.dat')
    #dns_cap=np.genfromtxt(filepath+'GWcap_DNS_maingrid.dat')

    nsbh_esc=np.genfromtxt(filepath+'Esc_NSBH_maingrid_v3.dat')
    nsbh_in=np.genfromtxt(filepath+'Incluster_NSBH_maingrid_v3.dat')
    #nsbh_cap=np.genfromtxt(filepath+'GWcap_NSBH_maingrid.dat')

    data_objects=[dns_esc, nsbh_esc, dns_in, nsbh_in]

    for i in range(len(data_objects)):
        B0=data_objects[i][:,-6]; B1=data_objects[i][:,-5]
        P0=data_objects[i][:,-4]; P1=data_objects[i][:,-3]
        modtype=data_objects[i][:,14]
        if i<2:
            tinspiral=data_objects[i][:,3]
            tmerge=np.add(data_objects[i][:,2],data_objects[i][:,3])
        else:
            tinspiral=np.subtract(data_objects[i][:,2], data_objects[i][:,10])
            #print tinspiral
            tmerge=data_objects[i][:,2]


        #print B0

        nbmsp=0; nbpsr=0; nsmsp=0; nspsr=0; npsrmsp=0
        npsr0=[-100]*len(B0); npsr1=[-100]*len(B0)

        for j in range(len(B0)):
            if B0[j]!=0 and B0[j]!=-100 and tmerge[j]<=14000. and modtype[j]==1:
                #B0_new, P0_new=B0[j], P0[j]
                B0_new, P0_new=psrev.single_psr_evolv(B0[j], P0[j], tinspiral[j])
                psr0=unit_convert.psr_deathline(P0_new, B0_new) 
                if psr0=="yes" and P0_new<=0.03: npsr0[j]=1
                if psr0=="yes" and P0_new>0.03: npsr0[j]=2#; print tinspiral[j]

            if B1[j]!=0 and B1[j]!=-100 and tmerge[j]<=14000. and modtype[j]==1:
                #B1_new, P1_new=B1[j], P1[j]
                B1_new, P1_new=psrev.single_psr_evolv(B1[j], P1[j], tinspiral[j])
                psr1=unit_convert.psr_deathline(P1_new, B1_new) 
                if psr1=="yes" and P1[j]<=0.03: npsr1[j]=1
                if psr1=="yes" and P1[j]>0.03: npsr1[j]=2#; print tinspiral[j]



        for k in range(len(npsr0)):
            if npsr0[k]==1 and npsr1[k]==1: 
                nbmsp+=1#; print npsr0[k], npsr1[k]
            elif npsr0[k]==2 and npsr1[k]==2: 
                nbpsr+=1#; print npsr0[k], npsr1[k]
            elif npsr0[k]>=1 and npsr1[k]>=1: 
                npsrmsp+=1#; print npsr0[k], npsr1[k]
            elif npsr0[k]==1 or npsr1[k]==1: 
                nsmsp+=1#; print npsr0[k], npsr1[k]
                #B0_new, P0_new=psrev.single_psr_evolv(B0[k], P0[k], tinspiral[k])
                #B1_new, P1_new=psrev.single_psr_evolv(B1[k], P1[k], tinspiral[k])
                #print unit_convert.psr_deathline(P1[k], B1[k])
                #print unit_convert.psr_deathline(P1_new, B1_new)
                #print B0_new, P0_new, B0[k], P0[k]
                #print B1_new, P1_new, B1[k], P1[k]
            elif npsr0[k]==2 or npsr1[k]==2: 
                nspsr+=1#; print npsr0[k], npsr1[k]
            else:
                pass


        print('nbmsp, nbpsr, npsrmsp, nsmsp, nspsr', nbmsp, nbpsr, npsrmsp, nsmsp, nspsr)


##Find the number of systems at redshift smaller than the input value
def find_num_atredshift(filepath, redshift):
    ##Read in the files
    dns_esc=np.genfromtxt(filepath+'Esc_DNS_maingrid_v3.dat')
    dns_in=np.genfromtxt(filepath+'Incluster_DNS_maingrid_v3.dat')
    #dns_cap=np.genfromtxt(filepath+'GWcap_DNS_maingrid.dat')

    nsbh_esc=np.genfromtxt(filepath+'Esc_NSBH_maingrid_v3.dat')
    nsbh_in=np.genfromtxt(filepath+'Incluster_NSBH_maingrid_v3.dat')
    #nsbh_cap=np.genfromtxt(filepath+'GWcap_NSBH_maingrid.dat')

    data_objects=[dns_esc, nsbh_esc, dns_in, nsbh_in]

    for i in range(len(data_objects)):
        modtype=data_objects[i][:,14]
        priflag=data_objects[i][:,11]
        if i<2:
            tinspiral=data_objects[i][:,3]
            tmerge=np.add(data_objects[i][:,2],data_objects[i][:,3])
        else:
            tinspiral=np.subtract(data_objects[i][:,2], data_objects[i][:,10])
            tmerge=data_objects[i][:,2]
        
        npri=0; ntot=0
        for j in range(len(tmerge)):
            tcut=unit_convert.redshifttot(redshift)*1000.   ##in Myr
            if tcut<tmerge[j]<=14000. and modtype[j]==1:
                ntot+=1
                if priflag[j]==1: npri+=1

        print('ntot, npri:', ntot, npri)



def get_id_formationchannel(pathlist, sourcepath):
    sourcedir=np.genfromtxt(pathlist, dtype=str)
    filepath=sourcedir[:,0]; status=sourcedir[:,1]

    escdns=np.genfromtxt(sourcepath+'Esc_DNS_maingrid_v3.dat')
    indns=np.genfromtxt(sourcepath+'Incluster_DNS_maingrid_v3.dat')
    gwdns=np.genfromtxt(sourcepath+'GWcap_DNS_maingrid.dat')
    escnsbh=np.genfromtxt(sourcepath+'Esc_NSBH_maingrid_v3.dat')
    innsbh=np.genfromtxt(sourcepath+'Incluster_NSBH_maingrid_v3.dat')
    gwnsbh=np.genfromtxt(sourcepath+'GWcap_NSBH_maingrid.dat')

    escall=[escdns, escnsbh]
    inall=[indns, innsbh]
    gwall=[gwdns, gwnsbh]

    filename=['DNS', 'NSBH']

    #esc_fc0=[[],[]]; esc_fc1=[[],[]]
    #in_fc0=[[],[]]; in_fc1=[[],[]]
    #gw_fc0=[[],[]]; gw_fc1=[[],[]]

    for i in range(2):
        esc_id0=escall[i][:,6]; esc_id1=escall[i][:,7]
        esc_model=escall[i][:,0]; esc_mrgtime=np.add(escall[i][:,2],escall[i][:,3])

        in_id0=inall[i][:,4]; in_id1=inall[i][:,5]
        in_model=inall[i][:,0]; in_codetime=inall[i][:,1]; in_mrgtime=inall[i][:,2]  ##In Myr

        gw_model=gwall[i][:,0]; gw_mrgtime=gwall[i][:,2]; gw_codetime=gwall[i][:,1]
        gw_k0=gwall[i][:,15]; gw_k1=gwall[i][:,16]; gw_k2=gwall[i][:,17]; gw_k3=gwall[i][:,18]
        for x in range(len(gw_k0)):
            gwk=np.array([gw_k0[x], gw_k1[x], gw_k2[x], gw_k3[x]])
            index_gw=np.where(gwk==13)[0]; index_gw_bh=np.where(gwk==14)[0]
            if i == 0:
                gw_id0=gwall[i][:,5+index_gw[0]]; gw_id1=gwall[i][:,5+index_gw[1]]
            else:
                gw_id0=gwall[i][:,5+index_gw[0]]; gw_id1=gwall[i][:,5+index_gw_bh[0]]


        f=open(sourcepath+filename[i]+'_fc.dat', 'a+')
        f.write('#1.Model 2.Type 3.ID0 4.ID1 5.FC0 6.FC1 7.Merger_flag\n')
        ###Type: 'ej' means ejected binaries, 'in' means binaries found in the semerge file, and 'gw' means binaries found in the collision file.

        for j in range(len(esc_model)):
            modelno=int(esc_model[j])
            filestr=filepath[modelno]+'initial'
            with open(filestr+'.esc.dat', 'r') as fesc:
                next(fesc)
                for line in fesc:
                    dataesc=line.split()
                    if int(dataesc[14])==1 and int(dataesc[17])==int(esc_id0[j]) and int(dataesc[18])==int(esc_id1[j]):
                            esc_fc0=int(dataesc[47]); esc_fc1=int(dataesc[48])
                            if esc_mrgtime[j]<=14000.: 
                                merger_flag='yes'
                            else: 
                                merger_flag='no'
                            typeflag='ej'
            f.write('%d %s %d %d %d %d %s\n'%(esc_model[j], typeflag, esc_id0[j], esc_id1[j], esc_fc0, esc_fc1, merger_flag))
            print(j)

        for k in range(len(in_model)):
            in_fc0=-100; in_fc1=-100
            modelno=int(in_model[k])
            filestr=filepath[modelno]+'initial'
            snaps=dyn.get_snapshots(filestr)

            if in_mrgtime[k]<=14000.: 
                merger_flag='yes'
            else: 
                merger_flag='no'

            typeflag='in'

            for y in range(len(snaps)-1, -1, -1):
                t_snap=dyn.get_time(snaps[y])
                if t_snap<in_codetime[k]:
                    print(snaps[y])
                    with gzip.open(snaps[y], 'r') as fsnap:
                        next(fsnap); next(fsnap)
                        for line in fsnap:
                            datasnap=line.split()
                            if int(datasnap[7])==1 and int(datasnap[10])==int(in_id0[k]) and int(datasnap[11])==int(in_id1[k]):
                                in_fc0=int(datasnap[49]); in_fc1=int(datasnap[50])
                                break                  

                            if int(datasnap[7])==1 and int(datasnap[10])==int(in_id1[k]) and int(datasnap[11])==int(in_id0[k]):
                                in_fc0=int(datasnap[50]); in_fc1=int(datasnap[49])
                                break

                            if int(datasnap[7])==1 and (int(datasnap[10])==int(in_id0[k]) or int(datasnap[11])==int(in_id0[k])):
                                if int(datasnap[10])==int(in_id0[k]):
                                    in_fc0=int(datasnap[49])
                                if int(datasnap[11])==int(in_id0[k]):
                                    in_fc0=int(datasnap[50])

                            if int(datasnap[7])==1 and (int(datasnap[10])==int(in_id1[k]) or int(datasnap[11])==int(in_id1[k])):
                                if int(datasnap[10])==int(in_id1[k]):
                                    in_fc1=int(datasnap[49])
                                if int(datasnap[11])==int(in_id1[k]):
                                    in_fc1=int(datasnap[50])

                            if int(datasnap[7])!=1 and int(datasnap[0])==int(in_id0[k]):
                                in_fc0=int(datasnap[61])
                                
                            if int(datasnap[7])!=1 and int(datasnap[0])==int(in_id1[k]):
                                in_fc1=int(datasnap[61])

                    break
                print(y)

            f.write('%d %s %d %d %d %d %s\n'%(in_model[k], typeflag, in_id0[k], in_id1[k], in_fc0, in_fc1, merger_flag))
            print(k)



        for l in range(len(gw_model)):
            gw_fc0=-100; gw_fc1=-100
            modelno=int(gw_model[l])
            filestr=filepath[modelno]+'initial'
            snaps=dyn.get_snapshots(filestr)

            if gw_mrgtime[l]<=14000.: 
                merger_flag='yes'
            else: 
                merger_flag='no'
            
            typeflag='gw'

            for z in range(len(snaps)-1, -1, -1):
                t_snap=dyn.get_time(snaps[z])
                if t_snap<gw_codetime[l]:
                    print(snaps[z])   
                    with gzip.open(snaps[z], 'r') as fsnap:
                        next(fsnap); next(fsnap)
                        for line in fsnap:
                            datasnap=line.split()
                            if int(datasnap[7])==1 and int(datasnap[10])==int(gw_id0[l]) and int(datasnap[11])==int(gw_id1[l]):
                                gw_fc0=int(datasnap[49]); gw_fc1=int(datasnap[50])
                                break

                            if int(datasnap[7])==1 and int(datasnap[10])==int(gw_id1[l]) and int(datasnap[11])==int(gw_id0[l]):
                                gw_fc0=int(datasnap[50]); gw_fc1=int(datasnap[49])
                                break

                            if int(datasnap[7])==1 and (int(datasnap[10])==int(gw_id0[l]) or int(datasnap[11])==int(gw_id0[l])):
                                if int(datasnap[10])==int(gw_id0[l]):
                                    gw_fc0=int(datasnap[49])
                                if int(datasnap[11])==int(gw_id0[l]):
                                    gw_fc0=int(datasnap[50])

                            if int(datasnap[7])==1 and (int(datasnap[10])==int(gw_id1[l]) or int(datasnap[11])==int(gw_id1[l])):
                                if int(datasnap[10])==int(gw_id1[l]):
                                    gw_fc1=int(datasnap[49])
                                if int(datasnap[11])==int(gw_id1[l]):
                                    gw_fc1=int(datasnap[50])
                
                            if int(datasnap[7])!=1 and int(datasnap[0])==int(gw_id0[l]):
                                gw_fc0=int(datasnap[61])
                                
                            if int(datasnap[7])!=1 and int(datasnap[0])==int(gw_id1[l]):
                                gw_fc1=int(datasnap[61])

                    break
                print(z)

            f.write('%d %s %d %d %d %d %s\n'%(gw_model[l], typeflag, gw_id0[l], gw_id1[l], gw_fc0, gw_fc1, merger_flag))
            print(l)                 

        f.close()


##formation channel fraction analysis
def fc_analysis(sourcefile, pathlist):
    datafc=np.genfromtxt(sourcefile, dtype=str)
    models=datafc[:,0]; tp=datafc[:,1]; fc0=datafc[:,4]; fc1=datafc[:,5]; mrgflag=datafc[:,6]
    fmab0=datafc[:,7]; fmab1=datafc[:,8]

    sourcedir=np.genfromtxt(pathlist, dtype=str)
    paths=sourcedir[:,0]; status=sourcedir[:,1]
    
    ntot=0; ntot_one_ecsn=0;  ntot_two_ecsn=0; ntot_zero_ecsn=0 
    for i in range(len(models)):
        if int(status[int(models[i])])==1 and mrgflag[i]=='yes':
            ntot+=1; print(models[i])

            if int(fc0[i])>4 and int(fc1[i])>4: 
                ntot_two_ecsn+=1

            elif int(fc0[i])>4 or int(fc1[i])>4:
                if int(fc0[i])==0 and int(fmab0[i])==7002:
                    ntot_two_ecsn+=1
                elif int(fc1[i])==0 and int(fmab1[i])==7002:
                    ntot_two_ecsn+=1
                else:
                    ntot_one_ecsn+=1

            elif int(fc0[i])==4 and int(fc1[i])==4:
                ntot_zero_ecsn+=1

            else:
                if (int(fc0[i])==0 and int(fmab0[i])==7002) and (int(fc1[i])==0 and int(fmab1[i])==7002):
                    ntot_two_ecsn+=1
                elif (int(fc0[i])==0 and int(fmab0[i])==7002) or (int(fc1[i])==0 and int(fmab1[i])==7002):
                    ntot_one_ecsn+=1
                else:
                    ntot_zero_ecsn+=1

    print(ntot, ntot_two_ecsn, ntot_one_ecsn, ntot_zero_ecsn)


##Find the disrupted time of the disrupted clusters
def find_disrupted_time(pathlist):
    sourcedir=np.genfromtxt(pathlist, dtype=str)
    filepath=sourcedir[:,0]; status=sourcedir[:,1]

    for i in range(len(filepath)):
        if int(status[i])!=1:
            snaps=dyn.get_snapshots(filepath[i]+'initial')
            lastsnap=snaps[-1]
            time=dyn.get_time(lastsnap)
            t_conv=dyn.conv('t', filepath[i]+'initial.conv.sh')
            print(time*t_conv/1000.)


##Find the NS-MS/Giant binaries in the esc files
def find_nsbin_escfile(pathlist):
    sourcedir=np.genfromtxt(pathlist, dtype=str)
    filepath=sourcedir[:,0]; status=sourcedir[:,1]

    fesc_nsbin=open(savepath+'/newruns/finaldata/nsbin_esc_extra.dat', 'a+')
    fesc_nsbin.write('#1.Model 2.Status 3.Tesc(Myr) 4.Tphyf(Myr) 5.M0 6.M1 7.ID0 8.ID1 9.K0 10.K1 11.Period(days) 12.Ecc 13.Metallicity 14.Primordial? 15.Dyn_eject? 16.Dyn_any?\n')
    
    for i in range(len(filepath)):
        n, zm, r_g, r_v = unit_convert.find_init_conditions(filepath[i])
        n_nsbin=0

        filestr=filepath[i]+'initial'
        t_conv=dyn.conv('t', filestr+'.conv.sh')

        ####For NS binaries that may become NSNS/NSBH at late times####

        escfile=filestr+'.esc.dat'
        with open(escfile, 'r') as fesc:
            next(fesc)
            for line in fesc:
                dataesc=line.split()
                if int(dataesc[14])==1:
                    if (int(dataesc[22])==13 and 0<int(dataesc[23])<10 and float(dataesc[16])>=7.) or (int(dataesc[23])==13 and 0<int(dataesc[22])<10 and float(dataesc[15])>=7.):
                        n_nsbin+=1

                        tesc=float(dataesc[1])*t_conv  ##In Myr
                        model_esc=i
                        timeesc=float(dataesc[1]); m0=float(dataesc[15]); m1=float(dataesc[16]); id0=int(dataesc[17]); id1=int(dataesc[18]); a=float(dataesc[19]); ecc=float(dataesc[20]); k0=int(dataesc[22]); k1=int(dataesc[23])

                        period=unit_convert.au_to_period(a, m0, m1)
                        tphyf=14000.-tesc

                        tf_esc, sno_esc, pribin, dyn_eject, dyn_any=find_formationtime_randombin([id0, id1], filestr)

                        fesc_nsbin.write('%d %d %f %f %f %f %d %d %d %d %f %f %f %d %d %d\n'%(model_esc, int(status[i]), tesc, tphyf, m0, m1, id0, id1, k0, k1, period, ecc, zm, pribin, dyn_eject, dyn_any))


        escfile2=filestr+'2.esc.dat'
        if os.path.isfile(escfile2) and os.path.getsize(escfile2) > 0:
            with open(escfile2, 'r') as fesc:
                for line in fesc:
                    dataesc=line.split()
                    if int(dataesc[14])==1:
                        if (int(dataesc[22])==13 and 0<int(dataesc[23])<10 and float(dataesc[16])>=7.) or (int(dataesc[23])==13 and 0<int(dataesc[22])<10 and float(dataesc[15])>=7.):
                            n_nsbin+=1

                            tesc=float(dataesc[1])*t_conv  ##In Myr
                            model_esc=i
                            timeesc=float(dataesc[1]); m0=float(dataesc[15]); m1=float(dataesc[16]); id0=int(dataesc[17]); id1=int(dataesc[18]); a=float(dataesc[19]); ecc=float(dataesc[20]); k0=int(dataesc[22]); k1=int(dataesc[23])

                            period=unit_convert.au_to_period(a, m0, m1)
                            tphyf=14000.-tesc

                            tf_esc, sno_esc, pribin, dyn_eject, dyn_any=find_formationtime_randombin([id0, id1], filestr)

                            fesc_nsbin.write('%d %d %f %f %f %f %d %d %d %d %f %f %f %d %d %d\n'%(model_esc, int(status[i]), tesc, tphyf, m0, m1, id0, id1, k0, k1, period, ecc, zm, pribin, dyn_eject, dyn_any))

        print('escaped:', i, status[i], n_nsbin)

    fesc_nsbin.close()


##Find the number of DNSs and NS-BHs before 6 Gyr (z>1) in snapshots
def find_ndns_snapshots(pathlist):
    sourcedir=np.genfromtxt(pathlist, dtype=str)
    paths=sourcedir[:,0]; status=sourcedir[:,1]

    for i in range(len(paths)):
        if status[i]=='1':
            filestr=paths[i]+'initial'

            t_conv=dyn.conv('t', filestr+'.conv.sh')
            datans=np.genfromtxt(filestr+'.ns.dat')
            time=datans[:,0]; Ndns=datans[:,7]; Nnsbh=datans[:,8]

            ntot_dns=0; ntot_nsbh=0
            for j in range(len(time)):
                if time[j]*t_conv<6000.:
                    ntot_dns+=Ndns[j]; ntot_nsbh+=Nnsbh[j]

            print(i, ntot_dns, ntot_nsbh)


##Find second generation merger of a light BH (from DNS merger) with a massive BH (~15-25 Msun)
def find_2gen_merger(sourcepath):
    escdns=np.genfromtxt(sourcepath+'Esc_DNS_maingrid_v3.dat')
    indns=np.genfromtxt(sourcepath+'Incluster_DNS_maingrid_v3.dat')
    gwdns=np.genfromtxt(sourcepath+'GWcap_DNS_maingrid.dat')

    escnsbh=np.genfromtxt(sourcepath+'Esc_NSBH_maingrid_v3.dat')
    innsbh=np.genfromtxt(sourcepath+'Incluster_NSBH_maingrid_v3.dat')
    gwnsbh=np.genfromtxt(sourcepath+'GWcap_NSBH_maingrid.dat')

    escbbh=np.genfromtxt(sourcepath+'Esc_BBH_maingrid.dat')
    inbbh=np.genfromtxt(sourcepath+'Incluster_BBH_maingrid.dat')
    gwbbh=np.genfromtxt(sourcepath+'GWcap_BBH_maingrid.dat')


    escall=[escdns, escnsbh, escbbh]
    inall=[indns, innsbh, inbbh]
    gwall=[gwdns, gwnsbh, gwbbh]

    IDM=list(inall[0][:,3])+list(gwall[0][:,4])
    #print(IDM)

    ID_esc_nsbh=escall[1][:,6]+escall[1][:,7]
    ID_in_nsbh=inall[1][:,4]+inall[1][:,5]
    #print(ID_in_nsbh)
    ID_gw_nsbh=gwall[1][:,5]+gwall[1][:,6]+gwall[1][:,7]+gwall[1][:,8]

    ID_esc_bbh=escall[2][:,6]+escall[2][:,7]
    ID_in_bbh=inall[2][:,4]+inall[2][:,5]
    print(ID_in_bbh)
    ID_gw_bbh=gwall[2][:,5]+gwall[2][:,6]+gwall[2][:,7]+gwall[2][:,8]


    for i in range(len(IDM)):
        if IDM[i] in ID_in_bbh or IDM[i] in ID_gw_bbh or IDM[i] in ID_esc_bbh:
            print(IDM[i])

    for j in range(len(IDM)):
        if IDM[j] in ID_in_nsbh or IDM[j] in ID_gw_nsbh or IDM[j] in ID_esc_nsbh:
            print(IDM[j])










            




