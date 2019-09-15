import numpy as np 
import os,sys
import scripts
import scripts1
import scripts2
import scripts3
import dynamics as dyn
import unit_convert as uc
import ecc_calc as gwcalc
import LISA_calculations as lisa

##Find the number of NS-MS collisions in the models
def get_NS_collision(pathlist, start, end):
    sourcedir=np.genfromtxt(pathlist, dtype=str)
    filepaths=sourcedir[:,0]; status=sourcedir[:,1]


    model=[]; model_status=[]; mm=[]; m0=[]; m1=[]; ktypem=[]; ktype0=[]; ktype1=[]; timem=[]; idm=[]; rm=[]
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
            if int(line[2])==2:  ##Single-single star collision
                if int(line[11])==13 or int(line[12])==13:
                    model.append(i); model_status.append(int(status[i])); timem.append(t_conv*float(line[0]))
                    mm.append(float(line[4])); m0.append(float(line[6])); m1.append(float(line[8]))
                    ktypem.append(int(line[10])); ktype0.append(int(line[11])); ktype1.append(int(line[12]))
                    idm.append(int(line[3])); rm.append(float(line[9])*l_conv)

        print i


    np.savetxt('/projects/b1095/syr904/projects/PULSAR2/newruns/tidal_capture/NSMS.dat', np.c_[model, timem, idm, rm, mm, m0, m1, ktypem, ktype0, ktype1, model_status], fmt='%d %f %d %f %f %f %f %d %d %d %d', header='1.Model 2.Time(Myr) 3.IDcoll 4.Radius(pc) 5.Mcoll 6.M0 7.M1 8.kcoll 9.k0 10.k1 11.model_status', comments='#', delimiter='')


##Total numbers of different collision systems
def get_num_NSMS(nscollfile):
    data=np.genfromtxt(nscollfile)
    status=data[:,10]; k0=data[:,8]; k1=data[:,9]

    count_ms=0; count_giant=0; count_wd=0; count_ns=0; count_bh=0
    for i in range(len(status)):
        if status[i]==1:
            if k0[i]<=1 or k1[i]<=1: count_ms+=1
            elif k0[i]<10 or k1[i]<10: count_giant+=1
            elif k0[i]<13 or k1[i]<13: count_wd+=1
            elif k0[i]==13 and k1[i]==13: count_ns+=1
            else: 
                count_bh+=1
                print i


    print count_ms, count_giant, count_wd, count_ns, count_bh



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


        print i

    print np.sum(count_allbinint[0]), np.sum(count_alldns[0]), np.sum(count_allmerge[0])
    print np.sum(count_allbinint[1]), np.sum(count_alldns[1]), np.sum(count_allmerge[1])
    print np.sum(count_allbinint[2]), np.sum(count_alldns[2]), np.sum(count_allmerge[2])

    return sma, m0_sma, m1_sma



            


