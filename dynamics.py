import numpy as np
from glob import glob
import collections
from collections import Counter
import os,sys
import subprocess
import matplotlib
import matplotlib.lines as mlines
#matplotlib.use('PDF')
import matplotlib.pyplot as plt
#import seaborn as sns
import gzip
import math
import re
import history_cmc as hcmc
import conversions

savepath='/projects/b1095/syr904/projects/SGRB'


def conv_dict(): return {'l':15, 't':19, 'm':7}    #?


def conv(unit,filepath):   # Returns the unit conversion multiplier given a simulation's *.conv.sh file and a unit (either 'l' or 't' or 'm')
    dict = conv_dict()
    from re import findall
    with open(filepath,'r') as f:
        head = [next(f) for x in range(24)]
    return float(head[dict[unit]].strip().split('=')[-1])
    #return float(findall('\d+[\.]?\d*',head[dict[unit]])[0])


def get_time(filepath):      # Returns the cluster's age for a given snapshot
    import gzip
    from re import findall
    with gzip.open(filepath,'r') as f: contents = f.readline()
    if not findall(b'\d+[\.]?\d*',contents):        # Returns time = 0 for snapshot files without a time header
        print('snapshot empty'); return float(0)
    else: return float(findall(b'\d+[\.]?\d*',contents)[0])


###Find the snapshots while excluding the 2Dproj snapshots
def get_snapshots(filestring):
    snaps=np.sort(glob(filestring+'.snap*.dat.gz'))
    #print snaps
    snap2d=np.sort(glob(filestring+'.snap*.2Dproj.dat.gz'))

    index=[]
    for m in range(len(snaps)):
        for n in range(len(snap2d)):
            if snaps[m]==snap2d[n]: index.append(m)

    snaps = [x for y, x in enumerate(snaps) if y not in index]

    return snaps



##Find observed rc, rh in models at a random time
def find_obsrcrh(filestring, time):  ##if time==-100, find observed rc, rh in models in the last snapshot
    snapobs=np.sort(glob(filestring+'.snap*.cluster_params.dat'))

    if time>0:
        for i in range(len(snapobs)):
            try:
                dataobs=np.genfromtxt(snapobs[i])
                tobs=dataobs[0,0]
                if tobs>=time:
                    rc=dataobs[0, 9]; rhl=dataobs[0, 10];  t_Gyr=dataobs[0, 0]/1000.0; mass=dataobs[0, 4]
                    break
            except:
                print("Empty File", "snapno", i, "model=", filestring)
                continue
    else:
        try:
            snapobs_last=snapobs[-1]
            dataobs=np.genfromtxt(snapobs_last)
            rc=dataobs[0, 9]; rhl=dataobs[0, 10];  t_Gyr=dataobs[0, 0]/1000.0; mass=dataobs[0, 4]
        except:
            print("Empty File", "snapno=", snapobs_last, "model=", filestring)
            snapobs_last=snapobs[-2]
            dataobs=np.genfromtxt(snapobs_last)
            rc=dataobs[0, 9]; rhl=dataobs[0, 10];  t_Gyr=dataobs[0, 0]/1000.0; mass=dataobs[0, 4]

    return rc, rhl, t_Gyr, mass


##Find observed rc, rh in models in the last snapshot
def find_obsrcrh_last(filestring):
    snapobs=np.sort(glob(filestring+'.snap*.cluster_params.dat'))
    try:
        snapobs_last=snapobs[-1]
        dataobs=np.genfromtxt(snapobs_last)
        rc=dataobs[0, 9]; rhl=dataobs[0, 10];  t_Gyr=dataobs[0, 0]/1000.0; mass=dataobs[0, 4]
    except:
        print("Empty File", "snapno=", snapobs_last, "model=", filestring)
        snapobs_last=snapobs[-2]
        dataobs=np.genfromtxt(snapobs_last)
        rc=dataobs[0, 9]; rhl=dataobs[0, 10];  t_Gyr=dataobs[0, 0]/1000.0; mass=dataobs[0, 4]


    return rc, rhl, t_Gyr, mass


def find_obsrcrh_lastsnap_allmodels(pathlist, start, end):
    pref='initial'
    sourcedir=np.genfromtxt(pathlist, dtype='|S')
    RC=[]; RHL=[]; T=[]; NBH=[]; NTOT=[]; M=[]
    for i in range(start, end):
        filepath=sourcedir[i]
        pref='initial'
        filestr=filepath+'/'+pref
        snapobs=np.sort(glob(filestr+'.snap*.obs_params.dat'))
        lastsnapobs=snapobs[-1]
        Rc, Rhl, T_Gyr, m=find_obsrcrh(lastsnapobs)
        RC.append(Rc); RHL.append(Rhl); T.append(T_Gyr); M.append(m)
        Nbh, Ntot=find_NBH_NTOT_last(filestr)
        NBH.append(float(Nbh)); NTOT.append(float(Ntot))
        
        print(i)

    np.savetxt('/projects/b1011/syr904/projects/PULSAR/kickgrid_runs/paper/data/kickgrid_obsproperty_newmodel.dat', np.c_[T, NTOT, NBH, RC, RHL, M], fmt='%f %d %d %f %f %f', delimiter=' ', header='t_Gyr Ntot Nbh rc rhl M', comments='#')


##Find BH in the last timestep
def find_NBH_NTOT_last(filestring):
    filebh=filestring+'.bh.dat'
    filedyn=filestring+'.dyn.dat'
    with open(filebh, 'r') as fbh:
        for line in fbh:pass
        lastbh=line
    databh=lastbh.split()
    nbhlast=int(databh[2])

    with open(filedyn, 'r') as fdyn:
        for line in fdyn:pass
        lastdyn=line
    datadyn=lastdyn.split()
    ntotlast=int(datadyn[3])

    return nbhlast, ntotlast


##Find BH in random timestep
def find_NBH_NTOT(filestring, time, tconv):  ##For time==-100, find BH in the last timestep
    nbh=0; ntot=0; mass=0

    filebh=filestring+'.bh.dat'
    databh=np.genfromtxt(filebh)

    try:
        databh2 = np.genfromtxt(filestring+'2'+'.bh.dat')
        databh=np.concatenate((databh, databh2))
        #print databh
    except:
        databh=databh

    if time>0:
        for i in range(len(databh[:,1])):
            if databh[:,1][i]*tconv>=time:
                nbh=databh[:,2][i]
                break
    else:
        nbh=databh[:,2][-1]


    filedyn=filestring+'.dyn.dat'
    datadyn=np.genfromtxt(filedyn)

    try:
        datadyn2 = np.genfromtxt(filestring+'2'+'.dyn.dat')
        datadyn=np.concatenate((datadyn,datadyn2))
        #print datadyn
    except:
        datadyn=datadyn

    if time>0:
        for j in range(len(datadyn[:,0])):
            if datadyn[:,0][j]*tconv>=time:
                ntot=datadyn[:,3][j]; mass=datadyn[:,4][j] ##mass in code unit
                break
    else:
        ntot=datadyn[:,3][-1]; mass=datadyn[:,4][-1] ##mass in code unit

    return nbh, ntot, mass


##Find Mtot, rc, rh, rho0 in the last timestep in the dyn file
def find_rcrh_mtotrho0_last(filestring):
    filedyn=filestring+'.dyn.dat'
    with open(filedyn, 'r') as fdyn:
        for line in fdyn: pass
        lastdyn=line
    datadyn=lastdyn.split()
    mtot=float(datadyn[4]); rc=float(datadyn[7]); rh=float(datadyn[20]); rhoc=float(datadyn[21])

    return mtot, rc, rh, rhoc    ##in code unit


##Find Mtot, rc, rh, rho0 at a random timestep in the dyn file
def find_rcrh_mtotrho0(filestring, time, tconv):   ##if time==-100, find Mtot, rc, rh, rho0 in the last timestep in the dyn file
    filedyn=filestring+'.dyn.dat'
    filedyn2=filestring+'2.dyn.dat'

    check=0

    if time>0:
        if os.path.isfile(filedyn2) and os.path.getsize(filedyn2)>0:
            with open(filedyn2, 'r') as fdyn:
                next(fdyn)
                next(fdyn)
                for line in fdyn:
                    datadyn=line.split()
                    if float(datadyn[0])*tconv>=time: 
                        mtot=float(datadyn[4]); rc=float(datadyn[7]); rh=float(datadyn[20]); rhoc=float(datadyn[21])
                        check=1
                        break

        if check==0:
            with open(filedyn, 'r') as fdyn:
                next(fdyn)
                next(fdyn)
                for line in fdyn:
                    datadyn=line.split()
                    if float(datadyn[0])*tconv>=time: 
                        mtot=float(datadyn[4]); rc=float(datadyn[7]); rh=float(datadyn[20]); rhoc=float(datadyn[21])
                        break

    else:
        if os.path.isfile(filedyn2) and os.path.getsize(filedyn2)>0:
            with open(filedyn2, 'r') as fdyn:
                for line in fdyn: pass
                lastdyn=line
            datadyn=lastdyn.split()
            mtot=float(datadyn[4]); rc=float(datadyn[7]); rh=float(datadyn[20]); rhoc=float(datadyn[21])

        else:
            with open(filedyn, 'r') as fdyn:
                for line in fdyn: pass
                lastdyn=line
            datadyn=lastdyn.split()
            mtot=float(datadyn[4]); rc=float(datadyn[7]); rh=float(datadyn[20]); rhoc=float(datadyn[21])


    return mtot, rc, rh, rhoc    ##in code unit


##Find core mass in the last snapshot
def find_mc(filestring, r_c):   ##r_c in code units
    mc=0; mc_wobh=0
    snaps=np.sort(glob(filestring+'.snap*.dat.gz'))
    with gzip.open(snaps[-1], 'r') as flastsnap:
        for _ in xrange(2): next(flastsnap)
        for line in flastsnap:
            datalast=line.split()
            r=float(datalast[2]); m=float(datalast[1])
            binflag=int(datalast[7]); k=int(datalast[14]); k0=int(datalast[17]); k1=int(datalast[18]); m0=float(datalast[8]); m1=float(datalast[9])

            if r<=r_c: 
                mc+=m
                if binflag!=1 and k!=14:
                    mc_wobh+=m
                if binflag==1 and (k0!=14 and k1!=14):
                    mc_wobh+=m
                if binflag==1 and (k0!=14 and k1==14):
                    mc_wobh+=m0
                if binflag==1 and (k0==14 and k1!=14):
                    mc_wobh+=m1
            if r>r_c: break

    return mc, mc_wobh


##Find the number of NSs and NS binaries in the initial.ns.dat file at a random timestep
def find_Nns(filestring, time, tconv): ##if time==-100, find the number of NSs and NS binaries in the initial.ns.dat file at the last snapshot
    datans=np.genfromtxt(filestring+'.ns.dat')
    if time>0:
        for i in range(len(datans[:,0])):
            if datans[:,0][i]*tconv>=time:
                nns=datans[:,1][i]; ndns=datans[:,7][i]; nnsbh=datans[:,8][i]
                npsr=datans[:,5][i]; nmsp=datans[:,6][i]
                break
    else:
        nns=datans[:,1][-1]; ndns=datans[:,7][-1]; nnsbh=datans[:,8][-1]
        npsr=datans[:,5][-1]; nmsp=datans[:,6][-1]

    return nns, ndns, nnsbh, npsr, nmsp


##Find the number of NSs and NS binaries in the initial.ns.dat file at the last snapshot
def find_Nns_last(filestring):
    datans=np.genfromtxt(filestring+'.ns.dat')
    nns=datans[-1,1]; ndns=datans[-1,7]; nnsbh=datans[-1,8]

    return nns, ndns, nnsbh


##Find Nbh, Mtot, Nns, Ndns for a given time, thetime in Myr
def find_clusterparameter_allmodel(pathlist, start, end, thetime):
    model=[]; NBH=[]; MTOT=[]; RC=[]; RH=[]; NNS=[]; NDNS=[]; NNSBH=[]; st=[]
    RC_obs=[]; RHL_obs=[]; NPSR=[]; NMSP=[]
    sourcedir=np.genfromtxt(pathlist, dtype=str)
    filepath=sourcedir[:,0]#; status=sourcedir[:,1]
    status=[0,0,0,0,0,0,0,0,0,0]
    #print filepath

    for i in range(start, end):
        #filestr=sourcedir[i]+'initial'
        filestr=filepath[i]+'initial'
        t_conv=conv('t', filestr+'.conv.sh')
        l_conv=conv('l', filestr+'.conv.sh')
        m_conv=conv('m', filestr+'.conv.sh')

        #filedyn=filestr+'.dyn.dat'
        #with open(filedyn, 'r') as fdyn:
        #    for line in fdyn: pass
        #    lastdyn=line
        #datadyn=lastdyn.split()
        #t_last=float(datadyn[0])*t_conv


        ##If the cluster dissolved, get the numbers for the last snapshot; if not, get the numbers within the last 3Gyr
        if int(status[i])!=1:
            ##Numbers at the last snapshot
            Nbh, Ntot, Mtot=find_NBH_NTOT(filestr, -100, t_conv)
            Mtot, Rc, Rh, Rhoc=find_rcrh_mtotrho0(filestr, -100, t_conv)
            Rc_obs, Rhl_obs, t_Gyr_obs, Mass_obs=find_obsrcrh(filestr, -100)
            Nns, Ndns, Nnsbh, Npsr, Nmsp=find_Nns(filestr, -100, t_conv)
            st.append(int(status[i]))

        else:
            ##Numbers at a certain time
            Nbh, Ntot, Mtot=find_NBH_NTOT(filestr, thetime, t_conv)
            Mtot, Rc, Rh, Rhoc=find_rcrh_mtotrho0(filestr, thetime, t_conv)
            Nns, Ndns, Nnsbh, Npsr, Nmsp=find_Nns(filestr, thetime, t_conv)
            Rc_obs, Rhl_obs, t_Gyr_obs, Mass_obs=find_obsrcrh(filestr, thetime)
            st.append(int(status[i]))


        model.append(i); NBH.append(Nbh); MTOT.append(Mtot*m_conv); RC.append(Rc*l_conv); RH.append(Rh*l_conv); RC_obs.append(Rc_obs); RHL_obs.append(Rhl_obs)
        NNS.append(Nns); NDNS.append(Ndns); NNSBH.append(Nnsbh); NPSR.append(Npsr); NMSP.append(Nmsp)

        print(i)


    np.savetxt('/projects/b1095/syr904/cmc/47Tuc/rundir/clusterproperty_lastsnap.dat', np.c_[model, NBH, MTOT, RC, RH, RC_obs, RHL_obs, NNS, NDNS, NNSBH, NPSR, NMSP, st], fmt='%d %d %f %f %f %f %f %d %d %d %d %d %d', header='1.Model 2.Nbh 3.Mtot(Msun) 4.rc(pc) 5.rh(pc) 6.rc_obs(pc) 7.rhl_obs(pc) 8.Nns 9.Ndns 10.Nnsbh 11.Npsr 12.Nmsp 13.Dissolved?', delimiter='', comments='#')



##Find Nbh, Mtot, Nns, Ndns at the last snapshot
def find_clusterparameter_allmodel_last(pathlist, start, end):
    model=[]; NBH=[]; MTOT=[]; RC=[]; RH=[]; NNS=[]; NDNS=[]; NNSBH=[]; st=[]
    RC_obs=[]; RHL_obs=[]
    sourcedir=np.genfromtxt(pathlist, dtype=str)
    filepath=sourcedir[:,0]#; status=sourcedir[:,1]
    #print filepath

    for i in range(start, end):
        filestr=sourcedir[i]+'initial'
        #filestr=filepath[i]+'initial'
        t_conv=conv('t', filestr+'.conv.sh')
        l_conv=conv('l', filestr+'.conv.sh')
        m_conv=conv('m', filestr+'.conv.sh')

        filedyn=filestr+'.dyn.dat'
        with open(filedyn, 'r') as fdyn:
            for line in fdyn: pass
            lastdyn=line
        datadyn=lastdyn.split()
        t_last=float(datadyn[0])*t_conv


        ##Numbers at the last snapshot
        Nbh, Ntot=find_NBH_NTOT_last(filestr)
        Mtot, Rc, Rh, Rhoc=find_rcrh_mtotrho0_last(filestr)
        Rc_obs, Rhl_obs, t_Gyr_obs, Mass_obs=find_obsrcrh_last(filestr)
        Nns, Ndns, Nnsbh=find_Nns_last(filestr)


        model.append(i); NBH.append(Nbh); MTOT.append(Mtot*m_conv); RC.append(Rc*l_conv); RH.append(Rh*l_conv); RC_obs.append(Rc_obs); RHL_obs.append(Rhl_obs)
        NNS.append(Nns); NDNS.append(Ndns); NNSBH.append(Nnsbh)

        print(i)


    np.savetxt(savepath+'/newruns/finaldata/clusterproperty_nondissolved_last.dat', np.c_[model, NBH, MTOT, RC, RH, RC_obs, RHL_obs, NNS, NDNS, NNSBH], fmt='%d %d %f %f %f %f %f %d %d %d', header='1.Model 2.Nbh 3.Mtot(Msun) 4.rc(pc) 5.rh(pc) 6.rc_obs(pc) 7.rhl_obs(pc) 8.Nns 9.Ndns 10.Nnsbh', delimiter='', comments='#')



##Print out histories of interactions for interesting systems
def printout_interact_history(pathlist, sourcefile, filename):
    paths=np.genfromtxt(pathlist, dtype=str)
    status=paths[:,1]; paths=paths[:,0]
    data=np.genfromtxt(sourcefile, usecols=[0, 6, 7])
    models=data[:,0]; id0=data[:,1]; id1=data[:,2]

    for k in range(len(id0)):
        filepath=paths[int(models[k])]
        filestr=filepath+'initial'
        t_conv=conv('t', filestr+'.conv.sh')
        hdict=hcmc.history_maker([id0[k]], [1], 'initial', filepath, 1.)
        #print id0[k]
        find_binint(hdict, id0[k], id1[k], models[k], t_conv, filename)
        #print a_in, a_out

        print(k)

    

##Extract parameters from history dictionary
def find_binint(hdict, theid, comid, modelno, tconv, fname):
    ain=0; ein=0; aout=0; eout=0
    timeint=0; typeint=0; minp=0; mout=0
    for j in hdict[theid]['binint']['binint'].keys():
        for i in hdict[theid]['binint']['binint'][j]['interaction']['input'].keys():
            idlenin=len(hdict[theid]['binint']['binint'][j]['interaction']['input'][i]['ids'])
            if idlenin==2:
                #print str(hdict[theid]['binint']['binint'][j]['interaction']['input'][i]['ids'][0]), str(hdict[theid]['binint']['binint'][j]['interaction']['input'][i]['ids'][1])
                if str(hdict[theid]['binint']['binint'][j]['interaction']['input'][i]['ids'][0]).rfind(':')<=-1 and str(hdict[theid]['binint']['binint'][j]['interaction']['input'][i]['ids'][1]).rfind(':')<=-1:
                    if long(hdict[theid]['binint']['binint'][j]['interaction']['input'][i]['ids'][0])==theid and long(hdict[theid]['binint']['binint'][j]['interaction']['input'][i]['ids'][1])==comid:
                        ain=float(hdict[theid]['binint']['binint'][j]['interaction']['input'][i]['a'])
                        ein=float(hdict[theid]['binint']['binint'][j]['interaction']['input'][i]['e'])
                        minp=float(hdict[theid]['binint']['binint'][j]['interaction']['input'][i]['m'][0])
                        #print ain, ein, minp
                
                    if long(hdict[theid]['binint']['binint'][j]['interaction']['input'][i]['ids'][1])==theid and long(hdict[theid]['binint']['binint'][j]['interaction']['input'][i]['ids'][0])==comid:
                        ain=float(hdict[theid]['binint']['binint'][j]['interaction']['input'][i]['a'])
                        ein=float(hdict[theid]['binint']['binint'][j]['interaction']['input'][i]['e'])
                        minp=float(hdict[theid]['binint']['binint'][j]['interaction']['input'][i]['m'][1])
                        #print ain, ein, minp


        if ain!=0:
            for o in hdict[theid]['binint']['binint'][j]['interaction']['output'].keys():
                idlenout=len(hdict[theid]['binint']['binint'][j]['interaction']['output'][o]['ids'])
                if idlenout==2:
                    if str(hdict[theid]['binint']['binint'][j]['interaction']['output'][o]['ids'][0]).rfind(':')<=-1 and str(hdict[theid]['binint']['binint'][j]['interaction']['output'][o]['ids'][1]).rfind(':')<=-1:
                        if long(hdict[theid]['binint']['binint'][j]['interaction']['output'][o]['ids'][0])==theid and long(hdict[theid]['binint']['binint'][j]['interaction']['output'][o]['ids'][1])==comid:
                            aout=float(hdict[theid]['binint']['binint'][j]['interaction']['output'][o]['a'][0])
                            eout=float(hdict[theid]['binint']['binint'][j]['interaction']['output'][o]['e'][0])
                            mout=float(hdict[theid]['binint']['binint'][j]['interaction']['output'][o]['m'][0])
                            #print aout, eout, mout

                        if long(hdict[theid]['binint']['binint'][j]['interaction']['output'][o]['ids'][1])==theid and long(hdict[theid]['binint']['binint'][j]['interaction']['output'][o]['ids'][0])==comid:
                            aout=float(hdict[theid]['binint']['binint'][j]['interaction']['output'][o]['a'][0])
                            eout=float(hdict[theid]['binint']['binint'][j]['interaction']['output'][o]['e'][0])
                            mout=float(hdict[theid]['binint']['binint'][j]['interaction']['output'][o]['m'][1])
                            #print aout, eout, mout

                if idlenout==3:
                    if str(hdict[theid]['binint']['binint'][j]['interaction']['output'][o]['ids'][0]).rfind(':')<=-1 and str(hdict[theid]['binint']['binint'][j]['interaction']['output'][o]['ids'][1]).rfind(':')<=-1:
                        if long(hdict[theid]['binint']['binint'][j]['interaction']['output'][o]['ids'][0])==theid and long(hdict[theid]['binint']['binint'][j]['interaction']['output'][o]['ids'][1])==comid:
                            aout=float(hdict[theid]['binint']['binint'][j]['interaction']['output'][o]['a'][0])
                            eout=float(hdict[theid]['binint']['binint'][j]['interaction']['output'][o]['e'][0])
                            mout=float(hdict[theid]['binint']['binint'][j]['interaction']['output'][o]['m'][0])
                            print(aout, eout, mout)

                        if long(hdict[theid]['binint']['binint'][j]['interaction']['output'][o]['ids'][1])==theid and long(hdict[theid]['binint']['binint'][j]['interaction']['output'][o]['ids'][0])==comid:
                            aout=float(hdict[theid]['binint']['binint'][j]['interaction']['output'][o]['a'][0])
                            eout=float(hdict[theid]['binint']['binint'][j]['interaction']['output'][o]['e'][0])
                            mout=float(hdict[theid]['binint']['binint'][j]['interaction']['output'][o]['m'][1])
                            print(aout, eout, mout)

                if aout!=0:
                        typeint=hdict[theid]['binint']['binint'][j]['interaction']['type']['type']
                        timeint=float(hdict[theid]['binint']['binint'][j]['interaction']['type']['time'])
                        timeint=timeint*tconv


        if ain!=0 and aout!=0:
            fbinint=open('/projects/b1095/syr904/projects/SGRB/newruns/finaldata/interactions/'+fname, 'a+', 0)
            fbinint.write('%d %d %d %f %s %f %f %f %f %f %f\n'%(modelno, theid, comid, timeint, typeint, ain, ein, minp, aout, eout, mout))

    #print ain, aout
    #return timeint, typeint, ain, ein, minp, aout, eout, mout


##Find the mass within a radius of the cluster. 
##Input radius is in the unit of arcsec. Input R_sun is in the radius of kpc
def find_mass_inradius(modelpath, snap2d_no, radius_arcsec, R_sun_kpc):
    #data2d=np.genfromtxt(modelpath+'initial.snap'+snap2d_no+'.2Dproj.dat.gz')
    #r2d=data2d[:,0]; Mr=data2d[:,9]   ##r2d in pc and Mr in Msun

    rcut=conversions.arcsec_to_pc(radius_arcsec, R_sun_kpc)

    Mtot_cut=0; Mtot=0
    with gzip.open(modelpath+'initial.snap'+snap2d_no+'.2Dproj.dat.gz', 'r') as f2d:
        next(f2d)
        next(f2d)
        for line in f2d:
            data2d=line.split()
            Mtot+=float(data2d[9])

            if float(data2d[0])<=rcut:
                Mtot_cut+=float(data2d[9])

    print(Mtot_cut, Mtot)



