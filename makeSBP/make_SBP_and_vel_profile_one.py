import numpy as np
import pandas as pd
import re
#import finder
#import finder2
#import finder3
import os,sys
import subprocess
import gzip
from glob import glob
#import scripts
import matplotlib.pyplot as plt
#import history_cmc
#import blackhole
#import extract_observed_prop as OBS
import random
import scipy.optimize as opt
import unit_convert as uc

sys.path.insert(1, '/projects/b1095/syr904/MyCodes/cmctoolkit')
import cmctoolkit as cmct


def conv_dict(): return {'l':15, 't':19, 'm':7, 'mstar':11}    #?


def conv(unit,filepath):   # Returns the unit conversion multiplier given a simulation's *.conv.sh file and a unit (either 'l' or 't' or 'm' or 'mstar')
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


def read_units(filestring):
    """reads from the supplied conv file and stores the physical units"""
    s = filestring+'.conv.sh'
    f=open(s,'r')
    count=0
    while count<10:
        a=f.readline()
        b=a.split('=')
        if b[0]=='massunitcgs':
            massunitcgs = float(b[1])
            count+=1
        elif b[0]=='massunitmsun':
            massunitmsun = float(b[1])
            count+=1
        elif b[0]=='mstarunitcgs':
            mstarunitcgs = float(b[1])
            count+=1
        elif b[0]=='mstarunitmsun':
            mstarunitmsun = float(b[1])
            count+=1
        elif b[0]=='lengthunitcgs':
            lengthunitcgs = float(b[1])
            count+=1
        elif b[0]=='lengthunitparsec':
            lengthunitparsec = float(b[1])
            count+=1
        elif b[0]=='timeunitcgs':
            timeunitcgs = float(b[1])
            count+=1
        elif b[0]=='timeunitsmyr':
            timeunitsmyr = float(b[1])
            count+=1
        elif b[0]=='nbtimeunitcgs':
            nbtimeunitcgs = float(b[1])
            count+=1
        elif b[0]=='nbtimeunitsmyr':
            nbtimeunitsmyr = float(b[1])
            count+=1
    f.close()
    units = []
    unittype = [('m_cgs', float), ('m_msun', float), ('mstar_cgs', float), ('mstar_msun', float), ('l_cgs', float), ('l_pc', float), ('t_cgs', float),('t_myr', float), ('nbt_cgs', float), ('nbt_myr', float)]
    units.append((massunitcgs, massunitmsun, mstarunitcgs, mstarunitmsun, lengthunitcgs, lengthunitparsec, timeunitcgs, timeunitsmyr, nbtimeunitcgs, nbtimeunitsmyr))
    units = np.array(units, dtype=unittype)
    #print(units)
    #return (massunitcgs, massunitmsun, mstarunitcgs, mstarunitmsun, lengthunitcgs, lengthunitparsec, timeunitcgs, timeunitsmyr, nbtimeunitcgs, nbtimeunitsmyr)
    return units


def find_max_snap(path,string):
    snapstring = path+string
    snapfiles = get_snapshots(snapstring)

    snapno_max = len(snapfiles)-1
    print(snapno_max)
    if snapno_max < 10:
        snapno_max_str = '000'+str(snapno_max)
    if snapno_max < 100 and snapno_max >= 10:
        snapno_max_str = '00'+str(snapno_max)
    if snapno_max < 1000 and snapno_max >= 100:
        snapno_max_str = '0'+str(snapno_max)
    if snapno_max < 10000 and snapno_max >= 1000:
        snapno_max_str = str(snapno_max)
    if snapno_max >= 10000:
        snapno_max_str = str(snapno_max)
    return snapno_max_str



def find_snap_time_array(path,string):
    snapno_max = int(find_max_snap(path,string))
    #print(snapno_max)
    units=read_units(path+string)
    km = units[0]['l_cgs']*1.0e-5
    time_units = units[0]['t_myr']

    t_array = []
    s_array = []
    for j in range(0,snapno_max+1):
        if j < 10:
            snapno = '000'+str(j)
        if j < 100 and j >= 10:
            snapno = '00'+str(j)
        if j < 1000 and j >= 100:
            snapno = '0'+str(j)
        if j < 10000 and j >= 1000:
            snapno = str(j)
        if j >= 10000:
            snapno = str(j)

        time = get_time(path+string+'.snap'+snapno+'.dat.gz')
        time = time*time_units      # DEFINE THE TIME OF THE SNAPFILE
        t_array.append(time)
        s_array.append(snapno)

    #print(t_array, s_array)
    return t_array, s_array



def convert_to_3d(r, vr, vt, SEEDY=100):
    random.seed(SEEDY)
    #costheta = np.random.uniform(-1, 1)
    #sintheta = (1-costheta**2.)**0.5

    if np.shape(r)==():
        #print 'came here'
        sintheta = np.random.uniform(low=-1., high=1.)
        phi = np.random.uniform(low=0., high=2.*np.pi)
        anglev = np.random.uniform(low=0., high=2.*np.pi)
    else:
        #print 'came here too'
        r = np.array(r)
        vr = np.array(vr)
        vt = np.array(vt)
    
        sintheta = np.random.uniform(low=-1., high=1., size=len(r))
        phi = np.random.uniform(low=0., high=2.*np.pi, size=len(r))
        anglev = np.random.uniform(low=0., high=2.*np.pi, size=len(r))

    print(sintheta, phi, anglev)
        
    costheta = (1-sintheta**2.)**0.5
    #costheta = (1-sintheta**2.)**0.5
    #phi = np.random.uniform(low=0., high=2.*np.pi, size=len(r))
    
    rz = r*sintheta
    rx = r*costheta*np.cos(phi)
    ry = r*costheta*np.sin(phi)
    
    #anglev = np.random.uniform(low=0., high=2.*np.pi, size=len(r))
    magv = (vr*vr + vt*vt)**0.5 
    thetadot = np.cos(anglev) * vt/r
    phidot = np.sin(anglev)*vt/(r*costheta)
    
    #vx = vr * np.sin(np.arccos(costheta)) * np.cos(phi) + r * thetadot * costheta * np.cos(phi) - r * phidot * np.sin(np.arccos(costheta)) * np.sin(phi)
    vx = vr * costheta * np.cos(phi) - r * phidot * costheta * np.sin(phi) - r * thetadot * sintheta * np.cos(phi)
    #vy = vr * np.sin(np.arccos(costheta)) * np.sin(phi) + r * thetadot * costheta * np.sin(phi) + r * phidot * np.sin(np.arccos(costheta)) * np.cos(phi)
    vy = vr * costheta * np.sin(phi) + r * phidot * costheta * np.cos(phi) - r * thetadot * sintheta * np.sin(phi)
    #vz = vr * costheta - r * thetadot * np.sin(np.arccos(costheta))
    vz = vr * sintheta + r * thetadot * costheta


    r3d = np.array([rx, ry, rz])
    v3d = np.array([vx, vy, vz])

    return r3d, v3d


def project_and_radially_sort(r3d, PROJ=(0,1)):
    r2d = (r3d[PROJ[0],:]**2. + r3d[PROJ[1],:]**2.)**0.5
    ind = np.argsort(r2d)
    return r2d, ind


def make_2D_projection(filestring, snapno, units, writefilename, SEEDY=100, PROJ=(0,1)):
    #units = scripts.read_units(filestring)
    lpc = units[0]['l_pc']
    kms = 1e-5 * units[0]['l_cgs']/units[0]['nbt_cgs']
    t_myr = get_time(filestring+'.snap'+snapno+'.dat.gz')*units[0]['t_myr']

    #writefilename=filestring+'.snap'+snapno+'.2Dproj.dat'
    writefile=open(writefilename, 'w')
    writefile.write("#t=%g\n#1.r2D(pc) 2.Ltot(Lsun) 3.binflag 4.startype 5.L(Lsun) 6.startype0 7.startype1 8.L0(Lsun) 9.L1(Lsun) 10.Mtot(Msun) 11.M0(Msun) 12.M1(Msun) 13.id 14.id0 15.id1 16.rx(pc) 17.ry(pc) 18.rz(pc) 19.vx(km/s) 20.vy(km/s) 21.vz(km/s)\n" %(t_myr))

    #read the snapfile
    snapfile = filestring+'.snap'+snapno+'.dat.gz'
    colnos = (2, 7, 14, 15, 17, 18, 19, 20, 1, 8, 9, 3, 4, 0, 10, 11)
    #0-r, 1-binflag, 2-startype, 3-L, 4-startype0, 5-startype1, 6-L0, 7-L1, 8-Mtot, 9-M0, 10-M1, 11-vr, 12-vt, 13-id, 14-id0, 15-id1
    data = np.genfromtxt(snapfile, usecols=colnos)
    #data = np.genfromtxt(snapfile)
    r = data[:,0]*lpc
    vr = data[:,11]*kms
    vt = data[:,12]*kms
    r3d, v3d = convert_to_3d(r, vr, vt, SEEDY=SEEDY)
    r2d, ind = project_and_radially_sort(r3d, PROJ=PROJ)
    valid_line=0
    print('N:', len(ind))
    for i in range(len(ind)):
        try:
            for j in range(len(data[ind[i]])):
                if str(data[ind[i],j])=='nan' or str(data[ind[i],j])=='inf':
                    valid_line = 0
                    print(data[ind[i],:])
                    raise StopIteration()
                else:
                    valid_line = 1
        except StopIteration:
            pass
        if valid_line==1:
            if data[ind[i],1]==1.:
                Ltot = data[ind[i],6]+data[ind[i],7]
                Mtot = data[ind[i],9]+data[ind[i],10]
            else:
                Ltot = data[ind[i],3]
                Mtot = data[ind[i],8]
            writefile.write("%g %g %g %g %g %g %g %g %g %g %g %g %d %d %d %g %g %g %g %g %g\n" %(r2d[ind[i]], Ltot, data[ind[i],1], data[ind[i],2], Ltot, data[ind[i],4], data[ind[i],5], data[ind[i],6], data[ind[i],7], Mtot, data[ind[i],9], data[ind[i],10], data[ind[i],13], data[ind[i],14], data[ind[i],15], r3d[0,ind[i]], r3d[1,ind[i]], r3d[2,ind[i]], v3d[0,ind[i]], v3d[1,ind[i]], v3d[2,ind[i]], ))
    writefile.close()
    return r, vr, vt, r3d, v3d


def fit_king_cum_func(x, p0, p1):
    Lcum = np.pi * p0 * p1 * p1 * np.log( np.abs(1 + (x/p1)**2.))
    return Lcum


def fit_king_cum_curvefit(r2d, Lcum, p0guess=[1e5, 1.]):
    #print 'came here'
    p_opt, p_cov = opt.curve_fit(fit_king_cum_func, r2d, Lcum, p0=p0guess)
    #print 'came here too'
    return p_opt, p_cov



def get_obs_props(filestring, snapno, FAC=1.):
    filename = filestring+'.snap'+snapno+'.2Dproj.dat.gz'
    t = get_time(filestring+'.snap'+snapno+'.2Dproj.dat.gz')
    #with gzip.open(filename, 'r') as f:
    #    try:
    #        for line in f:
    #            if line.rfind("#t=")>-1:
    #                t = float(line.split('=')[1].split()[0])
    #                raise StopIteration()
    #    except StopIteration:
    #        pass
    data = np.genfromtxt(filename, usecols=(0,1,9,18,19,20,))
    print(len(data))
    Mtot = np.sum(data[:,2])
    Mave = np.mean(data[:,2])
    Lcum = np.cumsum(data[:,1])
    Ltot = Lcum[-1]
    halfL = 0.5*Ltot
    try:
        for i in range(len(data)-1):
            L1, L2 = Lcum[i], Lcum[i+1]
            r1, r2 = data[i,0], data[i+1,0]
            if L1<=halfL<=L2:
                rhl = r1 + (r2-r1)*(halfL-L1)/(L2-L1)
                drhl = 0.5*(r2-r1)
                raise StopIteration()
    except StopIteration:
        N = i
        print('found rhl:', rhl, drhl, i, N)

    #guess_sigma = np.sum(data[:N, 1])/ (np.pi*data[N-1,0]**2.)
    #print N
    guess_sigma = Lcum[N]/(np.pi*data[N-1,0]**2.)
    guess_rc = 0.5*rhl
    p0guess = [guess_sigma, guess_rc]
    lim = rhl*FAC
    rselect, Lcumselect, count = [], [], 0
    while count<len(data) and data[count,0]<=lim:
        rselect.append(data[count,0])
        Lcumselect.append(Lcum[count])
        count+=1
        #print count
    rselect = np.array(rselect)
    Lcumselect = np.array(Lcumselect)
    p_opt, p_cov = fit_king_cum_curvefit(rselect, Lcumselect, p0guess=p0guess)
    p_err = np.sqrt(np.diag(p_cov))
    sigmac, sigmacerr, rc, rcerr = np.abs(p_opt[0]), np.abs(p_err[0]), np.abs(p_opt[1]), np.abs(p_err[1])
    print('sigmac:', sigmac, 'sigmacerr:', sigmacerr, 'rc:', rc, 'rcerr:', rcerr)
    vselect1, vselect2, vselect3 = [], [], []
    count = 0
    #print count, data[count,0]
    while data[count,0]<=rc and data[count,0]<=rhl:
        vselect1.append(data[count,3])
        vselect2.append(data[count,4])
        vselect3.append(data[count,5])
        count+=1
        #print count
    vselect1, vselect2, vselect3 = np.array(vselect1), np.array(vselect2), np.array(vselect3)
    vsigma = [np.std(vselect1), np.std(vselect2), np.std(vselect3)]
    vsigma = np.array(vsigma) 
    vsigmac = np.mean(vsigma)
    dvsigmac = np.std(vsigma)

    props = {'t': t,
        'rhl': rhl,
        'drhl': drhl,
        'rc': rc,
        'drc': rcerr,
        'sigmac': sigmac,
        'dsigmac': sigmacerr,
        'Ltot': Ltot,
        'M/L': Mtot/Ltot,
        'Mave': Mave,
        'vsigmac_rv': vsigmac,
        'dvsigmac_rv': dvsigmac
        }

    return props



def get_sbp_from_2D_projection(filestring, snapno, BINNO=50, LCUT=15):
    filename = filestring+'.snap'+snapno+'.2Dproj.dat.gz'
    print(filename)
    projfile = gzip.open(filename, 'r')
    projfile.seek(0)
    line = projfile.readline()
    print(line.split('t=')[1].split())
    t_myr = float(line.split('t=')[1].split()[0])
    print(projfile)
    data = np.loadtxt(projfile)
    writefilename=filestring+'.snap'+snapno+'.2D_SBP.dat'
    writefile=open(writefilename, 'w')
    writefile.write("#t=%g\n#1.r2Dlow(pc) 2.r2Dmid(pc) 3.r2Dhigh 4.Sigma(L/pc^2) 5.Sigmaerr(L/pc^2) 6.Sigma_n(1/pc^2) 7.Sigma_nerr(1/pc^2)\n" %(t_myr))
    writefilename1=filestring+'.snap'+snapno+'.2D_SBPLcut'+str(LCUT)+'.dat'
    writefile1=open(writefilename1, 'w')
    writefile1.write("#t=%g\n#1.r2Dlow(pc) 2.r2Dmid(pc) 3.r2Dhigh 4.Sigma(L/pc^2) 5.Sigmaerr(L/pc^2) 6.Sigma_n(1/pc^2) 7.Sigma_nerr(1/pc^2)\n" %(t_myr))

    lr2d = np.log10(data[:,0])
    lbinsize = (lr2d[-1]-lr2d[0])/float(BINNO)
    print(lbinsize)
    n2d_prev = 0
    for i in range(1, BINNO+1):
        lsum, lsumerr, n2d, n2derr = 0., 0., 0., 0.
        lsumcut, lsumcuterr, n2dcut, n2dcuterr = 0., 0., 0., 0.
        lr_high, lr_low = lr2d[0]+i*lbinsize, lr2d[0]+(i-1)*lbinsize
        lr_mid = (lr_low+lr_high)/2.
        area = np.pi * ((10.**lr_high)**2. - (10.**lr_low)**2.)
        try:
            for j in range(int(n2d_prev), len(lr2d)):
                #print n2d_prev, lr_low, lr_high, lr2d[j]
                if lr2d[j]<lr_high and lr2d[j]>=lr_low:
                    lsum = lsum + data[j,1]
                    n2d = n2d + 1
                    if data[j,1]<LCUT:
                        lsumcut += data[j,1]
                        n2dcut += 1
                else:
                    raise StopIteration()
        except StopIteration:
            print('got value:n2d=', n2d, 'n2d_prev=', n2d_prev)
        n2d_prev += n2d
        if n2d>2:
            sbp, sbperr = lsum/area, lsum/float(n2d)**0.5/area  
            snp, snperr = n2d/area, float(n2d)**0.5/area
            writefile.write('%g %g %g %g %g %g %g\n' %(10**lr_low, 10**lr_mid, 10**lr_high, sbp, sbperr, snp, snperr))
            if n2dcut>2:
                ############3   
                sbpcut, sbpcuterr = lsumcut/area, lsumcut/float(n2d)**0.5/area  
                snpcut, snpcuterr = n2dcut/area, float(n2dcut)**0.5/area
                ################
                writefile1.write('%g %g %g %g %g %g %g\n' %(10**lr_low, 10**lr_mid, 10**lr_high, sbpcut, sbpcuterr, snpcut, snpcuterr))

    projfile.close()
    writefile.close()
    writefile1.close()
    return t_myr


def get_sbp_from_2D_projection_ncut(filestring, snapno, BINNO=50, LCUT=15, NCUT=0.9):
    filename = filestring+'.snap'+snapno+'.2Dproj.dat.gz'
    print(filename)
    projfile = gzip.open(filename, 'r')
    projfile.seek(0)
    line = projfile.readline()
    line = line.decode()
    #print(line.split('t=')[1].split())
    t_myr = float(line.split('t=')[1].split()[0])
    print(projfile)
    data = np.loadtxt(projfile)
    writefilename=filestring+'.snap'+snapno+'.2D_SBP_NCUT'+str(NCUT)+'.dat'
    writefile=open(writefilename, 'w')
    writefile.write("#t=%g\n#1.r2Dlow(pc) 2.r2Dmid(pc) 3.r2Dhigh 4.Sigma(L/pc^2) 5.Sigmaerr(L/pc^2) 6.Sigma_n(1/pc^2) 7.Sigma_nerr(1/pc^2) 8.N_density(1/pc^2) 9.N_densityerr(1/pc^2)\n" %(t_myr))

    writefilename1=filestring+'.snap'+snapno+'.2D_SBPLcut'+str(LCUT)+'_NCUT'+str(NCUT)+'.dat'
    writefile1=open(writefilename1, 'w')
    writefile1.write("#t=%g\n#1.r2Dlow(pc) 2.r2Dmid(pc) 3.r2Dhigh 4.Sigma(L/pc^2) 5.Sigmaerr(L/pc^2) 6.Sigma_n(1/pc^2) 7.Sigma_nerr(1/pc^2) 8.N_density(1/pc^2) 9.N_densityerr(1/pc^2)\n" %(t_myr))


    lr2d = np.log10(data[:,0])
    lbinsize = (lr2d[-1]-lr2d[0])/float(BINNO)
    print(lbinsize)
    mtot=data[:,9]; m0=data[:,10]; m1=data[:,11]; biflag=data[:,2]; ksin=data[:,3]; k0=data[:,5]; k1=data[:,6]
    n2d_prev = 0
    for i in range(1, BINNO+1):
        lsum, lsumerr, n2d, n2derr = 0., 0., 0., 0.
        lsumcut, lsumcuterr, n2dcut, n2dcuterr = 0., 0., 0., 0.
        lr_high, lr_low = lr2d[0]+i*lbinsize, lr2d[0]+(i-1)*lbinsize
        lr_mid = (lr_low+lr_high)/2.
        area = np.pi * ((10.**lr_high)**2. - (10.**lr_low)**2.)
        nd, ndcut = 0., 0.

        try:
            for j in range(int(n2d_prev), len(lr2d)):
                #print n2d_prev, lr_low, lr_high, lr2d[j]
                if lr2d[j]<lr_high and lr2d[j]>=lr_low:
                    lsum = lsum + data[j,1]
                    n2d = n2d + 1
                    if NCUT!=-1 and biflag[j]!=1 and ksin[j]<10 and mtot[j]>=NCUT:
                        nd += 1
                    if NCUT!=-1 and biflag[j]==1 and ((k0[j]<10 and m0[j]>=NCUT) or (k1[j]<10 and m1[j]>=NCUT)):
                        nd += 1
                    if data[j,1]<LCUT:
                        lsumcut += data[j,1]
                        n2dcut += 1
                        if NCUT!=-1 and biflag[j]!=1 and ksin[j]<10 and mtot[j]>=NCUT:
                            ndcut += 1 
                        if NCUT!=-1 and biflag[j]==1 and ((k0[j]<10 and m0[j]>=NCUT) or (k1[j]<10 and m1[j]>=NCUT)):
                            ndcut += 1
                                                                                                                                          
                else:
                    raise StopIteration()
        except StopIteration:
            print('got value:n2d=', n2d, 'n2d_prev=', n2d_prev, 'nd=', nd)
        n2d_prev += n2d
        if n2d>2:
            sbp, sbperr = lsum/area, lsum/float(n2d)**0.5/area  
            snp, snperr = n2d/area, float(n2d)**0.5/area
            sdp, sdperr = nd/area, float(nd)**0.5/area
            writefile.write('%g %g %g %g %g %g %g %g %g\n' %(10**lr_low, 10**lr_mid, 10**lr_high, sbp, sbperr, snp, snperr, sdp, sdperr))
            if n2dcut>2:
                ############3   
                sbpcut, sbpcuterr = lsumcut/area, lsumcut/float(n2d)**0.5/area  
                snpcut, snpcuterr = n2dcut/area, float(n2dcut)**0.5/area
                sdpcut, sdpcuterr = ndcut/area, float(ndcut)**0.5/area                              
                ################
                writefile1.write('%g %g %g %g %g %g %g %g %g\n' %(10**lr_low, 10**lr_mid, 10**lr_high, sbpcut, sbpcuterr, snpcut, snpcuterr, sdpcut, sdpcuterr))

    projfile.close()
    writefile.close()
    writefile1.close()
    return t_myr


def velocity_dispersion_snap(path,string,snapno,ALL=0,Starinbin=200,mcut=0):
    units=read_units(path+string)
    km = units[0]['l_cgs']*1.0e-5
    time = units[0]['nbt_cgs']
    f = open('RV_model.dat','w')
    f1 = open('RV_model_sigma.dat','w')
    if ALL == 1:
        f4 = open(path+string+'.snap'+snapno+'.vel_dispersion_tight.dat','w')
    if ALL == 0:
        f4 = open(path+string+'.snap'+snapno+'.vel_dispersion_giants_vr_pm_'+str(Starinbin)+'_'+str(mcut)+'.dat','w')
    f2 = open('RV_obs_sigma.dat','w')
    f3 = open('giants.dat','w')
    print('#0)r2D(pc) 1)sigma_v(1D; km/s) 2)sigma_v_err(km/s) 3)sigma_pmr(1D; km/s) 4)sigma_pmr_err(km/s) 5)sigma_pmt(1D; km/s) 6)sigma_pmt_err(km/s) 7)sigma_pm(km/s) 8)sigma_pm_err(km/s) 9)sigma_pm2(km/s) 10)sigma_pm2_err(km/s)', file = f4)
####################
    f55 = gzip.open(path+string+'.snap'+snapno+'.dat.gz','r')
    lines55 = f55.readlines()
    #data = np.genfromtxt(path+string+'.snap'+snapno+'.dat.gz')
#####################
###################
    Vr = []; Vpm_r = []; Vpm_t = []; Vpm = []; Vpm2 = []
    R = []
    bin_count = 0
    count = 0
    bin_array = []
    for i in range(2,len(lines55)):
        line55 = lines55[i]
        data = line55.split()
        for k in range(0,21):
            data[k] = np.float(data[k])         
        bin = 0
        if ALL == 1:
            if data[1] <= 1e6:  # Looks at all stars, not a cut.
            #if data[i,3] <= 1.8 and data[i,3] >= -1.2:  # asks if in V-band mag range given by Caretta et al. 2009.
                r_xyz = [data[15], data[16], data[17]]  #x,y,z components of r in pc
                v_xyz = [data[18], data[19], data[20]]  #x,y,z components of velocity in km/s
                count = count + 1
                if data[2] == 1:
                    bin_count = bin_count + 1
                    bin = 1
                bin_array.append(bin)
                r = np.sqrt(r_xyz[0]**2. + r_xyz[1]**2.)
                #R_temp = conversions.pc_to_arcsec(r,d_heliocentric)  #Converts radius in pc to arcsec
                R.append(r)
                Vr.append(data[20])     # Just use the z component as your LOS direction
        if ALL == 0:
            v_r = data[3]*km/time # converts from code units to km/s
            v_t = data[4]*km/time
            r_km = data[2]*km  #convert r from code units to km
            if data[7] == 1:
               #if data[i,17] >= 2 and data[i,17] <= 9:
               #    for k in range(0,70):
               #        r,v = convert_to_3d(r_km, v_r, v_t)
               #        Vr.append(v[0])    # Use this if you want 1d vel dispersion
               #        r = np.sqrt(r[1]**2. + r[2]**2.)/3.086e13  # convert r from km to pc
               #        R.append(r)
               #        count = count + 1
                
               #if data[i,18] >= 2 and data[i,18] <= 9:
               #    for k in range(0,70):
               #        r,v = convert_to_3d(r_km, v_r, v_t)
               #        Vr.append(v[0])    # Use this if you want 1d vel dispersion
               #        r = np.sqrt(r[1]**2. + r[2]**2.)/3.086e13  # convert r from km to pc
               #        R.append(r)
               #        count = count + 1     

                if mcut == 0 and (2 <= data[17] <= 9 or 2 <= data[18] <= 9):
                   #for k in range(0,70):
                    r,v = convert_to_3d(r_km, v_r, v_t)
                    Vr.append(v[0])    # Use this if you want 1d vel dispersion
                    Vpm_r.append(v[1]); Vpm_t.append(v[2])
                    Vpm.append(np.sqrt(v[1]**2 + v[2]**2))
                    angle_pm = np.random.uniform(0, np.pi)
                    Vpm2.append(v_t*np.cos(angle_pm))
                    r = np.sqrt(r[1]**2. + r[2]**2.)/3.086e13  # convert r from km to pc
                    R.append(r)
                    count = count + 1

                elif (data[17]<10 and data[8]>=mcut) or (data[18]<10 and data[9]>=mcut):
                    r,v = convert_to_3d(r_km, v_r, v_t)
                    Vr.append(v[0])    # Use this if you want 1d vel dispersion
                    Vpm_r.append(v[1]); Vpm_t.append(v[2])
                    Vpm.append(np.sqrt(v[1]**2 + v[2]**2))
                    angle_pm = np.random.uniform(0, np.pi)
                    Vpm2.append(v_t*np.cos(angle_pm))
                    r = np.sqrt(r[1]**2. + r[2]**2.)/3.086e13  # convert r from km to pc
                    R.append(r)
                    count = count + 1



            if data[7] != 1:
                if mcut == 0 and 2 <= data[14] <= 9:
                    for k in range(0,1):
                    #for k in range(0,25):
                        r,v = convert_to_3d(r_km, v_r, v_t)
                        Vr.append(v[0])    # Use this if you want 1d vel dispersion
                        Vpm_r.append(v[1]); Vpm_t.append(v[2])
                        Vpm.append(np.sqrt(v[1]**2 + v[2]**2))
                        angle_pm = np.random.uniform(0, np.pi)
                        Vpm2.append(v_t*np.cos(angle_pm))
                        r = np.sqrt(r[1]**2. + r[2]**2.)/3.086e13  # convert r from km to pc
                        R.append(r)
                        count = count + 1

                elif data[14]<10 and data[1]>=mcut:
                    r,v = convert_to_3d(r_km, v_r, v_t)
                    Vr.append(v[0])    # Use this if you want 1d vel dispersion
                    Vpm_r.append(v[1]); Vpm_t.append(v[2])
                    Vpm.append(np.sqrt(v[1]**2 + v[2]**2))
                    angle_pm = np.random.uniform(0, np.pi)
                    Vpm2.append(v_t*np.cos(angle_pm))
                    r = np.sqrt(r[1]**2. + r[2]**2.)/3.086e13  # convert r from km to pc
                    R.append(r)
                    count = count + 1

################################################
    mean_model = np.mean(Vr) #Find global mean of all model RVs
    mean_model_pmr = np.mean(Vpm_r)
    mean_model_pmt = np.mean(Vpm_t)
    mean_model_pm = np.mean(Vpm)
    mean_model_pm2 = np.mean(Vpm2)
    print('mean=',mean_model, mean_model_pm, mean_model_pm2)
    array = np.zeros((len(Vr),6))
    for i in range(0,len(Vr)):
        array[i,0] = R[i]
        array[i,1] = Vr[i]
        array[i,2] = Vpm_r[i]
        array[i,3] = Vpm_t[i]
        array[i,4] = Vpm[i]
        array[i,5] = Vpm2[i]
    array = array[array[:,0].argsort()]  # sorts each measurement in order of radial position
    for i in range(0,len(Vr)):
        print(array[i,0], array[i,1], array[i,2], array[i,3], array[i,4], array[i,5], file = f) #Print each radius/LOS velocity pair in order of radial position   

#####################################   
    sigma_vel_array = []
    sigma_pmr_array = []
    sigma_pmt_array = []
    sigma_pm_array = []
    sigma_pm2_array = []
    R_array = []
    #mean_array = []
    #flag = 0
    sum_vel = 0
    sum_pmr = 0
    sum_pmt = 0
    sum_pm = 0; sum_pm2 = 0
    print(snapno, file = f1)
        ## Define each bin as having 2000 stars. Every 2000 stars, start new bin

    vel_array = []  # Makes an array with velocites of all stars within each bin
    pmr_array = []; pmt_array = []; pm_array = []; pm2_array = []
    r_array = []
    #bin_count = 0
    total_count = 0
    for j in range(0,len(array)):   
        if total_count <= Starinbin:
            vel_array.append(array[j,1])
            pmr_array.append(array[j,2])
            pmt_array.append(array[j,3])
            pm_array.append(array[j,4])
            pm2_array.append(array[j,5])
            r_array.append(array[j,0])
            total_count = total_count + 1
        else:
            count = 0.
            for k in range(0,len(vel_array)):
                r = np.mean(r_array)
                sum_vel = sum_vel + (vel_array[k]-mean_model)**2.
                sum_pmr = sum_pmr + (pmr_array[k]-mean_model_pmr)**2.
                sum_pmt = sum_pmt + (pmt_array[k]-mean_model_pmt)**2.
                sum_pm = sum_pm + (pm_array[k]-mean_model_pm)**2.
                sum_pm2 = sum_pm2 + (pm2_array[k]-mean_model_pm2)**2.
                count = count + 1.
            sigma_vel = np.sqrt(sum_vel/count)
            error_vel = np.sqrt(sigma_vel**2.0/(2.*count))
            sigma_pmr = np.sqrt(sum_pmr/count)
            error_pmr = np.sqrt(sigma_pmr**2.0/(2.*count))
            sigma_pmt = np.sqrt(sum_pmt/count)
            error_pmt = np.sqrt(sigma_pmt**2.0/(2.*count))
            sigma_pm = np.sqrt(sum_pm/count)
            error_pm = np.sqrt(sigma_pm**2.0/(2.*count))
            sigma_pm2 = np.sqrt(sum_pm2/count)
            error_pm2 = np.sqrt(sigma_pm2**2.0/(2.*count))

            print(r, sigma_vel, error_vel, sigma_pmr, error_pmr, sigma_pmt, error_pmt, sigma_pm, error_pm, sigma_pm2, error_pm2, file = f1)
            print(r, sigma_vel, error_vel, sigma_pmr, error_pmr, sigma_pmt, error_pmt, sigma_pm, error_pm, sigma_pm2, error_pm2, file = f4)
            #print r, 'sigma=',sigma, 'error=',error, 'count=',count, 'N_binaries=',bin_count,'N_stars=',total_count,'binary fraction=',float(bin_count)/float(total_count)
            sigma_vel_array.append(sigma_vel)
            sigma_pmr_array.append(sigma_pmr)
            sigma_pmt_array.append(sigma_pmt)
            sigma_pm_array.append(sigma_pm)
            sigma_pm2_array.append(sigma_pm2)
            R_array.append(np.mean(r_array))
            sum_vel = 0; sum_pmr = 0; sum_pmt = 0; sum_pm = 0; sum_pm2 = 0
            #flag = 0
            #bin_count = 0
            total_count = 0
            vel_array = []; pmr_array = []; pmt_array = []; pm_array = []; pm2_array = []
            r_array = []
    print(' ', file = f1)



def velocity_dispersion_2dsnap(path, snapno, Starinbin=200, mcut=0):
    fvel = open(path+'initial.snap'+snapno+'.vel_dispersion_vr_pm_'+str(Starinbin)+'_'+str(mcut)+'.dat','w+')
    print('#0)r2D(pc) 1)sigma_v(1D; km/s) 2)sigma_v_err(km/s) 3)sigma_pmr(1D; km/s) 4)sigma_pmr_err(km/s) 5)sigma_pmt(1D; km/s) 6)sigma_pmt_err(km/s)', file = fvel)

    binflag = []; VX = []; VY = []; VZ = []
    rd = []
    ktype = []; k0 = []; k1 = []
    m = []; m0 = []; m1 = []
    snap2d = path+'initial.snap'+snapno+'.2Dproj.dat.gz'
    with gzip.open(snap2d, 'r') as f2d:
        next(f2d); next(f2d)
        for line in f2d:
            data = line.split()
            binflag.append(int(data[2]))
            VX.append(float(data[18])); VY.append(float(data[19])); VZ.append(float(data[20]))
            rd.append(float(data[0]))
            ktype.append(int(data[3])); k0.append(int(data[5])); k1.append(int(data[6]))
            m.append(float(data[9])); m0.append(float(data[10])); m1.append(float(data[11]))

    print('read snap data')
    
    vel_r = []; vel_pmr = []; vel_pmt = []
    rpc = []
    for xx in range(len(ktype)):
        if binflag[xx] == 1:
            if mcut == 0 and 2<=k0[xx]<=9 or 2<=k1[xx]<=9:
                rpc.append(rd[xx])
                vel_r.append(VX[xx]); vel_pmr.append(VY[xx]); vel_pmt.append(VZ[xx])

            elif (k0[xx]<10 and m0[xx]>=mcut) or (k1[xx]<10 and m1[xx]>=mcut):
                rpc.append(rd[xx])
                vel_r.append(VX[xx]); vel_pmr.append(VY[xx]); vel_pmt.append(VZ[xx])

        else:
            if mcut == 0 and 2<=ktype[xx]<=9:
                rpc.append(rd[xx])
                vel_r.append(VX[xx]); vel_pmr.append(VY[xx]); vel_pmt.append(VZ[xx])

            elif ktype[xx]<10 and m[xx]>=mcut:
                rpc.append(rd[xx])
                vel_r.append(VX[xx]); vel_pmr.append(VY[xx]); vel_pmt.append(VZ[xx])


    array = np.zeros((len(vel_r),4))
    for i in range(0,len(vel_r)):
        array[i,0] = rpc[i]
        array[i,1] = vel_r[i]
        array[i,2] = vel_pmr[i]
        array[i,3] = vel_pmt[i]
    array = array[array[:,0].argsort()]  # sorts each measurement in order of radial position

    mean_r = np.mean(vel_r); mean_pmr = np.mean(vel_pmr); mean_pmt = np.mean(vel_pmt)
    print(mean_r, mean_pmr, mean_pmt)


    sigma_r_array = []; sigma_pmr_array = []; sigma_pmt_array = []; R_array = []
    velr_array = []; velpmr_array = []; velpmt_array = []
    r_array = []
    sum_r = 0; sum_pmr = 0; sum_pmt = 0

    total_count = 0
    for ii in range(0, len(array)):
        if total_count <= Starinbin:
            velr_array.append(array[ii,1])
            velpmr_array.append(array[ii,2])
            velpmt_array.append(array[ii,3])
            r_array.append(array[ii,0])
            total_count = total_count + 1
        else:
            count = 0.
            for k in range(0,len(velr_array)):
                r = np.mean(r_array)
                sum_r = sum_r + (velr_array[k]-mean_r)**2.
                sum_pmr = sum_pmr + (velpmr_array[k]-mean_pmr)**2.
                sum_pmt = sum_pmt + (velpmt_array[k]-mean_pmt)**2.
                count = count + 1.
            sigma_r = np.sqrt(sum_r/count)
            error_r = np.sqrt(sigma_r**2.0/(2.*count))
            sigma_pmr = np.sqrt(sum_pmr/count)
            error_pmr = np.sqrt(sigma_pmr**2.0/(2.*count))
            sigma_pmt = np.sqrt(sum_pmt/count)
            error_pmt = np.sqrt(sigma_pmt**2.0/(2.*count))

            print(r, sigma_r, error_r, sigma_pmr, error_pmr, sigma_pmt, error_pmt, file = fvel)
        
            sigma_r_array.append(sigma_r)
            sigma_pmr_array.append(sigma_pmr)
            sigma_pmt_array.append(sigma_pmt)
            R_array.append(np.mean(r_array))

            sum_r = 0; sum_pmr = 0; sum_pmt = 0
            total_count = 0
            velr_array = []; velpmr_array = []; velpmt_array = []
            r_array = []

    #return sigma_r_array, sigma_pmr_array, sigma_pmt_array, R_array


def main():
    sourcepath = '/projects/b1095/syr904/cmc/47Tuc/rundir/47Tuc/MOCHA47Tuc_150maxmass_rv1.8/'
    N=2000000
    Z=0.0038
    rv=1.8
    rg=7.4

    string = 'initial'
    units=read_units(sourcepath+string)
    km = units[0]['l_cgs']*1.0e-5
    time_units = units[0]['t_myr']
    m_cgs = units[0]['l_cgs']*1.0e-2
    kg = units[0]['m_cgs']*1.0e-3
    time_cgs = units[0]['t_cgs']
    nbt = units[0]['nbt_cgs']
    M_total_units = units[0]['m_msun']
    pc = units[0]['l_pc']

    print(sourcepath)
    time_array, snap_array = find_snap_time_array(sourcepath,string)
    print(snap_array)

    Delta = -3  #default -5000 for only making the last snapshot

    for k in range(len(time_array)-1,-1,Delta):
        time_snap = time_array[k]
        no_snap = snap_array[k]
        print('time=', time_snap, 'snapno=', no_snap)
        if time_snap > 14000.: 
            continue
        if time_snap < 10000.:
            print('stop!')
            break
        try:
            f = open(sourcepath+'initial.snap'+no_snap+'.vel_dispersion_vr_pm_700_0.75.dat','r')
            print(no_snap, 'is done')
            continue
        except:
            try:
                print('start')
                filename = sourcepath+'initial.snap'+no_snap+'.2Dproj.dat'
                make_2D_projection(sourcepath+'initial', no_snap, units, filename)
                os.system('gzip '+sourcepath+'initial.snap'+no_snap+'.2Dproj.dat')
                print('made 2D projection')

                ##### make cluster params file
                f2 = open(sourcepath+'initial.snap'+no_snap+'.cluster_params.dat','w')
                f5 = gzip.open(sourcepath+'initial.snap'+no_snap+'.dat.gz','r')
                lines5 = f5.readlines()
                props = get_obs_props(sourcepath+'initial', no_snap, FAC=1.)
                print('props=', props)
                rc = props['rc']
                print('rc=', rc)
    
                #Initialization
                count_obj=0; count=0; BHnonBH=0; BBH=0; NSnonNS=0; BNS=0    
                P=0; MS=0; G=0; WD=0; NS=0; BH=0
    
                for j in range(2,len(lines5)):
                    line5 = lines5[j]
                    data = line5.split()
                    #print data
                    #for j in range(0,len(data)-2):
                    #    if data[j] != 'na':
                    #        data[j] = np.float(data[j])   # Convert strings to floats
                    r = float(data[2])*pc
                    if r <= rc:
                        count_obj += 1
                    if int(data[7]) == 1:
                        if r <= rc:
                            count += 2
                        k1 = int(data[17])
                        k2 = int(data[18])
                        M1 = float(data[8])
                        M2 = float(data[9])
                        ID1 = int(data[10])
                        ID2 = int(data[11])
                        a = float(data[12])
                        e = float(data[13])
                        if k1 == 14 and k2 != 14:
                            BHnonBH += 1
                        if k2 == 14 and k1 != 14:
                            BHnonBH += 1
                        if k1 == 14 and k2 == 14:
                            BBH += 1
                        if k1 == 13 and k2 < 13:
                            NSnonNS += 1
                        if k2 == 13 and k1 < 13:
                            NSnonNS += 1
                        if k1 == 13 and k2 == 13:
                            BNS += 1
                    else:
                        if r <= rc:
                            count += 1
                        M = float(data[1])
                        k = int(data[14])
                        if k == 0 and M == 0.001:
                            P += 1
                        if k <= 1 and M > 0.01:
                            MS += 1
                        if k>= 10 and k <= 12:
                            WD += 1
                        if k == 13:
                            NS += 1
                        if k == 14:
                            BH += 1
                        if k >= 2 and k <= 9:
                            G += 1

                print('end of loop')
                number_density = count/(rc**3.)
                number_density2 = count_obj/(rc**3.)
                ##print 'number_density', number_density
    
                print("#time, number_density, props['Ltot'], props['M/L'], props['Mave'], props['drc'], props['drhl'], props['dsigmac'], props['dvsigmac_rv'], props['rc'], props['rhl'], props['sigmac'],  props['vsigmac_rv'], number_density2", file = f2)
                print(time_snap, number_density, props['Ltot'], props['M/L'], props['Mave'], props['drc'], props['drhl'], props['dsigmac'], props['dvsigmac_rv'], props['rc'], props['rhl'], props['sigmac'],  props['vsigmac_rv'], number_density2, N, Z, rv, rg, file = f2)
                print(time_snap, number_density, props['Ltot'], props['M/L'], props['Mave'], props['drc'], props['drhl'], props['dsigmac'], props['dvsigmac_rv'], props['rc'], props['rhl'], props['sigmac'],  props['vsigmac_rv'], number_density2, N, Z, rv, rg, file = f2)
                f2.close()
                f5.close()
                print('cluster_params done')    
        
                ###############
                print('made params file')
                #get_sbp_from_2D_projection(path+string, snapno)
                get_sbp_from_2D_projection_ncut(sourcepath+string, no_snap, BINNO=50, LCUT=12, NCUT=0.85)
                #get_sbp_from_2D_projection_ncut(sourcepath+string, no_snap, BINNO=50, LCUT=12, NCUT=0.75)
                print(no_snap, 'made SBP')

                velocity_dispersion_2dsnap(sourcepath,no_snap, Starinbin=700, mcut = 0.85)
                #velocity_dispersion_2dsnap(sourcepath,no_snap, Starinbin=500, mcut = 0)
                print('Made vel dispersion for',no_snap)
            except:
                print(no_snap, 'failed')



def read_keys(thekey):
    return re.findall(r'\d+\.\d+|\d+', thekey)


##Make all relevent profiles for the new CMC-COSMIC models using hdf5
def make_SBP_VDP_NDP_hdf5_smooth(modelpath, dgc):   ##dgc in kpc
    snap_h5 = modelpath + 'initial.snapshots.h5'

    t_conv = conv('t', modelpath+'initial.conv.sh')

    with pd.HDFStore(snap_h5) as snaps:
        snap_keys = snaps.keys()

    snap_no = []; snap_time = []
    for ii in range(len(snap_keys)):
        temp_key = read_keys(snap_keys[ii])
        temp_no = int(temp_key[0])
        temp_time = float(temp_key[1])

        snap_no.append(int(temp_key[0]))
        snap_time.append(float(temp_key[1]))

        #if temp_time*t_conv < 10000. or temp_time*t_conv >13000. or temp_no%2!=0:
        #    continue

        if temp_time*t_conv < 7000. or temp_time*t_conv >9000. or temp_no%2!=0:
            continue

        else:
            temp_snap = cmct.Snapshot(fname=snap_h5, snapshot_name=snap_keys[ii], 
                        conv=modelpath+'initial.conv.sh', 
                        dist=dgc, # distance to cluster in kpc
                        z=0.0038)


            ##Make surface brightness profile
            temp_snap.add_photometry('/projects/b1095/syr904/MyCodes/cmctoolkit/filt_index.txt')
            v_bincenter, v_profile = temp_snap.make_smoothed_brightness_profile('V', bins=80,
                                                min_mass=None, max_mass=None,
                                                max_lum=12, fluxdict=None,
                                                startypes=np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9]),
                                                min_logr=-2.0)


            # Make velocity dispersion profiles
            star_velbin_center, star_veldisp_profile, star_e_veldisp_profile = temp_snap.make_smoothed_veldisp_profile(bins=80,
                                                min_mass=0.85,
                                                max_mass=None,
                                                dmax=None,
                                                fluxdict=None,
                                                startypes=np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9]),
                                                min_logr=-1.5)

            star_velbin_arcsec = uc.pc2arcsec(dgc, star_velbin_center)


            ##Make number density profile
            star_numbin_center, star_profile, star_e_profile = temp_snap.make_smoothed_number_profile(bins=80,
                                                min_mass=0.85,
                                                max_mass=None,
                                                fluxdict=None,
                                                startypes=np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9]),
                                                min_logr=-2.0)

            star_numbin_arcsec = uc.pc2arcsec(dgc, star_numbin_center)
            star_profile_arcsec = star_profile/(uc.pc2arcsec(dgc, 1.)**2)
            star_e_profile_arcsec = star_e_profile/(uc.pc2arcsec(dgc, 1.)**2)

            fsbp = modelpath+'SBP'+str(temp_no)+'_LCUT12.txt'
            fvdp = modelpath+'VDP'+str(temp_no)+'_MCUT0d85.txt'
            fndp = modelpath+'NDP'+str(temp_no)+'_MCUT0d85.txt'

            np.savetxt(fsbp, np.c_[v_bincenter, v_profile], fmt = '%f %f', header = '1.r(arcsec) 2.miu_v(mag/arcsec^2) t='+str(temp_time), comments = '#')
            np.savetxt(fvdp, np.c_[star_velbin_arcsec, star_veldisp_profile, star_e_veldisp_profile], fmt = '%f %f %f', header = '1.r(arcsec) 2.sigma_v(km/s) 3.sigma_v_err(km/s) t='+str(temp_time), comments = '#')
            np.savetxt(fndp, np.c_[star_numbin_arcsec, star_profile_arcsec, star_e_profile_arcsec], fmt = '%f %f %f', header = '1.r(arcsec) 2.num(1/arcsec^2) 3.num_err(1/arcsec^2) t='+str(temp_time), comments = '#')

            print(snap_keys[ii])


    snap_no_sort, snap_time_sort = (list(t) for t in zip(*sorted(zip(snap_no, snap_time))))

    np.savetxt(modelpath+'snap_keys.txt', np.c_[snap_no_sort, snap_time_sort], fmt = '%d %f', header = '1.snap_no 2.snap_time', comments = '#')




##Make multiple different projections
def make_different_projection_snap(modelpath, snap_no):
    string = 'initial'
    units=read_units(modelpath+string)

    for ii in range(20):
        randomseed = ii*10+10
        r1 = random.choice([0,1,2])
        if r1 == 0: r2 = 1
        elif r1 == 1: r2 = 2
        else: r2 = 0
        projection = (int(r1), int(r2))

        print(projection)

        filename = modelpath+'SNAP'+snap_no+'_'+str(ii)+'.2Dproj.dat'
        make_2D_projection(modelpath+string, '0'+snap_no, units, filename, SEEDY=randomseed, PROJ=projection)
        os.system('gzip '+filename)

        print(ii)




