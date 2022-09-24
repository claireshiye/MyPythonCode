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

def read_keys(thekey):
    return re.findall(r'\d+\.\d+|\d+', thekey)


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



def convert_to_3d(r, vr, vt, SEEDY=10):
    #np.random.seed(SEEDY)
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


def make_2D_projection(modelpath, thekey, units, writefilename, dist_sun, themetal, SEEDY=10, PROJ=(0,1)):
    #units = scripts.read_units(filestring)

    snapno = read_keys(thekey)[0]; snaptime = float(read_keys(thekey)[1])
    print(type(snapno))

    lpc = units[0]['l_pc']
    kms = 1e-5 * units[0]['l_cgs']/units[0]['nbt_cgs']
    t_myr = snaptime*units[0]['t_myr']

    #read the snapfile
    snapfile = cmct.Snapshot(fname=modelpath+'initial.snapshots.h5', 
                             snapshot_name=thekey, conv=modelpath+'initial.conv.sh', 
                             dist=dist_sun, # distance to cluster in kpc
                             z=themetal)

    print('read_snapfile')


    #writefilename=modelpath+'initial.snap'+snapno+'.2Dproj.dat'
    writefile=open(writefilename, 'w+')
    writefile.write("#t=%g\n#1.r2D(pc) 2.Ltot(Lsun) 3.binflag 4.startype 5.L(Lsun) 6.startype0 7.startype1 8.L0(Lsun) 9.L1(Lsun) 10.Mtot(Msun) 11.M0(Msun) 12.M1(Msun) 13.id 14.id0 15.id1 16.rx(pc) 17.ry(pc) 18.rz(pc) 19.vx(km/s) 20.vy(km/s) 21.vz(km/s)\n" %(t_myr))

    print('write file')

    colnos = (2, 7, 14, 15, 17, 18, 19, 20, 1, 8, 9, 3, 4, 0, 10, 11)
    #0-r, 1-binflag, 2-startype, 3-L, 4-startype0, 5-startype1, 6-L0, 7-L1, 8-Mtot, 9-M0, 10-M1, 11-vr, 12-vt, 13-id, 14-id0, 15-id1

    #data = np.genfromtxt(snapfile)
    r = snapfile.data['r']*lpc
    vr = snapfile.data['vr']*kms
    vt = snapfile.data['vt']*kms
    r3d, v3d = convert_to_3d(r, vr, vt, SEEDY=SEEDY)
    r2d, ind = project_and_radially_sort(r3d, PROJ=PROJ)
    print('N:', len(ind))

    #valid_line=0
    valid_line=1
    for i in range(len(ind)):
        #try:
        #    for j in range(len(data[ind[i]])):
        #        if str(data[ind[i],j])=='nan' or str(data[ind[i],j])=='inf':
        #            valid_line = 0
        #            print(data[ind[i],:])
        #            raise StopIteration()
        #        else:
        #            valid_line = 1
        #except StopIteration:
        #    pass
        if valid_line==1:
            if snapfile.data['binflag'][ind[i]]==1.:
                Ltot = snapfile.data['bin_star_lum0_LSUN'][ind[i]] + snapfile.data['bin_star_lum1_LSUN'][ind[i]]
                Mtot = snapfile.data['m_MSUN'][ind[i]]
            else:
                Ltot = snapfile.data['luminosity_LSUN'][ind[i]]
                Mtot = snapfile.data['m_MSUN'][ind[i]]
            writefile.write("%g %g %g %g %g %g %g %g %g %g %g %g %d %d %d %g %g %g %g %g %g\n" %(r2d[ind[i]], Ltot, snapfile.data['binflag'][ind[i]], snapfile.data['startype'][ind[i]], Ltot, snapfile.data['bin_startype0'][ind[i]], snapfile.data['bin_startype1'][ind[i]], snapfile.data['bin_star_lum0_LSUN'][ind[i]], snapfile.data['bin_star_lum1_LSUN'][ind[i]], Mtot, snapfile.data['m0_MSUN'][ind[i]], snapfile.data['m1_MSUN'][ind[i]], snapfile.data['id'][ind[i]], snapfile.data['id0'][ind[i]], snapfile.data['id1'][ind[i]], r3d[0,ind[i]], r3d[1,ind[i]], r3d[2,ind[i]], v3d[0,ind[i]], v3d[1,ind[i]], v3d[2,ind[i]]))
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


##Something is wrong but i don't know what with this function. The velocity dispersions are not isotropic in the output.
def velocity_dispersion_old(path,string,thekey,Starinbin=200,mcut=0):
    snapno = read_keys(thekey)[0]; snaptime = float(read_keys(thekey)[1])

    units=read_units(path+string)
    lpc = units[0]['l_pc']
    kms = 1e-5 * units[0]['l_cgs']/units[0]['nbt_cgs']

    f4 = open(path+string+'.snap'+snapno+'.vel_dispersion_vr_pm_'+str(Starinbin)+'_'+str(mcut)+'.dat','w+')

    print('#0)r2D(pc) 1)sigma_v(1D; km/s) 2)sigma_v_err(km/s) 3)sigma_pmr(1D; km/s) 4)sigma_pmr_err(km/s) 5)sigma_pmt(1D; km/s) 6)sigma_pmt_err(km/s) 7)sigma_pm(km/s) 8)sigma_pm_err(km/s) 9)sigma_pm2(km/s) 10)sigma_pm2_err(km/s)', file = f4)
####################
    snapfile = cmct.Snapshot(fname=path+'initial.snapshots.h5', 
                             snapshot_name=thekey, conv=path+'initial.conv.sh', 
                             dist=4.52, # distance to cluster in kpc
                             z=0.0038)

    print('read snap hdf5')
#####################
###################
    Vr = []; Vpm_r = []; Vpm_t = []; Vpm = []
    R = []
    count = 0
    print(len(snapfile.data['r']))

    v_r = snapfile.data['vr']*kms # converts from code units to km/s
    v_t = snapfile.data['vt']*kms
    r_pc = snapfile.data['r']*lpc  #convert r from code units to pc

    r_all, v_all = convert_to_3d(r_pc, v_r, v_t)

    binflag = snapfile.data['binflag']
    ktype = snapfile.data['startype']; k0 = snapfile.data['bin_startype0']; k1 = snapfile.data['bin_startype1']
    m = snapfile.data['m_MSUN']; m0 = snapfile.data['m0_MSUN']; m1 = snapfile.data['m1_MSUN']

    #print(v_all)

    for i in range(len(r_all[0])):       
        if binflag[i] == 1:
            if mcut == 0 and (2 <= k0[i] <= 9 or 2 <= k1[i] <= 9):
               #for k in range(0,70):
                Vr.append(v_all[0,i])    # Use this if you want 1d vel dispersion
                Vpm_r.append(v_all[1,i]); Vpm_t.append(v_all[2,i])
                Vpm.append(np.sqrt(v_all[1,i]**2 + v_all[2,i]**2))
                rd = np.sqrt(r_all[1,i]**2. + r_all[2,i]**2.)
                R.append(rd)
                count = count + 1

            elif (k0[i]<10 and m0[i]>=mcut) or (k1[i]<10 and m1[i]>=mcut):
                #print(snapfile.data['bin_startype0'][i], snapfile.data['m0_MSUN'][i], snapfile.data['binflag'][i])
                Vr.append(v_all[0,i])    # Use this if you want 1d vel dispersion
                Vpm_r.append(v_all[1,i]); Vpm_t.append(v_all[2,i])
                Vpm.append(np.sqrt(v_all[1,i]**2 + v_all[2,i]**2))
                rd = np.sqrt(r_all[1,i]**2. + r_all[2,i]**2.)
                R.append(rd)
                count = count + 1


        else:
            if mcut == 0 and 2 <= ktype[i] <= 9:
                Vr.append(v_all[0,i])    # Use this if you want 1d vel dispersion
                Vpm_r.append(v_all[1,i]); Vpm_t.append(v_all[2,i])
                Vpm.append(np.sqrt(v_all[1,i]**2 + v_all[2,i]**2))
                rd = np.sqrt(r_all[1,i]**2. + r_all[2,i]**2.)
                R.append(rd)
                count = count + 1

            elif ktype[i]<10 and m[i]>=mcut:
                #print(snapfile.data['startype'][i], snapfile.data['m_MSUN'][i], snapfile.data['binflag'][i])
                Vr.append(v_all[0,i])    # Use this if you want 1d vel dispersion
                Vpm_r.append(v_all[1,i]); Vpm_t.append(v_all[2,i])
                Vpm.append(np.sqrt(v_all[1,i]**2 + v_all[2,i]**2))
                rd = np.sqrt(r_all[1,i]**2. + r_all[2,i]**2.)
                R.append(rd)
                count = count + 1

    print(Vr)
    print('get vr vt done')
################################################
    mean_model = np.mean(Vr) #Find global mean of all model RVs
    mean_model_pmr = np.mean(Vpm_r)
    mean_model_pmt = np.mean(Vpm_t)
    mean_model_pm = np.mean(Vpm)
    print('mean=',mean_model, mean_model_pm)
    array = np.zeros((len(Vr),5))
    for i in range(0,len(Vr)):
        array[i,0] = R[i]
        array[i,1] = Vr[i]
        array[i,2] = Vpm_r[i]
        array[i,3] = Vpm_t[i]
        array[i,4] = Vpm[i]
    array = array[array[:,0].argsort()]  # sorts each measurement in order of radial position

#####################################   
    sigma_vel_array = []
    sigma_pmr_array = []
    sigma_pmt_array = []
    sigma_pm_array = []
    R_array = []
    #mean_array = []
    #flag = 0
    sum_vel = 0
    sum_pmr = 0
    sum_pmt = 0
    sum_pm = 0

    vel_array = []  # Makes an array with velocites of all stars within each bin
    pmr_array = []; pmt_array = []; pm_array = []
    r_array = []
    #bin_count = 0
    total_count = 0
    for j in range(0,len(array)):   
        if total_count <= Starinbin:
            vel_array.append(array[j,1])
            pmr_array.append(array[j,2])
            pmt_array.append(array[j,3])
            pm_array.append(array[j,4])
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
                count = count + 1.
            sigma_vel = np.sqrt(sum_vel/count)
            error_vel = np.sqrt(sigma_vel**2.0/(2.*count))
            sigma_pmr = np.sqrt(sum_pmr/count)
            error_pmr = np.sqrt(sigma_pmr**2.0/(2.*count))
            sigma_pmt = np.sqrt(sum_pmt/count)
            error_pmt = np.sqrt(sigma_pmt**2.0/(2.*count))
            sigma_pm = np.sqrt(sum_pm/count)
            error_pm = np.sqrt(sigma_pm**2.0/(2.*count))

            print(r, sigma_vel, error_vel, sigma_pmr, error_pmr, sigma_pmt, error_pmt, sigma_pm, error_pm, file = f4)
            #print r, 'sigma=',sigma, 'error=',error, 'count=',count, 'N_binaries=',bin_count,'N_stars=',total_count,'binary fraction=',float(bin_count)/float(total_count)
            sigma_vel_array.append(sigma_vel)
            sigma_pmr_array.append(sigma_pmr)
            sigma_pmt_array.append(sigma_pmt)
            sigma_pm_array.append(sigma_pm)
            R_array.append(np.mean(r_array))
            sum_vel = 0; sum_pmr = 0; sum_pmt = 0; sum_pm = 0
            #flag = 0
            #bin_count = 0
            total_count = 0
            vel_array = []; pmr_array = []; pmt_array = []; pm_array = []
            r_array = []

    #return sigma_vel_array, sigma_pmr_array, sigma_pmt_array, R_array


def velocity_dispersion_hdf5(modelpath, thekey, Starinbin=200, mcut=0, hdf5flag = 0):
    snapno = read_keys(thekey)[0]; snaptime = float(read_keys(thekey)[1])

    fvel = open(modelpath+'initial.snap'+snapno+'.vel_dispersion_vr_pm_'+str(Starinbin)+'_'+str(mcut)+'.dat','w+')
    print('#0)r2D(pc) 1)sigma_v(1D; km/s) 2)sigma_v_err(km/s) 3)sigma_pmr(1D; km/s) 4)sigma_pmr_err(km/s) 5)sigma_pmt(1D; km/s) 6)sigma_pmt_err(km/s)', file = fvel)

    if hdf5flag == 1:
        snap_h5 = modelpath+'initial.snapshots.h5'
        snap = cmct.Snapshot(fname=snap_h5, snapshot_name=thekey, conv=modelpath+'initial.conv.sh', 
                            dist=4.52, # distance to cluster in kpc
                            z=0.0038)
        snap.make_2d_projection(seed=8675309)
        VX = snap.data['vx[KM/S]']; VY = snap.data['vy[KM/S]']; VZ = snap.data['vz[KM/S]']
        rd = snap.data['d[PC]']
        ktype = snap.data['startype']; k0 = snap.data['bin_startype0']; k1 = snap.data['bin_startype1']
        m = snap.data['m_MSUN']; m0 = snap.data['m0_MSUN']; m1 = snap.data['m1_MSUN']
        binflag = snap.data['binflag']
        print('read snap data')

    if hdf5flag == 0:
        binflag = []; VX = []; VY = []; VZ = []
        rd = []
        ktype = []; k0 = []; k1 = []
        m = []; m0 = []; m1 = []
        snap2d = modelpath+'initial.snap'+snapno+'.2Dproj.dat.gz'
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



def velocity_dispersion_2dsnap(path, snapno, Starinbin=200, mcut=0):
    fvel = open(path+'.snap'+snapno+'.vel_dispersion_vr_pm_'+str(Starinbin)+'_'+str(mcut)+'.dat','w+')
    print('#0)r2D(pc) 1)sigma_v(1D; km/s) 2)sigma_v_err(km/s) 3)sigma_pmr(1D; km/s) 4)sigma_pmr_err(km/s) 5)sigma_pmt(1D; km/s) 6)sigma_pmt_err(km/s)', file = fvel)

    binflag = []; VX = []; VY = []; VZ = []
    rd = []
    ktype = []; k0 = []; k1 = []
    m = []; m0 = []; m1 = []
    snap2d = path+'.snap'+snapno+'.2Dproj.dat.gz'
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
    


def main(sourcepath, N, Z, rv, rg, thedist, tlimlow, tlimhigh, deltastep):  ##tlimlow, tlimhigh in Myr
    #N=2950000
    #Z=0.0038
    #rv=3.5
    #rg=7.4

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

    prefix = 'initial'
    snap_h5 = prefix+'.snapshots.h5'
    #with pd.HDFStore(sourcepath+snap_h5) as snap_hdf:
    #    snap_keys = snap_hdf.keys()

    #time_array = []; snap_array = []
    #for ii in range(2, len(snap_keys)):
    #    temp_no, temp_t = read_keys(snap_keys[ii])
    #    time_array.append(time_units*float(temp_t))
    #    snap_array.append(temp_no)
    #print(time_array)

    #time_sort, keys_sort = (list(t) for t in zip(*sorted(zip(time_array,snap_keys))))
    #time_sort, no_sort = (list(t) for t in zip(*sorted(zip(time_array,snap_array))))

    data_snapkey = np.genfromtxt(sourcepath+'snap_keys.txt', dtype = str)
    temp_t = data_snapkey[:,1]
    time_sort = np.array([float(s) for s in data_snapkey[:,1]])*time_units
    no_sort = data_snapkey[:,0]
    print(no_sort, time_sort)
   
    print(sourcepath)
    Delta = deltastep  #default -5000 for only making the last snapshot

    for k in range(len(time_sort)-1,-1,Delta):
        time = time_sort[k]
        no_snap = no_sort[k]
        #key_snap = keys_sort[k]
        key_snap = '/'+no_snap+'(t='+temp_t[k]+')'
        print('time=', time, 'snapno=', no_snap, 'key_snap=', key_snap)
        if time > tlimhigh: 
            continue
        if time < tlimlow:
            print('stop!')
            break
        try:
            f = open(sourcepath+prefix+'.snap'+no_snap+'.cluster_params.dat','r')
            #f = open(sourcepath+'initial.snap'+no_snap+'.vel_dispersion_vr_pm_50_0.dat','r')
            #f = open(sourcepath+'initial.snap'+no_snap+'.test','r')
            print(no_snap, 'is done')
            continue
        except:
            try:
                make_2D_projection(sourcepath, key_snap, units, sourcepath+prefix+'.snap'+no_snap+'.2Dproj.dat', thedist, Z)
                os.system('gzip '+sourcepath+prefix+'.snap'+no_snap+'.2Dproj.dat')
                print('made 2D projection')

                ###### make cluster params file
                f2 = open(sourcepath+prefix+'.snap'+no_snap+'.cluster_params.dat','w')
                props = get_obs_props(sourcepath+prefix, no_snap, FAC=1.)
                print('props=', props)
                rc = props['rc']
                print('rc=', rc)
    
                ##Initialization
                count_obj=0; count=0
                #BHnonBH=0; BBH=0; NSnonNS=0; BNS=0    
                #P=0; MS=0; G=0; WD=0; NS=0; BH=0
                
                thesnap = cmct.Snapshot(fname=sourcepath+snap_h5, snapshot_name=key_snap, 
                    conv=sourcepath+'initial.conv.sh', 
                    dist=thedist, # distance to cluster in kpc
                    z=Z)
                rsnap = thesnap.data['r']; binflag = thesnap.data['binflag']
                mass = thesnap.data['m_MSUN']
                mbin0 = thesnap.data['m0_MSUN']; mbin1 = thesnap.data['m1_MSUN']
                ksin = thesnap.data['startype']
                kbin0 = thesnap.data['bin_startype0']; kbin1 = thesnap.data['bin_startype1']

                for j in range(2,len(rsnap)):
                    r = float(rsnap[j])*pc
                    if r <= rc:
                        count_obj += 1
                    if int(binflag[j]) == 1:
                        if r <= rc:
                            count += 2
                        #k1 = int(kbin0[j])
                        #k2 = int(kbin1[j])
                        #M1 = float(mbin0[j])
                        #M2 = float(mbin1[j])
                        #if k1 == 14 and k2 != 14:
                        #    BHnonBH += 1
                        #if k2 == 14 and k1 != 14:
                        #    BHnonBH += 1
                        #if k1 == 14 and k2 == 14:
                        #    BBH += 1
                        #if k1 == 13 and k2 < 13:
                        #    NSnonNS += 1
                        #if k2 == 13 and k1 < 13:
                        #    NSnonNS += 1
                        #if k1 == 13 and k2 == 13:
                        #    BNS += 1
                    else:
                        if r <= rc:
                            count += 1
                        #M = float(mass[j])
                        #kstar = int(ksin[j])
                        #if kstar == 0 and M == 0.001:
                        #    P += 1
                        #if kstar <= 1 and M > 0.01:
                        #    MS += 1
                        #if kstar>= 10 and k <= 12:
                        #    WD += 1
                        #if kstar == 13:
                        #    NS += 1
                        #if kstar == 14:
                        #    BH += 1
                        #if kstar >= 2 and k <= 9:
                        #    G += 1

                print('end of loop')
                number_density = count/(rc**3.)
                number_density2 = count_obj/(rc**3.)
                #print 'number_density', number_density
    
                print("#time, number_density, props['Ltot'], props['M/L'], props['Mave'], props['drc'], props['drhl'], props['dsigmac'], props['dvsigmac_rv'], props['rc'], props['rhl'], props['sigmac'],  props['vsigmac_rv'], number_density2", file = f2)
                #print(time, number_density, props['Ltot'], props['M/L'], props['Mave'], props['drc'], props['drhl'], props['dsigmac'], props['dvsigmac_rv'], props['rc'], props['rhl'], props['sigmac'],  props['vsigmac_rv'], number_density2, MS, G, BH+BHnonBH+2.*BBH, BHnonBH, BBH, NS+NSnonNS+2.*BNS, NSnonNS, BNS, N, Z, rv, rg, file = f2)
                #print(time, number_density, props['Ltot'], props['M/L'], props['Mave'], props['drc'], props['drhl'], props['dsigmac'], props['dvsigmac_rv'], props['rc'], props['rhl'], props['sigmac'],  props['vsigmac_rv'], number_density2, MS, G, BH+BHnonBH+2.*BBH, BHnonBH, BBH, NS+NSnonNS+2.*BNS, NSnonNS, BNS, N, Z, rv, rg, file = f2)
                print(time, number_density, props['Ltot'], props['M/L'], props['Mave'], props['drc'], props['drhl'], props['dsigmac'], props['dvsigmac_rv'], props['rc'], props['rhl'], props['sigmac'],  props['vsigmac_rv'], number_density2, N, Z, rv, rg, file = f2)
                print(time, number_density, props['Ltot'], props['M/L'], props['Mave'], props['drc'], props['drhl'], props['dsigmac'], props['dvsigmac_rv'], props['rc'], props['rhl'], props['sigmac'],  props['vsigmac_rv'], number_density2, N, Z, rv, rg, file = f2)

                f2.close()
                #f5.close()
                print('cluster_params done')    
        
                ###############
                print('made params file')
                ##get_sbp_from_2D_projection(path+string, snapno)
                #get_sbp_from_2D_projection_ncut(sourcepath+string, no_snap, BINNO=50, LCUT=12, NCUT=0.85)
                get_sbp_from_2D_projection_ncut(sourcepath+prefix, no_snap, BINNO=50, LCUT=12, NCUT=-1)
                print(key_snap, 'made SBP')
                ##velocity_dispersion_hdf5(sourcepath, key_snap, Starinbin=700, mcut = 0.85)
                ##velocity_dispersion_hdf5(sourcepath, key_snap, Starinbin=500, mcut = 0)
                #velocity_dispersion_2dsnap(sourcepath, no_snap, Starinbin=700, mcut = 0.85)
                #velocity_dispersion_2dsnap(sourcepath, no_snap, Starinbin=500, mcut = 0)
                velocity_dispersion_2dsnap(sourcepath+prefix, no_snap, Starinbin=700, mcut = 0)
                print('Made vel dispersion for',no_snap)

            except:
                print(key_snap, 'failed')
        #except:
        #    continue



##Make all relevent profiles for the new CMC-COSMIC models using hdf5
def make_SBP_VDP_NDP_hdf5_smooth(modelpath, dgc):   ##dgc in kpc
    snap_h5 = modelpath + 'initial.snapshots.h5'

    t_conv = conv('t', modelpath+'initial.conv.sh')

    #with pd.HDFStore(snap_h5) as snaps:
    #    snap_keys = snaps.keys()

    #snap_no = []; snap_time = []
    all_keys = np.genfromtxt(modelpath+'snap_keys.txt', dtype=str)

    string = 'initial'
    units=read_units(modelpath+string)

    for ii in range(len(all_keys[:,0])):
        #temp_key = read_keys(snap_keys[ii])
        #temp_no = int(temp_key[0])
        #temp_time = float(temp_key[1])

        #snap_no.append(int(temp_key[0]))
        #snap_time.append(float(temp_key[1]))

        #if temp_time*t_conv < 10000. or temp_time*t_conv >13000. or temp_no%2!=0:
        #    continue
         
        snap_key = '/'+all_keys[:,0][ii]+'(t='+all_keys[:,1][ii]+')'
        if float(all_keys[:,1][ii])*t_conv < 8000. or float(all_keys[:,1][ii])*t_conv > 13000. or int(all_keys[:,0][ii])%2!=0:
            continue

        else:
            temp_snap = cmct.Snapshot(fname=snap_h5, snapshot_name=snap_key, 
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

            fsbp = modelpath+'SBP'+all_keys[:,0][ii]+'_LCUT12_smooth.txt'
            fvdp = modelpath+'VDP'+all_keys[:,0][ii]+'_MCUT0d85_smooth.txt'
            fndp = modelpath+'NDP'+all_keys[:,0][ii]+'_MCUT0d85_smooth.txt'

            np.savetxt(fsbp, np.c_[v_bincenter, v_profile], fmt = '%f %f', header = '1.r(arcsec) 2.mu_v(mag/arcsec^2) t='+all_keys[:,1][ii], comments = '#')
            np.savetxt(fvdp, np.c_[star_velbin_arcsec, star_veldisp_profile, star_e_veldisp_profile], fmt = '%f %f %f', header = '1.r(arcsec) 2.sigma_v(km/s) 3.sigma_v_err(km/s) t='+all_keys[:,1][ii], comments = '#')
            np.savetxt(fndp, np.c_[star_numbin_arcsec, star_profile_arcsec, star_e_profile_arcsec], fmt = '%f %f %f', header = '1.r(arcsec) 2.num(1/arcsec^2) 3.num_err(1/arcsec^2) t='+all_keys[:,1][ii], comments = '#')

            print(snap_key)


    #snap_no_sort, snap_time_sort = (list(t) for t in zip(*sorted(zip(snap_no, snap_time))))

    #np.savetxt(modelpath+'snap_keys.txt', np.c_[snap_no_sort, snap_time_sort], fmt = '%d %.8f', header = '1.snap_no 2.snap_time', comments = '#')



##Make multiple different projections
def make_different_projection_snap(modelpath, snapno):
    all_keys = np.genfromtxt(modelpath+'snap_keys.txt', dtype=str)
    snaptime = all_keys[:,1][int(snapno)]
    snapkey = '/'+snapno+'(t='+snaptime+')'

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

        filename = modelpath+'SNAP'+snapno+'_'+str(ii)+'.2Dproj.dat'
        make_2D_projection(modelpath, snapkey, units, filename, SEEDY=randomseed, PROJ=projection)
        os.system('gzip '+filename)

        print(ii)











