from numpy import *
import gzip
import math
import constants
import scripts1
import scripts2
import scripts3
import history_cmc
import subprocess

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
    units = array(units, dtype=unittype)
    #return (massunitcgs, massunitmsun, mstarunitcgs, mstarunitmsun, lengthunitcgs, lengthunitparsec, timeunitcgs, timeunitsmyr, nbtimeunitcgs, nbtimeunitsmyr)
    return units

def find_correlations(string,snapno, binary, z, mcut):
    """finds anything that may be interesting to see correlations:
        string: file string
        snapno: snapno as a string
        binary: whether there were primordial binaries or not.  0: not, 1: yes"""

    snap = string+'.snap'+snapno+'.dat.gz'
    dyn = string+'.dyn.dat'
    from glob import glob
    out = string+'.out.dat'
    bin = string+'.bin.dat'

    units = read_units(string)
    convert_v = units[0]['l_cgs']/units[0]['nbt_cgs']/1e5
    convert_rho = units[0]['m_msun']/(units[0]['l_pc'])**3.
    t_myr = find_t_myr(string, snapno)
    m_guess = find_MS_turnoff(t_myr)
    m_to = find_MS_TO(t_myr, z, m_guess)

    f_snap = gzip.open(snap,'rb')
    line_snap = f_snap.readline()
    a_snap = line_snap.split()
    b_snap = a_snap[1]
    c_snap = b_snap.split('=')
    t_snap = float(c_snap[1])
    #print t_snap
    
    dyndata = loadtxt(dyn)
    for i in range(len(dyndata)):
        if (dyndata[i,0])==t_snap:
            m = dyndata[i,4]
            rc = dyndata[i,7]
            rh = dyndata[i,20]
            rho0 = dyndata[i,21]
            #lo behold: time unit in central v calculation is in nbody units
            v0_rms = dyndata[i,23]
            #print 'v=',v0_rms,convert_v    
            m_msun,rc_pc,rh_pc,rho0_msun_pc3,v0_rms_km_s = m*units[0]['m_msun'],rc*units[0]['l_pc'],rh*units[0]['l_pc'],rho0*convert_rho,v0_rms*convert_v
    
    try:
        bindata = loadtxt(bin)
        f_bc_array=[]
        f_b_array=[]
        rhb_pc_array=[]
        rhs_pc_array=[]
        #print 'opened bindata'
        for i in range(len(bindata)):
            if (bindata[i,0])-t_snap < 0.00001 and (bindata[i,0])-t_snap > -0.00001:
                #print 'found time'
                f_bc_array.append(bindata[i,10])
                f_b_array.append(bindata[i,11])
                rhb_pc_array.append(bindata[i,5]*units[0]['l_pc'])
                rhs_pc_array.append(bindata[i,4]*units[0]['l_pc'])
            f_bc=mean(f_bc_array)
            f_b=mean(f_b_array)
            rhb_pc=mean(rhb_pc_array)
            rhs_pc=mean(rhs_pc_array)
    except IOError:
        f_bc = 0.0
        f_b = 0.0
        rhb_pc = 0.0
        rhs_pc = rh_pc    
    
    BSS_count_core = find_BSS_core(string, snapno, rc, z, mcut)
    BSS_count = find_BSS(string,snapno, z, mcut)
    total_BSS = BSS_count[0]+BSS_count[1]
    #print 'look here', rc, rc_pc
    
    total_BSS_core = BSS_count_core[0]+BSS_count_core[1]
    if (BSS_count[0]+BSS_count[1])>0:
        history = history_cmc.history_maker(BSS_count[3], BSS_count[4], string, binary)
        #print history
        branching_ratio = history_cmc.branching_ratio(history, binary)
        ##print branching_ratio
        BSS_binint = total_BSS*branching_ratio['binint']
        BSS_pure_binint = total_BSS*branching_ratio['pure_binint']
        BSS_coll = total_BSS*branching_ratio['coll']
        BSS_pure_coll = total_BSS*branching_ratio['pure_coll']
        BSS_merger = total_BSS*branching_ratio['merger']
        BSS_pure_merger = total_BSS*branching_ratio['pure_merger']
        BSS_pure_MTB = total_BSS*branching_ratio['pure_mtb']
    else:
        history={}
        BSS_binint=0
        BSS_pure_binint=0
        BSS_coll=0
        BSS_pure_coll=0
        BSS_merger=0
        BSS_pure_merger=0
        BSS_pure_MTB=0

    #print "got BSS"
    
    startype_count = find_startype(string,snapno,4.,'HB' )
    startype_core_count = find_startype_core(string, snapno, rc, 4., 'HB')
    #print "got startype"
    total_startype = (startype_count[0]+startype_count[1])
    total_startype_core = (startype_core_count[0] + startype_core_count[1])

    try:
        outfile=open(out,'r')
    except IOError:
        #print 'could not find out.dat\n%s\n' %(string)
        temp1=string.split('/')
        temp2=''
        temp2+='/'
        for i in range(1,len(temp1)-1):
            temp2+=temp1[i]+'/'
        temp2+='*.o*'
        out=glob(temp2)[0]
        outfile = open(out,'r')
    line='initial'
    tcoll_array=[]
    while len(line)>0:
        line=outfile.readline()
        if line.rfind('Tcoll')>-1:
            tcoll=line.split()
            tcoll
            if float(tcoll[2])>(0.001*t_myr - 1.) and float(tcoll[2])<=(0.001*t_myr) and tcoll[6]!='nan':
                tcoll_array.append(float(tcoll[6]))
    mean_gamma=(mean(tcoll_array))**(-1.)

    snapdata = loadtxt(snap)
    ave_mc, ave_rc, count = 0., 0., 0
    try:
        for i in range(len(snapdata)):
            if snapdata[i,7]==0 and snapdata[i,2]<=rc:
                ave_mc += snapdata[i,1]
                ave_rc += snapdata[i,16]
                count += 1
            elif snapdata[i,2]>rc:
                raise StopIteration()
    except StopIteration:
        #print 'obtained core ave mass and stellar rad'
        pass
    ave_mc = (ave_mc/count)
    ave_rc = (ave_rc/count)

            

    dtype = [('m',float), ('rc', float), ('rh',float), ('rho',float), ('gamma',float), ('v0',float), ('t',float), ('mto', float), 
            ('fbc',float), ('fb',float), ('rhb',float), ('rhs',float), ('bss_sing',int), ('bss_bin',int), ('bssc_sing',int),
            ('bssc_bin',int), ('bss_binint',float), ('bss_pure_binint',float), ('bss_coll',float), ('bss_pure_coll',float), 
            ('bss_merger',float), ('bss_pure_merger',float), ('bss_pure_mtb', float), ('startype_sing',int), ('startype_bin',int), ('startypec_sing',int), 
            ('startypec_bin',int), ('ave_mc', float), ('ave_stellar_rad', float)]
    dummy=[]
    dummy.append( (m_msun,rc_pc,rh_pc,rho0_msun_pc3,mean_gamma,v0_rms_km_s,t_myr,m_to,f_bc,f_b,rhb_pc,rhs_pc,BSS_count[0],BSS_count[1],BSS_count_core[0], BSS_count_core[1], BSS_binint, BSS_pure_binint, BSS_coll, BSS_pure_coll, BSS_merger, BSS_pure_merger, BSS_pure_MTB, startype_count[0], startype_count[1], startype_core_count[0], startype_core_count[1], ave_mc, ave_rc) )
    dummy=array(dummy,dtype=dtype)
    

    #for i in range(len(dummy)):
    #    #print "%.3f " %(dummy[i],)
#
#    #print "\n"

    return dummy
    
    

def find_t_myr(filestring, snapno):
    """goes in the given snapshot and finds the physical time corresponding to that snap"""
    snapfile = filestring+'.snap'+snapno+'.dat.gz'

    f=gzip.open(snapfile,'r')
    line=f.readline()
    a=line.split()
    b=a[1]
    print(b)
    c=b.split('=')
    t=float(c[1])
    d=read_units(filestring)
    t_myr=t*d['t_myr'][0]
    f.close()
    return (t_myr)

def find_MS_turnoff(t):
    """given the time in Myr it finds the MS turn-off mass in Solar masses.  Very simple now.  Need to make the MS lifetime formula better. """
    t_yr = t*10**6
    lm = (9.921 - log10(t_yr))/3.6648
    m = 10**lm
    return(m)    

def find_t_ms(z, m):
    eta = log10(z/0.02)
    a1 = 1.593890e3+2.053038e3*eta+1.231226e3*eta**2.+2.327785e2*eta**3.
    a2 = 2.706708e3+ 1.483131e3*eta+ 5.772723e2*eta**2.+ 7.411230e1*eta**3.
    a3 = 1.466143e2 - 1.048442e2*eta - 6.795374e1*eta**2. - 1.391127e1*eta**3.
    a4 = 4.141960e-2 + 4.564888e-2*eta + 2.958542e-2*eta**2 + 5.571483e-3*eta**3.
    a5 = 3.426349e-1
    a6 = 1.949814e1 + 1.758178*eta - 6.008212*eta**2. - 4.470533*eta**3.
    a7 = 4.903830
    a8 = 5.212154e-2 + 3.166411e-2*eta - 2.750074e-3*eta**2. - 2.271549e-3*eta**3.
    a9 = 1.312179 - 3.294936e-1*eta + 9.231860e-2*eta**2. + 2.610989e-2*eta**3.
    a10 = 8.073972e-1


    m_hook = 1.0185 + 0.16015*eta + 0.0892*eta**2.
    m_HeF = 1.995 + 0.25*eta + 0.087*eta**2.
    m_FGB = 13.048*(z/0.02)**0.06/(1+0.0012*(0.02/z)**1.27)

    t_BGB = (a1+a2*m**4.+a3*m**5.5+m**7.)/(a4*m**2.+a5*m**7.)
    x = max([0.95,min([0.95-0.03*(eta+0.30103)]),0.99])
    mu = max(0.5, 1.0-0.01* max(a6/(m**a7) , a8+a9/m**a10))
    t_hook = mu*t_BGB

    t_MS = max(t_hook, x*t_BGB)

    return (t_MS)   ## t in Myr

def find_MS_TO(t, z, mguess):  ## t in myr
    tguess = find_t_ms(z, mguess)
    ##print mguess, tguess, (t-tguess)/t

    while (t-tguess)/t > 0.00005:
        mguess -= 0.00001
        tguess = find_t_ms(z, mguess)
        ##print mguess, tguess, (t-tguess)/t
    mto = mguess
    return mto
    


def find_startype(file_string,snapno,startype,startypename):
    """Input:
        file_string and the snapno to get the right snap file
        startype integer and the startype name, e.g., to find the giants it will be 3., 'giants'
        returns the startype number of stars in singlets and in binaries."""
    snapfile=file_string+'.snap'+snapno+'.dat.gz'
    convfile=file_string+'.conv.sh'
    wfile1=file_string+'.snap'+snapno+'.'+startypename+'.dat'

    t_yr = (find_t_myr(file_string, snapno))*10**6

    f=gzip.open(snapfile,'rb')
    f1 = open(wfile1, 'w')
    startype_sing_count = 0
    startype_bin_count = 0

    ##print header
    f1.write("#1:id #2:m[MSUN] #3:r #4:vr #5:vt #6:E #7:J #8:binflag #9:m0[MSUN] #10:m1[MSUN] #11:id0 #12:id1 #13:a[AU] #14:e #15:startype #16:luminosity[LSUN] #17:radius[RSUN]  #18:bin_startype0 #19:bin_startype1 #20:bin_star_lum0[LSUN] #21:bin_star_lum1[LSUN] #22:bin_star_radius0[RSUN] #23:bin_star_radius1[RSUN]\n")
    #f2.write("#1:id #2:m[MSUN] #3:r #4:vr #5:vt #6:E #7:J #8:binflag #9:m0[MSUN] #10:m1[MSUN] #11:id0 #12:id1 #13:a[AU] #14:e #15:startype #16:luminosity[LSUN] #17:radius[RSUN]  #18:bin_startype0 #19:bin_startype1 #20:bin_star_lum0[LSUN] #21:bin_star_lum1[LSUN] #22:bin_star_radius0[RSUN] #23:bin_star_radius1[RSUN]\n")
    for line in f:    
        a=line.split()
        if a[0]!='#' and a[0]!='#1:id':
            if a[7]=='0' and float(a[14])==startype:
                for i in range(len(a)):
                    f1.write("%s " % (a[i]))
                startype_sing_count += 1
                f1.write("\n")
                
            
            if a[7]!='0' and a[8]!='0' and a[9]!='0': 
                if float(a[17])==startype:
                    for i in range(len(a)):
                        f1.write("%s " % (a[i]))
                    startype_bin_count += 1
                    f1.write("\n")
                elif float(a[18])==startype:
                    for i in range(len(a)):
                        f1.write("%s " % (a[i]))
                    startype_bin_count += 1
                    f1.write("\n")
    f1.close()
    return (startype_sing_count, startype_bin_count,t_yr/10**6)


def find_startype_core(file_string,snapno, rc, startype, startypename):
    """takes the snapshot file, the conv.sh file, the file name to write the BSS in singles and the file name to write the BSS in binaries at that snapshot.  Returns the number of BSS in single within the core and the number of BSS in binaries within the core."""
    snapfile=file_string+'.snap'+snapno+'.dat.gz'
    convfile=file_string+'.conv.sh'
    
    t_yr = (find_t_myr(file_string, snapno))*10**6
    t_gyr = '%.1f' %(t_yr/10.**9.)

    #wfile1=file_string+'.snap'+snapno+'.'+str(t_gyr)+'.core.single_BSS.dat'
    #wfile2=file_string+'.snap'+snapno+'.'+str(t_gyr)+'.core.binary_BSS.dat'

    #m = find_MS_TO(float(t_yr))
    #m_cut = (1.+0.1)*m
    #print "looking at a snap at t = %f Gyr" % (t_yr/10**9)

    f=gzip.open(snapfile,'rb')
    f.seek(0)
    #f1 = open(wfile1,'w')
    #f2 = open(wfile2,'w')
    startype_sing_count = 0
    startype_bin_count = 0

    ##print header
    #f1.write("#1:id #2:m[MSUN] #3:r #4:vr #5:vt #6:E #7:J #8:binflag #9:m0[MSUN] #10:m1[MSUN] #11:id0 #12:id1 #13:a[AU] #14:e #15:startype #16:luminosity[LSUN] #17:radius[RSUN]  #18:bin_startype0 #19:bin_startype1 #20:bin_star_lum0[LSUN] #21:bin_star_lum1[LSUN] #22:bin_star_radius0[RSUN] #23:bin_star_radius1[RSUN]\n")
    #f2.write("#1:id #2:m[MSUN] #3:r #4:vr #5:vt #6:E #7:J #8:binflag #9:m0[MSUN] #10:m1[MSUN] #11:id0 #12:id1 #13:a[AU] #14:e #15:startype #16:luminosity[LSUN] #17:radius[RSUN]  #18:bin_startype0 #19:bin_startype1 #20:bin_star_lum0[LSUN] #21:bin_star_lum1[LSUN] #22:bin_star_radius0[RSUN] #23:bin_star_radius1[RSUN]\n")
    startype_ids = []
    startype_positions = []
    for line in f:    
        a=line.split()
        if a[0]!='#' and a[0]!='#1:id':
            if float(a[2]) <= float(rc):
                if a[7]=='0' and float(a[14])==3.:
                    #for i in range(len(a)):
                        #f1.write("%s " % (a[i]))
                    startype_sing_count += 1
                    #f1.write("\n")
                    
                
                if a[7]!='0' and a[8]!='0' and a[9]!='0': 
                    if float(a[17])==3.:
                        #for i in range(len(a)):
                            #f2.write("%s " % (a[i]))
                        startype_bin_count += 1
                        #f2.write("\n")
                    if float(a[18])==3.:
                        #for i in range(len(a)):
                            #f2.write("%s " % (a[i]))
                        startype_bin_count += 1
                        #f2.write("\n")
    f.close()
    #f1.close()
    #f2.close()
    return (startype_sing_count, startype_bin_count,t_yr/10**6)


    

def find_BSS(file_string,snapno, z, mcut):
    """takes the snapshot file, the conv.sh file, the file name to write the BSS in singles and the file name to write the BSS in binaries at that snapshot.  Returns the number of BSS in single and the number of BSS in binaries."""
    from glob import glob
    binintstring=file_string+'.binint.log'
    #binint=glob('*.binint.log')
    binint=glob(binintstring)
    binary=1
    if len(binint)==0:
        binary=0
    #print 'binary=',binary, len(binint), binint, file_string, binintstring
    snapfile=file_string+'.snap'+snapno+'.dat.gz'
    units = read_units(file_string)
    
    t_Myr = (find_t_myr(file_string, snapno))
    t_gyr = '%.1f' %(t_Myr/10.**3.)

    wfile1=file_string+'.snap'+snapno+'.'+str(t_gyr)+'.mcut'+str(mcut)+'.single_BSS.dat'
    wfile2=file_string+'.snap'+snapno+'.'+str(t_gyr)+'.mcut'+str(mcut)+'.binary_BSS.dat'

    mguess = find_MS_turnoff(float(t_Myr))
    m = find_MS_TO(float(t_Myr), z, mguess)
    m_cut = mcut*m
    #print "looking at a snap at t = %s Gyr, MS turnoff is %f MSun" % (t_gyr, m)

    f=gzip.open(snapfile,'rb')
    f.seek(0)
    f1 = open(wfile1,'w')
    f2 = open(wfile2,'w')
    bss_sing_count = 0
    bss_bin_count = 0

    #print 'header'
    #f1.write("#1:id #2:m[MSUN] #3:r #4:vr #5:vt #6:E #7:J #8:binflag #9:m0[MSUN] #10:m1[MSUN] #11:id0 #12:id1 #13:a[AU] #14:e #15:startype #16:luminosity[LSUN] #17:radius[RSUN]  #18:bin_startype0 #19:bin_startype1 #20:bin_star_lum0[LSUN] #21:bin_star_lum1[LSUN] #22:bin_star_radius0[RSUN] #23:bin_star_radius1[RSUN] 24.star_phi 25.binint? 26.pure_binint? 27.coll? 28.pure_coll? 29.se? 30.pure_se? 31.Leff(LSUN) 32.Teff(K)\n")
    f1.write("#1:id #2:m[MSUN] #3:r #4:vr #5:vt #6:E #7:J #8:binflag #9:m0[MSUN] #10:m1[MSUN] #11:id0 #12:id1 #13:a[AU] #14:e #15:startype #16:luminosity[LSUN] #17:radius[RSUN]  #18:bin_startype0 #19:bin_startype1 #20:bin_star_lum0[LSUN] #21:bin_star_lum1[LSUN] #22:bin_star_radius0[RSUN] #23:bin_star_radius1[RSUN] 24.star_phi 25.Leff(LSUN) 26.Teff(K)\n")
    #f2.write("#1:id #2:m[MSUN] #3:r #4:vr #5:vt #6:E #7:J #8:binflag #9:m0[MSUN] #10:m1[MSUN] #11:id0 #12:id1 #13:a[AU] #14:e #15:startype #16:luminosity[LSUN] #17:radius[RSUN]  #18:bin_startype0 #19:bin_startype1 #20:bin_star_lum0[LSUN] #21:bin_star_lum1[LSUN] #22:bin_star_radius0[RSUN] #23:bin_star_radius1[RSUN] 24.star_phi 25.binint? 26.pure_binint? 27.coll? 28.pure_coll? 29.se? 30.pure_se? 31.Leff(LSUN) 32.Teff(K)\n")
    f2.write("#1:id #2:m[MSUN] #3:r #4:vr #5:vt #6:E #7:J #8:binflag #9:m0[MSUN] #10:m1[MSUN] #11:id0 #12:id1 #13:a[AU] #14:e #15:startype #16:luminosity[LSUN] #17:radius[RSUN]  #18:bin_startype0 #19:bin_startype1 #20:bin_star_lum0[LSUN] #21:bin_star_lum1[LSUN] #22:bin_star_radius0[RSUN] #23:bin_star_radius1[RSUN] 24.star_phi 25.Leff(LSUN) 26.Teff(K)\n")

    BSS_ids = []
    BSS_positions = []
    for line in f:
        binint, pure_binint, coll, pure_coll, se, pure_se = 0, 0, 0, 0, 0, 0
        coll_time = []
        a=line.split()
        if a[0]!='#' and a[0]!='#1:id':
            if float(a[1])>m_cut and a[7]=='0' and float(a[14])<2.:
                for i in range(len(a)):
                    f1.write("%s " % (a[i]))
                bss_sing_count += 1
                dummy_id = int(a[0])
                BSS_ids.append(dummy_id)
                BSS_positions.append(float(a[2]))
                #f1.write("\n")
                #now get info about the BSS
                #hist=history_cmc.history_maker([dummy_id], [float(a[2])], file_string, binary)
                #if binary==1:
                #    #print 'binary', binary
                #    if len(hist[dummy_id]['binint']['binint'].keys())>0:
                #        binint = 1
                #        if len(hist[dummy_id]['coll'].keys())==0 and len(hist[dummy_id]['se'].keys())==0:
                #            pure_binint = 1
                #if len(hist[dummy_id]['coll'].keys())>0:
                #    coll = 1
                #    for i in hist[dummy_id]['coll'].keys():
                #        #print 'key', i
                #        coll_time.append(hist[dummy_id]['coll'][i]['coll_params']['time'])
                #        coll_time.sort()
                #        last_coll=coll_time[-1]
                #    if binary==1:
                #        #print 'binary', binary
                #        if len(hist[dummy_id]['binint']['binint'].keys())==0 and len(hist[dummy_id]['se'].keys())==0:
                #            pure_coll = 1
                #if len(hist[dummy_id]['se'].keys())>0:
                #    se = 1
                #    if binary==1 and len(hist[dummy_id]['binint'].keys())==0 and len(hist[dummy_id]['coll'].keys())==0:
                #        pure_se = 1

                #f1.write("%d %d %d %d %d %d " %(binint, pure_binint, coll, pure_coll, se, pure_se))
                Teff = find_T(float(a[15]), float(a[16]))
                Leff = float(a[15])
                f1.write("%f %f\n" %(Leff, Teff))


                
            
            if a[7]!='0' and a[8]!='0' and a[9]!='0': 
                if (float(a[8]))>m_cut and float(a[17])<2.:
                    for i in range(len(a)):
                        f2.write("%s " % (a[i]))
                    bss_bin_count += 1
                    dummy_id=int(a[10])
                    BSS_ids.append(dummy_id)
                    BSS_positions.append(float(a[2]))

                    #now get info about the BSS
                    #hist=history_cmc.history_maker([dummy_id], [float(a[2])], file_string, binary)
                    #if binary==1 and len(hist[dummy_id]['binint']['binint'].keys())>0:
                    #    binint = 1
                    #    if len(hist[dummy_id]['coll'].keys())==0 and len(hist[dummy_id]['se'].keys())==0:
                    #        pure_binint = 1
                    #if len(hist[dummy_id]['coll'].keys())>0:
                    #    coll = 1
                    #    for i in hist[dummy_id]['coll'].keys():
                    #        coll_time.append(hist[dummy_id]['coll'][i]['coll_params']['time'])
                    #        coll_time.sort()
                    #        last_coll=coll_time[-1]
                    #    if binary==1:
                    #        if binary==1 and len(hist[dummy_id]['binint'].keys())==0 and len(hist[dummy_id]['se'].keys())==0:
                    #            pure_coll = 1
                    #if len(hist[dummy_id]['se'].keys())>0:
                    #    se = 1
                    #    if binary==1 and len(hist[dummy_id]['binint'].keys())==0 and len(hist[dummy_id]['coll'].keys())==0:
                    #        pure_se = 1

                    #f2.write("%d %d %d %d %d %d " %(binint, pure_binint, coll, pure_coll, se, pure_se))
                    L0,R0,L1,R1=float(a[19]), float(a[21]), float(a[20]), float(a[22])
                    T0=find_T(L0, R0)
                    T1=find_T(L1, R1)
                    Leff=L0+L1
                    Teff=(L0*T0+L1*T1)/Leff
                    f2.write("%f %f\n" %(Leff, Teff))


                    #f2.write("\n")
                if (float(a[9]))>m_cut and float(a[18])<2.:
                    for i in range(len(a)):
                        f2.write("%s " % (a[i]))
                    bss_bin_count += 1
                    dummy_id = int(a[11])
                    BSS_ids.append(dummy_id)
                    BSS_positions.append(float(a[2]))
                    
                    #now get info about the BSS
                    #hist=history_cmc.history_maker([dummy_id], [float(a[2])], file_string, binary)
                    #if binary==1 and len(hist[dummy_id]['binint']['binint'].keys())>0:
                    #    binint = 1
                    #    if len(hist[dummy_id]['coll'].keys())==0 and len(hist[dummy_id]['se'].keys())==0:
                    #        pure_binint = 1
                    #if len(hist[dummy_id]['coll'].keys())>0:
                    #    coll = 1
                    #    for i in hist[dummy_id]['coll'].keys():
                    #        coll_time.append(hist[dummy_id]['coll'][i]['coll_params']['time'])
                    #        coll_time.sort()
                    #        last_coll=coll_time[-1]
                    #    if binary==1:
                    #        if binary==1 and len(hist[dummy_id]['binint'].keys())==0 and len(hist[dummy_id]['se'].keys())==0:
                    #            pure_coll = 1
                    #if len(hist[dummy_id]['se'].keys())>0:
                    #    se = 1
                    #    if binary==1 and len(hist[dummy_id]['binint'].keys())==0 and len(hist[dummy_id]['coll'].keys())==0:
                    #        pure_se = 1

                    #f2.write("%d %d %d %d %d %d " %(binint, pure_binint, coll, pure_coll, se, pure_se))

                    L0,R0,L1,R1=float(a[19]), float(a[21]), float(a[20]), float(a[22])
                    T0=find_T(L0, R0)
                    T1=find_T(L1, R1)
                    Leff=L0+L1
                    Teff=(L0*T0+L1*T1)/Leff
                    f2.write("%f %f\n" %(Leff, Teff))

                    

                    #f2.write("\n")
    f.close()
    f1.close()
    f2.close()
    return (bss_sing_count, bss_bin_count,t_Myr, BSS_ids, BSS_positions, m)



def find_BSS_core(file_string,snapno, rc, z, mcut):
    """takes the snapshot file, the conv.sh file, the file name to write the BSS in singles and the file name to write the BSS in binaries at that snapshot.  Returns the number of BSS in single within the core and the number of BSS in binaries within the core."""
    snapfile=file_string+'.snap'+snapno+'.dat.gz'
    convfile=file_string+'.conv.sh'
    
    t_Myr = find_t_myr(file_string, snapno)
    t_gyr = '%.1f' %(t_Myr/10.**3.)

    #wfile1=file_string+'.snap'+snapno+'.'+str(t_gyr)+'.core.single_BSS.dat'
    #wfile2=file_string+'.snap'+snapno+'.'+str(t_gyr)+'.core.binary_BSS.dat'

    mguess = find_MS_turnoff(float(t_Myr))
    m = find_MS_TO(t_Myr, z, mguess)
    m_cut = mcut*m
    #print "looking at a snap at t = %s Gyr, MS turnoff is %f MSun" % (t_gyr, m)

    f=gzip.open(snapfile,'rb')
    f.seek(0)
    #f1 = open(wfile1,'w')
    #f2 = open(wfile2,'w')
    bss_sing_count = 0
    bss_bin_count = 0

    ##print header
    #f1.write("#1:id #2:m[MSUN] #3:r #4:vr #5:vt #6:E #7:J #8:binflag #9:m0[MSUN] #10:m1[MSUN] #11:id0 #12:id1 #13:a[AU] #14:e #15:startype #16:luminosity[LSUN] #17:radius[RSUN]  #18:bin_startype0 #19:bin_startype1 #20:bin_star_lum0[LSUN] #21:bin_star_lum1[LSUN] #22:bin_star_radius0[RSUN] #23:bin_star_radius1[RSUN]\n")
    #f2.write("#1:id #2:m[MSUN] #3:r #4:vr #5:vt #6:E #7:J #8:binflag #9:m0[MSUN] #10:m1[MSUN] #11:id0 #12:id1 #13:a[AU] #14:e #15:startype #16:luminosity[LSUN] #17:radius[RSUN]  #18:bin_startype0 #19:bin_startype1 #20:bin_star_lum0[LSUN] #21:bin_star_lum1[LSUN] #22:bin_star_radius0[RSUN] #23:bin_star_radius1[RSUN]\n")
    BSS_ids = []
    BSS_positions = []
    for line in f:    
        a=line.split()
        if a[0]!='#' and a[0]!='#1:id':
            if float(a[2]) <= float(rc):
                if float(a[1])>m_cut and a[7]=='0' and float(a[14])<2.:
                    #for i in range(len(a)):
                        #f1.write("%s " % (a[i]))
                    bss_sing_count += 1
                    BSS_ids.append(int(a[0]))
                    BSS_positions.append(float(a[2]))
                    #f1.write("\n")
                    
                
                if a[7]!='0' and a[8]!='0' and a[9]!='0': 
                    if (float(a[8]))>m_cut and float(a[17])<2.:
                        #for i in range(len(a)):
                            #f2.write("%s " % (a[i]))
                        bss_bin_count += 1
                        BSS_ids.append(int(a[10]))
                        BSS_positions.append(float(a[2]))
                        #f2.write("\n")
                    if (float(a[9]))>m_cut and float(a[18])<2.:
                        #for i in range(len(a)):
                            #f2.write("%s " % (a[i]))
                        bss_bin_count += 1
                        BSS_ids.append(int(a[11]))
                        BSS_positions.append(float(a[2]))
                        #f2.write("\n")
    f.close()
    #f1.close()
    #f2.close()
    return (bss_sing_count, bss_bin_count,t_Myr, BSS_ids, BSS_positions)



def BS_r_cum(singlefile, binaryfile, giantfile, rc, writefile, nbin):
    sing_data=loadtxt(singlefile)
    bin_data=loadtxt(binaryfile)
    giant_data=loadtxt(giantfile)
    f=open(writefile,'w')
    
    #initialize
    sing_roverrc=zeros(len(sing_data))
    bin_roverrc=zeros(len(bin_data))
    both_roverrc=zeros(len(bin_data)+len(sing_data))
    g_roverrc=zeros(len(giant_data))
    
    #getting the scaled radial distribution length scaled to the core radius 
    for i in range(len(sing_data)):
        sing_roverrc[i]=sing_data[i,2]/rc
    for i in range(len(bin_data)):
        bin_roverrc[i]=bin_data[i,2]/rc
    for i in range(len(giant_data)):
        g_roverrc[i]=giant_data[i,2]/rc
    for i in range(len(sing_data)):
        both_roverrc[i]=sing_roverrc[i]
    for i in range(len(bin_data)):
        both_roverrc[i+len(sing_data)]=bin_roverrc[i]


    nbinmem = len(g_roverrc)/nbin
    sing_cum, bin_cum, all_cum = 0., 0., 0.
    for i in range(1,nbin):
        temp_sing, temp_bin, temp_all = 0., 0., 0.
        rmin, rmax = g_roverrc[(i-1)*nbinmem], g_roverrc[i*nbinmem]
        #print rmin, rmax, (rmin+rmax)/2.
        for j in range(len(sing_roverrc)):
            if sing_roverrc[j] < rmax and sing_roverrc[j]>=rmin:
                temp_sing += 1
        for j in range(len(bin_roverrc)):
            if bin_roverrc[j] < rmax and bin_roverrc[j]>=rmin:
                temp_bin += 1
        for j in range(len(both_roverrc)):
            if both_roverrc[j] < rmax and both_roverrc[j]>=rmin:
                temp_all += 1
        r = (rmin+rmax)/2.
        temp_sing, dtemp_sing = temp_sing/float(nbinmem), temp_sing**0.5/float(nbinmem)
        temp_bin, dtemp_bin = temp_bin/float(nbinmem), temp_bin**0.5/float(nbinmem)
        temp_all, dtemp_all  = temp_all/float(nbinmem), temp_all**0.5/float(nbinmem)
        sing_cum += temp_sing
        bin_cum += temp_bin
        all_cum += temp_all

        f.write("%f %f %f %f %f %f %f %f %f %f\n" %(r, temp_sing, dtemp_sing, temp_bin, dtemp_bin, temp_all, dtemp_all, sing_cum, bin_cum, all_cum) )
        ##print r, temp_sing, temp_bin, temp_all, sing_cum, bin_cum, all_cum
    f.close()

    
                
def BS_r(singlefile, binaryfile, giantfile, rc, writefile):
    """gives the normalized radial distribution of the BSS population
    to get this there is some preprocessing needs to be done
    1. run hrdiag on your favorite snapshot
    2. from hrdiag file filter out the singles with the BSS cut as you like in one file
    3. from hrdiag file filter out the binaries with the BSS cut as you like in another file
    4. from hrdiag file filter out the Giants 
    5. feed these 3 files and the value of rc at that time from *.dyn.dat"""
    sing_data=loadtxt(singlefile)
    bin_data=loadtxt(binaryfile)
    giant_data=loadtxt(giantfile)
    f=open(writefile,'w')
    
    #initialize
    sing_roverrc=zeros(len(sing_data))
    bin_roverrc=zeros(len(bin_data))
    both_roverrc=zeros(len(bin_data)+len(sing_data))
    g_roverrc=zeros(len(giant_data))
    
    #getting the scaled radial distribution length scaled to the core radius 
    for i in range(len(sing_data)):
        sing_roverrc[i]=sing_data[i,8]/rc
    for i in range(len(bin_data)):
        bin_roverrc[i]=bin_data[i,8]/rc
    for i in range(len(giant_data)):
        g_roverrc[i]=giant_data[i,8]/rc
    for i in range(len(sing_data)):
        both_roverrc[i]=sing_roverrc[i]
    for i in range(len(bin_data)):
        both_roverrc[i+len(sing_data)]=bin_roverrc[i]

    ##print sing_roverrc, bin_roverrc, both_roverrc
    ##print len(sing_roverrc), len(bin_roverrc), len(both_roverrc)

    
    
    #get histogram for the scaled radii
    amin, amax = min(g_roverrc[:]), max(g_roverrc[:])
    hd_both=histogram(both_roverrc,bins=50, range=(amin, amax))
    hd_sing=histogram(sing_roverrc,bins=50, range=(amin, amax))
    hd_bin=histogram(bin_roverrc,bins=50, range=(amin, amax))
    hd_g=histogram(g_roverrc,bins=50, range=(amin, amax))

    #print hd_both, hd_sing, hd_bin
    #print len(hd_both[0]), len(hd_sing[0]), len(hd_bin[0])

    #normalize the numbers with the giant population
    normalized_both=zeros((len(hd_both[0]),3))
    normalized_sing=zeros((len(hd_sing[0]),3))
    normalized_bin=zeros((len(hd_sing[0]),3))
    ##print normalized_both, normalized_sing, normalized_bin
    ##print len(normalized_both), len(normalized_sing), len(normalized_bin)
    for i in range(1,len(hd_both[1])):
        count=0
        for j in range(len(g_roverrc)):
            if g_roverrc[j]<hd_both[1][i] and g_roverrc[j]>hd_both[1][i-1]:
                count += 1
        #print count
        normalized_both[i-1]=[ float(hd_both[0][i-1])/count , (hd_both[0][i-1]**(-0.5) + count**-0.5)*float(hd_both[0][i-1])/count , hd_both[1][i-1] ]
        f.write ("%f %f %f\n" % (normalized_both[i-1,2],normalized_both[i-1,0],normalized_both[i-1,1]))
    
    for i in range(1,len(hd_sing[1])):
        count=0
        for j in range(len(g_roverrc)):
            if g_roverrc[j]<hd_sing[1][i] and g_roverrc[j]>hd_sing[1][i-1]:
                count += 1
        #print count
        normalized_sing[i-1]=[ float(hd_sing[0][i-1])/count , (hd_sing[0][i-1]**(-0.5) + count**-0.5)*float(hd_sing[0][i-1])/count , hd_sing[1][i-1] ]
    
    for i in range(1,len(hd_bin[1])):
        count=0
        for j in range(len(g_roverrc)):
            if g_roverrc[j]<hd_bin[1][i] and g_roverrc[j]>hd_bin[1][i-1]:
                count += 1
        #print count
        normalized_bin[i-1]=[ float(hd_bin[0][i-1])/count , (hd_bin[0][i-1]**(-0.5) + count**-0.5)*float(hd_bin[0][i-1])/count, hd_bin[1][i-1] ]
    
    ##print normalized_both, normalized_sing, normalized_bin
    ##print len(normalized_both), len(normalized_sing), len(normalized_bin)
    

    #now plot
    import gracePlot
    gpl=gracePlot.gracePlot()
    gpl.hold()
    gpl.plot(normalized_both[:,2],normalized_both[:,0],normalized_both[:,1])

    f.close()

    #gpl.plot(normalized_sing[:,2],normalized_sing[:,0],normalized_sing[:,1])
    #gpl.plot(normalized_bin[:,2],normalized_bin[:,0],normalized_bin[:,1])    


#
#    return (normalized_both, normalized_sing, normalized_bin)
        

def find_problematic(id,filename):
    """takes the id and then returns if that id was somewhere in the given file"""    
    f=open(filename,'r')
    f.seek(0)
    try:
        while f:
            line=f.readline()
            if line.rfind(id)>-1:
                ##print "%s\n" % (line)
                problem=line
            elif line=='':
                raise StopIteration()
    except StopIteration:
        f.seek(0)
    f.close()
    return (problem)

def c_w0(c):
    """taking the concentration parameter gives the w0 value within range c (0.5,3.5)"""
    w0 = -2.47882 + 9.53235*c - 0.843304*c**2 - 2.34742*c**3 + 1.22119*c**4 - 0.167361*c**5
    return w0 

def w0_c(w0):
    """taking the w0 gives the concentration parameter within range w0 (2,15)"""
    #c=0.247266 + 0.2079*w0 - 0.057241*w0**2 + 0.0145234*w0**3 - 0.00120862*w0**4 -0.0000330051*w0**5
    c1 = 1.34813 - 1.56801*w0 + 1.00755*w0**2. - 0.29303*w0**3. + 0.0438717*w0**4. - 0.00274824*w0**5. - 6.67751e-05*w0**6. + 1.98723e-05*w0**7. - 1.06025e-06*w0**8 + 1.91532e-08*w0**9.
    c2 = -4.5853 + 9.0324 * w0 - 7.0123 * w0**2 + 3.0991 * w0**3 - 0.84764 * w0**4 + 0.15002 * w0**5 - 0.017426 * w0**6 + 0.0013168 * w0**7 - 6.227e-05 * w0**8 + 1.6728e-06 * w0**9 - 1.948e-08 * w0**10
    return c1, c2

def metallicity(m,s):
    """give the metalicity and tell how to convert
    options are:
        ztofe/h
        fe/htoz"""

    if s=='ztofe/h':
        m1 = log10(m/0.02)
    elif s=='fe/htoz':
        m1 = 0.02*10**m
    return m1

def find_rc(file_string,snapno):
    """finds r_c at this snapfile"""
    snapfile=file_string+'.snap'+snapno+'.dat.gz'
    dynfile=file_string+'.dyn.dat'

    snap=gzip.open(snapfile,'rb')
    line=snap.readline()
    data = line.split()
    data1=data[1].split('=')
    time=data1[1]
    
    dyndata=loadtxt(dynfile)
    for i in range(len(dyndata)):
        if dyndata[i,0]==float(time):
            rc = dyndata[i,7]

    return rc
    


def filter_single(filter,T):
    data=loadtxt(filter)
    sum=0.0
    normalize=0.0
    const=2*constants.h*constants.c**2
    bin_min = data[0,0] - (data[1,0] - data[0,0])
    for i in range (len(data)):
        sum+= const / (data[i,0]*constants.Angstrom)**5 / (exp(constants.c * constants.h/(data[i,0]*constants.Angstrom)/constants.k / T) -1) * (data[i,0]-bin_min)*constants.Angstrom * data[i,1] 
        normalize+= (data[i,0]-bin_min)*constants.Angstrom * data[i,1]
        bin_min = data[i,0]
    normalized = (sum/normalize)/4/pi/85000./constants.PC
    mag = 1.2 -2.5*log10(normalized)
    return (const,normalized,normalize,sum,mag)

def find_T(L,R):
    """L in LSUN and R in RSUN: returns T in K"""
    try:
        if R>10.**-100:
            T= (L*constants.Lsun / (4*pi*(R*constants.Rsun)**2) / constants.sigma)**0.25
        else:
            T=0.001
    except OverflowError:
        #print 'problem in radius'
        T= 0.001

    return (T)

def find_g(M,R):
    """M in MSun and R in RSun: returns g in CGS"""
    g = constants.G*M*constants.Msun/(R*constants.Rsun)**2.
    return(g)


def filter_singles(filterdata,L,R):
    #filterdata = loadtxt(filter)
    sum = 0.0
    sun_sum = 0.0 
    normalize = 0.0 
    const = 2*constants.h*constants.c**2
    T = (L*constants.Lsun / (4*pi*(R*constants.Rsun)**2) / constants.sigma)**0.25
    
    bin_min = filterdata[0,0] - (filterdata[1,0] - filterdata[0,0])
    
    for i in range(len(filterdata)):
        a1 = const/(filterdata[i,0]*constants.Angstrom)**5
        a2 = ( exp(constants.c * constants.h / filterdata[i,0]/constants.Angstrom/constants.k/T) -1 )**-1.
        sun_a2 = ( exp(constants.c * constants.h / filterdata[i,0]/constants.Angstrom/constants.k/constants.Tsun) -1 )**-1.
        a3 = (filterdata[i,0] - bin_min)*constants.Angstrom

        sum += a1*a2*a3*filterdata[i,1]
        sun_sum += a1*sun_a2*a3*filterdata[i,1]
        normalize += a3*filterdata[i,1]    
        
        #sum += const / (filterdata[i,0]*constants.Angstrom)**5. / (exp(constants.c * constants.h / constants.k / filterdata[i,0]*constants.Angstrom / T) - 1) * (filterdata[i,0] - bin_min)*constants.Angstrom * filterdata[i,1]
    
        #normalize += (filterdata[i,0] - bin_min)*constants.Angstrom * filterdata[i,1]        
        
        bin_min = filterdata[i,0]
        ##print (sum, normalize, bin_min)
    ##print (sum,sun_sum, T)
    #normalized = (sum/normalize) * (4*pi*R**2) * pi 
    #sun_normalized = (sun_sum/normalize) * (4*pi*R**2) * pi
    #fraction = normalized/L/constants.Lsun
    #mag = -2.5*log10(normalized/sun_normalized)
    lum = sum * 4*pi*(R*constants.Rsun)**2 * pi
    absolute_mag = -2.5*log10(lum) + 92.0218
    #mag = absolute_mag + 5*log10(dist_mod) - 5
    return (absolute_mag)


def filter_binaries(filterdata,L1,R1,L2,R2):
    sum = 0.0
    const = 2*constants.h*constants.c**2
    T1 = (L1*constants.Lsun / (4*pi*R1**2*constants.Rsun**2) / constants.sigma)**0.25
    T2 = (L2*constants.Lsun / (4*pi*R2**2*constants.Rsun**2) / constants.sigma)**0.25
    bin_min = filterdata[0,0] - (filterdata[1,0] - filterdata[0,0])
    for i in range(len(filterdata)):
        a1 = const / (filterdata[i,0]*constants.Angstrom)**5
        a2_1 = 4*pi*(R1*constants.Rsun)**2 * pi * ( exp(constants.c * constants.h / filterdata[i,0]/constants.Angstrom/constants.k/T1) -1 )**-1.
        a2_2 = 4*pi*(R2*constants.Rsun)**2 * pi * ( exp(constants.c * constants.h / filterdata[i,0]/constants.Angstrom/constants.k/T2) -1 )**-1.
        a3 = (filterdata[i,0] - bin_min)*constants.Angstrom

        sum += a1*a2_1*a3*filterdata[i,1]
        sum += a1*a2_2*a3*filterdata[i,1]
        
        bin_min = filterdata[i,0]
    
    lum = sum 
    absolute_mag = -2.5*log10(lum) + 92.0218
    #mag = absolute_mag + 5*log10(dist_mod) - 5
    return (absolute_mag)



def hrdiag(snap,filter1,filter2,writefile):
    #print 'started'
    snap=loadtxt(snap)
    #print 'snap done'
    f1=loadtxt(filter1)
    #print 'filter1 done'
    f2=loadtxt(filter2)
    #print 'filter2 done'
    na=-100
    fwrite=open(writefile,'w')
    #print 'open file done'
    fwrite.write("#1.binflag 2.startype0 3.startype1 4.id 5.id0 6.id1 7.m0[MSUN] 8.m1[MSUN] 9.r 10.rad0[RSUN] 11.rad1[RSUN] 12.L0[erg/s] 13.L1[erg/s] 14.B-I 15.B\n")
    #print 'header done'

    for i in range(len(snap)):
        if snap[i,7]==0.:
            magb=filter_singles(f1,snap[i,15],snap[i,16])
            magi=filter_singles(f2,snap[i,15],snap[i,16])
            fwrite.write("%d %d %d %ld %ld %ld %g %g %g %g %g %g %g %g %g\n" % ( snap[i,7],snap[i,14],na,snap[i,0],na,na,snap[i,1],na,snap[i,2],snap[i,16],na,snap[i,15],na,magb-magi,magb ))
            #print (i,snap[i,7])
        elif snap[i,7]==1.:
            magb=filter_binaries(f1,snap[i,19],snap[i,21],snap[i,20],snap[i,22])
            magi=filter_binaries(f2,snap[i,19],snap[i,21],snap[i,20],snap[i,22])
            fwrite.write("%d %d %d %ld %ld %ld %g %g %g %g %g %g %g %g %g\n" % ( snap[i,7],snap[i,17],snap[i,18],snap[i,0],snap[i,10],snap[i,11],snap[i,8],snap[i,9],snap[i,2],snap[i,21],snap[i,22],snap[i,19],snap[i,20],magb-magi,magb ))
            #print (i,snap[i,7])

    fwrite.close()


def hrdiag_L_T(filestring, snapno):
    """Give it a snap file and it gives you L and T_eff for all stars
    for binaries the temperature is a luminosity averaged temperature"""
    #print 'started' 
    snapfile=filestring+'.snap'+snapno+'.dat.gz'
    snap=genfromtxt(snapfile)
    #print 'snap done'
    writefile=filestring+'.snap'+snapno+'.hrdiag.dat'
    fwrite=open(writefile,'w')
    na=-100
    fwrite.write("#1.binflag 2.startype0 3.startype1 4.id 5.id0 6.id1 7.m0[MSUN] 8.m1[MSUN] 9.r 10.rad0[RSUN] 11.rad1[RSUN] 12.L0[LSUN] 13.L1[LSUN] 14.T0(K) 15.T1(K) 16.log10(T_eff/K) 17.log10(Leff/LSUN)\n")  #T_eff and T0 are the same for singles. For binaries Teff is L averaged T
    #print 'header done'

    for i in range(len(snap)):
        if snap[i,7]!=1:
            T0=find_T(snap[i,15],snap[i,16])
            fwrite.write("%d %d %d %ld %ld %ld %g %g %g %g %g %g %g %g %d %g %g\n" % ( snap[i,7],snap[i,14],na,snap[i,0],na,na,snap[i,1],na,snap[i,2],snap[i,16],na,snap[i,15],na,T0,na,log10(T0),log10(snap[i,15])))
            ##print (i,snap[i,7])
        if snap[i,7]==1:
            T0=find_T(snap[i,19],snap[i,21])
            T1=find_T(snap[i,20],snap[i,22])
            Teff=(snap[i,19]*T0 + snap[i,20]*T1)/(snap[i,19]+snap[i,20])
            fwrite.write("%d %d %d %ld %ld %ld %g %g %g %g %g %g %g %g %d %g %g\n" % ( snap[i,7],snap[i,17],snap[i,18],snap[i,0],snap[i,10],snap[i,11],snap[i,8],snap[i,9],snap[i,2],snap[i,21],snap[i,22],snap[i,19],snap[i,20],T0,T1,log10(Teff),log10(snap[i,19]+snap[i,20])))
            ##print (i,snap[i,7])
    fwrite.close()


def prune_hrdiag_L_T(filestring, snapno, smoothlength):
    hrdfile=filestring+'.snap'+snapno+'.hrdiag.dat'
    writefile=filestring+'.snap'+snapno+'.hrdiag_pruned.dat'
    f=open(writefile, 'w')

    data=loadtxt(hrdfile)
    for i in range(len(data)):
        if data[i,0]==0:  #single
            if data[i,1]<2. and i%smoothlength==0:  #prune some single MSSs
                for j in range(len(data[i,:])):
                    f.write("%f " %(data[i,j],))
                f.write("\n")
            elif data[i,1]>9 and data[i,1]<13 and i%100==0:
                for j in range(len(data[i,:])):
                    f.write("%f " %(data[i,j],))
                f.write("\n")
            elif (data[i,1]>2 and data[i,1]<9):  #all other stars don't prune
                for j in range(len(data[i,:])):
                    f.write("%f " %(data[i,j],))
                f.write("\n")
        elif data[i,0]==1:   #binaries
            if data[i,1]<2. and data[i,2]<2. and i%100==0:   #prune some MS binaries
                for j in range(len(data[i,:])):
                    f.write("%f " %(data[i,j],))
                f.write("\n")
            elif data[i,1]>2 or data[i,2]>2:                        ##print everything else
                for j in range(len(data[i,:])):
                    f.write("%f " %(data[i,j],))
                f.write("\n")
    f.close()    


def color_comp(snap,writefile,filter1,filter2,filter3,filter4):
    #print 'started'
    snap=loadtxt(snap)
    #print 'snap done'
    f1=loadtxt(filter1)
    #print 'filter1 done'
    f2=loadtxt(filter2)
    #print 'filter2 done'
    f3=loadtxt(filter3)
    #print 'filter3 done'
    f4=loadtxt(filter4)
    #print 'filter4 done'
    na=-100
    fwrite=open(writefile,'w')
    #print 'open file done'
    fwrite.write("#1.binflag 2.startype0 3.startype1 4.id 5.id0 6.id1 7.m0[MSUN] 8.m1[MSUN] 9.r 10.rad0[RSUN] 11.rad1[RSUN] 12.L0[erg/s] 13.L1[erg/s] 14.g1(cm/s^2) 15.g2(cm/s^2) 16.T1 (K) 17.T2(K) 18.f435w 19.f555w 20.f606w 21.f814w\n")
    #print 'header done'

    for i in range(len(snap)):
        if snap[i,7]==0.:
            g1=constants.G*snap[i,1]*constants.Msun/(snap[i,16]*constants.Rsun)**2.
            T1=find_T(snap[i,15],snap[i,16])
            magb=filter_singles(f1,snap[i,15],snap[i,16])
            magv=filter_singles(f2,snap[i,15],snap[i,16])
            magr=filter_singles(f3,snap[i,15],snap[i,16])
            magi=filter_singles(f4,snap[i,15],snap[i,16])
            fwrite.write("%d %d %d %ld %ld %ld %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n" % ( snap[i,7],snap[i,14],na,snap[i,0],na,na,snap[i,1],na,snap[i,2],snap[i,16],na,snap[i,15],na,g1,na,T1,na,magb,magv,magr,magi ))
            #print (i,snap[i,7],len(snap))
        elif snap[i,7]==1.:
            g1=constants.G*snap[i,8]*constants.Msun/(snap[i,21]*constants.Rsun)**2.
            g2=constants.G*snap[i,9]*constants.Msun/(snap[i,22]*constants.Rsun)**2.
            T1=find_T(snap[i,19],snap[i,21])
            T2=find_T(snap[i,20],snap[i,22])
            magb=filter_binaries(f1,snap[i,19],snap[i,21],snap[i,20],snap[i,22])
            magv=filter_binaries(f2,snap[i,19],snap[i,21],snap[i,20],snap[i,22])
            magr=filter_binaries(f3,snap[i,19],snap[i,21],snap[i,20],snap[i,22])
            magi=filter_binaries(f4,snap[i,19],snap[i,21],snap[i,20],snap[i,22])
            fwrite.write("%d %d %d %ld %ld %ld %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n" % ( snap[i,7],snap[i,17],snap[i,18],snap[i,0],snap[i,10],snap[i,11],snap[i,8],snap[i,9],snap[i,2],snap[i,21],snap[i,22],snap[i,19],snap[i,20],g1,g2,T1,T2,magb,magv,magr,magi ))
            #print (i,snap[i,7],len(snap))

    fwrite.close()

def find_rtidal(mc,vg=220.,rg=8000.):
    """takes 
    the cluster mass (mc) in solar mass, 
    the galactic circular velocity (vg) in km/s, 
    and 
    the galactocentric distance (rg) in pc
    returns: the tidal radius of the cluster in pc"""
    r_tidal = (constants.G * mc*constants.Msun / 2 / (vg*constants.km)**2.)**(1./3.) * (rg*constants.PC)**(2./3.) / constants.PC
    return r_tidal


def find_t_vs_rgb(filestring, z):
    from glob import glob
    files = glob('*.dat.gz')
    dynfile = filestring+'.dyn.dat'
    dyndata = loadtxt(dynfile)
    convfile = filestring+'.conv.sh'
    units = read_units(convfile)
    wfile = filestring+'.t_rgb_bss.dat'
    wfile = open(wfile, 'w')
    wfile.write("#1.t(Myr) 2.SBSS 3.BBSS 4.tot_BSS 5.SRGB 6.BRGB 7.tot_RGB 8.SBSS_c 9.BBSS_c 10.tot_BSS_c\n")
    for i in range(len(files)):
        filestring = files[i].split('.snap')[0]
        snapno = files[i].split('.snap')[1].split('.')[0]
        snapfile = gzip.open(str(files[i]),'rb')
        t_snap = float(snapfile.readline().split('=')[1].split('[')[0])
        t_Myr = t_snap*units[7]
        m_guess = find_MS_turnoff(t_Myr)
        m_to = find_MS_TO(t_Myr, z, mguess)
        m_cut = 1.1*m_to
        snapdata = loadtxt(snapfile)

        for j in range(len(dyndata)):
            if t_snap == dyndata[j,0]:
                rc_snap = dyndata[j,7]
                rh_snap = dyndata[j,20]

        #get bss
        bss_sing_count = 0
        bss_bin_count = 0
        bss_sing_count_core = 0
        bss_bin_count_core = 0
        giants_sing_count = 0
        giants_bin_count = 0
        for j in range(len(snapdata)):
            if snapdata[j,0]!='#' and snapdata[j,0]!='#1:id':
                if snapdata[j,1]>m_cut and snapdata[j,7]==0 and snapdata[j,14]<2.:
                    bss_sing_count += 1
                    if snapdata[j,2]<=rc_snap:
                        bss_sing_count_core += 1
                if snapdata[j,7]==0 and snapdata[j,14]==3:
                    giants_sing_count += 1
            
                if snapdata[j,7]!=0. and snapdata[j,8]!=0. and snapdata[j,9]!=0.: 
                    if snapdata[j,8]>m_cut and snapdata[j,17]<2.:
                        bss_bin_count += 1
                        if snapdata[j,2]<=rc_snap:
                            bss_bin_count_core += 1
                    if snapdata[j,9]>m_cut and snapdata[j,18]<2.:
                        bss_bin_count += 1
                        if snapdata[j,2]<=rc_snap:
                            bss_bin_count_core += 1
                    if snapdata[j,17]==3 or snapdata[j,18]==3:
                        giants_bin_count += 1


        #print t_snap*units[7]/10.**3., 'Gyr'
        tot_bss_snap = bss_sing_count + bss_bin_count 
        tot_bss_core_snap = bss_sing_count_core + bss_bin_count_core 
        tot_giants_snap = giants_sing_count + giants_bin_count 

        dummy = ( t_snap*units[7], bss_sing_count, bss_bin_count, tot_bss_snap, giants_sing_count, giants_bin_count, tot_giants_snap, bss_sing_count_core, bss_bin_count_core, tot_bss_core_snap )
        for j in range(len(dummy)):
            wfile.write("%f " %(dummy[j]))
        wfile.write("\n")

        snapfile.close()

    wfile.close()



def two_point_correlation(a1, a2):
    """given two arrays of numbers of equal dimension gives the 2-point correlation coefficient.  """

    a1_bar = mean(a1)
    a2_bar = mean(a2)
    a1_sd = std(a1)
    a2_sd = std(a2)

    temp_sum=0
    try:
        if len(a1)!=len(a2):
            raise StopIteration()
        else:
            for i in range(len(a1)):
                temp_sum += (a1[i]-a1_bar)*(a2[i]-a2_bar)/(a1_sd*a2_sd)
    except StopIteration:
        #print "bad arrays: arrays have different dimensions"
        pass
        
    cor = float(temp_sum)/float(len(a1))

    return cor

def weighted_pairwise_correlation(a1, a2, weight1, weight2):
    a1_bar = average(a1, weights=weight1)
    a2_bar = average(a2, weights=weight2)
    #calculate sd's
    a1_sd = weighted_sd(a1, weight1)
    a2_sd = weighted_sd(a2, weight2)
    #calculate covariance
    a1a2_cov = weighted_covariance(a1, a2, weight1, weight2)
    #print a1_bar, a2_bar, a1_sd, a2_sd, a1a2_cov
    wpc = a1a2_cov/(a1_sd*a2_sd)

    return wpc
    

def weighted_sd(a, weight):
    weight = array(weight)
    
    try:
        if len(a)!=len(weight):
            raise StopIteration()
        else:
            #normalize weights
            sqsum = 0.
            w = weight*len(weight)/sum(weight)
            a_bar = average(a, weights=w)
            for i in range(len(a)):
                sqsum += w[i]*( a[i]-a_bar )**2.
            a_sd = (sqsum/len(a))**0.5
            #print a_sd
    except StopIteration:
        #print "bad input: arrays have different dimensions"
        pass
    
    return a_sd

def weighted_covariance(a1, a2, weight1, weight2):
    weight1 = array(weight1)
    weight2 = array(weight2)
    try:
        if len(a1)!=len(weight1) or len(a2)!=len(weight2) or len(a1)!=len(a2):
            raise StopIteration()
        else:
            a1_bar = average(a1, weights=weight1)
            a2_bar = average(a2, weights=weight2)
            weight = weight1+weight2
            #normalize weights
            w = weight*len(weight)/sum(weight)
            #print w
            temp_sum = 0.
            for i in range(len(a1)):
                #print temp_sum, a1[i], a1_bar, a2[i], a2_bar
                temp_sum += w[i] * ( a1[i] - a1_bar ) * ( a2[i] - a2_bar )
            temp_cov = temp_sum/len(a1)
    except StopIteration:
        #print "bad input: arrays have different dimensions"
        pass
    return temp_cov




def tidal_mass_loss(filestring):
    """creates file with tidal mass loss and stellar evolution mass loss
    takes the filestring and creates a file"""

    dynfile=filestring+'.dyn.dat'
    outfile=filestring+'.mass_loss.dat'
    convfile=filestring+'.conv.sh'

    data=loadtxt(dynfile)
    units=read_units(convfile)
    f=open(outfile, 'w')

    sum_Mse=0.
    f.write("#1.t(Myr) 2.M/M(0) 3.M_SE/M(0) 4.dM_SE/M(0)\n")
    for i in range(len(data)):
        sum_Mse += data[i,25]/units[1]
        f.write("%f %f %f %f\n" %(data[i,0]*units[7], data[i,4], data[0,4]-sum_Mse, data[i,25]/units[1]))
    
    f.close()

def mass_vs_tdiss(filestring):

    dynfile=filestring+'.dyn.dat'
    convfile=filestring+'.conv.sh'

    units=read_units(convfile)
    data=loadtxt(dynfile, usecols=(0,4))
    tdiss=[-100, -100, -100, -100, -100, -100, -100, -100]
    count=0
    for j in (0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.05):
        success=0
        for i in range(len(data)):
            if data[i,1]<j and success==0:
                t1=data[i-1,0]
                t2=data[i,0]
                m1=data[i-1,1]
                m2=data[i,1]
                t=t1 - (j-m1)*(t2-t1)/(m2-m1)
                tdiss[count]=t*units[7]
                success=1
                count += 1
    #print "%f %f %f %f %f %f %f %f %f" %(units[1], tdiss[0], tdiss[1], tdiss[2], tdiss[3], tdiss[4], tdiss[5], tdiss[6], tdiss[7])


def sdss_rdist_cumulative(hrdiagfile, mcut, nbin):
    """take the file created by hrdiag_L_T and create a cumulative 
    radial distribution normalized with the RGBs"""
    data = loadtxt(hrdiagfile)
    bss_sing_r = []
    bss_bin_r = []
    bss_all_r = []
    rgb_r = []
    for i in range(len(data)):
        if data[i,0]==0:  #singles
            if data[i,6]>mcut and data[i,1]<2:  #BSSs
                bss_all_r.append(data[i,8])
                bss_sing_r.append(data[i,8])
        if data[i,0]==1: #binaries
            if data[i,6]>mcut and data[i,1]<2:  #BSSs
                bss_all_r.append(data[i,8])
                bss_bin_r.append(data[i,8])
            if data[i,7]>mcut and data[i,2]<2:  #BSSs
                bss_all_r.append(data[i,8])
                bss_bin_r.append(data[i,8])
        if data[i,1]==4 or data[i,2]==4:
            rgb_r.append(data[i,8])

    maxr_g, minr_g = max(rgb_r), min(rgb_r)
    maxr_bss = max(bss_all_r)
    if maxr_g >= maxr_bss:
        maxr = maxr_g
    else:
        maxr = maxr_bss
    minr = minr_g

    #print len(bss_all_r), len(bss_sing_r), len(bss_bin_r), len(rgb_r)
    
    cum_bss_all = cumulate_array(bss_all_r, minr, maxr, nbin)
    cum_bss_sing = cumulate_array(bss_sing_r, minr, maxr, nbin)
    cum_bss_bin = cumulate_array(bss_bin_r, minr, maxr, nbin)
    cum_rgb = cumulate_array(rgb_r, minr, maxr, nbin)

    binsize = (maxr-minr)/nbin
    for i in range(len(cum_rgb)):
        print(cum_rgb[i][0], cum_bss_all[i][1], cum_bss_sing[i][1], cum_bss_bin[i][1], cum_rgb[i][1])
        

def cumulate_array(a, minvalue, maxvalue, nbin):
    """Cumulates from array a from minvalue with steps of 
    binsize=(maxvalue-minvalue)/nbin and #prints out the values"""
    binsize = (maxvalue-minvalue)/nbin
    c_array = []
    for i in range(nbin+1):
        c = 0
        a_temp = minvalue+(i+1)*binsize
        for j in range(len(a)):
            if a[j] < a_temp:
                c += 1
        c_array.append([a_temp, c])
    return c_array
        

def nbss_mcut(hrdiagfile, mcut):
    data = loadtxt(hrdiagfile)
    for j in range(100):
        n, ns, nbin = 0, 0, 0
        for i in range(len(data)):
            if data[i,0]==0 and data[i,1]<2 and data[i,6]>mcut:
                n += 1
                ns += 1
            elif data[i,1]==1:
                if data[i,1]<2 and data[i,6]>mcut:
                    n += 1
                    nbin += 1
                elif data[i,2]<2 and data[i,7]>mcut:
                    n += 1
                    nbin += 1
            print(mcut, n, ns, nbin)
            mcut += 0.01


def make_3D(r):
    """takes some quantity, like the radial position of a star and creates the 3D position 
    by randomly creating the angles \theta and \phi.  Assumes spherical symmetry.  """
    import random
    costheta = random.uniform(-1, 1)
    sintheta = (1-costheta**2.)**0.5
    phi = random.uniform(0, 4*pi)
    rz = r*costheta
    rx = r*sintheta*cos(phi)
    ry = r*sintheta*sin(phi)

    return (rx, ry, rz)

def surface_density_profile_L(filestring, snapno, seedy, proj, nbins, bintype):
    """takes the filestring for the run, snapno of interest, seed the random generator, 
    projections like (0, 1, 2) gives 3D projection, (0, 1) gives 2d suppressing z
    no of bins wanted and the binning type, 0:equal member in each bin 1:bins are equidistant 
    in log10(r2D)"""
    #first learn the physical units
    units = read_units(filestring)
    lpc = units[0]['l_pc']
    
    writefile=filestring+'.snap'+snapno+'2D_sbp.dat'
    writefile=open(writefile, 'w')

    #read the snapfile
    snapfile = filestring+'.snap'+snapno+'.dat.gz'
    colnos = (2, 7, 15, 19, 20, 14, 17, 18)
    data = loadtxt(snapfile, usecols=colnos)

    dtype = [('r', float), ('x', float), ('y', float), ('z', float), ('r2D', float), ('logr2D', float), ('binflag', int), ('L', float), ('type0', int), ('bintype0', int), ('bintype1', int)]
    a = []
    import random
    random.seed(seedy)
    valid_line = 1
    for i in range(len(data)):
        try:
            for j in range(len(data[i])):
                if str(data[i,j])=='nan' or str(data[i,j])=='inf':
                    valid_line = 0
                    raise StopIteration()
                else:
                    valid_line = 1
        except StopIteration:
            #print "there was a nan or inf"
            pass

        if valid_line == 1:
            rs = make_3D(data[i,0])
            r2d = 0
            for j in proj:
                r2d += rs[j]**2. 
            r2d = r2d**0.5

            if data[i,1]==0:        #singles
                a.append((data[i,0]*lpc, rs[0]*lpc, rs[1]*lpc, rs[2]*lpc, r2d*lpc, log10(r2d*lpc), data[i,1], data[i,2], data[i,5], data[i,6], data[i,7]))
            elif data[i,1]==1:        #binaries
                a.append((data[i,0]*lpc, rs[0]*lpc, rs[1]*lpc, rs[2]*lpc, r2d*lpc, log10(r2d*lpc), data[i,1], data[i,3]+data[i,4], data[i,5], data[i,6], data[i,7]))

    a = array(a, dtype=dtype)
    unsorted_a = a
    a = sort(a, order='logr2D')
    #print a[:]['logr2D']
    #for equal member binning
    if bintype==0:
        binsize = len(a)/nbins  
        for i in range(binsize, len(a), binsize):
            if a[i]['type0']!=3 and a[i]['bintyppe0']!=3 and a[i]['bintype1']!=3:
                lsum = sum(a[i-binsize:i]['L'])
            ##print lsum, i, binsize, i-binsize
            #for j in range(i-binsize, i):
            #    #print a[j]['L']
            area = 4.*pi*(a[i]['r2D']**2. - a[i-binsize]['r2D']**2.)
            lsum, lsum_err = lsum/area, lsum/float(binsize**0.5)/area
            n2d, n2d_err = binsize/area, float(binsize**0.5)/area
            r2d_av = (a[i]['r2D']+a[i-binsize]['r2D'])/2.
            writefile.write("%f %f %f %f %f\n" %(r2d_av, lsum, lsum_err, n2d, n2d_err))
            ##print r2d_av, lsum, lsum_err, n2d, n2d_err


    #for binning equidistant in logr2D
    elif bintype==1:
        lbinsize = abs((a[-1]['logr2D']-a[0]['logr2D'])/nbins)
        n2d_prev=0    
        for i in range(1, nbins):
            lsum, lsum_err, rsum, n2d, n2d_err = 0., 0., 0., 0., 0.
            lsum_cut, lsum_cut_err, n2d_cut, n2d_cut_err = 0., 0., 0., 0.
            lr_upperbound, lr_lowerbound = a[0]['logr2D']+i*lbinsize, a[0]['logr2D']+(i-1)*lbinsize
            lr2D_av = (lr_upperbound+lr_lowerbound)/2.
            area = pi*( (10.**lr_upperbound)**2. - (10**lr_lowerbound)**2. )
            try:
                for j in range(int(n2d_prev), len(a)):
                    if a[j]['logr2D']<lr_upperbound and a[j]['logr2D']>=lr_lowerbound:
                        lsum += a[j]['L'] 
                        n2d += 1
                        if a[j]['L']<=20.:
                            lsum_cut += a[j]['L']
                            n2d_cut += 1
                    else:
                        raise StopIteration()
            except StopIteration:
                print(n2d, n2d_prev, 'got the values')
            n2d_prev += n2d
            ##print lr2D_av, n2d
            if n2d>1:    #just to avoid the points with 100% poisson error
                lsum, lsum_err = lsum/area, lsum/float(n2d)**(0.5)/area
                n2d, n2d_err = n2d/area, float(n2d)**0.5/area
                if n2d_cut>2:
                    lsum_cut, lsum_cut_err = lsum_cut/area, lsum_cut/float(n2d_cut)**(0.5)/area
                    n2d_cut, n2d_cut_err = n2d_cut/area, float(n2d_cut)**0.5/area
                
                
                writefile.write("%f %f %f %f %f %f %f %f %f\n" %(10**lr2D_av, lsum, lsum_err, n2d, n2d_err, lsum_cut, lsum_cut_err, n2d_cut, n2d_cut_err))
                ##print 10**lr2D_av, lsum, lsum_err, n2d, n2d_err


    writefile.close()

    return a



def surface_density_profile_L_MS_giants(filestring, snapno, seedy, proj, nbins, bintype, MS_cut, giant_cut):
    """takes the filestring for the run, snapno of interest, seed the random generator, 
    projections like (0, 1, 2) gives 3D projection, (0, 1) gives 2d suppressing z
    no of bins wanted and the binning type, 0:equal member in each bin 1:bins are equidistant 
    in log10(r2D)"""
    #first learn the physical units
    units = read_units(filestring)
    lpc = units[0]['l_pc']
    
    writefile=filestring+'.snap'+snapno+'2D_sbp_MS_cut_'+str(MS_cut)+'giants_cut_'+str(giant_cut)+'.dat'
    writefile=open(writefile, 'w')
    writefile.write("#1.r2D(pc) 2.L_dens(Lsun/pc^2) 3.L_dens_err(Lsun/pc^2) 4.N_dens(1/pc^2) 5.N_densn_err(1/pc^2) 6.L_dens_cut(Lsun/pc^2) 7.L_dens_cut_err(Lsun/pc^2) 8.N_dens_cut(1/pc^2) 9.N_dens_cut_err(1/pc^2)\n")

    #read the snapfile
    snapfile = filestring+'.snap'+snapno+'.dat.gz'
    colnos = (2, 7, 15, 19, 20, 14, 17, 18)
    data = loadtxt(snapfile, usecols=colnos)

    dtype = [('r', float), ('x', float), ('y', float), ('z', float), ('r2D', float), ('logr2D', float), ('binflag', int), ('L', float), ('type0', int), ('bintype0', int), ('bintype1', int)]
    a = []
    import random
    random.seed(seedy)
    valid_line = 1
    for i in range(len(data)):
        try:
            for j in range(len(data[i])):
                if str(data[i,j])=='nan' or str(data[i,j])=='inf':
                    valid_line = 0
                    raise StopIteration()
                else:
                    valid_line = 1
        except StopIteration:
            #print "there was a nan or inf"
            pass

        if valid_line == 1:
            rs = make_3D(data[i,0])
            r2d = 0
            for j in proj:
                r2d += rs[j]**2. 
            r2d = r2d**0.5

            if data[i,1]==0:        #singles
                a.append((data[i,0]*lpc, rs[0]*lpc, rs[1]*lpc, rs[2]*lpc, r2d*lpc, log10(r2d*lpc), data[i,1], data[i,2], data[i,5], data[i,6], data[i,7]))
            elif data[i,1]==1:        #binaries
                a.append((data[i,0]*lpc, rs[0]*lpc, rs[1]*lpc, rs[2]*lpc, r2d*lpc, log10(r2d*lpc), data[i,1], data[i,3]+data[i,4], data[i,5], data[i,6], data[i,7]))

    a = array(a, dtype=dtype)
    unsorted_a = a
    a = sort(a, order='logr2D')
    #print a[:]['logr2D']
    #for equal member binning
    if bintype==0:
        binsize = len(a)/nbins  
        for i in range(binsize, len(a), binsize):
            if a[i]['type0']!=3 and a[i]['bintyppe0']!=3 and a[i]['bintype1']!=3:
                lsum = sum(a[i-binsize:i]['L'])
            ##print lsum, i, binsize, i-binsize
            #for j in range(i-binsize, i):
            #    #print a[j]['L']
            area = 4.*pi*(a[i]['r2D']**2. - a[i-binsize]['r2D']**2.)
            lsum, lsum_err = lsum/area, lsum/float(binsize**0.5)/area
            n2d, n2d_err = binsize/area, float(binsize**0.5)/area
            r2d_av = (a[i]['r2D']+a[i-binsize]['r2D'])/2.
            writefile.write("%f %f %f %f %f\n" %(r2d_av, lsum, lsum_err, n2d, n2d_err))
            ##print r2d_av, lsum, lsum_err, n2d, n2d_err


    #for binning equidistant in logr2D
    elif bintype==1:
        lbinsize = abs((a[-1]['logr2D']-a[0]['logr2D'])/nbins)
        n2d_prev=0    
        for i in range(1, nbins):
            lsum, lsum_err, rsum, n2d, n2d_err = 0., 0., 0., 0., 0.
            lsum_cut, lsum_cut_err, n2d_cut, n2d_cut_err = 0., 0., 0., 0.
            lr_upperbound, lr_lowerbound = a[0]['logr2D']+i*lbinsize, a[0]['logr2D']+(i-1)*lbinsize
            lr2D_av = (lr_upperbound+lr_lowerbound)/2.
            area = pi*( (10.**lr_upperbound)**2. - (10**lr_lowerbound)**2. )
            try:
                for j in range(int(n2d_prev), len(a)):
                    #choose only stars that are above 0.2 Msun and less than 20 Lsun. Exclude compact objects
                    #single, not compact object, above MS_cut luminosity, below giant cut luminosity                                        ##print 'came here'
                        #Now bin them efficiently
                    if a[j]['logr2D']<lr_upperbound and a[j]['logr2D']>=lr_lowerbound:
                        #print 'came here'
                        #if (a[j]['type0']<10 and a[j]['type0']>-1) or (a[j]['bintype0']>-1 and (a[j]['bintype0']<10 or a[j]['bintype1']<10)) :

                        lsum += a[j]['L'] 
                        n2d += 1
                        if a[j]['L']>MS_cut and a[j]['L']<giant_cut:
                            lsum_cut += a[j]['L']
                            n2d_cut += 1
                    else:
                        raise StopIteration()
            except StopIteration:
                print(n2d, n2d_prev, 'got the values')
            n2d_prev += n2d
            ##print lr2D_av, n2d
            if n2d>1:    #just to avoid the points with 100% poisson error
                lsum, lsum_err = lsum/area, lsum/float(n2d)**(0.5)/area
                n2d, n2d_err = n2d/area, float(n2d)**0.5/area
                if n2d_cut>2:
                    lsum_cut, lsum_cut_err = lsum_cut/area, lsum_cut/float(n2d_cut)**(0.5)/area
                    n2d_cut, n2d_cut_err = n2d_cut/area, float(n2d_cut)**0.5/area
                
                
                writefile.write("%f %f %f %f %f %f %f %f %f\n" %(10**lr2D_av, lsum, lsum_err, n2d, n2d_err, lsum_cut, lsum_cut_err, n2d_cut, n2d_cut_err))
                #print 10**lr2D_av, lsum, lsum_err, n2d, n2d_err


    writefile.close()

    return a









def get_rc_obs(filestring, snapno, skipno, maxno, sampleno):
    """Need to run surface_density_profile_L first.  This script uses the surface brightness profile 
    files written as a part of the surface_density_profile_L script.  
    Makes the max L from the SBP by avaeraging over the first maxno points
    Then it runs a running avaerage of sampleno over all the subsequent points to find where for the 
    first time the running average reaches 0.5*max L
    skipno denotes the number of data points to skip for the max L calculation.  this is useful since 
    the first few datapoints in the SBP is plain noise
    This is the rc obs"""

    sbpfile = filestring+'.snap'+snapno+'2D_sbp.dat'
    #print sbpfile
    sbpdata = loadtxt(sbpfile)

    maxl = mean( sbpdata[skipno:skipno+maxno, 1] )
    maxl_cut = mean( sbpdata[skipno:skipno+maxno, 5] )
    maxn = mean( sbpdata[skipno:skipno+maxno, 3] )
    maxn_cut = mean( sbpdata[skipno:skipno+maxno, 7] )

    #print maxl, maxl_cut, maxn, maxn_cut

    iskip = skipno+maxno
    success_l, success_l_cut, success_n, success_n_cut = 0, 0, 0, 0
    
    try:
        for i in range(len(sbpdata)):
            temp_l = mean( sbpdata[iskip+i:iskip+i+sampleno, 1] )
            temp_l_cut = mean( sbpdata[iskip+i:iskip+i+sampleno, 5] )
            temp_n = mean( sbpdata[iskip+i:iskip+i+sampleno, 3] )
            temp_n_cut = mean( sbpdata[iskip+i:iskip+i+sampleno, 7] )
            #print maxl, maxl_cut, maxn, maxn_cut, temp_l, temp_l_cut, temp_n, temp_n_cut
            if temp_l<=0.5*maxl and success_l==0:
                rc_l, rc_l_err = (sbpdata[iskip+i, 0] + sbpdata[iskip+i+sampleno,0])/2., (sbpdata[iskip+i+sampleno,0] - sbpdata[iskip+i, 0])/2.
                success_l = 1
                #print 'rc_l, rc_l_err', rc_l, rc_l_err
            if temp_l_cut<=0.5*maxl_cut and success_l_cut==0:
                rc_l_cut, rc_l_cut_err = (sbpdata[iskip+i, 0] + sbpdata[iskip+i+sampleno,0])/2., (sbpdata[iskip+i+sampleno,0] - sbpdata[iskip+i, 0])/2.
                success_l_cut = 1
                #print 'rc_l_cut, rc_l_cut_err', rc_l_cut, rc_l_cut_err
            if temp_n<=0.5*maxn and success_n==0:
                rc_n, rc_n_err = (sbpdata[iskip+i, 0] + sbpdata[iskip+i+sampleno,0])/2., (sbpdata[iskip+i+sampleno,0] - sbpdata[iskip+i, 0])/2.
                success_n = 1
                #print 'rc_n, rc_n_err', rc_n, rc_n_err
            if temp_n_cut<=0.5*maxn_cut and success_n_cut==0:
                rc_n_cut, rc_n_cut_err = (sbpdata[iskip+i, 0] + sbpdata[iskip+i+sampleno,0])/2., (sbpdata[iskip+i+sampleno,0] - sbpdata[iskip+i, 0])/2.
                success_n_cut = 1
                #print 'rc_n_cut, rc_n_cut_err', rc_n_cut, rc_n_cut_err
            if success_l == 1 and success_l_cut == 1 and success_n == 1 and success_n_cut == 1:
                raise StopIteration()
    except StopIteration:
        pass


    return rc_l, rc_l_err, rc_l_cut, rc_l_cut_err, rc_n, rc_n_err, rc_n_cut, rc_n_cut_err


def get_half_light_obs(filestring, snapno, seedy, proj, nbins, bintype):
    a = surface_density_profile_L(filestring, snapno, seedy, proj, nbins, bintype)
    
    total_l = 0.
    total_l_cut=0.
    for i in range(len(a)):
        if str(a[i]['L'])!='nan' and str(a[i]['logr2D'])!='nan':
            total_l += a[i]['L']
            if a[i]['L']<20.:
                total_l_cut += a[i]['L']
    half_light = total_l*0.5
    half_light_cut = total_l_cut*0.5
    #print total_l, total_l_cut, half_light, half_light_cut

    light, light_cut = 0., 0.
    for i in range(len(a)):
        if str(a[i]['L'])!='nan' and str(a[i]['logr2D'])!='nan':
            light_lower = light
            light = light + a[i]['L']
            light_upper = light
            if light_lower <= half_light and light_upper > half_light:
                #print 'checking for rh'
                rh_lower, rh_upper = a[i]['r2D'], a[i+1]['r2D']
                #print half_light, rh_lower, rh_upper
                rh = rh_lower + (half_light - light_lower) * (rh_upper - rh_lower)/(light_upper-light_lower)
                rh_err = (rh_upper - rh_lower)/2.
            #now find rh with a cutoff L<20LSun
            if a[i]['L']<20.:
                light_cut_lower = light_cut
                light_cut = light_cut+a[i]['L']
                light_cut_upper = light_cut
                if light_cut_lower <= half_light_cut and light_cut_upper > half_light_cut:
                    #print 'checking for rh_cut' 
                    rh_lower_cut, rh_upper_cut = a[i]['r2D'], a[i+1]['r2D']
                    #print half_light_cut, rh_lower_cut, rh_upper_cut
                    rh_cut = rh_lower_cut + (half_light_cut - light_cut_lower) * (rh_upper_cut - rh_lower_cut)/(light_cut_upper-light_cut_lower)
                    rh_cut_err = (rh_upper_cut - rh_lower_cut)/2.

    return rh, rh_err, rh_cut, rh_cut_err




def get_rcobs_rhobs(filestring, snapno, seedy, proj, nbins, bintype):
    """takes the filestring for the run, snapno of interest, seed the random generator, 
    projections like (0, 1, 2) gives 3D projection, (0, 1) gives 2d suppressing z
    no of bins wanted and the binning type, 0:equal member in each bin 1:bins are equidistant 
    in log10(r2D)"""
    #first learn the physical units
    units = read_units(filestring)
    lpc = units[0]['l_pc']

    #read the snapfile
    snapfile = filestring+'.snap'+snapno+'.dat.gz'
    colnos = (2, 7, 15, 19, 20, 14, 17, 18)
    data = loadtxt(snapfile, usecols=colnos)

    dtype = [('r', float), ('x', float), ('y', float), ('z', float), ('r2D', float), ('logr2D', float), ('binflag', int), ('L', float), ('type0', int), ('bintype0', int), ('bintype1', int)]
    a = []
    import random
    random.seed(seedy)
    for i in range(len(data)):
        if str(data[i,0])!='nan' and str(data[i,1])!='nan' and str(data[i,2])!='nan' and str(data[i,3])!='nan' and str(data[i,4])!='nan' and str(data[i,5])!='nan' and str(data[i,6])!='nan' and str(data[i,7])!='nan' and str(data[i,0])!='inf' and str(data[i,1])!='inf' and str(data[i,2])!='inf' and str(data[i,3])!='inf' and str(data[i,4])!='inf' and str(data[i,5])!='inf' and str(data[i,6])!='inf' and str(data[i,7])!='inf':
            rs = make_3D(data[i,0])
            r2d = 0
            for j in proj:
                r2d += rs[j]**2. 
            r2d = r2d**0.5

            if data[i,1]==0:        #singles
                a.append((data[i,0]*lpc, rs[0]*lpc, rs[1]*lpc, rs[2]*lpc, r2d*lpc, log10(r2d*lpc), data[i,1], data[i,2], data[i,5], data[i,6], data[i,7]))
            elif data[i,1]==1:        #binaries
                a.append((data[i,0]*lpc, rs[0]*lpc, rs[1]*lpc, rs[2]*lpc, r2d*lpc, log10(r2d*lpc), data[i,1], data[i,3]+data[i,4], data[i,5], data[i,6], data[i,7]))

    a = array(a, dtype=dtype)
    a = sort(a, order='logr2D')
    #print a[:]['logr2D']
    #print a[:]['L']

    total_l = sum(a[:]['L'])
    half_light = total_l*0.5

    light = 0.
    for i in range(len(a)):
        light_lower = light
        light = light + a[i]['L']
        light_upper = light
        if light_lower <= half_light and light_upper > half_light:
            #print 'checking for rh'
            rh_lower, rh_upper = a[i]['r2D'], a[i+1]['r2D']
            #print half_light, rh_lower, rh_upper
            rh = rh_lower + (half_light - light_lower) * (rh_upper - rh_lower)/(light_upper-light_lower)
    
    #finding rcobs
    obsparam = zeros((nbins, 5))
    #for equal member binning
    if bintype==0:
        binsize = len(a)/nbins  
        for i in range(binsize, len(a), binsize):
            if a[i]['type0']!=3 and a[i]['bintyppe0']!=3 and a[i]['bintype1']!=3:
                lsum = sum(a[i-binsize:i]['L'])
            ##print lsum, i, binsize, i-binsize
            #for j in range(i-binsize, i):
            #    #print a[j]['L']
            area = 4.*pi*(a[i]['r2D']**2. - a[i-binsize]['r2D']**2.)
            lsum, lsum_err = lsum/area, lsum/float(binsize**0.5)/area
            n2d, n2d_err = binsize/area, float(binsize**0.5)/area
            r2d_av = (a[i]['r2D']+a[i-binsize]['r2D'])/2.
            #writefile.write("%f %f %f %f %f\n" %(r2d_av, lsum, lsum_err, n2d, n2d_err))
            ##print r2d_av, lsum, lsum_err, n2d, n2d_err


    #for binning equidistant in logr2D
    elif bintype==1:
        lbinsize = abs((a[-1]['logr2D']-a[0]['logr2D'])/nbins)
        n2d_prev=0
        indarray = 0
        for i in range(1, nbins):
            lsum, lsum_err, rsum, n2d, n2d_err = 0., 0., 0., 0., 0.
            lsum_cut, lsum_cut_err, n2d_cut, n2d_cut_err = 0., 0., 0., 0.
            lr_upperbound, lr_lowerbound = a[0]['logr2D']+i*lbinsize, a[0]['logr2D']+(i-1)*lbinsize
            lr2D_av = (lr_upperbound+lr_lowerbound)/2.
            area = pi*( (10.**lr_upperbound)**2. - (10**lr_lowerbound)**2. )
            try:
                for j in range(n2d_prev, len(a)):
                    if a[j]['logr2D']<lr_upperbound and a[j]['logr2D']>=lr_lowerbound:
                        lsum += a[j]['L'] 
                        n2d += 1
                        if a[j]['L']<=20.:
                            lsum_cut += a[j]['L']
                            n2d_cut += 1
                    else:
                        raise StopIteration()
            except StopIteration:
                print(n2d, n2d_prev, 'got the values')
            n2d_prev += n2d
            ##print lr2D_av, n2d
            if n2d>1:    #just to avoid the points with 100% poisson error
                lsum, lsum_err = lsum/area, lsum/float(n2d)**(0.5)/area
                lsum_cut, lsum_cut_err = lsum_cut/area, lsum_cut/float(n2d_cut)**(0.5)/area
                n2d, n2d_err = n2d/area, float(n2d)**0.5/area
                n2d_cut, n2d_cut_err = n2d_cut/area, float(n2d_cut)**0.5/area
                #writefile.write("%f %f %f %f %f %f %f %f %f\n" %(10**lr2D_av, lsum, lsum_err, n2d, n2d_err, lsum_cut, lsum_cut_err, n2d_cut, n2d_cut_err))
                obsparam[indarray,0],obsparam[indarray,1],obsparam[indarray,2],obsparam[indarray,3], obsparam[indarray,4] = 10**lr2D_av, lsum, lsum_err, lsum_cut, lsum_cut_err
                indarray += 1
                ##print 10**lr2D_av, lsum, lsum_err, n2d, n2d_err

        #print obsparam
        maxlsum, maxlsumcut = mean(obsparam[:10,1]), mean(obsparam[:10,3])
        #print maxlsum, maxlsumcut
        
        success_rc, k = 0, 1
        while success_rc==0:
            templ = mean(obsparam[k:k+15, 1])
            if templ <= 0.5*maxlsum:
                rcobs = obsparam[k,0]
                #print "rcobs found", obsparam[k,0]
                success_rc=1
            else:
                #print "rcobs not found", obsparam[k,0]
                k+=1

        success_rc_cut, k = 0, 1
        while success_rc_cut==0:
            templ = mean(obsparam[k:k+15, 3])
            if templ <= 0.5*maxlsumcut:
                rcobscut = obsparam[k,0]
                #print "rcobscut found", obsparam[k,0]
                success_rc_cut=1
            else:
                #print "rcobscut not found", obsparam[k,0]
                k+=1


    return rcobs, rcobscut, rh

    


def king_fit(r, n, rc_guess, rt_guess):
    """fits a single mass King profile given the r array, the n array and 
    the initial guesses for rc and rt"""
    rc = []
    rt = []
    for i in range(100):
        rc.append(rc_guess*(1.-50*0.01)+0.01*rc_guess*i)
        rt.append(rt_guess*(1.-50*0.01)+0.01*rt_guess*i)
    a = []
    dtype = [('rc', float), ('rt', float), ('norm', float), ('sqsum', float)]
    for i in range(len(rc)):
        for j in range(len(rt)):
            norm = n[len(n)/2]/( ( 1/ (1+ (r[len(r)/2]/rc[i])**2.)**0.5 ) - ( 1/ (1+ (rt[j]/rc[i])**2.)**0.5 ) )**2.
            sqsum = 0.
            for k in range(len(n)):
                n_k = norm*( 1/( 1+ (r[k]/rc[i])**2. )**0.5 - 1/( 1+ (rt[j]/rc[i])**2. )**0.5 )**2.
                sqsum += (n[k] - n_k)**2.
            a.append( (rc[i], rt[j], norm, sqsum) ) 
    a = array(a, dtype=dtype)
    a = sort(a, order='sqsum')

    return a

def smooth_data(filename, writefile, smoothtype, smoothlength):
    """smoothens data with the smoothing length smoothlength and writes in writefile
    smoothing length can be either number of lines or intervals of time: type=n or t, respectively."""
    data=loadtxt(filename)
    writefile=open(writefile, 'w')
    for j in range(len(data[0,:])):
        writefile.write("%f " %(data[0,j]))
    writefile.write("\n")
    counter=[]

    if smoothtype=='n':
        for i in range(smoothlength, len(data)):
                if i%smoothlength==0:
                    for j in range(len(data[i,:])):
                            temp=mean(data[i-smoothlength:i,j])
                            writefile.write("%f " %(temp))
                    writefile.write("\n")
    if smoothtype=='t':
        interval = data[-1,0]/smoothlength
        #print interval
        count = 1
        counter = []
        for i in range(1,len(data)):
            #print data[i,0], count*interval, count, data[i-1,0], data[i,0]>=count*interval and data[i-1,0]<count*interval
            if data[i,0]>=count*interval and data[i-1,0]<count*interval:
                count += 1
                counter.append(i)
        for i in range(1,len(counter)):
            for j in range(len(data[counter[i],:])):
                temp=mean( data[counter[i-1]:counter[i], j] )
                writefile.write("%f " %(temp))
            writefile.write("\n")
    for j in range(len(data[-1,:])):
        writefile.write("%f " %(data[-1,j]))
    writefile.write("\n")


    writefile.close()

    return smoothtype, counter


def prune_data(filename, writefilename, prunetype, prunelength):
    """prunens data with the pruning length prunelength and writes in writefile
    pruneing length can be either number of lines or intervals of time: type=n or t, respectively."""
    f = open(filename, 'r')
    writefile=open(writefilename, 'w')
    f.seek(0)
    for line in f:
        if line.rfind('#')>-1: #Commented lines
            writefile.write("%s" %(line))
    f.close()
            
    data=loadtxt(filename)
    for j in range(len(data[0,:])):
        writefile.write("%f " %(data[0,j]))
    writefile.write("\n")
    counter=[]
    if prunetype=='n':
        for i in range(prunelength, len(data)):
                if i%prunelength==0:
                    for j in range(len(data[i,:])):
                            temp=data[i,j]
                            writefile.write("%f " %(temp))
                    writefile.write("\n")
    if prunetype=='t':
        interval = data[-1,0]/prunelength
        #print interval
        count = 1
        counter = []
        for i in range(1,len(data)):
            #print data[i,0], count*interval, count, data[i-1,0], data[i,0]>=count*interval and data[i-1,0]<count*interval
            if data[i,0]>=count*interval and data[i-1,0]<count*interval:
                count += 1
                counter.append(i)
        for i in range(1,len(counter)):
            for j in range(len(data[counter[i],:])):
                temp=data[counter[i], j] 
                writefile.write("%f " %(temp))
            writefile.write("\n")
    for j in range(len(data[-1,:])):
        writefile.write("%f " %(data[-1,j]))
    writefile.write("\n")


    writefile.close()

    return prunetype, counter


def get_roche_r(q, a):
    R = a*0.49*q**(2./3.) / (0.6*q**(2./3.) + log(1+q**(1./3.)))

    return R

def get_period(a,m1,m2): 
    """Given a, m1, and m2 of a binary calculate keplerian period.  All input in CGS.  Output in days.  """
    p = (4.*pi**2.*(a**3.)/constants.G/(m1+m2))**0.5 /3600./24.    #p in days

    return p

def get_semimajor_axis(p,m1,m2): 
    """Given period, m1, and m2 of a binary calculate keplerian period.  All input in CGS.  Output in AU.  """
    a = p**(2./3.) * ( constants.G*(m1+m2)/4./pi**2. )**(1./3.) / constants.AU 
    return a

def plot_P_e(bsfilename):
    data = loadtxt(bsfilename)
    p = []
    for i in range(len(data)):
        p.append(get_period(data[i,12]*constants.AU, data[i,8]*constants.Msun, data[i,9]*constants.Msun))
    
    
    #separate
    msms=[]
    msg=[]
    for i in range(len(data)):
        if data[i,17]<2 and data[i,18]<2:
            msms.append((p[i],data[i,13]))
        else:
            msg.append((p[i],data[i,13]))




    dtype = [('p', float), ('e', float)] 
    msms = array(msms, dtype=dtype)
    msg = array(msg, dtype=dtype)

    #print msms[:]['p'], msms[:]['e']
    #now plot
    import gracePlot
    gpl=gracePlot.gracePlot()
    gpl.hold()
    gpl.plot(msms[:]['p'], msms[:]['e'], symbols=1)
    gpl.plot(msg[:]['p'], msg[:]['e'], symbols=1)

    return len(msms), len(msg)


def plot_hrdiag(filename, colnos):
    data=loadtxt(filename, usecols=colnos)
    single=[]
    binary=[]
    dtype=[('B-V', float), ('B', float)]
    for i in range(len(data)):
        if data[i,0]==1:
            binary.append((data[i,1]-data[i,2], data[i,1]))
        else:
            single.append((data[i,1]-data[i,2], data[i,1]))
    single=array(single,dtype=dtype)
    binary=array(binary,dtype=dtype)

    import gracePlot
    gpl=gracePlot.gracePlot()
    gpl.hold()
    gpl.plot(binary[:]['B-V'], binary[:]['B'], symbols=1)
    gpl.plot(single[:]['B-V'], single[:]['B'], symbols=1)

    


def get_fb_bss(filestring):
    from glob import glob
    name=filestring+'.snap*.dat.gz'
    filenames = glob(name)[-3:]
    t, bss_tot, fb_bss = [], [], []
    for i in range(len(filenames)):
        temp=filenames[i]
        n1=temp.split('.snap')[0]
        n2=temp.split('.snap')[1].split('.dat')[0]
        bss = find_BSS(n1, n2)
        t.append(bss[2]) 
        bss_tot.append((bss[0]+bss[1])) 
        fb_bss.append(float(bss[1])/float(bss[0]+bss[1]))
    t_myr, dt_myr = mean(t), std(t)
    tot_bss, dtot_bss = mean(bss_tot), std(bss_tot)
    bss_fb, dbss_fb = mean(fb_bss), std(fb_bss)
    return (t_myr, dt_myr, tot_bss, dtot_bss, bss_fb, dbss_fb)

def m_P_correlation(bsfilenames):
    m=[]
    dtype=[('m', float), ('type', int), ('P', float)]

    for k in range(len(bsfilenames)):
        data = loadtxt(bsfilenames[k])
        p, m0, m1, type0, type1 = [], [], [], [], []
        for i in range(len(data)):
            if data[i,25]==0:
                p.append(get_period(data[i,12]*constants.AU, data[i,8]*constants.Msun, data[i,9]*constants.Msun))
                m0.append(data[i,8])
                m1.append(data[i,9])
                type0.append(data[i,17])
                type1.append(data[i,18])

        for i in range(len(m0)):
            if type0[i]<2 and type1[i]<2:
                if m0[i]>m1[i]:
                    m.append((m0[i], type0[i], p[i]))
                else:
                    m.append((m1[i], type1[i], p[i]))    #if MSS-MSS take the more massive to be BSS mass
            elif type0[i]<2 and type1[i]>1:
                m.append((m0[i], type0[i], p[i]))
            elif type0[i]>1 and type1[i]<2:
                m.append((m1[i], type1[i], p[i]))
    
    m = array(m, dtype=dtype)
    cor = two_point_correlation(m[:]['m'], m[:]['P'])

    #now plot
    import gracePlot
    gpl=gracePlot.gracePlot()
    gpl.hold()
    gpl.plot(m[:]['m'], m[:]['P'], symbols=1)


    return cor, m 
    

def distribution(dataarray, number, param):
    """takes a one D data array and creates its distribution function.  
    return: a two dim array with the PDF value and the midpoint of the 
    bin.  
        Equal number of members per bin: param=0-> number=#members in each bin
        Equal bin width: param=1 -> number=#of bins """


    
def projection(r, projecttuple, seedy, samplesize):
    """takes a vector of any dimension and projects it in any other dimension"""
    
    import random
    random.seed(seedy)
    proj_r_list = []
    for i in range(samplesize):
        costheta = random.uniform(-1, 1)
        sintheta = (1-costheta**2.)**0.5
        phi = random.uniform(0, 4*pi)
        r_vec = []
        r_vec.append(r*costheta)
        r_vec.append(r*sintheta*cos(phi))
        r_vec.append(r*sintheta*sin(phi))
        tmp = 0.
        for j in projecttuple:
            tmp += r_vec[j]**2.
        tmp = tmp**0.5
        ##print tmp
        proj_r_list.append(tmp)
    r_projected = mean(proj_r_list)
    return (r_projected)



def create_peter(filestring, snapno):
    """from a snapshot file it extracts star id, stellartype, surface g, T_eff, L, stellar mass, 
    stellar radius
    format is:
    binflag, id, id0, id1, type0, type1, m0, m1, R0, R1, L0, L1, g0, g1, T0, T1, total_L, T_eff
    for single stars (binid = 0):
    id1, m1, R1, L1, g1, T1 are all set to zero.  total_L = L0, T_eff = T0
    for binary stars (binid = 1):
    total_L = L0 + L1, T_eff = (L0*T0 + L1*T1)/(L0+L1)"""

    snapfile = filestring+'.snap'+snapno+'.dat.gz'
    writefile = filestring+'.snap'+snapno+'.extracted_params.dat'
    writefile = open(writefile, 'w')

    units = read_units(filestring)

    data = loadtxt(snapfile)

    #note time of snap
    f_snap = gzip.open(snapfile,'rb')
    line_snap = f_snap.readline()
    a_snap = line_snap.split()[1].split('=')[1]
    t_myr = float(a_snap)*units['t_myr']
    #b_snap = a_snap[1]
    #c_snap = b_snap.split('=')
    #t_snap = float(c_snap[1])
    ##print t_snap


    #initialize
    binflag, id, id0, id1, type0, type1 = 0, 0, 0, 0, 0, 0
    m0, m1, R0, R1, L0, L1, g0, g1, T0, T1, total_L, T_eff = 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. 

    ##print header
    writefile.write("#t = %f Myr\n" %(t_myr,))
    writefile.write("#1.binflag 2.id 3.id0 4.id1 5.type0 6.type1 7.m0(MSUN) 8.m1(MSUN) 9.R0(RSUN) 10.R1(RSUN) 11.L0(LSUN) 12.L1(LSUN) 13.g0(g/cm^2) 14.g1(g/cm^2) 15.T0(K) 16.T1(K) 17.total_L(LSUN) 18.T_eff(K)\n")

    for i in range(len(data)):
        if data[i,7]==0:  #singles
            binflag = 0
            id, id0, id1 = data[i,0], data[i,10], data[i,11]
            type0, type1 = data[i,14], 0
            m0, m1 = data[i,1], 0.
            R0, R1 = data[i,16], 0.
            L0, L1 = data[i,15], 0.
            g0, g1 = find_g(m0, R0), 0.
            T0, T1 = find_T(L0, R0), 0.
            total_L = (L0+L1)
            T_eff = (L0*T0 + L1*T1)/(L0+L1)
        elif data[i,7]==1:
            binflag = 1
            id, id0, id1 = data[i,0], data[i,10], data[i,11]
            type0, type1 = data[i,17], data[i,18]
            m0, m1 = data[i,8], data[i,9]
            R0, R1 = data[i,21], data[i,22]
            L0, L1 = data[i,19], data[i,20]
            g0, g1 = find_g(m0, R0), find_g(m1, R1)
            T0, T1 = find_T(L0, R0), find_T(L1, R1)
            total_L = (L0+L1)
            T_eff = (L0*T0 + L1*T1)/(L0+L1)
        #now #print
        params = (binflag, id, id0, id1, type0, type1, m0, m1, R0, R1, L0, L1, g0, g1, T0, T1, total_L, T_eff)
        for j in range(len(params)):
            writefile.write("%f " %(params[j],))
        writefile.write("\n")

    writefile.close()
    f_snap.close()
    
    return t_myr
        



def lognormal(mode, amin, amax, seedy, samplesize, nbin):
    
    
    import random
    random.seed(seedy)
    a_good = []
    b=zeros(samplesize)
    success = 0
    while success<=samplesize:
        X = random.uniform(0,1)
        a = (amin + X*(amax-amin)) 
        fa = 0.33/( (a/mode)**0.33 + (mode/a)**0.33 ) / a
        b = random.uniform(0,1) * 3.
        if b <= fa:
            success += 1
            a_good.append(a)
    #print a_good
    
    dist, lower_edges = histogram(a_good[:],bins=nbin, range=(amin, amax), normed=True )
    p, mid_a = [], []

    #norm, binwidth = 0., (amax-amin)/nbin 
    #for i in range(len(dist)):
    #    norm += dist[i]*binwidth
    for i in range(len(dist)):
        binwidth = lower_edges[i+1] - lower_edges[i]
        p.append(dist[i]/binwidth)
        mid_a.append((lower_edges[i]+lower_edges[i+1])/2.)
    
    #now check whether the above distribution is indeed following f(a)
    temp_norm = 0.
    fa_list = []
    for i in range(nbin):
        temp_a = mid_a[i]
        temp_fa = 0.33/( (temp_a/mode)**0.33 + (mode/temp_a)**0.33 ) / a
    
        temp_norm += temp_fa*binwidth
        fa_list.append(temp_fa)
    for i in range(len(fa_list)):
        fa_list[i] = fa_list[i]/temp_norm
    
    import gracePlot
    gpl=gracePlot.gracePlot()
    gpl.hold()
    gpl.plot(mid_a[:], p[:], symbols=1)
    gpl.plot(mid_a[:], fa_list[:], symbols=1)

    return (a_good, p, mid_a)


def distribution_rejection(ntry, seedy, arange, sigma, mean, samplesize):
    import random
    random.seed(seedy)
    #first define the function
    
    result=[]
    step = (arange[1]-arange[0])/ntry
    for i in range(ntry):
        a=arange[0]+i*step
        fa =  1./(a*sigma*(2*pi)**0.5) * exp (- (log(a) - mean)**2. / (2*sigma**2.) )
        #fa = 0.33/( (a/0.1)**0.33 + (0.1/a)**0.33 ) / a
        result.append((a,fa))
    result = array(result)
    #print 'max, min', max(result[:,1]), min(result[:,1])

    a_good, success = [], 0
    while success<=samplesize:
        X = random.uniform(0,1)
        a = (arange[0] + X*(arange[1]-arange[0])) 
        fa = 1./(a*sigma*(2*pi)**0.5) * exp (- (log(a) - mean)**2. / (2*sigma**2.) )
        b = random.uniform(0,1) * 3.
        if b <= fa:
            success += 1
            a_good.append(a)
    #print a_good
    #now plot
    import gracePlot
    gpl=gracePlot.gracePlot()
    gpl.hold()
    gpl.plot(result[:,0], result[:,1], symbols=1)

    return result


def lognormal_hurley_way(mode, amin, amax, seedy, samplesize, eta):
    
    
    import random
    random.seed(seedy)
    a_good = []
    success = 0
    temp = 0.5*(eta**0.33 + 1./(eta**0.33))
    k = arccos(temp**-1.)

    while success<=samplesize:
        W = -1. +random.uniform(0,1) * 2.
        temp1 = 1./( cos(k*W) ) + tan( k*W )
        #print temp1, W, k
        a = (temp1**(1./0.33)) * mode
        if a>amin and a<=amax:
            success += 1
            a_good.append(a)
    return (a_good)


def double_degenerate(directorystring, string):
    from glob import glob
    import gzip
    import history_cmc as hc

    filestring = directorystring+'/'+string+'*.dat.gz'
    files = glob(filestring)
    dynfile = directorystring+'/'+string+'.dyn.dat'
    dyndata = loadtxt(dynfile)
    writefile=directorystring+'/'+string+'.double_degenerate1.dat'
    f1=open(writefile,'w')

    for i in range(len(files)):
        #first find time
        f=gzip.open(files[i], 'rb')
        line = f.readline()
        t = float(line.split()[1].split('=')[1])
        #then find rc for that time
        try:
            for j in range(len(dyndata)):
                if dyndata[j,0]==t:
                    rcnb = dyndata[j,24]
                    raise StopIteration()
                elif dyndata[j,0]<t and dyndata[j+1,0]>t:
                    rcnb = (dyndata[j,24]+dyndata[j+1,24])/2.
                    raise StopIteration()
                else:
                    rcnb = -100.
        except StopIteration:
            pass

        #print files[i], 'rc=', rcnb, 't=', t
        #now find double degenerate binaries in the core
        data = loadtxt(files[i])
        nbc, nbcdd, primordial_nbc_dd = 0, 0, 0

        for j in range(len(data)):
            if data[j,2]<rcnb:
                if data[j,7]==1:
                    nbc += 1
                    if data[j,17]>9. and data[j,18]>9. and data[j,17]<14. and data[j,18]<14.:
                        nbcdd += 1
                        id0, id1, pos = data[j,10], data[j,11], data[j,2]
                        history = hc.history_maker([id0, id1], [pos, pos], 'test1', 1)
                        primordial = []
                        for k in history.keys():
                            if len(history[k]['binint']['binint'].keys())==0 and len(history[k]['coll'].keys())==0:
                                primordial.append(1)
                            else:
                                primordial.append(0)
                        if primordial[0]==1 and primordial[1]==1:
                            primordial_nbc_dd += 1
        
        f1.write("%g %g %g %g %g\n" %(t, rcnb, nbc, nbcdd, primordial_nbc_dd))
    


def calculate_least_sqaure(array1, array2):
    """calculates the 1D square differences and finds the rms difference between array1 and array2
    len(array1) must be = len(array2)"""
    sq_sum, not_good = 0., 0
    mean1, mean2 = mean(array1), mean(array2)
    n = len(array1)

    if len(array1) is not len(array2):
        not_good = 1

    if not_good==0:
        for i in range(len(array1)):
            sq_sum = sq_sum + (array1[i] - array2[i])**2.
    
    #adjust due to difference in means
    #sq_sum = sq_sum - n*(mean1 - mean2)**2.

    ls = sq_sum/float(n)

    return ls

def find_trh_Djorgovski(filestring, snapno):
    import glob
    import numpy as np
    a = find_sorted_2D_light(filestring, snapno, 100, (1,2,))

    total_l = 0.
    total_l_cut=0.
    for i in range(len(a)):
        if str(a[i]['L'])!='nan' and str(a[i]['logr2D'])!='nan':
            total_l += a[i]['L']
            if a[i]['L']<20.:
                total_l_cut += a[i]['L']
    half_light = total_l*0.5
    half_light_cut = total_l_cut*0.5
    #print total_l, total_l_cut, half_light, half_light_cut
    M_tot = 2*total_l
    M_mean = 1./3.
    N_tot = M_tot/M_mean


    light, light_cut = 0., 0.
    DONE1, DONE2 = 0, 0
    try:
        for i in range(len(a)):
            if str(a[i]['L'])!='nan' and str(a[i]['logr2D'])!='nan':
                light_lower = light
                light = light + a[i]['L']
                light_upper = light
                if light_lower <= half_light and light_upper > half_light:
                    #print 'checking for rh'
                    rh_lower, rh_upper = a[i]['r2D'], a[i+1]['r2D']
                    #print half_light, rh_lower, rh_upper
                    rh = rh_lower + (half_light - light_lower) * (rh_upper - rh_lower)/(light_upper-light_lower)
                    rh_err = (rh_upper - rh_lower)/2.
                    DONE1 = 1
                #now find rh with a cutoff L<20LSun
                if a[i]['L']<20.:
                    light_cut_lower = light_cut
                    light_cut = light_cut+a[i]['L']
                    light_cut_upper = light_cut
                    if light_cut_lower <= half_light_cut and light_cut_upper > half_light_cut:
                        #print 'checking for rh_cut' 
                        rh_lower_cut, rh_upper_cut = a[i]['r2D'], a[i+1]['r2D']
                        #print half_light_cut, rh_lower_cut, rh_upper_cut
                        rh_cut = rh_lower_cut + (half_light_cut - light_cut_lower) * (rh_upper_cut - rh_lower_cut)/(light_cut_upper-light_cut_lower)
                        rh_cut_err = (rh_upper_cut - rh_lower_cut)/2.
                        DONE2 = 1
            if DONE1 == 1 and DONE2 == 1:
                #print 'Found rh and rhcut'
                raise StopIteration()
    except StopIteration:
        pass


    trh = 2.055 * 1e6 * M_tot**(0.5) * rh_cut**1.5 / (M_mean * log(0.4*N_tot)) 


    return trh, rh, rh_err, rh_cut, rh_cut_err 




def find_sorted_2D_light(filestring, snapno, seedy, proj):
    #first learn the physical units
    units = read_units(filestring)
    lpc = units[0]['l_pc']
    
    #read the snapfile
    snapfile = filestring+'.snap'+snapno+'.dat.gz'
    colnos = (2, 7, 15, 19, 20, 14, 17, 18)
    data = loadtxt(snapfile, usecols=colnos)

    dtype = [('r', float), ('x', float), ('y', float), ('z', float), ('r2D', float), ('logr2D', float), ('binflag', int), ('L', float), ('type0', int), ('bintype0', int), ('bintype1', int)]
    a = []
    import random
    random.seed(seedy)
    valid_line = 1
    for i in range(len(data)):
        try:
            for j in range(len(data[i])):
                if str(data[i,j])=='nan' or str(data[i,j])=='inf':
                    valid_line = 0
                    raise StopIteration()
                else:
                    valid_line = 1
        except StopIteration:
            #print "there was a nan or inf"
            pass

        if valid_line == 1:
            rs = make_3D(data[i,0])
            r2d = 0
            for j in proj:
                r2d += rs[j]**2. 
            r2d = r2d**0.5

            if data[i,1]==0:        #singles
                a.append((data[i,0]*lpc, rs[0]*lpc, rs[1]*lpc, rs[2]*lpc, r2d*lpc, log10(r2d*lpc), data[i,1], data[i,2], data[i,5], data[i,6], data[i,7]))
            elif data[i,1]==1:        #binaries
                a.append((data[i,0]*lpc, rs[0]*lpc, rs[1]*lpc, rs[2]*lpc, r2d*lpc, log10(r2d*lpc), data[i,1], data[i,3]+data[i,4], data[i,5], data[i,6], data[i,7]))

    a = array(a, dtype=dtype)
    unsorted_a = a
    a = sort(a, order='logr2D')


    return a





def find_trh(filename, writefile, filestring):
    units = read_units(filestring)
    scalet = units['t_myr'][0]
    writefile = open(writefile, 'w')
    writefile.write("#1.t (Myr) 2.trh(Myr)\n")
    f=open(filename, 'r')
    f.seek(0)
    for line in f:
        trhpoint, tpoint = line.rfind('trh='), line.rfind('TotalTime=')
        if trhpoint>-1:
            trh = line.split('=')[1].split()[0]
            #print 'trh', float(trh)*scalet, '\n'
            writefile.write("%f\n" %(float(trh)*scalet,))
            
        if tpoint>-1:
            t = line.split('=')[2].split()[0]
            #print 't', float(t)*scalet, '\n'
            writefile.write("%f " %(float(t)*scalet,))

    f.close()
    writefile.close()


def count_massive_star(mass, nbodydatafile):
    data = loadtxt(nbodydatafile)
    mmin, mmax, mmean = min(data[:,0]), max(data[:,0]), mean(data[:,0])
    scale = 0.1/mmin

    #print mmin*scale, mmean*scale, mmax*scale

    n = 0
    for i in range(len(data)):
        if data[i,0]*scale >= mass:
            n += 1
    return n


def rhoc_for_single_ms(filestring, writefile):
    dynfile = filestring+'.dyn.dat'
    snapstring = filestring+'.snap*.dat.gz'
    from glob import glob
    snaps = glob(snapstring)
    #print snaps

    units = read_units(filestring)

    writefile=open(writefile, 'w')
    writefile.write("#1.t(Myr) 2.single_ms_rhoc(MSun/pc^3) 3.binary_ms_rhoc(MSun/pc^3) 4.all_rhoc(MSun/pc^3)\n")

    for i in range(len(snaps)):
        f_snap = gzip.open(snaps[i],'rb')
        line_snap = f_snap.readline()
        t_snap = float(line_snap.split()[1].split('=')[1]) 
        #print 't_snap', t_snap
        dyndata = loadtxt(dynfile)
        #print dyndata[-1]
        try:
            for j in range(len(dyndata)-1, -1, -1):   #start from the last line, it will be quicker
                if dyndata[j,0]==t_snap:
                    rc = dyndata[j,7]
                    rhoc = dyndata[j,21]*units[0]['m_msun']/units[0]['l_pc']
                    raise StopIteration()
        except StopIteration:
            print('found rc', rc, 'at time', t_snap)
        
        snapdata=loadtxt(snaps[i])
        #print 'done loading snap'
        single_ms_mass = 0.
        binary_ms_mass = 0.
        for j in range(len(snapdata)):
            binflag = snapdata[j,7]
            startype, startype1, startype2 = snapdata[j,14], snapdata[j,17], snapdata[j,18]
            mass, mass1, mass2 = snapdata[j,1], snapdata[j,8], snapdata[j,9]
            r = snapdata[j,2]
            if binflag==0 and startype<2 and r<rc:
                single_ms_mass += mass
                #print single_ms_mass
            if binflag==1 and (startype1<2 or startype2<2) and r<rc:
                binary_ms_mass += (mass1 + mass2)
                #print binary_ms_mass
        single_ms_rhoc = single_ms_mass/rc**3./units[0]['l_pc']
        binary_ms_rhoc = binary_ms_mass/rc**3./units[0]['l_pc']
        #print t_snap, rc, rhoc, single_ms_rhoc, binary_ms_rhoc 
        writefile.write("%f %f %f %f\n" %(t_snap*units[0]['t_myr'], single_ms_rhoc, binary_ms_rhoc, rhoc))

        f_snap.close()

    writefile.close()




def copy_last_snap(basename, DESTBASENAME='/Volumes/MyBook_raid2/for_alison_more/', MIN_T_MYR=9000.):
    """
        Copy snapshots if the time corresponding to the snap is >= MIN_T_MYR.
    """
    import create_plot as cplot
    from glob import glob
    import create_table as ctable

    dirs = cplot.get_the_dirs(basename)
    
    for i in range(len(dirs)):
        temp_dir=dirs[i]
        destdir=dirs[i].split('/')[-1]
        filestring=dirs[i]+'/initial*.gz'
        files=glob(filestring)
        if len(files)>0:
            #print len(files)
            dest = DESTBASENAME+'/'+destdir
            t_string=dirs[i]+'/initial'
            lastsnapno = files[-1].split('snap')[1].split('.')[0]
            last_t_myr = find_t_myr(t_string,lastsnapno)
            #print dirs[i], last_t_myr
            if last_t_myr>=MIN_T_MYR:
                #print dest
                dataout, dataerr = subprocess.Popen([r"mkdir",dest], stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()
                #print dataout, dataerr, dest
                fstring=dest+'/info.txt'
                #print 'came here', fstring
                writefile = open(fstring,'w')
                #print 'came here too'
                writefile.write("#1.t(myr) 2.z 3.m_TO(Msun) 4.snapno\n")
                for i in range(len(files)):
                    snapno = files[i].split('snap')[1].split('.')[0]
                    t_myr = find_t_myr(t_string,snapno)
                    #print snapno, t_myr    
                    if t_myr >= MIN_T_MYR:
                        #subprocess.Popen([r"mkdir",dest])
                        dataout, dataerr = subprocess.Popen([r"cp",files[i],dest],stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()
                        m_guess = find_MS_turnoff(t_myr)
                        m_to = find_MS_TO(t_myr, 0.001, m_guess)
                        writefile.write("%f %f %f %s\n" %(t_myr,0.001,m_to, snapno))
                writefile.close()


def L_to_apparentmag(Lstar, dstar):
    """L: Luminosity of the star in units of Lsun
    d: Distance of the star from the Earth in units of pc
    returns magstar: apparent visual magnitude of the star.  """
    
    magstar = constants.magsun - 2.5 * log10(Lstar * ( constants.AU/(dstar*constants.PC) )**2. )

    return magstar

def find_M_over_L(filename):
    data = loadtxt(filename)
    mtot = sum(data[:,1])
    ltot = 0.
    for i in xrange(len(data)):
        if data[i,7] == 0.:  #single
            ltot += data[i,15]
        elif data[i,7] == 1.:  #binaries
            ltot += data[i,19]+data[i,20]
        else:
            print('not single, not binary. bad outcome')
    return mtot, ltot, mtot/ltot

