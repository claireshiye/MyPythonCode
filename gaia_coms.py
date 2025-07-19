from glob import glob #for looping over files
import gzip # for opening files
import pandas as pd # for manipulating data
from astropy.constants import G #import gravitational constant
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
import os
import dynamics as dyn
import unit_convert as uc
import scripts

from scipy.interpolate import interp1d
from scipy import stats
import scipy.integrate as integrate

yearsc=31557600.
twopi=6.283185307179586
Kconst=9.87*10**-48 ##yr/G^2
Gconst=6.674*10**-8 ##cm3*g-1*s-2
Gconst_sun = 4.30091*10**-3 ##pc*M_sun**-1*(km/s)2
clight=3*10**10 ##cm/s
Msun=2*10**33 ##gram
Rsun=6.957e+10 ##cm
AU=1.496*10**13  ##cm
AU_Rsun=214.93946938362 ##AU to R_sun
PC=3.086*10**18  ##cm
PC_Rsun = 44334448.0068964 ##pc to R_sun

sourcefile = np.genfromtxt('/projects/b1095/syr904/projects/GAIA_COMS/path_allfinished_newruns_maingrid.dat',
             dtype='str')
paths = sourcefile[:,0]

from cosmic.sample.initialbinarytable import InitialBinaryTable
from cosmic.evolve import Evolve

def has_length(x):
    return isinstance(x, (list, tuple, dict, str, np.ndarray, pd.Series))


def evolve_binaries_to_z0(ktype1,ktype2,m1,m2,semimajor,ecc,ospin1,ospin2,B1,B2,mass0_1,mass0_2,epoch1,epoch2,rad0,rad1,lum0,lum1,massc0,massc1,radc0,radc1,menv0,menv1,renv0,renv1,tms0,tms1,bacc0,bacc1,tacc0,tacc1,disruption_time,formation_time,metal):

    num_to_evolve = len(ktype1)
    #num_to_evolve = 1
    print('num_to_evolve in COSMIC', num_to_evolve)

    binary_set = InitialBinaryTable.InitialBinaries(m1=np.ones(num_to_evolve), 
                                                    m2=np.ones(num_to_evolve),
                                                    porb=1e4*np.ones(num_to_evolve), 
                                                    ecc=np.zeros(num_to_evolve), 
                                                    tphysf=np.ones(num_to_evolve), 
                                                    kstar1=np.zeros(num_to_evolve), 
                                                    kstar2=np.zeros(num_to_evolve), 
                                                    metallicity=metal*np.ones(num_to_evolve))

    
    BSEDict = {'xi': 1.0, 'bhflag': 1, 'neta': 0.5, 'windflag': 2, 'wdflag': 1, 'alpha1': 1.0, 
               'pts1': 0.05, 'pts3': 0.02, 'pts2': 0.01, 'epsnov': 0.001, 'hewind': 0.5, 
               'ck': 1000, 'bwind': 0.0, 'lambdaf': 0.0, 'mxns': 2.5, 'beta': -1.0, 
               'tflag': 1, 'acc2': 1.5, 'grflag' : 1, 'remnantflag': 3, 'ceflag': 0, 
               'eddfac': 1.0, 'ifflag': 0, 'bconst': 3000, 'sigma': 265.0, 'gamma': -2.0, 
               'pisn': 45.0, 'natal_kick_array' : [[-100.0,-100.0,-100.0,-100.0,0.0], [-100.0,-100.0,-100.0,-100.0,0.0]], 
               'bhsigmafrac' : 1.0, 'polar_kick_angle' : 90, 
               'qcrit_array' : [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0], 
               'cekickflag' : 2, 'cehestarflag' : 0, 'cemergeflag' : 0, 'ecsn' : 2.25, 
               'ecsn_mlow' : 1.6, 'aic' : 1, 'ussn' : 0, 'sigmadiv' :-20.0, 'qcflag' : 1, 'eddlimflag' : 0, 
               'fprimc_array' : [2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0], 
               'bhspinflag' : 0, 'bhspinmag' : 0.0, 'rejuv_fac' : 1.0, 'rejuvflag' : 0, 'htpmb' : 1, 
               'ST_cr' : 1, 'ST_tide' : 1, 'bdecayfac' : 1, 'rembar_massloss' : 0.5, 'kickflag' : -1, 
               'zsun' : 0.02, 'bhms_coll_flag' : 0, 'don_lim' : -1, 'acc_lim' : -1,'rtmsflag' : 0, 'wd_mass_lim' : 1}
    
    
    ## Evolve for 1 Myr year to set initC
    bpp, bcm, initC, kick_info = Evolve.evolve(initialbinarytable=binary_set, BSEDict=BSEDict)
    #print('initialization')
    #print(initC)
    
    ## Now replace initC with the variables from CMC (36 COSMIC params)
    #InitC.columns:['kstar_1', 'kstar_2', 'mass_1', 'mass_2', 'porb', 'ecc', 'metallicity','binfrac', 'tphysf', 'mass0_1', 'mass0_2', 'rad_1', 'rad_2', 'lum_1','lum_2', 'massc_1', 'massc_2', 'radc_1', 'radc_2', 'menv_1', 'menv_2','renv_1', 'renv_2', 'omega_spin_1', 'omega_spin_2', 'B_1', 'B_2','bacc_1', 'bacc_2', 'tacc_1', 'tacc_2', 'epoch_1', 'epoch_2', 'tms_1','tms_2', 'bhspin_1', 'bhspin_2', 'tphys', 
    #'neta', 'bwind', 'hewind','alpha1', 'lambdaf', 'ceflag', 'tflag', 'ifflag', 'wdflag', 'pisn','rtmsflag', 'bhflag', 'remnantflag', 'grflag', 'bhms_coll_flag','wd_mass_lim', 'cekickflag', 'cemergeflag', 'cehestarflag', 'mxns','pts1', 'pts2', 'pts3', 'ecsn', 'ecsn_mlow', 'aic', 'ussn', 'sigma','sigmadiv', 'bhsigmafrac', 'polar_kick_angle', 'beta', 'xi', 'acc2','epsnov', 'eddfac', 'gamma', 'don_lim', 'acc_lim', 'bdecayfac','bconst', 'ck', 'windflag', 'qcflag', 'eddlimflag', 'dtp', 'randomseed','bhspinflag', 'bhspinmag', 'rejuv_fac', 'rejuvflag', 'htpmb', 'ST_cr','ST_tide', 'rembar_massloss', 'zsun', 'kickflag', 'bin_num','natal_kick_1', 'phi_1', 'theta_1', 'mean_anomaly_1', 'randomseed_1','natal_kick_2', 'phi_2', 'theta_2', 'mean_anomaly_2', 'randomseed_2','qcrit_0', 'qcrit_1', 'qcrit_2', 'qcrit_3', 'qcrit_4', 'qcrit_5','qcrit_6', 'qcrit_7', 'qcrit_8', 'qcrit_9', 'qcrit_10', 'qcrit_11','qcrit_12', 'qcrit_13', 'qcrit_14', 'qcrit_15', 'fprimc_0', 'fprimc_1','fprimc_2', 'fprimc_3', 'fprimc_4', 'fprimc_5', 'fprimc_6', 'fprimc_7','fprimc_8', 'fprimc_9', 'fprimc_10', 'fprimc_11', 'fprimc_12','fprimc_13', 'fprimc_14', 'fprimc_15']
    
    initC["kstar_1"] = ktype1
    initC["kstar_2"] = ktype2
    initC["mass_1"] = m1
    initC["mass_2"] = m2
    initC["ecc"] = ecc
    initC["porb"] = uc.au_to_period(np.array(semimajor), np.array(m1), np.array(m2))
    initC["tphysf"] = 13780.0-np.array(formation_time)-np.array(disruption_time)  #maximum evolution time in Myr
    initC["mass0_1"] = mass0_1
    initC["mass0_2"] = mass0_2
    initC["rad_1"] = rad0
    initC["rad_2"] = rad1
    initC["lum_1"] = lum0
    initC["lum_2"] = lum1
    initC["massc_1"] = massc0
    initC["massc_2"] = massc1
    initC["radc_1"] = radc0
    initC["radc_2"] = radc1
    initC["menv_1"] = menv0
    initC["menv_2"] = menv1
    initC["renv_1"] = renv0
    initC["renv_2"] = renv1
    initC["omega_spin_1"] = ospin1
    initC["omega_spin_2"] = ospin2
    initC["B_1"] = B1
    initC["B_2"] = B2
    initC["bacc_1"] = bacc0
    initC["bacc_2"] = bacc1
    initC["tacc_1"] = tacc0
    initC["tacc_2"] = tacc1
    initC["epoch_1"] = epoch1
    initC["epoch_2"] = epoch2
    initC["tms_1"] = tms0
    initC["tms_2"] = tms1
    initC["tphys"] = disruption_time ##evolution time in Myr
 
    ##Attach two new columns
    initC["sep"] = np.array(semimajor)*215.032
    initC["dtp"] = initC["tphysf"]-initC["tphys"]# -100 ## 2 fixes what i assume is some weird rounding error
    #initC["dtp"] = 0
    ## Not sure this is needed anymore but check with Katie!
    #print(initC[["tphys","tphysf","dtp"]])
    ## evolve to z=0

    #print("before min",initC["dtp"].min(),"max",initC["dtp"].max())
    #print(initC.iloc[np.argmin(initC["dtp"])][["mass_1","tphys","tphysf"]])
    #print(initC.iloc[np.argmax(initC["dtp"])][["mass_1","tphys","tphysf"]])

    #print(semimajor)
    #print('replace properties')
    #print(initC)

    bpp, bcm, initC, kick_info = Evolve.evolve(initialbinarytable=initC)
    #print("after min",initC["dtp"].min(),"max",initC["dtp"].max())
    #print(bcm)

    mass_1 = -100.*np.ones_like(ktype1)
    lumin_1 = -100.*np.ones_like(ktype1)
    teff_1 = -100.*np.ones_like(ktype1)
    kstar_1 = -100.*np.ones_like(ktype1)
    B_1 = -100.*np.ones_like(ktype1)
    ospin_1 = -100.*np.ones_like(ktype1)
    RRLO_1 = -100.*np.ones_like(ktype1)
    mass_2 = -100.*np.ones_like(ktype1)
    lumin_2 = -100.*np.ones_like(ktype1)
    teff_2 = -100.*np.ones_like(ktype1)
    kstar_2 = -100.*np.ones_like(ktype1)
    B_2 = -100.*np.ones_like(ktype1)
    ospin_2 = -100.*np.ones_like(ktype1)
    RRLO_2 = -100.*np.ones_like(ktype1)
    porb = -100.*np.ones_like(ktype1)
    ecc = -100.*np.ones_like(ktype1)
    bin_state =  -100.*np.ones_like(ktype1)

    # If one fails and we didn't write the final snapshot, this will force it to skip them
    idx = []
    for pp in range(num_to_evolve):
        if has_length(bcm.loc[pp]['tphys']):
            idx.append(pp)
            mass_1[pp] = bcm.loc[pp]['mass_1'].iloc[1]
            lumin_1[pp] = bcm.loc[pp]['lum_1'].iloc[1]
            teff_1[pp] = bcm.loc[pp]['teff_1'].iloc[1]
            kstar_1[pp] = bcm.loc[pp]['kstar_1'].iloc[1]
            B_1[pp] = bcm.loc[pp]['B_1'].iloc[1]
            ospin_1[pp] = bcm.loc[pp]['omega_spin_1'].iloc[1]
            RRLO_1[pp] = bcm.loc[pp]['RRLO_1'].iloc[1]

            mass_2[pp] = bcm.loc[pp]['mass_2'].iloc[1]
            lumin_2[pp] = bcm.loc[pp]['lum_2'].iloc[1]
            teff_2[pp] = bcm.loc[pp]['teff_2'].iloc[1]
            kstar_2[pp] = bcm.loc[pp]['kstar_2'].iloc[1]
            B_2[pp] = bcm.loc[pp]['B_2'].iloc[1]
            ospin_2[pp] = bcm.loc[pp]['omega_spin_2'].iloc[1]
            RRLO_2[pp] = bcm.loc[pp]['RRLO_2'].iloc[1]

            porb[pp] = bcm.loc[pp]['porb'].iloc[1]
            ecc[pp] = bcm.loc[pp]['ecc'].iloc[1]
            bin_state[pp] = bcm.loc[pp]['bin_state'].iloc[1]


    #idx = np.arange(len(mass_1))[bcm[1::2]["bin_num"].to_numpy()]    

    #(mass_1[idx],lumin_1[idx],teff_1[idx],kstar_1[idx],B_1[idx],ospin_1[idx],RRLO_1[idx],mass_2[idx],lumin_2[idx],teff_2[idx],kstar_2[idx],B_2[idx],ospin_2[idx],RRLO_2[idx],porb[idx],ecc[idx],bin_state[idx]) = bcm[1::2][["mass_1","lum_1","teff_1","kstar_1","B_1","omega_spin_1","RRLO_1","mass_2","lum_2","teff_2","kstar_2","B_2","omega_spin_2","RRLO_2","porb","ecc","bin_state"]].to_numpy().T

    #print(bcm[1::2][["mass_1","kstar_1", "mass_2","kstar_2","porb","ecc","bin_state"]].to_numpy().T)


    print('COSMIC DONE')

    return mass_1,lumin_1,teff_1,kstar_1,B_1,ospin_1,RRLO_1,mass_2,lumin_2,teff_2,kstar_2,B_2,ospin_2,RRLO_2,porb,ecc,bin_state#, bcm


def conversions(file_path):
    """
    Create dictionary of conversion from code units to physical units.
    
    Parameters: 
    
    file_path: string
        the file path of the cluster to get conversions for
        
    Returns: 
    
    variables: dictionary?
        the conversion factors between code units and physical units
    """
    
    variables = {}
    with open(file_path, 'r') as file:
        for line in file:
            # Strip leading and trailing whitespace
            line = line.strip()
            # Skip comments and empty lines
            if line.startswith('#') or not line:
                continue
            # Split line by '=' into variable name and value
            if '=' in line:
                var_name, var_value = line.split('=', 1)
                var_name = var_name.strip()
                var_value = var_value.strip()
                # Try to convert numeric values to float
                try:
                    var_value = float(var_value)
                except ValueError:
                    pass
                variables[var_name] = var_value
    return variables


def ids_primordial(firstsnap):
    id0_pri = []; id1_pri = []
    #####snap#####
    with gzip.open(firstsnap, 'r') as ffirst:
        next(ffirst)
        next(ffirst)
        for line in ffirst:
            datafirst=line.split()
            if int(datafirst[7])==1:
                id0_pri.append(int(datafirst[10])); id1_pri.append(int(datafirst[11])) 

    return id0_pri, id1_pri


def extract_NSMSs():
    ###############################ESCFILE#####################################
    ### initialize NSMS catalog ###
    NSMS = pd.DataFrame(columns = ['path', 'id0','id1', 'm0[MSUN]', 'm1[MSUN]', 'tb', 'e', 'bin_startype0', 'bin_startype1', 'tcount', 't','B0','B1', 'P0', 'P1'])
    
    #directory = '/projects/b1091/CMC_Grid_March2019/rundir' #define the directory that the clusters are kept in 
    
    #/projects/b1091/CMC_Grid_March2019/rundir/rv0.5/rg2/z0.0002/2e5/initial.esc.dat
    ###
    
    ### Get ejected NS-MS binaries ###
    ### columns for initial.esc.dat files 
    esc_columns = ['tcount', 't', 'm', 'r', 'vr', 'vt', 'r_peri', 'r_apo', 'Rtidal', 'phi_rtidal', 'phi_zero', 'E', 'J', 'id', 'binflag', 'm0[MSUN]', 'm1[MSUN]',\
                   'id0', 'id1', 'a', 'e', 'startype', 'bin_startype0', 'bin_startype1', 'rad0', 'rad1', 'tb', 'lum0', 'lum1', 'massc0', 'massc1', 'radc0', 'radc1', \
                   'menv0', 'menv1', 'renv0', 'renv1', 'tms0', 'tms1', 'dmdt0', 'dmdt1', 'radrol0', 'radrol1', 'ospin0', 'ospin1', 'B0', 'B1', 'formation0', 'formation1', \
                   'bacc0', 'bacc1', 'tacc0', 'tacc1', 'mass0_0', 'mass0_1', 'epoch0', 'epoch1', 'bhspin', 'bhspin1', 'bhspin2', 'ospin', 'B', 'formation']
    
    #for rv_folder in glob(f'{directory}/*'): #for each possible rv
    #    for rg_folder in glob(f'{rv_folder}/*'): #for each possible rg
    #        for z_folder in glob(f'{rg_folder}/*'): # for each possible metallicity 
    #            for n_folder in glob(f'{z_folder}/*'): # for each possible particle numbers
        
    for xx in range(0,144):
    
        #filepath = (f'{n_folder}/initial.esc.dat') #pick the file for the escaped objects
        filepath = paths[xx]+'initial.esc.dat' #pick the file for the escaped objects
    
        print('Opening file', filepath)
    
        t_conv = dyn.conv('t', paths[xx]+'initial.conv.sh')
    
        #with open(f'{filepath}', 'r') as file: # open the latest file for the current cluster in the loop
        with open(filepath, 'r') as file: # open the latest file for the current cluster in the loop
            df = pd.read_csv(file, delim_whitespace=True, skiprows=1, names=esc_columns) #save it to a dataframe
    
            df = df.apply(pd.to_numeric, errors='coerce')
    
            df1 = df[(df["bin_startype0"] == 13) & ((df["bin_startype1"] == 1) | (df["bin_startype1"] == 0))]
            df2 = df[(df["bin_startype1"] == 13) & ((df["bin_startype0"] == 1) | (df["bin_startype0"] == 0))] 
            df = pd.concat([df1,df2], ignore_index=True)
    
            df['P0'] = 2*np.pi*3.154e7/df['ospin0']
            df['P1'] = 2*np.pi*3.154e7/df['ospin1']
    
            #df['t'] = df['t']*t_conv
            
            #df['path'] = f'{n_folder}'
            df['path'] = paths[xx]
    
    
            columns_to_append = ['path','id0','id1', 'm0[MSUN]', 'm1[MSUN]', 'tb', 'e', 'bin_startype0', 'bin_startype1', 'tcount', 't','B0','B1', 'P0', 'P1'] #MSUN, MSUN, days 
    
            # Create new DataFrames with only the selected columns
            df = df[columns_to_append]
    
            # Append rows of df2_selected to df1_selected
            NSMS = pd.concat([NSMS, df], ignore_index=True)
    
        if os.path.isfile(paths[xx]+'initial2.esc.dat') and os.path.getsize(paths[xx]+'initial2.esc.dat') > 0:
            with open(paths[xx]+'initial2.esc.dat', 'r') as file: # open the latest file for the current cluster in the loop
                df = pd.read_csv(file, delim_whitespace=True, skiprows=1, names=esc_columns) #save it to a dataframe
        
                df1 = df[(df["bin_startype0"] == 13) & ((df["bin_startype1"] == 1) | (df["bin_startype1"] == 0))]
                df2 = df[(df["bin_startype1"] == 13) & ((df["bin_startype0"] == 1) | (df["bin_startype0"] == 0))] 
                df = pd.concat([df1,df2], ignore_index=True)
        
                df['P0'] = 2*np.pi*3.154e7/df['ospin0']
                df['P1'] = 2*np.pi*3.154e7/df['ospin1']
        
                #df['t'] = df['t']*t_conv
                
                #df['path'] = f'{n_folder}'
                df['path'] = paths[xx]
        
        
                columns_to_append = ['path','id0','id1', 'm0[MSUN]', 'm1[MSUN]', 'tb', 'e', 'bin_startype0', 'bin_startype1', 'tcount', 't','B0','B1', 'P0', 'P1'] #MSUN, MSUN, days 
        
                # Create new DataFrames with only the selected columns
                df = df[columns_to_append]
        
                # Append rows of df2_selected to df1_selected
                NSMS = pd.concat([NSMS, df], ignore_index=True)
    
    
    ### Save data to a csv file
    savepath = '/projects/b1095/syr904/projects/GAIA_COMS' #specify path for catalog to be saved to 
    NSMS.to_csv(savepath+'/NSMS_catalog_escaped.csv')
    
    print('ESCFILE done')
    ###
    
    ###############################PULSARFILE#####################################
    ### initialize NSMS catalog ###
    #NSMS = pd.DataFrame(columns = ['path', 'id0','id1', 'm0[MSUN]', 'm1[MSUN]', 'tb', 'e', 'bin_startype0', 'bin_startype1', 'tcount', 't','B0','B1', 'P0', 'P1'])
    #        
    #### Get retained NS-MS binaries ###
    #morepulsar_columns = ['tcount', 'TotalTime', 'binflag', 'id0', 'id1', 'm0[MSUN]', 'm1[MSUN]', 'B0', 'B1', 'P0', 'P1', 'bin_startype0', 'bin_startype1', \
    #    'a', 'e', 'radrol0', 'radrol1', 'dmdt0', 'dmdt1', 'r', 'vr', 'vt', 'bacc0', 'bacc1', 'tacc0', 'tacc1']
    #
    ##for rv_folder in glob(f'{directory}/*'): #for each possible rv
    ##    for rg_folder in glob(f'{rv_folder}/*'): #for each possible rg
    ##        for z_folder in glob(f'{rg_folder}/*'): # for each possible metallicity 
    ##            for n_folder in glob(f'{z_folder}/*'): # for each possible particle numbers
    #for yy in range(0,144):
    #    #filepath = (f'{n_folder}/initial.morepulsars.dat') 
    #    filepath = paths[yy]+'initial.morepulsars.dat'
    #
    #    print("Opening file", filepath)
    #
    #    t_conv = dyn.conv('t', paths[yy]+'initial.conv.sh')
    #
    #    #with open(f'{filepath}', 'r') as file: # open the latest file for the current cluster in the loop
    #    with open(filepath, 'r') as file: # open the latest file for the current cluster in the loop
    #        df = pd.read_csv(file, delim_whitespace=True, skiprows=1, names=morepulsar_columns) #save it to a dataframe
    #        #print(df[df['binflag']==1]['bin_startype1'])
    #
    #        df1 = df[(df["bin_startype0"] == 13) & ((df["bin_startype1"] == 1) | (df["bin_startype1"] == 0))]
    #        df2 = df[(df["bin_startype1"] == 13) & ((df["bin_startype0"] == 1) | (df["bin_startype0"] == 0))] 
    #        df = pd.concat([df1,df2], ignore_index=True)
    #
    #        Gunit = G.to(u.AU**3*u.Msun**-1*u.day**-2)
    #
    #        df['tb'] = np.sqrt(4*np.pi**2 / (Gunit*(df['m0[MSUN]']+df['m1[MSUN]'])) * df['a']**3)
    #        df['t'] = df['TotalTime']*t_conv
    #        
    #        #df['path'] = f'{n_folder}'
    #        df['path'] = paths[yy]
    #
    #        columns_to_append = ['path', 'id0','id1','m0[MSUN]', 'm1[MSUN]', 'tb', 'e', 'bin_startype0', 'bin_startype1', 'tcount', 't', 'B0','B1', 'P0', 'P1'] 
    #
    #        # Create new DataFrames with only the selected columns
    #        df = df[columns_to_append]
    #
    #        # Append rows of df2_selected to df1_selected
    #        NSMS = pd.concat([NSMS, df], ignore_index=True)
    #
    #    if os.path.isfile(paths[yy]+'initial2.morepulsars.dat') and os.path.getsize(paths[yy]+'initial2.morepulsars.dat') > 0:
    #        with open(paths[yy]+'initial2.morepulsars.dat', 'r') as file: # open the latest file for the current cluster in the loop
    #            df = pd.read_csv(file, delim_whitespace=True, skiprows=1, names=morepulsar_columns) #save it to a dataframe
    #            #print(df[df['binflag']==1]['bin_startype1'])
    #    
    #            df1 = df[(df["bin_startype0"] == 13) & ((df["bin_startype1"] == 1) | (df["bin_startype1"] == 0))]
    #            df2 = df[(df["bin_startype1"] == 13) & ((df["bin_startype0"] == 1) | (df["bin_startype0"] == 0))] 
    #            df = pd.concat([df1,df2], ignore_index=True)
    #    
    #            Gunit = G.to(u.AU**3*u.Msun**-1*u.day**-2)
    #    
    #            df['tb'] = np.sqrt(4*np.pi**2 / (Gunit*(df['m0[MSUN]']+df['m1[MSUN]'])) * df['a']**3)
    #            df['t'] = df['TotalTime']*t_conv
    #            
    #            #df['path'] = f'{n_folder}'
    #            df['path'] = paths[yy]
    #    
    #            columns_to_append = ['path', 'id0','id1','m0[MSUN]', 'm1[MSUN]', 'tb', 'e', 'bin_startype0', 'bin_startype1', 'tcount', 't', 'B0','B1', 'P0', 'P1'] 
    #    
    #            # Create new DataFrames with only the selected columns
    #            df = df[columns_to_append]
    #    
    #            # Append rows of df2_selected to df1_selected
    #            NSMS = pd.concat([NSMS, df], ignore_index=True)
    #
    #
    ####
    #
    #### Drop double counted clusters ###
    #NSMS = NSMS.drop_duplicates(subset=['id0', 'id1'],keep='last')
    #
    #### Save data to a csv file
    #savepath = '/projects/b1095/syr904/projects/GAIA_COMS' #specify path for catalog to be saved to 
    #NSMS.to_csv(savepath+'/NSMS_catalog_incluster.csv')
    #
    #print('PULSARFILE done')


def extract_escaped_COMSs(comsfile):
    data_dict = {"kstar_1":[],"kstar_2":[],"mass_1":[],"mass_2":[],"semimajor":[],"ecc":[],"ospin_1":[],"ospin_2":[],"B_1":[],"B_2":[],"mass0_1":[],"mass0_2":[],"epoch_1":[],"epoch_2":[],"rad_1":[],"rad_2":[],"lum_1":[],"lum_2":[],"massc_1":[],"massc_2":[],"radc_1":[],"radc_2":[],"menv_1":[],"menv_2":[],"renv_1":[],"renv_2":[],"tms_1":[],"tms_2":[],"bacc_1":[],"bacc_2":[], "tacc_1":[],"tacc_2":[],"disruption_time":[],"formation_time":[],"id1":[], "id2":[], "Clus_no":[], "CO_flag":[], "Primodial_flag":[], "Mass(CMC)":[], "rv(CMC)":[], "rg(CMC)":[], "z(CMC)":[], "T[Myr](CMC)":[]}
    #fwrite = open('/projects/b1095/syr904/projects/GAIA_COMS/coms_properties_presentday_'+mapflag+'_'+vfile+'.dat', 'w+')
    #fwrite.write('#1.Clus_No 2.Tdisrupt[Myr](galaxy) 3.ID0 4.ID1 5.M0[Msun] 6.M1[Msun] 7.K0 8.K1 9.Porb[Days] 10.ECC 11.Lum[Lsun] 12.Teff[K] 13.CO_flag 14.Primodial_flag 15.B[G] 16.Pspin[sec] 17.Mass(CMC) 18.rv(CMC) 19.rg(CMC) 20.z(CMC) 21.T[Myr](CMC) 22.Tform[Myr](galaxy)\n')
     
    data_coms = np.genfromtxt(comsfile, dtype='str', usecols=(0,1,16,17,18,19,20,21))
    clus_no = np.unique(data_coms[:,0]); tdisrupt = data_coms[:,1].astype(float)
    clus_n = data_coms[:,2]; clus_rv = data_coms[:,3]; clus_rg = data_coms[:,4]; clus_z = data_coms[:,5]
    clus_t = data_coms[:,6].astype(float)
    tform = data_coms[:,7].astype(float)

    print('Start extracting escaped binaries')

    ###############################ESCFILE#####################################
    directory = '/projects/b1091/CMC_Grid_March2019/rundir' #define the directory that the clusters are kept in
    #/projects/b1091/CMC_Grid_March2019/rundir/rv0.5/rg2/z0.0002/2e5/initial.esc.dat 
        
    ### Get ejected binaries ###
    ### columns for initial.esc.dat files 
    #esc_columns = ['tcount', 't', 'm', 'r', 'vr', 'vt', 'r_peri', 'r_apo', 'Rtidal', 'phi_rtidal', 'phi_zero', 'E', 'J', 'id', 'binflag', 'm0[MSUN]', 'm1[MSUN]',\
    #               'id0', 'id1', 'a', 'e', 'startype', 'bin_startype0', 'bin_startype1', 'rad0', 'rad1', 'tb', 'lum0', 'lum1', 'massc0', 'massc1', 'radc0', 'radc1', \
    #               'menv0', 'menv1', 'renv0', 'renv1', 'tms0', 'tms1', 'dmdt0', 'dmdt1', 'radrol0', 'radrol1', 'ospin0', 'ospin1', 'B0', 'B1', 'formation0', 'formation1', \
    #               'bacc0', 'bacc1', 'tacc0', 'tacc1', 'mass0_0', 'mass0_1', 'epoch0', 'epoch1', 'bhspin', 'bhspin1', 'bhspin2', 'ospin', 'B', 'formation']


    xx = 0
    while xx < len(clus_no):
        #if clus_no[xx]!='10021':
        #    xx+=1
        #    continue

        #check = 0
        for yy in range(len(clus_n)):

            if data_coms[:,0][yy]==clus_no[xx]:
                #check=1

                filedir = directory+'/rv'+clus_rv[yy]+'/rg'+clus_rg[yy]+'/z'+clus_z[yy]+'/'+clus_n[yy]
                filepath = filedir+'/initial.esc.dat' #pick the file for the escaped objects
    
                #print('Opening file', filepath)
    
                t_conv = dyn.conv('t', filedir+'/initial.conv.sh')

                snapfirst = filedir+'/initial.snap0000.dat.gz'
                id0_snap0, id1_snap0 = ids_primordial(snapfirst)
    
                with open(filepath, 'r') as file: # open the latest file for the current cluster in the loop
                    next(file)
                    for line in file:
                        data = line.split()

                        if int(data[14])==1 and float(data[1])*t_conv<=clus_t[yy]:
                            if ((float(data[15])>1. or float(data[16])>1.) and (int(data[22])<=1 or int(data[23])<=1)):
                                primordial = 0
                                for hh in range(len(id0_snap0)):
                                    if (int(data[17]) == int(id0_snap0[hh]) and int(data[18]) == int(id1_snap0[hh])) or (int(data[17]) == int(id1_snap0[hh]) and int(data[18]) == int(id0_snap0[hh])):
                                        primordial = 1
                                        break
                
                                #tage =  float(data[1])*t_conv - float(data[56]) ##this is the epoch; age = tphys - epoch
                                #if tage < 0 :
                                #    print(int(data[17]), int(data[18]), filepath)


                                data_dict['mass_1'].append(float(data[15]))
                                data_dict['mass_2'].append(float(data[16]))
                                data_dict['semimajor'].append(float(data[19]))  ##in AU
                                data_dict['ecc'].append(float(data[20])) 
                                data_dict['kstar_1'].append(int(data[22]))
                                data_dict['kstar_2'].append(int(data[23])) 
                                data_dict['rad_1'].append(float(data[24])) 
                                data_dict['rad_2'].append(float(data[25])) 
                                data_dict['lum_1'].append(float(data[27])) 
                                data_dict['lum_2'].append(float(data[28]))
                                data_dict['massc_1'].append(float(data[29])) 
                                data_dict['massc_2'].append(float(data[30]))  
                                data_dict['radc_1'].append(float(data[31])) 
                                data_dict['radc_2'].append(float(data[32])) 
                                data_dict['menv_1'].append(float(data[33])) 
                                data_dict['menv_2'].append(float(data[34]))
                                data_dict['renv_1'].append(float(data[35])) 
                                data_dict['renv_2'].append(float(data[36])) 
                                data_dict['tms_1'].append(float(data[37])) 
                                data_dict['tms_2'].append(float(data[38]))
                                #data_dict['RRLO_1'].append(float(data[41])) 
                                #data_dict['RRLO_2'].append(float(data[42]))
                                data_dict['ospin_1'].append(float(data[43])) 
                                data_dict['ospin_2'].append(float(data[44])) 
                                data_dict['B_1'].append(float(data[45])) 
                                data_dict['B_2'].append(float(data[46])) 
                                data_dict['bacc_1'].append(float(data[49])) 
                                data_dict['bacc_2'].append(float(data[50])) 
                                data_dict['tacc_1'].append(float(data[51])) 
                                data_dict['tacc_2'].append(float(data[52]))
                                data_dict['mass0_1'].append(float(data[53])) 
                                data_dict['mass0_2'].append(float(data[54]))  
                                data_dict['epoch_1'].append(float(data[55])) 
                                data_dict['epoch_2'].append(float(data[56]))
                                
                                data_dict["disruption_time"].append(tdisrupt[yy])
                                data_dict["formation_time"].append(tform[yy])
                                data_dict["id1"].append(int(data[17]))
                                data_dict["id2"] .append(int(data[18]))
                                data_dict["Clus_no"].append(clus_no[xx])
                                data_dict["Primodial_flag"].append(primordial)
                                data_dict["Mass(CMC)"].append(clus_n[yy])
                                data_dict["rv(CMC)"].append(clus_rv[yy])
                                data_dict["rg(CMC)"].append(clus_rg[yy]) 
                                data_dict["z(CMC)"].append(clus_z[yy])
                                data_dict["T[Myr](CMC)"].append(float(data[1])*t_conv)


                if os.path.isfile(filedir+'/initial2.esc.dat') and os.path.getsize(filedir+'/initial2.esc.dat') > 0:
                    with open(filedir+'/initial2.esc.dat', 'r') as file: # open the latest file for the current cluster in the loop
                        next(file)
                        for line in file:
                            data = line.split()
        
                            if int(data[14])==1 and float(data[1])*t_conv<=clus_t[yy]:
                                if ((float(data[15])>1. or float(data[16])>1.) and (int(data[22])<=1 or int(data[23])<=1)):
                                    primordial = 0
                                    for hh in range(len(id0_snap0)):
                                        if (int(data[17]) == int(id0_snap0[hh]) and int(data[18]) == int(id1_snap0[hh])) or (int(data[17]) == int(id1_snap0[hh]) and int(data[18]) == int(id0_snap0[hh])):
                                            primordial = 1
                                            break
                        
                                    #tage =  float(data[1])*t_conv - float(data[56]) ##this is the epoch; age = tphys - epoch
                                    #if tage < 0 :
                                    #    print(int(data[17]), int(data[18]), filepath)

                                    data_dict['mass_1'].append(float(data[15]))
                                    data_dict['mass_2'].append(float(data[16]))
                                    data_dict['semimajor'].append(float(data[19]))  ##in AU
                                    data_dict['ecc'].append(float(data[20])) 
                                    data_dict['kstar_1'].append(int(data[22]))
                                    data_dict['kstar_2'].append(int(data[23])) 
                                    data_dict['rad_1'].append(float(data[24])) 
                                    data_dict['rad_2'].append(float(data[25])) 
                                    data_dict['lum_1'].append(float(data[27])) 
                                    data_dict['lum_2'].append(float(data[28]))
                                    data_dict['massc_1'].append(float(data[29])) 
                                    data_dict['massc_2'].append(float(data[30]))  
                                    data_dict['radc_1'].append(float(data[31])) 
                                    data_dict['radc_2'].append(float(data[32])) 
                                    data_dict['menv_1'].append(float(data[33])) 
                                    data_dict['menv_2'].append(float(data[34]))
                                    data_dict['renv_1'].append(float(data[35])) 
                                    data_dict['renv_2'].append(float(data[36])) 
                                    data_dict['tms_1'].append(float(data[37])) 
                                    data_dict['tms_2'].append(float(data[38]))
                                    #data_dict['RRLO_1'].append(float(data[41])) 
                                    #data_dict['RRLO_2'].append(float(data[42]))
                                    data_dict['ospin_1'].append(float(data[43])) 
                                    data_dict['ospin_2'].append(float(data[44])) 
                                    data_dict['B_1'].append(float(data[45])) 
                                    data_dict['B_2'].append(float(data[46])) 
                                    data_dict['bacc_1'].append(float(data[49])) 
                                    data_dict['bacc_2'].append(float(data[50])) 
                                    data_dict['tacc_1'].append(float(data[51])) 
                                    data_dict['tacc_2'].append(float(data[52]))
                                    data_dict['mass0_1'].append(float(data[53])) 
                                    data_dict['mass0_2'].append(float(data[54]))  
                                    data_dict['epoch_1'].append(float(data[55])) 
                                    data_dict['epoch_2'].append(float(data[56]))
                                    
                                    data_dict["disruption_time"].append(tdisrupt[yy])
                                    data_dict["formation_time"].append(tform[yy])
                                    data_dict["id1"].append(int(data[17]))
                                    data_dict["id2"] .append(int(data[18]))
                                    data_dict["Clus_no"].append(clus_no[xx])
                                    data_dict["Primodial_flag"].append(primordial)
                                    data_dict["Mass(CMC)"].append(clus_n[yy])
                                    data_dict["rv(CMC)"].append(clus_rv[yy])
                                    data_dict["rg(CMC)"].append(clus_rg[yy]) 
                                    data_dict["z(CMC)"].append(clus_z[yy])
                                    data_dict["T[Myr](CMC)"].append(float(data[1])*t_conv)
        

            
                break
        
        #if clus_no[xx]=='10021':
        #    break
        xx+=1
        #break


    num_to_evolve = len(data_dict['mass_1'])
    print('num_to_evolve', num_to_evolve)
    print('Extract ESC done')

    return data_dict
    


def print_N_NSMS_psrfile(pathlist, start, end, readflag):
    if readflag == 1:
        sourcedir=np.genfromtxt(pathlist, dtype='str')
        status=sourcedir[:,1]; 
        sourcedir=sourcedir[:,0]
    else:
        sourcedir=pathlist
        status = [1]
    
    savepath = '/projects/b1095/syr904/projects/GAIA_COMS/'
    fh = open(savepath+'NSMS_allruns.dat', 'w+') ###changed file name
    fh.write('#1.Model 2.Totaltime[Myr] 3.N_NSMS\n') ### changed to only total time and NS-MS count
    for i in range(start, end):
        filepath = sourcedir[i]

        t_conv = dyn.conv('t', filepath+'initial.conv.sh')

        psrfiles = glob(filepath+'*.morepulsars.dat')
        allt = []
        for xx in range(len(psrfiles)):
            with open(psrfiles[xx], 'r') as fpsr:
                next(fpsr)
                for line in fpsr:
                    data = line.split()
                    allt.append(float(data[1]))
                    break
        #print(allt)
        if len(allt)<1:
            continue
        psrfiles, allt = (list(t) for t in zip(*sorted(zip(psrfiles, allt))))
        print(psrfiles, allt)


        print(sourcedir[i])
        for j in range(len(psrfiles)):
            N_NSMS=0 ###deleted other NS systems
            t_old = allt[j]
            with open(psrfiles[j], 'r') as fpsr:
                next(fpsr)
                for line in fpsr:
                    datapsr=line.split()
                    #print(datapsr)
                    t_curr = float(datapsr[1])
                    if t_curr != t_old:
                        T = t_old*t_conv
                        fh.write('%s %f %d\n'%(filepath,T,N_NSMS)) ###deleted other NS systems
                        N_NSMS=0 ###deleted other NS systems
                        t_old = t_curr

                    if int(datapsr[2])==1:
                        if int(datapsr[11])==13:
                            if int(datapsr[12])<2: N_NSMS+=1

                        if int(datapsr[12])==13:
                            if int(datapsr[11])<2: N_NSMS+=1

            T = t_old*t_conv              
            fh.write('%s %f %d \n'%(filepath, T, N_NSMS))
            #print(j)

    fh.close()


##Extract and grouping the number of CO--MS from the catalog models
def extract_n_coms(coflag):
    pathlist = np.genfromtxt('/projects/b1095/syr904/projects/GAIA_COMS/path_allfinished_newruns_maingrid.dat', dtype=str)
    paths = pathlist[:,0]; status = pathlist[:,1]

    bin_size = 400

    n_model_mass = [0,0,0,0]; n_model_rv = [0,0,0,0]; n_model_z = [0,0,0]; n_model_rg = [0,0,0]
    for ii in range(len(paths)):

        ##Initial Conditions
        s=paths[ii].split('/')
        n_star=float(s[-2])
        z=float(s[-3][1:])
        rg=int(s[-4][2:])
        rv=float(s[-5][2:])

        #if rg>2:
        #    continue
        
        if n_star==200000.:# and status[ii]=='1': 
            n_model_mass[0]+=1
        if n_star==400000.:# and status[ii]=='1': 
            n_model_mass[1]+=1
        if n_star==800000.:# and status[ii]=='1': 
            n_model_mass[2]+=1
        if n_star==1600000.:# and status[ii]=='1': 
            n_model_mass[3]+=1
            
        if rv==4.:# and status[ii]=='1': 
            n_model_rv[0]+=1
        if rv==2.:# and status[ii]=='1': 
            n_model_rv[1]+=1
        if rv==1.:# and status[ii]=='1': 
            n_model_rv[2]+=1
        if rv==0.5:# and status[ii]=='1': 
            n_model_rv[3]+=1
            
            
        if z==0.0002:# and status[ii]=='1': 
            n_model_z[0]+=1
        if z==0.002:# and status[ii]=='1': 
            n_model_z[1]+=1
        if z==0.02:# and status[ii]=='1': 
            n_model_z[2]+=1
            
            
        if rg==2:# and status[ii]=='1': 
            n_model_rg[0]+=1
        if rg==8:# and status[ii]=='1': 
            n_model_rg[1]+=1
        if rg==20:# and status[ii]=='1': 
            n_model_rg[2]+=1
        

    ##Grouping models        
    n_coms_rv = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_coms_rv_average = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_coms_rv_average_std = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_coms_rv_median = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]

    n_coms_mass = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_coms_mass_average = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_coms_mass_average_std = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_coms_mass_median = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]

    n_coms_z = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_coms_z_average = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_coms_z_average_std = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_coms_z_median = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]


    n_coms_rg = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_coms_rg_average = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_coms_rg_average_std = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_coms_rg_median = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    
    ncoms_scatter_n2e5 = [[] for _ in range(bin_size)]
    ncoms_scatter_n4e5 = [[] for _ in range(bin_size)]
    ncoms_scatter_n8e5 = [[] for _ in range(bin_size)]
    ncoms_scatter_n16e5 = [[] for _ in range(bin_size)]

    ncoms_scatter_rv4 = [[] for _ in range(bin_size)]
    ncoms_scatter_rv2 = [[] for _ in range(bin_size)]
    ncoms_scatter_rv1 = [[] for _ in range(bin_size)]
    ncoms_scatter_rv05 = [[] for _ in range(bin_size)]

    ncoms_scatter_z00002 = [[] for _ in range(bin_size)]
    ncoms_scatter_z0002 = [[] for _ in range(bin_size)]
    ncoms_scatter_z002= [[] for _ in range(bin_size)]

    ncoms_scatter_rg2 = [[] for _ in range(bin_size)]
    ncoms_scatter_rg8 = [[] for _ in range(bin_size)]
    ncoms_scatter_rg20= [[] for _ in range(bin_size)]


    t_all = np.linspace(0, 13000., bin_size+1)
    for kk in range(0, 144):
        print(paths[kk])

        ##Initial Conditions
        s=paths[kk].split('/')
        n_star=float(s[-2])
        z=float(s[-3][1:])
        rg=int(s[-4][2:])
        rv=float(s[-5][2:])

        #if rg>2:
        #    continue
        
        t_conv = dyn.conv('t', paths[kk]+'initial.conv.sh')
        if coflag=='BH':
            dataco = np.genfromtxt(paths[kk]+'initial.bh.dat')
            times = np.array(dataco[:,1])*t_conv
            n_coms = dataco[:,10]; n_copms = dataco[:,11]

        if coflag=='NS':
            dataco = np.genfromtxt(paths[kk]+'initial.ns.dat')
            times = np.array(dataco[:,0])*t_conv
            n_coms = dataco[:,10]; n_copms = dataco[:,11]
        
        ##Interpolate the number of CO data
        f = interp1d(times, n_coms, kind='nearest')
        t_interpld = np.linspace(0, np.max(times), 3*bin_size)
        n_coms_new = f(t_interpld)
        #print(n_coms_new)
    
        n_mass = [[],[],[],[]]; n_rv = [[],[],[],[]]; n_z = [[],[],[]]; n_rg = [[],[],[]]
        for jj in range(len(t_all)-1):
            #print(jj)
            n_mass_temp = [0,0,0,0]; n_rv_temp = [0,0,0,0]; n_z_temp = [0,0,0]; n_rg_temp = [0,0,0]
            count_mass = [0,0,0,0]; count_rv = [0,0,0,0]; count_z = [0,0,0]; count_rg = [0,0,0]
            
            ##Group by initial mass
            if n_star==200000.:# and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_mass_temp[0]+=n_coms_new[i]
                        count_mass[0]+=1  ##multiple time steps may belong to the same bin
        
            if n_star==400000.:# and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_mass_temp[1]+=n_coms_new[i]
                        count_mass[1]+=1
        
            if n_star==800000.:# and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_mass_temp[2]+=n_coms_new[i]
                        count_mass[2]+=1
        
            if n_star==1600000.:# and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_mass_temp[3]+=n_coms_new[i]
                        count_mass[3]+=1
            
            ##Group by initial rv   
            if rv==4.:# and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_rv_temp[0]+=n_coms_new[i]
                        count_rv[0]+=1
        
            if rv==2.:# and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_rv_temp[1]+=n_coms_new[i]
                        count_rv[1]+=1
        
            if rv==1.:# and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_rv_temp[2]+=n_coms_new[i]
                        count_rv[2]+=1
        
            if rv==0.5:# and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_rv_temp[3]+=n_coms_new[i]
                        count_rv[3]+=1
                    
        
            ##Group by metallicity
            if z==0.0002:# and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_z_temp[0]+=n_coms_new[i]
                        count_z[0]+=1
        
            if z==0.002:# and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_z_temp[1]+=n_coms_new[i]
                        count_z[1]+=1
        
            if z==0.02:# and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_z_temp[2]+=n_coms_new[i]
                        count_z[2]+=1
                    
                    
            ##Group by galactocentric distance
            if rg==2:# and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_rg_temp[0]+=n_coms_new[i]
                        count_rg[0]+=1
        
            if rg==8:# and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_rg_temp[1]+=n_coms_new[i]
                        count_rg[1]+=1
        
            if rg==20:# and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_rg_temp[2]+=n_coms_new[i]
                        count_rg[2]+=1
    
            #print(count_rv[0])
        
            for x in range(4):
                if count_rv[x]!=0:
                    n_rv_temp[x] = n_rv_temp[x]/count_rv[x]
                if count_mass[x]!=0:
                    n_mass_temp[x] = n_mass_temp[x]/count_mass[x]
                    
                n_rv[x].append(n_rv_temp[x])
                n_mass[x].append(n_mass_temp[x])
        
            for x in range(3):
                if count_z[x]!=0:
                    n_z_temp[x] = n_z_temp[x]/count_z[x]
                
                n_z[x].append(n_z_temp[x])     
                
                if count_rg[x]!=0:
                    n_rg_temp[x] = n_rg_temp[x]/count_rg[x]
                
                n_rg[x].append(n_rg_temp[x])
                            
            
        for y in range(4):
            n_coms_rv[y] = n_coms_rv[y]+np.array(n_rv[y])
            n_coms_mass[y] = n_coms_mass[y]+np.array(n_mass[y])
            n_coms_rv_average[y] = n_coms_rv_average[y] + np.array(n_rv[y])/n_model_rv[y]
            n_coms_mass_average[y] = n_coms_mass_average[y] + np.array(n_mass[y])/n_model_mass[y]
            
        for y in range(3):
            n_coms_z[y] = n_coms_z[y]+np.array(n_z[y])
            n_coms_z_average[y] = n_coms_z_average[y] + np.array(n_z[y])/n_model_z[y]
            n_coms_rg[y] = n_coms_rg[y]+np.array(n_rg[y])
            n_coms_rg_average[y] = n_coms_rg_average[y] + np.array(n_rg[y])/n_model_rg[y]


        ##Group by initial mass
        if n_star==200000.:# and status[kk]=='1':
            #print(len(n_mass[0]))
            ncoms_scatter_n2e5 = np.hstack((ncoms_scatter_n2e5, np.split(np.array(n_mass[0]),len(n_mass[0]))))
            #print(ncoms_scatter_n2e5)
        
        if n_star==400000.:# and status[kk]=='1':
            ncoms_scatter_n4e5 = np.hstack((ncoms_scatter_n4e5, np.split(np.array(n_mass[1]), len(n_mass[1]))))
        
        if n_star==800000.:# and status[kk]=='1':
            ncoms_scatter_n8e5 = np.hstack((ncoms_scatter_n8e5, np.split(np.array(n_mass[2]), len(n_mass[2]))))

        if n_star==1600000.:# and status[kk]=='1':
            ncoms_scatter_n16e5 = np.hstack((ncoms_scatter_n16e5, np.split(np.array(n_mass[3]), len(n_mass[3]))))
            

        ##Group by initial rv   
        if rv==4.:# and status[kk]=='1':
            ncoms_scatter_rv4 = np.hstack((ncoms_scatter_rv4, np.split(np.array(n_rv[0]), len(n_rv[0]))))
        
        if rv==2.:# and status[kk]=='1':
            ncoms_scatter_rv2 = np.hstack((ncoms_scatter_rv2, np.split(np.array(n_rv[1]), len(n_rv[1]))))
        
        if rv==1.:# and status[kk]=='1':
            ncoms_scatter_rv1 = np.hstack((ncoms_scatter_rv1, np.split(np.array(n_rv[2]), len(n_rv[2]))))
        
        if rv==0.5:# and status[kk]=='1':
            ncoms_scatter_rv05 = np.hstack((ncoms_scatter_rv05, np.split(np.array(n_rv[3]), len(n_rv[3]))))
                    
    
        ##Group by metallicity
        if z==0.0002:# and status[kk]=='1':
            ncoms_scatter_z00002 = np.hstack((ncoms_scatter_z00002, np.split(np.array(n_z[0]), len(n_z[0]))))
        
        if z==0.002:# and status[kk]=='1':
            ncoms_scatter_z0002 = np.hstack((ncoms_scatter_z0002, np.split(np.array(n_z[1]), len(n_z[1]))))
        
        if z==0.02:# and status[kk]=='1':
            ncoms_scatter_z002 = np.hstack((ncoms_scatter_z002, np.split(np.array(n_z[2]), len(n_z[2]))))
                    
                    
        ##Group by galactocentric distance
        if rg==2:# and status[kk]=='1':
            ncoms_scatter_rg2 = np.hstack((ncoms_scatter_rg2, np.split(np.array(n_rg[0]), len(n_rg[0]))))
        
        if rg==8:# and status[kk]=='1':
            ncoms_scatter_rg8 = np.hstack((ncoms_scatter_rg8, np.split(np.array(n_rg[1]), len(n_rg[1]))))
        
        if rg==20:# and status[kk]=='1':
            ncoms_scatter_rg20 = np.hstack((ncoms_scatter_rg20, np.split(np.array(n_rg[2]), len(n_rg[2]))))

    for ii in range(4):
        for xx in range(bin_size):
            if ii == 0:
                n_coms_mass_average_std[ii][xx]+=np.std(ncoms_scatter_n2e5[xx])
                n_coms_rv_average_std[ii][xx]+=np.std(ncoms_scatter_rv4[xx])

                n_coms_mass_median[ii][xx]+=np.median(ncoms_scatter_n2e5[xx])
                n_coms_rv_median[ii][xx]+=np.median(ncoms_scatter_rv4[xx])
            if ii == 1:
                n_coms_mass_average_std[ii][xx]+=np.std(ncoms_scatter_n4e5[xx])
                n_coms_rv_average_std[ii][xx]+=np.std(ncoms_scatter_rv2[xx])

                n_coms_mass_median[ii][xx]+=np.median(ncoms_scatter_n4e5[xx])
                n_coms_rv_median[ii][xx]+=np.median(ncoms_scatter_rv2[xx])
            if ii == 2:
                n_coms_mass_average_std[ii][xx]+=np.std(ncoms_scatter_n8e5[xx])
                n_coms_rv_average_std[ii][xx]+=np.std(ncoms_scatter_rv1[xx])

                n_coms_mass_median[ii][xx]+=np.median(ncoms_scatter_n8e5[xx])
                n_coms_rv_median[ii][xx]+=np.median(ncoms_scatter_rv1[xx])
            if ii == 3:
                n_coms_mass_average_std[ii][xx]+=np.std(ncoms_scatter_n16e5[xx])
                n_coms_rv_average_std[ii][xx]+=np.std(ncoms_scatter_rv05[xx])

                n_coms_mass_median[ii][xx]+=np.median(ncoms_scatter_n16e5[xx])
                n_coms_rv_median[ii][xx]+=np.median(ncoms_scatter_rv05[xx])

    for ii in range(3):
        for xx in range(bin_size):
            if ii == 0:
                n_coms_z_average_std[ii][xx]+=np.std(ncoms_scatter_z00002[xx])
                n_coms_rg_average_std[ii][xx]+=np.std(ncoms_scatter_rg2[xx])

                n_coms_z_median[ii][xx]+=np.median(ncoms_scatter_z00002[xx])
                n_coms_rg_median[ii][xx]+=np.median(ncoms_scatter_rg2[xx])
            if ii == 1:
                n_coms_z_average_std[ii][xx]+=np.std(ncoms_scatter_z0002[xx])
                n_coms_rg_average_std[ii][xx]+=np.std(ncoms_scatter_rg8[xx])

                n_coms_z_median[ii][xx]+=np.median(ncoms_scatter_z0002[xx])
                n_coms_rg_median[ii][xx]+=np.median(ncoms_scatter_rg8[xx])
            if ii == 2:
                n_coms_z_average_std[ii][xx]+=np.std(ncoms_scatter_z002[xx])
                n_coms_rg_average_std[ii][xx]+=np.std(ncoms_scatter_rg20[xx])

                n_coms_z_median[ii][xx]+=np.median(ncoms_scatter_z002[xx])
                n_coms_rg_median[ii][xx]+=np.median(ncoms_scatter_rg20[xx])
    
    print(n_coms_mass_average_std[0])

    for z in range(4):
        n_coms_rv[z] = np.insert(n_coms_rv[z], 0, 0.); n_coms_rv_average[z] = np.insert(n_coms_rv_average[z], 0, 0.)
        n_coms_rv_average_std[z] = np.insert(n_coms_rv_average_std[z], 0, 0.)
        n_coms_rv_median[z] = np.insert(n_coms_rv_median[z], 0, 0.)

        n_coms_mass[z] = np.insert(n_coms_mass[z], 0, 0.); n_coms_mass_average[z] = np.insert(n_coms_mass_average[z], 0, 0.)
        n_coms_mass_average_std[z] = np.insert(n_coms_mass_average_std[z], 0, 0.)
        n_coms_mass_median[z] = np.insert(n_coms_mass_median[z], 0, 0.)

    for z in range(3):
        n_coms_z[z] = np.insert(n_coms_z[z], 0, 0.); n_coms_z_average[z] = np.insert(n_coms_z_average[z], 0, 0.)
        n_coms_z_average_std[z] = np.insert(n_coms_z_average_std[z], 0, 0.)
        n_coms_z_median[z] = np.insert(n_coms_z_median[z], 0, 0.)


        n_coms_rg[z] = np.insert(n_coms_rg[z], 0, 0.); n_coms_rg_average[z] = np.insert(n_coms_rg_average[z], 0, 0.)
        n_coms_rg_average_std[z] = np.insert(n_coms_rg_average_std[z], 0, 0.)
        n_coms_rg_median[z] = np.insert(n_coms_rg_median[z], 0, 0.)

    if coflag=='NS':
        filenames = ['nnsms_mass_age_all.dat', 'nnsms_rv_age_all.dat', 'nnsms_z_age_all.dat', 'nnsms_rg_age_all.dat']
    if coflag=='BH':
        filenames = ['nbhms_mass_age_all.dat', 'nbhms_rv_age_all.dat', 'nbhms_z_age_all.dat', 'nbhms_rg_age_all.dat']

    np.savetxt('/projects/b1095/syr904/projects/GAIA_COMS/'+filenames[0], np.c_[t_all, n_coms_mass[0], n_coms_mass[1], n_coms_mass[2], n_coms_mass[3], n_coms_mass_average[0], n_coms_mass_average[1], n_coms_mass_average[2], n_coms_mass_average[3], n_coms_mass_average_std[0], n_coms_mass_average_std[1], n_coms_mass_average_std[2], n_coms_mass_average_std[3], n_coms_mass_median[0], n_coms_mass_median[1], n_coms_mass_median[2], n_coms_mass_median[3]], fmt = '%f %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g', header = '1.time(Myr) 2.N_2e5 3.N_4e5 4.N_8e5 5.N_16e5 6.N_2e5_ave 7.N_4e5_ave 8.N_8e5_ave 9.N_16e5_ave 10.N_2e5_ave_std 11.N_4e5_ave_std 12.N_8e5_ave_std 13.N_16e5_ave_std 14.N_2e5_med 15.N_4e5_med 16.N_8e5_med 17.N_16e5_med', comments = '#', delimiter = ' ')

    np.savetxt('/projects/b1095/syr904/projects/GAIA_COMS/'+filenames[1], np.c_[t_all, n_coms_rv[0], n_coms_rv[1], n_coms_rv[2], n_coms_rv[3], n_coms_rv_average[0], n_coms_rv_average[1], n_coms_rv_average[2], n_coms_rv_average[3], n_coms_rv_average_std[0], n_coms_rv_average_std[1], n_coms_rv_average_std[2], n_coms_rv_average_std[3], n_coms_rv_median[0], n_coms_rv_median[1], n_coms_rv_median[2], n_coms_rv_median[3]], fmt = '%f %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g', header = '1.time(Myr) 2.rv_4 3.rv_2 4.rv_1 5.rv_0.5 6.rv_4_ave 7.rv_2_ave 8.rv_1_ave 9.rv_0.5_ave 10.rv_4_ave_std 11.rv_2_ave_std 12.rv_1_ave_std 13.rv_0.5_ave_std 14.rv_4_med 15.rv_2_med 16.rv_1_med 17.rv_0.5_med', comments = '#', delimiter = ' ')
    
    np.savetxt('/projects/b1095/syr904/projects/GAIA_COMS/'+filenames[2], np.c_[t_all, n_coms_z[0], n_coms_z[1], n_coms_z[2], n_coms_z_average[0], n_coms_z_average[1], n_coms_z_average[2], n_coms_z_average_std[0], n_coms_z_average_std[1], n_coms_z_average_std[2], n_coms_z_median[0], n_coms_z_median[1], n_coms_z_median[2]], fmt = '%f %g %g %g %g %g %g %g %g %g %g %g %g', header = '1.time(Myr) 2.z_0.0002 3.z_0.002 4.z_0.02 5.z_0.0002_ave 6.z_0.002_ave 7.z_0.02_ave 8.z_0.0002_ave_std 9.z_0.002_ave_std 10.z_0.02_ave_std 11.z_0.0002_med 12.z_0.002_med 13.z_0.02_med', comments = '#', delimiter = ' ')

    np.savetxt('/projects/b1095/syr904/projects/GAIA_COMS/'+filenames[3], np.c_[t_all, n_coms_rg[0], n_coms_rg[1], n_coms_rg[2], n_coms_rg_average[0], n_coms_rg_average[1], n_coms_rg_average[2], n_coms_rg_average_std[0], n_coms_rg_average_std[1], n_coms_rg_average_std[2],n_coms_rg_median[0], n_coms_rg_median[1], n_coms_rg_median[2] ], fmt = '%f %g %g %g %g %g %g %g %g %g %g %g %g', header = '1.time(Myr) 2.rg_2 3.rg_8 4.rg_20 5.rg_2_ave 6.rg_8_ave 7.rg_20_ave 8.rg_2_ave_std 9.rg_8_ave_std 10.rg_20_ave_std 11.rg_2_med 12.rg_8_med 13.rg_20_med', comments = '#', delimiter = ' ')




def mapping_cmc_galaxy_nsms(mapflag):
    pathlist = np.genfromtxt('/projects/b1095/syr904/projects/GAIA_COMS/path_allfinished_newruns_maingrid.dat', dtype=str)
    paths = pathlist[:,0]; status = pathlist[:,1]

    ###Extract cluster parameters from galaxy simulation
    galaxy_csv = '/projects/b1095/syr904/projects/GAIA_COMS/nnsms_simulations.csv'
    df = pd.read_csv(galaxy_csv, 
        usecols = ['disruption_timescale', 'formation_time', 'feh', 'cluster_radius_after', 'N_NSMS_rounded', 'CMC_mass', 'flag_disruption', 'cluster_radius_initial'])
    print(len(df.index))

    df_selected = df[(df['flag_disruption']==1) & (df['cluster_radius_after'] !=-10)] # & (df['N_NSMS_rounded']>0)
    ##Note that by including the condition df['cluster_radius_after'] !=-10 we miss 81 disrupted clusters.
    print(len(df_selected.index))
    clus_no = np.array(df_selected.index)
    print(clus_no, len(clus_no))

    clus_mass = np.array(df_selected['CMC_mass']).astype('object')
    clus_mass_str = clus_mass
    clus_feh = np.array(df_selected['feh']); clus_z = (uc.metallicity(clus_feh, 'fe/htoz'))
    clus_rg_disrupt = np.array(df_selected['cluster_radius_after'])
    clus_rg_initial = np.array(df_selected['cluster_radius_initial'])
    t_disrupt = np.array(df_selected['disruption_timescale'])*1000# - np.array(df_selected['formation_time']))*1000.  ##in Myr
    t_disrupt_low = t_disrupt - 0.1*t_disrupt
    t_disrupt_up =  t_disrupt + 0.1*t_disrupt
    print(clus_mass, clus_z, clus_rg_initial)

    ###Randomly sample rv since the galaxy simulation does not vary rv
    rv = np.array(['0.5', '1', '2', '4'])
    rv_sample = np.random.choice(rv, len(clus_mass))

    ###Mapping the galaxy simulated clusters to CMC initial conditions
    clus_mass_str[clus_mass_str==200000]='2e5'; clus_mass_str[clus_mass_str==400000]='4e5'
    clus_mass_str[clus_mass_str==800000]='8e5'; clus_mass_str[clus_mass_str==1600000]='1.6e6'

    if mapflag == 'initial_rg':
        clus_z_str = []; clus_rg_str = []

        for xx in range(len(clus_z)):
            if clus_z[xx]>0.0065:
                clus_z_str.append('0.02')
            elif clus_z[xx]>0.00065 and clus_z[xx]<=0.0065:
                clus_z_str.append('0.002')
            else:
                clus_z_str.append('0.0002')
        
        #clus_z[clus_z>0.0065]='0.02'; clus_z[(clus_z>0.00065) & (clus_z<=0.0065)]='0.002'; clus_z[clus_z<=0.00065]='0.0002'

        for yy in range(len(clus_rg_initial)):
            if clus_rg_initial[yy]>14:
                clus_rg_str.append('20')
            elif clus_rg_initial[yy]>5 and clus_rg_initial[yy]<=14:
                clus_rg_str.append('8')
            else:
                clus_rg_str.append('2')

    elif mapflag == 'disrupt_rg':
        clus_z_str = []; clus_rg_str = []

        for xx in range(len(clus_z)):
            if clus_z[xx]>0.0065:
                clus_z_str.append('0.02')
            elif clus_z[xx]>0.00065 and clus_z[xx]<=0.0065:
                clus_z_str.append('0.002')
            else:
                clus_z_str.append('0.0002')

        #clus_z[clus_z>0.0065]='0.02'; clus_z[(clus_z>0.00065) & (clus_z<=0.0065)]='0.002'; clus_z[clus_z<=0.00065]='0.0002'

        for yy in range(len(clus_rg_disrupt)):
            if clus_rg_disrupt[yy]>14:
                clus_rg_str.append('20')
            elif clus_rg_disrupt[yy]>5 and clus_rg_disrupt[yy]<=14:
                clus_rg_str.append('8')
            else:
                clus_rg_str.append('2')

    elif mapflag == 'nomap':
        ##### Also randomly sampling rg and z #####
    
        rg = np.array(['2', '8', '20']) 
        clus_rg_str = np.random.choice(rg, len(clus_mass))
        
        z = np.array(['0.02', '0.002', '0.0002'])
        clus_z_str = np.random.choice(z, len(clus_mass))


        #clus_rg[clus_rg>14]='20'; clus_rg[(clus_rg>5) & (clus_rg<=14)]='8'; clus_rg[clus_rg<=5]='2'

    print(clus_mass_str, clus_z_str, clus_rg_str)


    ###Extracting NS-MS binaries
    #K0 = []; K1 = []; M0 = []; M1 = []; B = []; P = []; SMA = []; ECC = []; ID0 = []; ID1 = []; Tmyr = []
    #K0_model = []; K1_model = []; M0_model = []; M1_model = []; B_model = []; P_model = []; SMA_model = []; ECC_model = []; ID0_model = []; ID1_model = []
    #Tmyr_model = []; Clus_No_model = []
    fwrite = open('/projects/b1095/syr904/projects/GAIA_COMS/nsms_properties_'+mapflag+'.dat', 'w+')
    fwrite.write('#1.Clus_No 2.T[Myr] 3.ID0 4.ID1 5.M0[Msun] 6.M1[Msun] 7.K0 8.K1 9.SMA[AU] 10.ECC 11.B[G] 12.P[sec] 13.Mass(CMC) 14.rv(CMC) 15.rg(CMC) 16.z(CMC)\n')
    for ii in range(len(clus_mass)):
        thepath = '/projects/b1091/CMC_Grid_March2019/rundir/rv'+rv_sample[ii]+'/rg'+clus_rg_str[ii]+'/z'+clus_z_str[ii]+'/'+clus_mass_str[ii]+'/'
        t_conv = dyn.conv('t', thepath+'initial.conv.sh')
        psrfile = thepath+'initial.morepulsars.dat'
        #print(clus_no[ii])
        
        #K0_model = []; K1_model = []; M0_model = []; M1_model = []; B_model = []; P_model = []; SMA_model = []; ECC_model = []; ID0_model = []; ID1_model = []
        #Tmyr_model = []
        thetime=14000.
        with open(psrfile, 'r') as fpsr:
            next(fpsr)
            for line in fpsr:
                data = line.split()
                #if float(data[1])*t_conv > t_disrupt_low[ii] and float(data[1])*t_conv < t_disrupt_up[ii]:

                if float(data[1])*t_conv>thetime:
                    break

                if float(data[1])*t_conv>=t_disrupt[ii]:
                    thetime = float(data[1])*t_conv
                    if int(data[2])==1:
                        if int(data[11])==13 and (int(data[12])==0 or int(data[12])==1):
                            #K0_model.append(int(data[11])); K1_model.append(int(data[12]))
                            #M0_model.append(float(data[5])); M1_model.append(float(data[6]))
                            #B_model.append(float(data[7])); P_model.append(float(data[9]))
                            #SMA_model.append(float(data[13])); ECC_model.append(float(data[14]))
                            #ID0_model.append(int(data[3])); ID1_model.append(int(data[4]))
                            #Tmyr_model.append(float(data[1])*t_conv)
                            #Clus_No_model.append(clus_no[ii])
                            fwrite.write('%d %f %d %d %f %f %d %d %f %f %e %f %s %s %s %s\n'%(clus_no[ii], float(data[1])*t_conv, int(data[3]), int(data[4]),
                                float(data[5]), float(data[6]), int(data[11]), int(data[12]), float(data[13]), float(data[14]), float(data[7]), float(data[9]), clus_mass_str[ii], rv_sample[ii], clus_rg_str[ii], clus_z_str[ii]))

                        if int(data[12])==13 and (int(data[11])==0 or int(data[11])==1):
                            #K0_model.append(int(data[12])); K1_model.append(int(data[11]))
                            #M0_model.append(float(data[6])); M1_model.append(float(data[5]))
                            #B_model.append(float(data[8])); P_model.append(float(data[10]))
                            #SMA_model.append(float(data[13])); ECC_model.append(float(data[14]))
                            #ID0_model.append(int(data[4])); ID1_model.append(int(data[3]))
                            #Tmyr_model.append(float(data[1])*t_conv)
                            #Clus_No_model.append(clus_no[ii])

                            fwrite.write('%d %f %d %d %f %f %d %d %f %f %e %f %s %s %s %s\n'%(clus_no[ii], float(data[1])*t_conv, int(data[4]), int(data[3]),
                                float(data[6]), float(data[5]), int(data[12]), int(data[11]), float(data[13]), float(data[14]), float(data[8]), float(data[10]), clus_mass_str[ii], rv_sample[ii], clus_rg_str[ii], clus_z_str[ii]))

    fwrite.close()

    print('DONE')



def mapping_cmc_galaxy_coms(mapflag):
    pathlist = np.genfromtxt('/projects/b1095/syr904/projects/GAIA_COMS/path_allfinished_newruns_maingrid.dat', dtype=str)
    paths = pathlist[:,0]; status = pathlist[:,1]

    ###Extract cluster parameters from galaxy simulation
    galaxy_csv = '/projects/b1095/syr904/projects/GAIA_COMS/nnsms_simulations.csv'
    df = pd.read_csv(galaxy_csv, 
        usecols = ['disruption_timescale', 'formation_time', 'feh', 'cluster_radius_after', 'CMC_mass', 'flag_disruption', 'cluster_radius_initial'])
    print('length df.index', len(df.index))

    df_selected = df[(df['flag_disruption']==1) & (df['cluster_radius_after'] !=-10)] # & (df['N_NSMS_rounded']>0)
    ##Note that by including the condition df['cluster_radius_after'] !=-10 we miss 81 disrupted clusters.
    print('disrupted length df.index', len(df_selected.index))
    clus_no = np.array(df_selected.index)
    print('clus_no and length', clus_no, len(clus_no))

    clus_mass = np.array(df_selected['CMC_mass']).astype('object')
    clus_mass_str = clus_mass
    clus_feh = np.array(df_selected['feh']); clus_z = (uc.metallicity(clus_feh, 'fe/htoz'))
    clus_rg_disrupt = np.array(df_selected['cluster_radius_after'])
    clus_rg_initial = np.array(df_selected['cluster_radius_initial'])
    t_disrupt = np.array(df_selected['disruption_timescale'])*1000# - np.array(df_selected['formation_time']))*1000.  ##in Myr
    t_formation = np.array(df_selected['formation_time'])*1000. ##in Myr
    t_disrupt_low = t_disrupt - 0.1*t_disrupt
    t_disrupt_up =  t_disrupt + 0.1*t_disrupt
    print('clus_mass, clus_z, clus_rg_initial', clus_mass, clus_z, clus_rg_initial)

    ###Randomly sample rv since the galaxy simulation does not vary rv
    rv = np.array(['0.5', '1', '2', '4'])
    rv_sample = np.random.choice(rv, len(clus_mass))

    ###Mapping the galaxy simulated clusters to CMC initial conditions
    clus_mass_str[clus_mass_str==200000]='2e5'; clus_mass_str[clus_mass_str==400000]='4e5'
    clus_mass_str[clus_mass_str==800000]='8e5'; clus_mass_str[clus_mass_str==1600000]='1.6e6'

    if mapflag == 'initial_rg':
        clus_z_str = []; clus_rg_str = []

        for xx in range(len(clus_z)):
            if clus_z[xx]>0.0065:
                clus_z_str.append('0.02')
            elif clus_z[xx]>0.00065 and clus_z[xx]<=0.0065:
                clus_z_str.append('0.002')
            else:
                clus_z_str.append('0.0002')
        
        for yy in range(len(clus_rg_initial)):
            if clus_rg_initial[yy]>14:
                clus_rg_str.append('20')
            elif clus_rg_initial[yy]>5 and clus_rg_initial[yy]<=14:
                clus_rg_str.append('8')
            else:
                clus_rg_str.append('2')

    elif mapflag == 'disrupt_rg':
        clus_z_str = []; clus_rg_str = []

        for xx in range(len(clus_z)):
            if clus_z[xx]>0.0065:
                clus_z_str.append('0.02')
            elif clus_z[xx]>0.00065 and clus_z[xx]<=0.0065:
                clus_z_str.append('0.002')
            else:
                clus_z_str.append('0.0002')

        for yy in range(len(clus_rg_disrupt)):
            if clus_rg_disrupt[yy]>14:
                clus_rg_str.append('20')
            elif clus_rg_disrupt[yy]>5 and clus_rg_disrupt[yy]<=14:
                clus_rg_str.append('8')
            else:
                clus_rg_str.append('2')

    elif mapflag == 'nomap':
        ##### Also randomly sampling rg and z #####
    
        rg = np.array(['2', '8', '20']) 
        clus_rg_str = np.random.choice(rg, len(clus_mass))
        
        z = np.array(['0.02', '0.002', '0.0002'])
        clus_z_str = np.random.choice(z, len(clus_mass))


    elif mapflag == 'metal_map':
        clus_z_str = []
        for xx in range(len(clus_z)):
            if clus_z[xx]>0.0065:
                clus_z_str.append('0.02')
            elif clus_z[xx]>0.00065 and clus_z[xx]<=0.0065:
                clus_z_str.append('0.002')
            else:
                clus_z_str.append('0.0002')

        ##### Also randomly sampling rg #####
    
        rg = np.array(['2', '8', '20']) 
        clus_rg_str = np.random.choice(rg, len(clus_mass))
        


    #print(clus_mass_str, clus_z_str, clus_rg_str)
    print('clus_mass length', len(clus_mass))


    ###Extracting CO-MS binaries
    ###Inputs needed for COSMIC
    #cosmic_names_list = ["ktype1","ktype2","m1","m2","semimajor","ecc","ospin1","ospin2", "B1", "B2", "mass0_1","mass0_2","epoch1","epoch2","rad0","rad1","lum0","lum1","massc0","massc1","radc0","radc1","menv0","menv1","renv0","renv1","tms0","tms1","tacc0","tacc1","disruption_time","formation_time"]
    #cmc_names_list = ["id1", "id2", "Clus_no","T[Myr](galaxy)", "CO_flag", "Primodial_flag", "Mass(CMC)", "rv(CMC)", "rg(CMC)", "z(CMC)", "T[Myr](CMC)"]

    data_dict = {"kstar_1":[],"kstar_2":[],"mass_1":[],"mass_2":[],"semimajor":[],"ecc":[],"ospin_1":[],"ospin_2":[],"B_1":[],"B_2":[],"mass0_1":[],"mass0_2":[],"epoch_1":[],"epoch_2":[],"rad_1":[],"rad_2":[],"lum_1":[],"lum_2":[],"massc_1":[],"massc_2":[],"radc_1":[],"radc_2":[],"menv_1":[],"menv_2":[],"renv_1":[],"renv_2":[],"tms_1":[],"tms_2":[],"bacc_1":[],"bacc_2":[], "tacc_1":[],"tacc_2":[],"disruption_time":[],"formation_time":[],"id1":[], "id2":[], "Clus_no":[], "CO_flag":[], "Primodial_flag":[], "Mass(CMC)":[], "rv(CMC)":[], "rg(CMC)":[], "z(CMC)":[], "T[Myr](CMC)":[]}


    ii=0; count=0
    #num_of_binary = 0
    while ii < len(clus_mass):
        #if clus_no[ii]<16575: 
        #    ii+=1
        #    continue

        thepath = '/projects/b1091/CMC_Grid_March2019/rundir/rv'+rv_sample[ii]+'/rg'+clus_rg_str[ii]+'/z'+clus_z_str[ii]+'/'+clus_mass_str[ii]+'/'
        t_conv = dyn.conv('t', thepath+'initial.conv.sh')
        snapshots = dyn.get_snapshots(thepath+'initial')
        #print(clus_no[ii])

        snapfirst = thepath+'initial.snap0000.dat.gz'
        id0_snap0, id1_snap0 = ids_primordial(snapfirst)

        check = 0
        for jj in range(len(snapshots)):
            snaptime = dyn.get_time(snapshots[jj])*t_conv
            if t_disrupt_low[ii] <= snaptime < t_disrupt_up[ii] or snaptime >= t_disrupt[ii]:
                check = 1
                with gzip.open(snapshots[jj]) as fsnap:
                    next(fsnap); next(fsnap)
                    for line in fsnap:
                        data = line.split()
                        if int(data[7])==1:
                            ###WD-MS
                            if ((int(data[17])==11 or int(data[17])==12) and (int(data[18])==0 or int(data[18])==1)) or ((int(data[18])==11 or int(data[18])==12) and (int(data[17])==0 or int(data[17])==1)):
                                primordial = 0
                                for hh in range(len(id0_snap0)):
                                    if (int(data[10]) == int(id0_snap0[hh]) and int(data[11]) == int(id1_snap0[hh])) or (int(data[10]) == int(id1_snap0[hh]) and int(data[11]) == int(id0_snap0[hh])):
                                        primordial = 1
                                        break

                                #num_of_binary+=1
                
                                #tage = snaptime - float(data[58]) ##this is the epoch; age = tphys - epoch
                                #if tage < 0 :
                                #    print(int(data[10]), int(data[11]), snapshots[jj])

                                ##34 COSMIC params
                                data_dict['mass_1'].append(float(data[8]))
                                data_dict['mass_2'].append(float(data[9]))
                                data_dict['semimajor'].append(float(data[12]))  ##in AU
                                data_dict['ecc'].append(float(data[13])) 
                                data_dict['kstar_1'].append(int(data[17]))
                                data_dict['kstar_2'].append(int(data[18])) 
                                data_dict['rad_1'].append(float(data[26])) 
                                data_dict['rad_2'].append(float(data[27])) 
                                data_dict['lum_1'].append(float(data[29])) 
                                data_dict['lum_2'].append(float(data[30]))
                                data_dict['massc_1'].append(float(data[31])) 
                                data_dict['massc_2'].append(float(data[32]))  
                                data_dict['radc_1'].append(float(data[33])) 
                                data_dict['radc_2'].append(float(data[34])) 
                                data_dict['menv_1'].append(float(data[35])) 
                                data_dict['menv_2'].append(float(data[36]))
                                data_dict['renv_1'].append(float(data[37])) 
                                data_dict['renv_2'].append(float(data[38])) 
                                data_dict['tms_1'].append(float(data[39])) 
                                data_dict['tms_2'].append(float(data[40]))
                                #data_dict['RRLO_1'].append(float(data[43])) 
                                #data_dict['RRLO_2'].append(float(data[44]))
                                data_dict['ospin_1'].append(float(data[45])) 
                                data_dict['ospin_2'].append(float(data[46])) 
                                data_dict['B_1'].append(float(data[47])) 
                                data_dict['B_2'].append(float(data[48])) 
                                data_dict['bacc_1'].append(float(data[51])) 
                                data_dict['bacc_2'].append(float(data[52])) 
                                data_dict['tacc_1'].append(float(data[53])) 
                                data_dict['tacc_2'].append(float(data[54]))
                                data_dict['mass0_1'].append(float(data[55])) 
                                data_dict['mass0_2'].append(float(data[56]))  
                                data_dict['epoch_1'].append(float(data[57])) 
                                data_dict['epoch_2'].append(float(data[58]))
                                
                                data_dict["disruption_time"].append(t_disrupt[ii])
                                data_dict["formation_time"].append(t_formation[ii])
                                data_dict["id1"].append(int(data[10]))
                                data_dict["id2"] .append(int(data[11]))
                                data_dict["Clus_no"].append(clus_no[ii])
                                data_dict["CO_flag"].append('WD')
                                data_dict["Primodial_flag"].append(primordial)
                                data_dict["Mass(CMC)"].append(clus_mass_str[ii])
                                data_dict["rv(CMC)"].append(rv_sample[ii])
                                data_dict["rg(CMC)"].append(clus_rg_str[ii]) 
                                data_dict["z(CMC)"].append(clus_z_str[ii])
                                data_dict["T[Myr](CMC)"].append(snaptime)


                            ###NS-MS
                            if (int(data[17])==13 and (int(data[18])==0 or int(data[18])==1)) or (int(data[18])==13 and (int(data[17])==0 or int(data[17])==1)):
                                primordial = 0
                                for hh in range(len(id0_snap0)):
                                    if (int(data[10]) == int(id0_snap0[hh]) and int(data[11]) == int(id1_snap0[hh])) or (int(data[10]) == int(id1_snap0[hh]) and int(data[11]) == int(id0_snap0[hh])):
                                        primordial = 1
                                        break

                                #num_of_binary+=1
                                
                                #Pspin = twopi*yearsc/float(data[45])

                                #tage = snaptime - float(data[58]) ##this is the epoch; age = tphys - epoch
                                #if tage < 0 :
                                #    print(int(data[10]), int(data[11]), snapshots[jj])

                                data_dict['mass_1'].append(float(data[8]))
                                data_dict['mass_2'].append(float(data[9]))
                                data_dict['semimajor'].append(float(data[12]))  ##in AU
                                data_dict['ecc'].append(float(data[13])) 
                                data_dict['kstar_1'].append(int(data[17]))
                                data_dict['kstar_2'].append(int(data[18])) 
                                data_dict['rad_1'].append(float(data[26])) 
                                data_dict['rad_2'].append(float(data[27])) 
                                data_dict['lum_1'].append(float(data[29])) 
                                data_dict['lum_2'].append(float(data[30]))
                                data_dict['massc_1'].append(float(data[31])) 
                                data_dict['massc_2'].append(float(data[32]))  
                                data_dict['radc_1'].append(float(data[33])) 
                                data_dict['radc_2'].append(float(data[34])) 
                                data_dict['menv_1'].append(float(data[35])) 
                                data_dict['menv_2'].append(float(data[36]))
                                data_dict['renv_1'].append(float(data[37])) 
                                data_dict['renv_2'].append(float(data[38])) 
                                data_dict['tms_1'].append(float(data[39])) 
                                data_dict['tms_2'].append(float(data[40]))
                                #data_dict['RRLO_1'].append(float(data[43])) 
                                #data_dict['RRLO_2'].append(float(data[44]))
                                data_dict['ospin_1'].append(float(data[45])) 
                                data_dict['ospin_2'].append(float(data[46])) 
                                data_dict['B_1'].append(float(data[47])) 
                                data_dict['B_2'].append(float(data[48])) 
                                data_dict['bacc_1'].append(float(data[51])) 
                                data_dict['bacc_2'].append(float(data[52])) 
                                data_dict['tacc_1'].append(float(data[53])) 
                                data_dict['tacc_2'].append(float(data[54]))
                                data_dict['mass0_1'].append(float(data[55])) 
                                data_dict['mass0_2'].append(float(data[56]))  
                                data_dict['epoch_1'].append(float(data[57])) 
                                data_dict['epoch_2'].append(float(data[58]))
                                
                                data_dict["disruption_time"].append(t_disrupt[ii])
                                data_dict["formation_time"].append(t_formation[ii])
                                data_dict["id1"].append(int(data[10]))
                                data_dict["id2"] .append(int(data[11]))
                                data_dict["Clus_no"].append(clus_no[ii])
                                data_dict["CO_flag"].append('NS')
                                data_dict["Primodial_flag"].append(primordial)
                                data_dict["Mass(CMC)"].append(clus_mass_str[ii])
                                data_dict["rv(CMC)"].append(rv_sample[ii])
                                data_dict["rg(CMC)"].append(clus_rg_str[ii]) 
                                data_dict["z(CMC)"].append(clus_z_str[ii])
                                data_dict["T[Myr](CMC)"].append(snaptime)



                            ###BH-MS
                            if (int(data[17])==14 and (int(data[18])==0 or int(data[18])==1)) or (int(data[18])==14 and (int(data[17])==0 or int(data[17])==1)):
                                primordial = 0
                                for hh in range(len(id0_snap0)):
                                    if (int(data[10]) == int(id0_snap0[hh]) and int(data[11]) == int(id1_snap0[hh])) or (int(data[10]) == int(id1_snap0[hh]) and int(data[11]) == int(id0_snap0[hh])):
                                        primordial = 1
                                        break

                                #num_of_binary+=1

                                #tage = snaptime - float(data[58]) ##this is the epoch; age = tphys - epoch
                                #if tage < 0 :
                                #    print(int(data[10]), int(data[11]), snapshots[jj])

                                data_dict['mass_1'].append(float(data[8]))
                                data_dict['mass_2'].append(float(data[9]))
                                data_dict['semimajor'].append(float(data[12]))  ##in AU
                                data_dict['ecc'].append(float(data[13])) 
                                data_dict['kstar_1'].append(int(data[17]))
                                data_dict['kstar_2'].append(int(data[18])) 
                                data_dict['rad_1'].append(float(data[26])) 
                                data_dict['rad_2'].append(float(data[27])) 
                                data_dict['lum_1'].append(float(data[29])) 
                                data_dict['lum_2'].append(float(data[30]))
                                data_dict['massc_1'].append(float(data[31])) 
                                data_dict['massc_2'].append(float(data[32]))  
                                data_dict['radc_1'].append(float(data[33])) 
                                data_dict['radc_2'].append(float(data[34])) 
                                data_dict['menv_1'].append(float(data[35])) 
                                data_dict['menv_2'].append(float(data[36]))
                                data_dict['renv_1'].append(float(data[37])) 
                                data_dict['renv_2'].append(float(data[38])) 
                                data_dict['tms_1'].append(float(data[39])) 
                                data_dict['tms_2'].append(float(data[40]))
                                #data_dict['RRLO_1'].append(float(data[43])) 
                                #data_dict['RRLO_2'].append(float(data[44]))
                                data_dict['ospin_1'].append(float(data[45])) 
                                data_dict['ospin_2'].append(float(data[46])) 
                                data_dict['B_1'].append(float(data[47])) 
                                data_dict['B_2'].append(float(data[48])) 
                                data_dict['bacc_1'].append(float(data[51])) 
                                data_dict['bacc_2'].append(float(data[52])) 
                                data_dict['tacc_1'].append(float(data[53])) 
                                data_dict['tacc_2'].append(float(data[54]))
                                data_dict['mass0_1'].append(float(data[55])) 
                                data_dict['mass0_2'].append(float(data[56]))  
                                data_dict['epoch_1'].append(float(data[57])) 
                                data_dict['epoch_2'].append(float(data[58]))
                                
                                data_dict["disruption_time"].append(t_disrupt[ii])
                                data_dict["formation_time"].append(t_formation[ii])
                                data_dict["id1"].append(int(data[10]))
                                data_dict["id2"] .append(int(data[11]))
                                data_dict["Clus_no"].append(clus_no[ii])
                                data_dict["CO_flag"].append('BH')
                                data_dict["Primodial_flag"].append(primordial)
                                data_dict["Mass(CMC)"].append(clus_mass_str[ii])
                                data_dict["rv(CMC)"].append(rv_sample[ii])
                                data_dict["rg(CMC)"].append(clus_rg_str[ii]) 
                                data_dict["z(CMC)"].append(clus_z_str[ii])
                                data_dict["T[Myr](CMC)"].append(snaptime)


                break

        
        if check==0:
            print(ii, thepath, t_disrupt[ii])
            rv_sample[ii] = np.random.choice(rv)
            count+=1

            if count>=4:
                ii+=1
                count=0

        else:
            ii+=1
            count=0

        #if num_of_binary>2:
        #    break
        #print(ii)

    print('Mapping and Extraction DONE')
    return data_dict


##Run the extraction CO-MS function and evolve the binaries to a Hubble time through COSMIC
def run_extract_to_cosmic(mapflag, vfile):
    binary_dict = mapping_cmc_galaxy_coms(mapflag)

    #print(binary_dict["kstar_1"],binary_dict["kstar_2"],binary_dict["mass_1"],binary_dict["mass_2"],binary_dict["semimajor"],binary_dict["ecc"],binary_dict["disruption_time"],binary_dict["formation_time"],[float(s) for s in binary_dict["z(CMC)"]])
    
    #for xx in range(2):
    #    print(binary_dict["kstar_1"][xx],binary_dict["kstar_2"][xx],binary_dict["mass_1"][xx],binary_dict["mass_2"][xx],binary_dict["semimajor"][xx],binary_dict["ecc"][xx],binary_dict["ospin_1"][xx],binary_dict["ospin_2"][xx],binary_dict["B_1"][xx],binary_dict["B_2"][xx],binary_dict["mass0_1"][xx],binary_dict["mass0_2"][xx],binary_dict["epoch_1"][xx],binary_dict["epoch_2"][xx],binary_dict["rad_1"][xx],binary_dict["rad_2"][xx],binary_dict["lum_1"][xx],binary_dict["lum_2"][xx],binary_dict["massc_1"][xx],binary_dict["massc_2"][xx],binary_dict["radc_1"][xx],binary_dict["radc_2"][xx],binary_dict["menv_1"][xx],binary_dict["menv_2"][xx],binary_dict["renv_1"][xx],binary_dict["renv_2"][xx],binary_dict["tms_1"][xx],binary_dict["tms_2"][xx],binary_dict["bacc_1"][xx],binary_dict["bacc_2"][xx],binary_dict["tacc_1"][xx],binary_dict["tacc_2"][xx],binary_dict["disruption_time"][xx],binary_dict["formation_time"][xx],binary_dict["z(CMC)"][xx])

    mass_1,lum_1,teff_1,kstar_1,B_1,ospin_1,RRLO_1,mass_2,lum_2,teff_2,kstar_2,B_2,ospin_2,RRLO_2,porb,ecc,bin_state = evolve_binaries_to_z0(binary_dict["kstar_1"],binary_dict["kstar_2"],binary_dict["mass_1"],binary_dict["mass_2"],binary_dict["semimajor"],binary_dict["ecc"],binary_dict["ospin_1"],binary_dict["ospin_2"],binary_dict["B_1"],binary_dict["B_2"],binary_dict["mass0_1"],binary_dict["mass0_2"],binary_dict["epoch_1"],binary_dict["epoch_2"],binary_dict["rad_1"],binary_dict["rad_2"],binary_dict["lum_1"],binary_dict["lum_2"],binary_dict["massc_1"],binary_dict["massc_2"],binary_dict["radc_1"],binary_dict["radc_2"],binary_dict["menv_1"],binary_dict["menv_2"],binary_dict["renv_1"],binary_dict["renv_2"],binary_dict["tms_1"],binary_dict["tms_2"],binary_dict["bacc_1"],binary_dict["bacc_2"],binary_dict["tacc_1"],binary_dict["tacc_2"],binary_dict["T[Myr](CMC)"],binary_dict["formation_time"],[float(s) for s in binary_dict["z(CMC)"]])
 
    #print(mass_1, mass_2, porb, ecc)


    #data_dict = {"kstar_1":[],"kstar_2":[],"mass_1":[],"mass_2":[],"semimajor":[],"ecc":[],"ospin_1":[],"ospin_2":[],"B_1":[],"B_2":[],"mass0_1":[],"mass0_2":[],"epoch_1":[],"epoch_2":[],"rad_1":[],"rad_2":[],"lum_1":[],"lum_2":[],"massc_1":[],"massc_2":[],"radc_1":[],"radc_2":[],"menv_1":[],"menv_2":[],"renv_1":[],"renv_2":[],"tms_1":[],"tms_2":[],"bacc_1":[],"bacc_2":[], "tacc_1":[],"tacc_2":[],"disruption_time":[],"formation_time":[],"id1":[], "id2":[], "Clus_no":[], "CO_flag":[], "Primodial_flag":[], "Mass(CMC)":[], "rv(CMC)":[], "rg(CMC)":[], "z(CMC)":[], "T[Myr](CMC)":[]}
    fwrite = open('/projects/b1095/syr904/projects/GAIA_COMS/coms_properties_presentday_'+mapflag+'_'+vfile+'.dat', 'w+')
    fwrite.write('#1.Clus_No 2.Tdisrupt[Myr](galaxy) 3.ID0 4.ID1 5.M0[Msun] 6.M1[Msun] 7.K0 8.K1 9.Porb[Days] 10.ECC 11.Lum[Lsun] 12.Teff[K] 13.CO_flag 14.Primodial_flag 15.B[G] 16.Pspin[sec] 17.Mass(CMC) 18.rv(CMC) 19.rg(CMC) 20.z(CMC) 21.T[Myr](CMC) 22.Tform[Myr](galaxy)\n')

    num_COMS = 0
    for xx in range(len(bin_state)):
        if bin_state[xx] == 0 and ecc[xx]!=-1. and mass_1[xx]!=-100.:  ##Make sure it's still a binary
            if kstar_2[xx]<=1 and RRLO_2[xx]<1:  ##Check for star type and Roche Lobe overflow
                num_COMS+=1
                Pspin1 = twopi*yearsc/ospin_1[xx]
                fwrite.write('%d %f %d %d %f %f %d %d %f %f %f %f %s %d %e %f %s %s %s %s %f %f\n'%(binary_dict['Clus_no'][xx],binary_dict['disruption_time'][xx],binary_dict['id1'][xx],binary_dict['id2'][xx],mass_1[xx],mass_2[xx],kstar_1[xx],kstar_2[xx],porb[xx],ecc[xx],lum_2[xx],teff_2[xx],binary_dict['CO_flag'][xx],binary_dict['Primodial_flag'][xx],B_1[xx],Pspin1,binary_dict['Mass(CMC)'][xx],binary_dict['rv(CMC)'][xx],binary_dict['rg(CMC)'][xx],binary_dict['z(CMC)'][xx],binary_dict['T[Myr](CMC)'][xx],binary_dict['formation_time'][xx]))

            if kstar_1[xx]<=1 and RRLO_1[xx]<1:  ##Check for star type and Roche Lobe overflow
                num_COMS+=1
                Pspin2 = twopi*yearsc/ospin_2[xx]
                fwrite.write('%d %f %d %d %f %f %d %d %f %f %f %f %s %d %e %f %s %s %s %s %f %f\n'%(binary_dict['Clus_no'][xx],binary_dict['disruption_time'][xx],binary_dict['id1'][xx],binary_dict['id2'][xx],mass_2[xx],mass_1[xx],kstar_2[xx],kstar_1[xx],porb[xx],ecc[xx],lum_1[xx],teff_1[xx],binary_dict['CO_flag'][xx],binary_dict['Primodial_flag'][xx],B_2[xx],Pspin2,binary_dict['Mass(CMC)'][xx],binary_dict['rv(CMC)'][xx],binary_dict['rg(CMC)'][xx],binary_dict['z(CMC)'][xx],binary_dict['T[Myr](CMC)'][xx],binary_dict['formation_time'][xx]))

    fwrite.close()
    print('num_COMS', num_COMS)
    print('Extract+COSMIC DONE')
            

def run_extract_escaped_evolve_cosmic(comsfile, mapflag, vfile):
    binary_dict = extract_escaped_COMSs(comsfile)

    #for xx in range(2):
    #    print(binary_dict["kstar_1"][xx],binary_dict["kstar_2"][xx],binary_dict["mass_1"][xx],binary_dict["mass_2"][xx],binary_dict["semimajor"][xx],binary_dict["ecc"][xx],binary_dict["ospin_1"][xx],binary_dict["ospin_2"][xx],binary_dict["B_1"][xx],binary_dict["B_2"][xx],binary_dict["mass0_1"][xx],binary_dict["mass0_2"][xx],binary_dict["epoch_1"][xx],binary_dict["epoch_2"][xx],binary_dict["rad_1"][xx],binary_dict["rad_2"][xx],binary_dict["lum_1"][xx],binary_dict["lum_2"][xx],binary_dict["massc_1"][xx],binary_dict["massc_2"][xx],binary_dict["radc_1"][xx],binary_dict["radc_2"][xx],binary_dict["menv_1"][xx],binary_dict["menv_2"][xx],binary_dict["renv_1"][xx],binary_dict["renv_2"][xx],binary_dict["tms_1"][xx],binary_dict["tms_2"][xx],binary_dict["bacc_1"][xx],binary_dict["bacc_2"][xx],binary_dict["tacc_1"][xx],binary_dict["tacc_2"][xx],binary_dict["disruption_time"][xx],binary_dict["formation_time"][xx],binary_dict["z(CMC)"][xx])

    mass_1,lum_1,teff_1,kstar_1,B_1,ospin_1,RRLO_1,mass_2,lum_2,teff_2,kstar_2,B_2,ospin_2,RRLO_2,porb,ecc,bin_state = evolve_binaries_to_z0(binary_dict["kstar_1"],binary_dict["kstar_2"],binary_dict["mass_1"],binary_dict["mass_2"],binary_dict["semimajor"],binary_dict["ecc"],binary_dict["ospin_1"],binary_dict["ospin_2"],binary_dict["B_1"],binary_dict["B_2"],binary_dict["mass0_1"],binary_dict["mass0_2"],binary_dict["epoch_1"],binary_dict["epoch_2"],binary_dict["rad_1"],binary_dict["rad_2"],binary_dict["lum_1"],binary_dict["lum_2"],binary_dict["massc_1"],binary_dict["massc_2"],binary_dict["radc_1"],binary_dict["radc_2"],binary_dict["menv_1"],binary_dict["menv_2"],binary_dict["renv_1"],binary_dict["renv_2"],binary_dict["tms_1"],binary_dict["tms_2"],binary_dict["bacc_1"],binary_dict["bacc_2"],binary_dict["tacc_1"],binary_dict["tacc_2"],binary_dict["T[Myr](CMC)"],binary_dict["formation_time"],[float(s) for s in binary_dict["z(CMC)"]])
 
    #print(mass_1, mass_2, porb, ecc)

    fwrite = open('/projects/b1095/syr904/projects/GAIA_COMS/coms_escaped_presentday_'+mapflag+'_'+vfile+'.dat', 'w+')
    fwrite.write('#1.Clus_No 2.Tdisrupt[Myr](galaxy) 3.ID0 4.ID1 5.M0[Msun] 6.M1[Msun] 7.K0 8.K1 9.Porb[Days] 10.ECC 11.Lum[Lsun] 12.Teff[K] 13.CO_flag 14.Primodial_flag 15.B[G] 16.Pspin[sec] 17.Mass(CMC) 18.rv(CMC) 19.rg(CMC) 20.z(CMC) 21.T[Myr](CMC) 22.Tform[Myr](galaxy)\n')

    num_escaped_COMS = 0
    for xx in range(len(bin_state)):
        if bin_state[xx] == 0 and ecc[xx]!=-1. and mass_1[xx]!=-100.:  ##Make sure it's still a binary
            if kstar_2[xx]<=1 and RRLO_2[xx]<1 and kstar_1[xx]>10:  ##Check for star type and Roche Lobe overflow
                num_escaped_COMS+=1

                Pspin1 = twopi*yearsc/ospin_1[xx]

                if kstar_1[xx]==14:coflag='BH'
                elif kstar_1[xx]==13:coflag='NS'
                else:coflag='WD'

                fwrite.write('%d %f %d %d %f %f %d %d %f %f %f %f %s %d %e %f %s %s %s %s %f %f\n'%(int(binary_dict['Clus_no'][xx]),binary_dict['disruption_time'][xx],binary_dict['id1'][xx],binary_dict['id2'][xx],mass_1[xx],mass_2[xx],kstar_1[xx],kstar_2[xx],porb[xx],ecc[xx],lum_2[xx],teff_2[xx],coflag,binary_dict['Primodial_flag'][xx],B_1[xx],Pspin1,binary_dict['Mass(CMC)'][xx],binary_dict['rv(CMC)'][xx],binary_dict['rg(CMC)'][xx],binary_dict['z(CMC)'][xx],binary_dict['T[Myr](CMC)'][xx],binary_dict['formation_time'][xx]))

            if kstar_1[xx]<=1 and RRLO_1[xx]<1 and kstar_2[xx]>10:  ##Check for star type and Roche Lobe overflow
                num_escaped_COMS+=1

                Pspin2 = twopi*yearsc/ospin_2[xx]

                if kstar_2[xx]==14:coflag='BH'
                elif kstar_2[xx]==13:coflag='NS'
                else:coflag='WD'

                fwrite.write('%d %f %d %d %f %f %d %d %f %f %f %f %s %d %e %f %s %s %s %s %f %f\n'%(int(binary_dict['Clus_no'][xx]),binary_dict['disruption_time'][xx],binary_dict['id1'][xx],binary_dict['id2'][xx],mass_2[xx],mass_1[xx],kstar_2[xx],kstar_1[xx],porb[xx],ecc[xx],lum_1[xx],teff_1[xx],coflag,binary_dict['Primodial_flag'][xx],B_2[xx],Pspin2,binary_dict['Mass(CMC)'][xx],binary_dict['rv(CMC)'][xx],binary_dict['rg(CMC)'][xx],binary_dict['z(CMC)'][xx],binary_dict['T[Myr](CMC)'][xx],binary_dict['formation_time'][xx]))

    fwrite.close()
    print('num_escaped_COMS', num_escaped_COMS)
    print('Extract+COSMIC DONE')

