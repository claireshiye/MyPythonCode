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

from scipy.interpolate import interp1d
from scipy import stats
import scipy.integrate as integrate

sourcefile = np.genfromtxt('/projects/b1095/syr904/projects/GAIA_NSMS/path_allfinished_newruns_maingrid.dat',
             dtype='str')
paths = sourcefile[:,0]

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
    savepath = '/projects/b1095/syr904/projects/GAIA_NSMS' #specify path for catalog to be saved to 
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
    #savepath = '/projects/b1095/syr904/projects/GAIA_NSMS' #specify path for catalog to be saved to 
    #NSMS.to_csv(savepath+'/NSMS_catalog_incluster.csv')
    #
    #print('PULSARFILE done')


def print_N_NSMS_psrfile(pathlist, start, end, readflag):
    if readflag == 1:
        sourcedir=np.genfromtxt(pathlist, dtype='str')
        status=sourcedir[:,1]; 
        sourcedir=sourcedir[:,0]
    else:
        sourcedir=pathlist
        status = [1]
    
    savepath = '/projects/b1095/syr904/projects/GAIA_NSMS/'
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


##Extract and grouping the number of NS--MS from the catalog models
def extract_n_nsms():
    pathlist = np.genfromtxt('/projects/b1095/syr904/projects/GAIA_NSMS/path_allfinished_newruns_maingrid.dat', dtype=str)
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
    n_nsms_rv = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_nsms_rv_average = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_nsms_rv_average_std = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_nsms_rv_median = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]

    n_nsms_mass = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_nsms_mass_average = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_nsms_mass_average_std = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_nsms_mass_median = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]

    n_nsms_z = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_nsms_z_average = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_nsms_z_average_std = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_nsms_z_median = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]


    n_nsms_rg = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_nsms_rg_average = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_nsms_rg_average_std = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    n_nsms_rg_median = [np.zeros(bin_size),np.zeros(bin_size),np.zeros(bin_size)]
    
    nnsms_scatter_n2e5 = [[] for _ in range(bin_size)]
    nnsms_scatter_n4e5 = [[] for _ in range(bin_size)]
    nnsms_scatter_n8e5 = [[] for _ in range(bin_size)]
    nnsms_scatter_n16e5 = [[] for _ in range(bin_size)]

    nnsms_scatter_rv4 = [[] for _ in range(bin_size)]
    nnsms_scatter_rv2 = [[] for _ in range(bin_size)]
    nnsms_scatter_rv1 = [[] for _ in range(bin_size)]
    nnsms_scatter_rv05 = [[] for _ in range(bin_size)]

    nnsms_scatter_z00002 = [[] for _ in range(bin_size)]
    nnsms_scatter_z0002 = [[] for _ in range(bin_size)]
    nnsms_scatter_z002= [[] for _ in range(bin_size)]

    nnsms_scatter_rg2 = [[] for _ in range(bin_size)]
    nnsms_scatter_rg8 = [[] for _ in range(bin_size)]
    nnsms_scatter_rg20= [[] for _ in range(bin_size)]


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
        datans = np.genfromtxt(paths[kk]+'initial.ns.dat')
        times = np.array(datans[:,0])*t_conv
        n_nsms = datans[:,10]; n_nspms = datans[:,11]
        
        ##Interpolate the number of NS data
        f = interp1d(times, n_nsms, kind='nearest')
        t_interpld = np.linspace(0, np.max(times), 3*bin_size)
        n_nsms_new = f(t_interpld)
        #print(n_nsms_new)
    
        n_mass = [[],[],[],[]]; n_rv = [[],[],[],[]]; n_z = [[],[],[]]; n_rg = [[],[],[]]
        for jj in range(len(t_all)-1):
            #print(jj)
            n_mass_temp = [0,0,0,0]; n_rv_temp = [0,0,0,0]; n_z_temp = [0,0,0]; n_rg_temp = [0,0,0]
            count_mass = [0,0,0,0]; count_rv = [0,0,0,0]; count_z = [0,0,0]; count_rg = [0,0,0]
            
            ##Group by initial mass
            if n_star==200000.:# and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_mass_temp[0]+=n_nsms_new[i]
                        count_mass[0]+=1  ##multiple time steps may belong to the same bin
        
            if n_star==400000.:# and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_mass_temp[1]+=n_nsms_new[i]
                        count_mass[1]+=1
        
            if n_star==800000.:# and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_mass_temp[2]+=n_nsms_new[i]
                        count_mass[2]+=1
        
            if n_star==1600000.:# and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_mass_temp[3]+=n_nsms_new[i]
                        count_mass[3]+=1
            
            ##Group by initial rv   
            if rv==4.:# and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_rv_temp[0]+=n_nsms_new[i]
                        count_rv[0]+=1
        
            if rv==2.:# and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_rv_temp[1]+=n_nsms_new[i]
                        count_rv[1]+=1
        
            if rv==1.:# and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_rv_temp[2]+=n_nsms_new[i]
                        count_rv[2]+=1
        
            if rv==0.5:# and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_rv_temp[3]+=n_nsms_new[i]
                        count_rv[3]+=1
                    
        
            ##Group by metallicity
            if z==0.0002:# and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_z_temp[0]+=n_nsms_new[i]
                        count_z[0]+=1
        
            if z==0.002:# and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_z_temp[1]+=n_nsms_new[i]
                        count_z[1]+=1
        
            if z==0.02:# and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_z_temp[2]+=n_nsms_new[i]
                        count_z[2]+=1
                    
                    
            ##Group by galactocentric distance
            if rg==2:# and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_rg_temp[0]+=n_nsms_new[i]
                        count_rg[0]+=1
        
            if rg==8:# and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_rg_temp[1]+=n_nsms_new[i]
                        count_rg[1]+=1
        
            if rg==20:# and status[kk]=='1':
                for i in range(len(t_interpld)):
                    if t_all[jj] <= t_interpld[i] < t_all[jj+1]:
                        n_rg_temp[2]+=n_nsms_new[i]
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
            n_nsms_rv[y] = n_nsms_rv[y]+np.array(n_rv[y])
            n_nsms_mass[y] = n_nsms_mass[y]+np.array(n_mass[y])
            n_nsms_rv_average[y] = n_nsms_rv_average[y] + np.array(n_rv[y])/n_model_rv[y]
            n_nsms_mass_average[y] = n_nsms_mass_average[y] + np.array(n_mass[y])/n_model_mass[y]
            
        for y in range(3):
            n_nsms_z[y] = n_nsms_z[y]+np.array(n_z[y])
            n_nsms_z_average[y] = n_nsms_z_average[y] + np.array(n_z[y])/n_model_z[y]
            n_nsms_rg[y] = n_nsms_rg[y]+np.array(n_rg[y])
            n_nsms_rg_average[y] = n_nsms_rg_average[y] + np.array(n_rg[y])/n_model_rg[y]


        ##Group by initial mass
        if n_star==200000.:# and status[kk]=='1':
            #print(len(n_mass[0]))
            nnsms_scatter_n2e5 = np.hstack((nnsms_scatter_n2e5, np.split(np.array(n_mass[0]),len(n_mass[0]))))
            #print(nnsms_scatter_n2e5)
        
        if n_star==400000.:# and status[kk]=='1':
            nnsms_scatter_n4e5 = np.hstack((nnsms_scatter_n4e5, np.split(np.array(n_mass[1]), len(n_mass[1]))))
        
        if n_star==800000.:# and status[kk]=='1':
            nnsms_scatter_n8e5 = np.hstack((nnsms_scatter_n8e5, np.split(np.array(n_mass[2]), len(n_mass[2]))))

        if n_star==1600000.:# and status[kk]=='1':
            nnsms_scatter_n16e5 = np.hstack((nnsms_scatter_n16e5, np.split(np.array(n_mass[3]), len(n_mass[3]))))
            

        ##Group by initial rv   
        if rv==4.:# and status[kk]=='1':
            nnsms_scatter_rv4 = np.hstack((nnsms_scatter_rv4, np.split(np.array(n_rv[0]), len(n_rv[0]))))
        
        if rv==2.:# and status[kk]=='1':
            nnsms_scatter_rv2 = np.hstack((nnsms_scatter_rv2, np.split(np.array(n_rv[1]), len(n_rv[1]))))
        
        if rv==1.:# and status[kk]=='1':
            nnsms_scatter_rv1 = np.hstack((nnsms_scatter_rv1, np.split(np.array(n_rv[2]), len(n_rv[2]))))
        
        if rv==0.5:# and status[kk]=='1':
            nnsms_scatter_rv05 = np.hstack((nnsms_scatter_rv05, np.split(np.array(n_rv[3]), len(n_rv[3]))))
                    
    
        ##Group by metallicity
        if z==0.0002:# and status[kk]=='1':
            nnsms_scatter_z00002 = np.hstack((nnsms_scatter_z00002, np.split(np.array(n_z[0]), len(n_z[0]))))
        
        if z==0.002:# and status[kk]=='1':
            nnsms_scatter_z0002 = np.hstack((nnsms_scatter_z0002, np.split(np.array(n_z[1]), len(n_z[1]))))
        
        if z==0.02:# and status[kk]=='1':
            nnsms_scatter_z002 = np.hstack((nnsms_scatter_z002, np.split(np.array(n_z[2]), len(n_z[2]))))
                    
                    
        ##Group by galactocentric distance
        if rg==2:# and status[kk]=='1':
            nnsms_scatter_rg2 = np.hstack((nnsms_scatter_rg2, np.split(np.array(n_rg[0]), len(n_rg[0]))))
        
        if rg==8:# and status[kk]=='1':
            nnsms_scatter_rg8 = np.hstack((nnsms_scatter_rg8, np.split(np.array(n_rg[1]), len(n_rg[1]))))
        
        if rg==20:# and status[kk]=='1':
            nnsms_scatter_rg20 = np.hstack((nnsms_scatter_rg20, np.split(np.array(n_rg[2]), len(n_rg[2]))))

    for ii in range(4):
        for xx in range(bin_size):
            if ii == 0:
                n_nsms_mass_average_std[ii][xx]+=np.std(nnsms_scatter_n2e5[xx])
                n_nsms_rv_average_std[ii][xx]+=np.std(nnsms_scatter_rv4[xx])

                n_nsms_mass_median[ii][xx]+=np.median(nnsms_scatter_n2e5[xx])
                n_nsms_rv_median[ii][xx]+=np.median(nnsms_scatter_rv4[xx])
            if ii == 1:
                n_nsms_mass_average_std[ii][xx]+=np.std(nnsms_scatter_n4e5[xx])
                n_nsms_rv_average_std[ii][xx]+=np.std(nnsms_scatter_rv2[xx])

                n_nsms_mass_median[ii][xx]+=np.median(nnsms_scatter_n4e5[xx])
                n_nsms_rv_median[ii][xx]+=np.median(nnsms_scatter_rv2[xx])
            if ii == 2:
                n_nsms_mass_average_std[ii][xx]+=np.std(nnsms_scatter_n8e5[xx])
                n_nsms_rv_average_std[ii][xx]+=np.std(nnsms_scatter_rv1[xx])

                n_nsms_mass_median[ii][xx]+=np.median(nnsms_scatter_n8e5[xx])
                n_nsms_rv_median[ii][xx]+=np.median(nnsms_scatter_rv1[xx])
            if ii == 3:
                n_nsms_mass_average_std[ii][xx]+=np.std(nnsms_scatter_n16e5[xx])
                n_nsms_rv_average_std[ii][xx]+=np.std(nnsms_scatter_rv05[xx])

                n_nsms_mass_median[ii][xx]+=np.median(nnsms_scatter_n16e5[xx])
                n_nsms_rv_median[ii][xx]+=np.median(nnsms_scatter_rv05[xx])

    for ii in range(3):
        for xx in range(bin_size):
            if ii == 0:
                n_nsms_z_average_std[ii][xx]+=np.std(nnsms_scatter_z00002[xx])
                n_nsms_rg_average_std[ii][xx]+=np.std(nnsms_scatter_rg2[xx])

                n_nsms_z_median[ii][xx]+=np.median(nnsms_scatter_z00002[xx])
                n_nsms_rg_median[ii][xx]+=np.median(nnsms_scatter_rg2[xx])
            if ii == 1:
                n_nsms_z_average_std[ii][xx]+=np.std(nnsms_scatter_z0002[xx])
                n_nsms_rg_average_std[ii][xx]+=np.std(nnsms_scatter_rg8[xx])

                n_nsms_z_median[ii][xx]+=np.median(nnsms_scatter_z0002[xx])
                n_nsms_rg_median[ii][xx]+=np.median(nnsms_scatter_rg8[xx])
            if ii == 2:
                n_nsms_z_average_std[ii][xx]+=np.std(nnsms_scatter_z002[xx])
                n_nsms_rg_average_std[ii][xx]+=np.std(nnsms_scatter_rg20[xx])

                n_nsms_z_median[ii][xx]+=np.median(nnsms_scatter_z002[xx])
                n_nsms_rg_median[ii][xx]+=np.median(nnsms_scatter_rg20[xx])
    
    print(n_nsms_mass_average_std[0])

    for z in range(4):
        n_nsms_rv[z] = np.insert(n_nsms_rv[z], 0, 0.); n_nsms_rv_average[z] = np.insert(n_nsms_rv_average[z], 0, 0.)
        n_nsms_rv_average_std[z] = np.insert(n_nsms_rv_average_std[z], 0, 0.)
        n_nsms_rv_median[z] = np.insert(n_nsms_rv_median[z], 0, 0.)

        n_nsms_mass[z] = np.insert(n_nsms_mass[z], 0, 0.); n_nsms_mass_average[z] = np.insert(n_nsms_mass_average[z], 0, 0.)
        n_nsms_mass_average_std[z] = np.insert(n_nsms_mass_average_std[z], 0, 0.)
        n_nsms_mass_median[z] = np.insert(n_nsms_mass_median[z], 0, 0.)

    for z in range(3):
        n_nsms_z[z] = np.insert(n_nsms_z[z], 0, 0.); n_nsms_z_average[z] = np.insert(n_nsms_z_average[z], 0, 0.)
        n_nsms_z_average_std[z] = np.insert(n_nsms_z_average_std[z], 0, 0.)
        n_nsms_z_median[z] = np.insert(n_nsms_z_median[z], 0, 0.)


        n_nsms_rg[z] = np.insert(n_nsms_rg[z], 0, 0.); n_nsms_rg_average[z] = np.insert(n_nsms_rg_average[z], 0, 0.)
        n_nsms_rg_average_std[z] = np.insert(n_nsms_rg_average_std[z], 0, 0.)
        n_nsms_rg_median[z] = np.insert(n_nsms_rg_median[z], 0, 0.)

    
    filenames = ['nnsms_mass_age_all.dat', 'nnsms_rv_age_all.dat', 'nnsms_z_age_all.dat', 'nnsms_rg_age_all.dat']
    np.savetxt('/projects/b1095/syr904/projects/GAIA_NSMS/'+filenames[0], np.c_[t_all, n_nsms_mass[0], n_nsms_mass[1], n_nsms_mass[2], n_nsms_mass[3], n_nsms_mass_average[0], n_nsms_mass_average[1], n_nsms_mass_average[2], n_nsms_mass_average[3], n_nsms_mass_average_std[0], n_nsms_mass_average_std[1], n_nsms_mass_average_std[2], n_nsms_mass_average_std[3], n_nsms_mass_median[0], n_nsms_mass_median[1], n_nsms_mass_median[2], n_nsms_mass_median[3]], fmt = '%f %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g', header = '1.time(Myr) 2.N_2e5 3.N_4e5 4.N_8e5 5.N_16e5 6.N_2e5_ave 7.N_4e5_ave 8.N_8e5_ave 9.N_16e5_ave 10.N_2e5_ave_std 11.N_4e5_ave_std 12.N_8e5_ave_std 13.N_16e5_ave_std 14.N_2e5_med 15.N_4e5_med 16.N_8e5_med 17.N_16e5_med', comments = '#', delimiter = ' ')
    np.savetxt('/projects/b1095/syr904/projects/GAIA_NSMS/'+filenames[1], np.c_[t_all, n_nsms_rv[0], n_nsms_rv[1], n_nsms_rv[2], n_nsms_rv[3], n_nsms_rv_average[0], n_nsms_rv_average[1], n_nsms_rv_average[2], n_nsms_rv_average[3], n_nsms_rv_average_std[0], n_nsms_rv_average_std[1], n_nsms_rv_average_std[2], n_nsms_rv_average_std[3], n_nsms_rv_median[0], n_nsms_rv_median[1], n_nsms_rv_median[2], n_nsms_rv_median[3]], fmt = '%f %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g', header = '1.time(Myr) 2.rv_4 3.rv_2 4.rv_1 5.rv_0.5 6.rv_4_ave 7.rv_2_ave 8.rv_1_ave 9.rv_0.5_ave 10.rv_4_ave_std 11.rv_2_ave_std 12.rv_1_ave_std 13.rv_0.5_ave_std 14.rv_4_med 15.rv_2_med 16.rv_1_med 17.rv_0.5_med', comments = '#', delimiter = ' ')
    
    np.savetxt('/projects/b1095/syr904/projects/GAIA_NSMS/'+filenames[2], np.c_[t_all, n_nsms_z[0], n_nsms_z[1], n_nsms_z[2], n_nsms_z_average[0], n_nsms_z_average[1], n_nsms_z_average[2], n_nsms_z_average_std[0], n_nsms_z_average_std[1], n_nsms_z_average_std[2], n_nsms_z_median[0], n_nsms_z_median[1], n_nsms_z_median[2]], fmt = '%f %g %g %g %g %g %g %g %g %g %g %g %g', header = '1.time(Myr) 2.z_0.0002 3.z_0.002 4.z_0.02 5.z_0.0002_ave 6.z_0.002_ave 7.z_0.02_ave 8.z_0.0002_ave_std 9.z_0.002_ave_std 10.z_0.02_ave_std 11.z_0.0002_med 12.z_0.002_med 13.z_0.02_med', comments = '#', delimiter = ' ')
    np.savetxt('/projects/b1095/syr904/projects/GAIA_NSMS/'+filenames[3], np.c_[t_all, n_nsms_rg[0], n_nsms_rg[1], n_nsms_rg[2], n_nsms_rg_average[0], n_nsms_rg_average[1], n_nsms_rg_average[2], n_nsms_rg_average_std[0], n_nsms_rg_average_std[1], n_nsms_rg_average_std[2],n_nsms_rg_median[0], n_nsms_rg_median[1], n_nsms_rg_median[2] ], fmt = '%f %g %g %g %g %g %g %g %g %g %g %g %g', header = '1.time(Myr) 2.rg_2 3.rg_8 4.rg_20 5.rg_2_ave 6.rg_8_ave 7.rg_20_ave 8.rg_2_ave_std 9.rg_8_ave_std 10.rg_20_ave_std 11.rg_2_med 12.rg_8_med 13.rg_20_med', comments = '#', delimiter = ' ')



def mapping_cmc_galaxy_nsms(mapflag):
    pathlist = np.genfromtxt('/projects/b1095/syr904/projects/GAIA_NSMS/path_allfinished_newruns_maingrid.dat', dtype=str)
    paths = pathlist[:,0]; status = pathlist[:,1]

    ###Extract cluster parameters from galaxy simulation
    galaxy_csv = '/projects/b1095/syr904/projects/GAIA_NSMS/nnsms_simulations.csv'
    df = pd.read_csv(galaxy_csv, 
        usecols = ['disruption_timescale', 'formation_time', 'feh', 'cluster_radius_after', 'N_NSMS_rounded', 'CMC_mass', 'flag_disruption', 'cluster_radius_initial'])
    print(len(df.index))

    df_selected = df[(df['flag_disruption']==1) & (df['cluster_radius_after'] !=-10)] # & (df['N_NSMS_rounded']>0)
    print(len(df_selected.index))
    clus_no = np.array(df_selected.index)
    print(clus_no, len(clus_no))

    clus_mass = np.array(df_selected['CMC_mass']).astype('object')
    clus_mass_str = clus_mass
    clus_feh = np.array(df_selected['feh']); clus_z = (uc.metallicity(clus_feh, 'fe/htoz'))
    clus_rg_disrupt = np.array(df_selected['cluster_radius_after'])
    clus_rg_initial = np.array(df_selected['cluster_radius_initial'])
    t_disrupt = (np.array(df_selected['disruption_timescale']) - np.array(df_selected['formation_time']))*1000.  ##in Myr
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
    fwrite = open('/projects/b1095/syr904/projects/GAIA_NSMS/nsms_properties_'+mapflag+'.dat', 'w+')
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





            