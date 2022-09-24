import numpy as np
import constants
import os
import sys

from cosmic.sample import InitialCMCTable

# Set vg and rg for tidal radius calculation
vg=220.
rg = 7400  ##pc

N = [2700000,3000000,3300000]
Rv = [3, 4, 5]
Alpha1 = [0.4, 0.5]

for ii in range(len(N)):
    for jj in range(len(Rv)):
        for kk in range(len(Alpha1)):
            if N[ii]==3000000 and Rv[jj]==4 and Alpha1[kk]==0.4:
                continue


            ## CALCULATE TIDAL RADIUS IN UNIT OF THE VIRIAL RADIUS
            Singles, Binaries = InitialCMCTable.sampler('cmc', binfrac_model=0.022, primary_model='custom', alphas = [-Alpha1[kk], -2.8], mcuts = [0.08, 0.8, 150.],
                ecc_model='thermal', porb_model='log_uniform', qmin=0.1,
                cluster_profile='elson', met=0.0038, size=N[ii], gamma=2.1,
                r_max=200, params='ElsonProfile.ini',seed=12345, virial_radius=Rv[jj],
                tidal_radius = 45)
            mass_init = Singles.mass_of_cluster
            print(mass_init)
            r_tidal = (constants.G * mass_init * constants.Msun / 2. / (vg*constants.km)**2.)**(1./3.) * (rg*constants.PC)**(2./3.) / constants.PC
            r_tidal_code = r_tidal/Rv[jj]


            ## CREATE FOLDER
            DIR = '/projects/b1095/syr904/cmc/CMC-COSMIC/master_tc_test/ver_0601/near_fit_models/47Tuc_elson_rv'+str(Rv[jj])+'_n'+str(N[ii])+'_a1'+str(Alpha1[kk])+'_tcon'
            os.system('mkdir '+DIR)

            ## CREATE GENERATION FILE
            fgen = open(DIR+'/generate_elson_profile.py','w')
            print("# Using COSMIC's initial condition generators", file = fgen)
            print("# This imports the version of COSMIC that was downloaded and installed with CMC", file = fgen)
            print("# See https://clustermontecarlo.github.io/CMC-COSMIC/initialconditions/index.html for more details", file = fgen)
            print("from cosmic.sample import InitialCMCTable", file = fgen)
            print("", file =  fgen)
            print("# Generate the Singles and Binaries Pandas tables.", file = fgen)
            print("# See (https://cosmic-popsynth.github.io/COSMIC/runpop/index.html) for more details on the options here", file = fgen)
            print("Singles, Binaries = InitialCMCTable.sampler('cmc', binfrac_model=0.022, primary_model='custom',alphas = ["+str(-Alpha1[kk])+", -2.8], mcuts = [0.08, 0.8, 150.], ecc_model='thermal', porb_model='log_uniform', qmin=0.1, cluster_profile='elson', met=0.0038, size="+str(N[ii])+", gamma=2.1, r_max=200, params='ElsonProfile.ini', seed=12345, virial_radius="+str(Rv[jj])+", tidal_radius="+str(r_tidal_code)+")", file = fgen)
            print("", file =  fgen)
            print("print(Singles.mass_of_cluster)", file = fgen)
            print("# Scale the Cluster to Henon units (G = 1, M_cluster = 1, E_cluster = -0.25)", file = fgen)
            print("# Note that this option is automatically done in InitialCMCTable.write if the cluster", file = fgen)
            print("# isn't already normalized, but we do it here explicitly.", file = fgen)
            print("#InitialCMCTable.ScaleToNBodyUnits(Singles,Binaries)", file = fgen)
            print("", file =  fgen)
            print("# Save them to an hdf5 file for CMC", file = fgen)
            print('InitialCMCTable.write(Singles, Binaries, filename="'+DIR+'/elson.hdf5")', file = fgen)
            fgen.close()


            ## CREATE SUBMIT FILE
            fsub = open(DIR+'/submit_job.sh','w')
            print("#!/bin/bash", file = fsub)
            print("", file = fsub)
            print("#SBATCH -J rv"+str(Rv[jj])+"_n"+str(N[ii])+"_a1"+str(Alpha1[kk]), file = fsub)
            print("#SBATCH --mail-user=shiye2015@u.northwestern.edu", file = fsub)
            print("#SBATCH --error=error.out", file = fsub)
            print("#SBATCH --output=output.out", file = fsub)
            print("#SBATCH --nodes=4", file = fsub)
            print("#SBATCH --ntasks-per-node=28", file = fsub)
            print("#SBATCH --mem=0", file = fsub)
            print("#SBATCH --time=720:00:00", file = fsub)
            print("#SBATCH --account=b1095", file = fsub)
            print("#SBATCH --partition=grail-std", file = fsub)
            print("", file = fsub)
            print("module purge all", file = fsub)
            print("module load cmake/3.15.4", file = fsub)
            print("module load hdf5/1.10.7-openmpi-4.0.5-intel-19.0.5.281", file = fsub)
            print("module load mpi/openmpi-4.0.5-intel-19.0.5.281", file = fsub)
            print("module load gsl/2.5-intel-19.0.5.281", file = fsub)
            print("#########module load list", file = fsub)
            print("", file = fsub)
            print("###mpirun -np X(=n_nodes*n_cores) <exe> > output", file = fsub)
            print("mpirun -np 112 ./cmc ElsonProfile.ini initial", file = fsub)


            ## COPY SUBMIT FILE AND INI FILE
            #os.system('scp /projects/b1095/syr904/cmc/CMC-COSMIC/master_tc_test/ver_0601/MOCHA47Tuc_elson_rv4_3e6_tcon/submit_job.sh '+DIR+'/.')
            os.system('scp /projects/b1095/syr904/cmc/CMC-COSMIC/master_tc_test/ver_0601/MOCHA47Tuc_elson_rv4_3e6_tcon/cmc '+DIR+'/.')
            os.system('scp /projects/b1095/syr904/cmc/CMC-COSMIC/master_tc_test/ver_0601/MOCHA47Tuc_elson_rv4_3e6_tcon/ElsonProfile.ini '+DIR+'/.')

            os.system('python '+DIR+'/generate_elson_profile.py')

