#!/bin/bash

#SBATCH -J fb10_nssnap
#SBATCH --mail-user=shiye2015@u.northwestern.edu
#SBATCH --error=pyerror.out
#SBATCH --output=pyoutput.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --mem=100G
#SBATCH --time=72:00:00
#SBATCH --account=b1094
#SBATCH --partition=ciera-std

#export PATH=$PATH:/projects/b1095/syr904/projects/GCE/catalog/

module load python/anaconda3.6


#python -c "import cluster_sampling_v1 as csv1; csv1.main(5000, 2500, 3000, 'cluster_sample_initial_M_RG_dissol0_fcl0.024_ffa0.012_xcut4.997654049537297', 'cluster_sample_disrupt_M_RG_dissol0_fcl0.024_ffa0.012_xcut4.997654049537297', '_v6')"
#python -c "import cluster_sampling_v1 as csv1; csv1.read_property_all(2000, 2500, 'cluster_sample_disrupt_M_RG_dissol1_fcl0.024_ffa0.012_xcut4.997654049537297', 'cluster_sample_property_M_RG_dissol1_fcl0.024_ffa0.012_xcut4.997654049537297', '_v5')"
#python cluster_sampling_v1.py
#python BSS_hdf5.py
#python -c "import GCE; GCE.print_esc_Nns()"

#python -c "import ns_hdf5 as nhf; nhf.print_Nns_snap('/projects/b1095/syr904/cmc/CMC-COSMIC/master_tc_test/ver_0601/MOCHA47Tuc_elson_rv3.5_3e6_tcon_fb10/')"

python compact_object_nave.py

#module load anaconda

