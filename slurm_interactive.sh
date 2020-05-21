#!/bin/bash

#SBATCH --job-name=interactive
#SBATCH --mail-user=shiye2015@u.northwestern.edu
#SBATCH --error=pyerror.out
#SBATCH --output=pyoutput.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=12:00:00
#SBATCH --account=b1095
#SBATCH --partition=grail-std
#SBATCH --mem-per-cpu=1G

module purge all
module load python
