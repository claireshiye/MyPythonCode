#!/bin/bash

#MSUB -N pyfindpulsar
#MSUB -M shiye2015@u.northwestern.edu
#MSUB -e pyerror.out
#MSUB -o pyoutput.out
#MSUB -l nodes=1:ppn=16
####MSUB -l naccesspolicy=singlejob
#MSUB -l walltime=23:0:0
#MSUB -A b1011
#MSUB -q ligo


cd /projects/b1011/syr904/cmc/cmc_pulsar_2/rundir/kick50_IMF20_fb50
python -c 'from spin import straightline; straightline(sourcedir="/projects/b1011/syr904/cmc/cmc_pulsar_2/rundir/kick50_IMF20_fb50")'

#module load anaconda

