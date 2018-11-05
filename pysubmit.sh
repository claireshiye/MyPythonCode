#!/bin/bash

#MSUB -N enc26
#MSUB -M shiye2015@u.northwestern.edu
#MSUB -e pyerror.out
#MSUB -o pyoutput.out
#MSUB -l nodes=2:ppn=6
####MSUB -l naccesspolicy=singlejob
#MSUB -l walltime=5:23:0:0
#MSUB -A b1011
#MSUB -q ligo


cd /projects/b1011/syr904/MyCodes/PythonCode
python -c 'from ns import get_nenc; get_nenc("/projects/b1011/syr904/projects/PULSAR/kickgrid_runs/paper/data/path_newmodel.dat", 25, 26)'


#module load anaconda

