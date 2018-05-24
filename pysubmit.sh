#!/bin/bash

#MSUB -N nskickgrid
#MSUB -M shiye2015@u.northwestern.edu
#MSUB -e pyerror.out
#MSUB -o pyoutput.out
#MSUB -l nodes=2:ppn=6
####MSUB -l naccesspolicy=singlejob
#MSUB -l walltime=2:23:0:0
#MSUB -A b1011
#MSUB -q ligo


cd /projects/b1011/syr904/MyCodes/PythonCode
#python -c 'from ns import get_allid_BP; get_allid_BP("/projects/b1011/sourav/new_runs/kick_grid/rv1/kickscale_1.0", "/projects/b1011/syr904/projects/PULSAR/kickgrid_runs/kickgrid_mspid.dat", "/projects/b1011/syr904/projects/PULSAR/kickgrid_runs/history/23/BP")'
python -c 'from ns import get_id_allmodel_position_10Gyr; get_id_allmodel_position_10Gyr("/projects/b1011/syr904/projects/PULSAR/kickgrid_runs/kickgrid_mspid.dat", "/projects/b1011/syr904/projects/PULSAR/kickgrid_runs/kickgrid_path.dat")'
#python -c 'from ns import get_normalpsr_ini; get_normalpsr_ini("/projects/b1011/syr904/projects/PULSAR/kickgrid_runs/kickgrid_path.dat", "/projects/b1011/syr904/projects/PULSAR/kickgrid_runs/history/normalpsr")'


#module load anaconda

