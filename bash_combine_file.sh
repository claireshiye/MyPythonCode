#!/usr/bin/bash

filename='/projects/b1095/syr904/cmc/cmc-mpi-tidalcapture/rvgrid/standard_models_tcon/path_rvgrid_tcon.dat'
while read line; do
# reading each line
echo $line
#cat $line'msp_last.dat' >> '/projects/b1095/syr904/cmc/cmc-mpi-tidalcapture/rvgrid/standard_models_tcon/msp_last_rvgrid_tcon.dat'
#cat $line'normalpsr_last.dat' >> '/projects/b1095/syr904/cmc/cmc-mpi-tidalcapture/rvgrid/standard_models_tcon/normalpsr_last_rvgrid_tcon.dat'
cat $line'NS_MS_last.dat' >> '/projects/b1095/syr904/cmc/cmc-mpi-tidalcapture/rvgrid/standard_models_tcon/NSMS_last_rvgrid_tcon.dat'
cat $line'NS_WD_last.dat' >> '/projects/b1095/syr904/cmc/cmc-mpi-tidalcapture/rvgrid/standard_models_tcon/NSWD_last_rvgrid_tcon.dat'
done < $filename
