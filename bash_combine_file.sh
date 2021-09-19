#!/usr/bin/bash

filepath='/projects/b1095/syr904/cmc/cmc-mpi-tidalcapture/rvgrid/standard_models_tcon/test18_mpi1/'
filename=$filepath'path_text18_tcon_fb20.dat'
while read line; do
# reading each line
echo $line
cat $line'msp_last.dat' >> $filepath'msp_last_fb20.dat'
cat $line'normalpsr_last.dat' >> $filepath'normalpsr_last_fb20.dat'
cat $line'NSMS_last.dat' >> $filepath'NSMS_last_fb20.dat'
cat $line'NSWD_last.dat' >> $filepath'NSWD_last_fb20.dat'
done < $filename
