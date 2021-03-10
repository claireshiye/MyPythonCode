#!/usr/bin/bash

filename='/projects/b1095/syr904/cmc/cmc-mpi-tidalcapture/rvgrid/path_rvgrid.dat'
while read line; do
# reading each line
echo $line
scp syr904@quest.it.northwestern.edu:$line"msp_last.dat" /Users/shiye/Documents/Research/tidal_capture/rvgrid
scp syr904@quest.it.northwestern.edu:$line"normalpsr_last.dat" /Users/shiye/Documents/Research/tidal_capture/rvgrid
done < $filename
