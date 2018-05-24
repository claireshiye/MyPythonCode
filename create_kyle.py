import numpy as np
import constants
import os
import sys

# Set vg and rg for tidal radius calculation
vg=220.
rg=8000.
BSE_BHFLAG = '4'

rv = 1
rv_str = str(int(rv))
n = 8
N = int(n*1.e5)
#N = int(i/2.*100000)
n_str = str(n)
N_str = str(N)
wo = 5.0
wo_str = str(int(wo))
#kickscale_array = [0.005,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.4,0.6,0.8,1.0]
kickscale_array = [0.17]
for j in range(0,1):#(8,12):
	kickscale = kickscale_array[j]
	kickscale_str = str(kickscale)
	BHSIGMAFRAC = kickscale_str
	binfrac = 0.05
	binfrac_str = str(binfrac)

	# CREATE FOLDER
	DIR = 'kickscale_'+kickscale_str
	os.system('mkdir '+DIR)
	# CREATE GENFITS.SH FILE
	f = open('genfits.sh','w')
	print>>f, '../../bin/cmc_mkking -s 100 -w '+wo_str+' -N '+N_str+' -o test1.fits'
	print>>f, '../../bin/cmc_setimf -i test1.fits -o test2.fits -R 100 -I 0 -m 0.08 -M 150 && rm test1.fits'
	print>>f, '../../bin/cmc_setunits -i test2.fits -o test3.fits -R '+rv_str+' -Z 0.001 && rm test2.fits'
	print>>f, '../../bin/cmc_setstellar -i test3.fits -o test4.fits && rm test3.fits'
	bin_num = str(binfrac*N)
	print>>f, '../../bin/cmc_addbinaries -i test4.fits -o input.fits -s 100 -N '+bin_num+' -l 0 -b 0 && rm test4.fits'
	f.close()
	# COPY SUBMIT FILE AND CMC FILE
	os.system('scp submit_job.sh '+DIR+'/.')
	os.system('scp cmc '+DIR+'/.')
	# RUN GENFITS.SH
	os.system('sh genfits.sh 1>genfits.out 2>genfits.err')
	# MOVE INPUT.FITS TO THE CORRECT DIRECTORY
	os.system('scp input.fits '+DIR+'/.')
	# READ IN THE INITIAL MASS
	f2 = open('genfits.err','r')
	lines = f2.readlines()
	line = lines[2]
	parsed = line.split('=')
	mc = float(parsed[1])

	# CALCULATE THE TIDAL RADIUS
	#r_tidal = (constants.G * mc*constants.Msun / 2. / (vg*constants.km)**2.)**(1./3.) * (rg*constants.PC)**(2./3.) / constants.PC
	r_tidal = 111.04
	r_tidal_str = str(r_tidal)

	# CREATE THE INPUT.CMC FILE
	f3 = open(DIR+'/input.cmc','w')
	print>>f3,'# required parameters'
	print>>f3,'GAMMA = 0.01'
	print>>f3,'INPUT_FILE = input.fits'
	print>>f3,'MASS_PC = 0.0001,0.0003,0.0005,0.0007,0.0009,0.001,0.003,0.005,0.007,0.009,0.01,0.03,0.05,0.07,0.09,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.99'
	print>>f3,'MASS_BINS = 0.1,1.0,10.0,100.0,1000.0 \n'
	print>>f3,'# additional parameters'
	print>>f3,'OVERWRITE_RVIR = '+rv_str
	print>>f3,'OVERWRITE_RTID = '+r_tidal_str+'/'+rv_str
	print>>f3,'OVERWRITE_Z = 0.001'
	print>>f3,'THETASEMAX = 1.4142'
	print>>f3,'BH_AVEKERNEL = 3'
	print>>f3,'T_MAX_PHYS = 12.0'
	print>>f3,'WRITE_PULSAR_INFO = 1'
	print>>f3,'T_MAX_COUNT = 1000000000'
	print>>f3,'MIN_CHUNK_SIZE = 40'
	print>>f3,'SS_COLLISION = 1'
	print>>f3,'STELLAR_EVOLUTION = 1'
	print>>f3,'BINBIN = 1'
	print>>f3,'BINSINGLE = 1'
	print>>f3,'THREEBODYBINARIES = 1'
	print>>f3,'MIN_BINARY_HARDNESS = 1.0'
	print>>f3,'ONLY_FORM_BH_THREEBODYBINARIES = 0'
	print>>f3,'IDUM = 1709217'
	print>>f3,'SNAPSHOTTING = 1'
	print>>f3,'SNAPSHOT_DELTACOUNT = 7000   #IF TOO MANY SNAPSHOTS INCREASE THIS'
	print>>f3,'WRITE_STELLAR_INFO = 0'
	print>>f3,'STOPATCORECOLLAPSE = 0'
	print>>f3,'WRITE_EXTRA_CORE_INFO = 1'
	print>>f3,'TIDAL_TREATMENT = 1'
	print>>f3,'WRITE_BH_INFO = 1'
	print>>f3,'BH_SNAPSHOTTING = 1'
	print>>f3,'BH_SNAPSHOT_DELTACOUNT = 5000  #IF TOO MANY SNAPSHOTS INCREASE THIS'
	print>>f3,'BSE_BHSIGMAFRAC = 1'
	print>>f3,'BSE_FBKICKSWITCH = 1'
	print>>f3,'TERMINAL_ENERGY_DISPLACEMENT = 100000 \n'
	print>>f3,'#BSE Parameters'
	print>>f3,'BSE_MXNS=2.5'
	print>>f3,'BSE_BHFLAG='+BSE_BHFLAG
	if BSE_BHFLAG == '4':
		print>>f3,'BSE_BHSIGMAFRAC='+BHSIGMAFRAC
	print>>f3,'BSE_NSFLAG=3  # Rapid SN model from Fryer 2012'
	print>>f3,'BSE_WINDFLAG=3 # Vink winds plus LBV winds for all stars past WD limit'
	print>>f3,'BSE_ALPHA1=1.0'
	print>>f3,'BSE_PPSN=1'
	print>>f3,'#'
	print>>f3,'#Checkpoint Parameters'
	print>>f3,'MAX_WCLOCK_TIME=3456000'
	print>>f3,'CHECKPOINT_INTERVAL=43200'
	print>>f3,'CHECKPOINTS_TO_KEEP=2'
	print>>f3,'BSE_BHSPINFLAG=0'
	print>>f3,'TERMINAL_ENERGY_DISPLACEMENT=50000'

	# DELETE THE GENFITS FILES
	os.system('rm genfits.sh genfits.out genfits.err input.fits')

