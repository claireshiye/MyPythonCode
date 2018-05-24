import numpy as np
import os
from glob import glob
import subprocess

def create_directories(dirname):
	dataout,dataerr=subprocess.Popen([r"mkdir","-p",dirname],stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()
	print "creating directory:", dirname
	print dataout, dataerr

def copy_files(sourcename, destname):
	dataout,dataerr=subprocess.Popen([r"cp",sourcename,destname],stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()
	print 'copying file:', sourcename, destname
	print dataout, dataerr


def modify_input(sourcename, destname, SIGMAFRAC=1.):
	writefile = open(destname, 'w')
	with open(sourcename, 'r') as f:
		for line in f:
			if line.rfind('BHSIGMAFRAC')>-1:
				writefile.write("BSE_BHSIGMAFRAC=%g\n" %(SIGMAFRAC))
			else:
				writefile.write("%s" %(line))
	writefile.close()

def modify_submit(sourcename, destname, SIGMAFRAC=1.):
	writefile = open(destname, 'w')
	with open(sourcename, 'r') as f:
		for line in f:
			if line.rfind("-N n2w5")>-1:
				writefile.write("#MSUB -N n8w5r2k%g\n" %(SIGMAFRAC))
			else:
				writefile.write("%s" %(line))
	writefile.close()


#kicks = [0.01, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1.]
kicks = [0.005, 0.02, 0.03, 0.04, 0.06, 0.07, 0.08, 0.09,]
basename = '/projects/b1011/sourav/new_runs/kick_grid'
fitsfilename = 'input_n8e5_f5_m0.08_150_w5_rg8.fits'
exfilename = 'cmc'
inputfilename = 'input.cmc'
submitfilename = 'submit_job.sh'

for i in range(len(kicks)):
	dirname = basename+'/kickscale_'+str(kicks[i])
	create_directories(dirname)
	src = basename+'/'+fitsfilename
	dest = dirname+'/'+fitsfilename
	copy_files(src, dest)
	src = basename+'/'+exfilename
	dest = dirname+'/'+exfilename
	copy_files(src,dest)
	src = basename+'/'+inputfilename
	dest = dirname+'/'+inputfilename
	modify_input(src, dest,kicks[i])
	src = basename+'/'+submitfilename
	dest = dirname+'/'+submitfilename
	modify_submit(src, dest,kicks[i])
	

	
		
