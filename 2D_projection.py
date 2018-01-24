import numpy as np
from glob import glob
import blackhole as bh
import extract_observed_prop as obs

#path=np.genfromtxt('/projects/b1011/sourav/new_runs/kick_grid/kick_grid_path.dat', dtype='|S')
path=['/projects/b1011/sourav/new_runs/kick_grid/rv1/kickscale_0.005']
#path=['/projects/b1011/sourav/new_runs/kick_grid/rv1/kickscale_0.01','/projects/b1011/sourav/new_runs/kick_grid/rv1/kickscale_0.05','/projects/b1011/sourav/new_runs/kick_grid/rv1/kickscale_0.1','/projects/b1011/sourav/new_runs/kick_grid/rv1/kickscale_0.2']

def conv_dict(): return {'l':15, 't':19}    #?


def conv(unit,filepath):   # Returns the unit conversion multiplier given a simulation's *.conv.sh file and a unit (either 'l' or 't')
    dict = conv_dict()
    from re import findall
    with open(filepath,'r') as f:
        head = [next(f) for x in xrange(24)]
    return float(findall('\d+[\.]?\d*',head[dict[unit]])[0])


def get_time(filepath):      # Returns the cluster's age for a given snapshot
    import gzip
    from re import findall
    with gzip.open(filepath,'r') as f: contents = f.readline()
    if not findall('\d+[\.]?\d*',contents):        # Returns time = 0 for snapshot files without a time header
        print 'snapshot empty'; return float(0)
    else: return float(findall('\d+[\.]?\d*',contents)[0])


def step_of_loop(maxsnapno):
    if maxsnapno < 100: delta=-1
    if maxsnapno>=100 and maxsnapno<300: delta=-2
    if maxsnapno>=300 and maxsnapno<600: delta=-4
    if maxsnapno>=600 and maxsnapno<900: delta=-6
    if maxsnapno>=900: delta=-8

    return delta


def print_obspara(filestring, snapno):
    #units = obs.read_units(filepath+'/'+'initial')
    f_obs_params = open(filestring+'.snap'+snapno+'.obs_params.dat','w')
    props = obs.get_obs_props(filestring, snapno, FAC=1.)
    print>>f_obs_params, '#0.Ltot, 1.M/L, 2.Mave, 3.drc, 4.drhl, 5.dsigmac, 6.dvsigmac_rv, 7.rc, 8.rhl, 9.sigmac, 10.t, 11.vsigmac_rv, 12.M_total, 13.N_BH, 14.N_BHsingle, 15.N_BHBH, 16.N_BHnBH, 17.N_BHNS, 18.N_BHWD, 19.N_BHMS, 20.N_BHG'
    print>>f_obs_params, props['Ltot'], props['M/L'], props['Mave'], props['drc'], props['drhl'], props['dsigmac'], props['dvsigmac_rv'], props['rc'], props['rhl'], props['sigmac'], props['t'], props['vsigmac_rv'], M_total, N_BH, N_BHsingle, N_BHBH, N_BHnBH, N_BHNS, N_BHWD, N_BHMS, N_BHG
    print>>f_obs_params, props['Ltot'], props['M/L'], props['Mave'], props['drc'], props['drhl'], props['dsigmac'], props['dvsigmac_rv'], props['rc'], props['rhl'], props['sigmac'], props['t'], props['vsigmac_rv'], M_total, N_BH, N_BHsingle, N_BHBH, N_BHnBH, N_BHNS, N_BHWD, N_BHMS, N_BHG



for k in range(len(path)):
	snaps=np.sort(glob(path[k]+'/'+'initial.snap*.dat.gz'))
    snapno_max=len(snaps)-1
    step=int(step_of_loop(snapno_max))
    filestring = path[k]+'/initial'
    units = obs.read_units(filestring)
	for n in range(len(snaps), 0, step):
		#t_conv=conv('t',path[k]+'/'+'initial.conv.sh')
    	#time=get_time(snaps[n])*t_conv
    	#if time >= 10000.0 and time <=12000.0:
    	snapno=str(n).zfill(4)
		obs.make_2D_projection(filestring, snapno, units, SEEDY=100, PROJ=(0,1))
        print_obspara(filestring, snapno)
		#bh.get_sbp_from_2D_projection(filestring, snapno, BINNO=50, LCUT=15)

    if n > 0:
        snapno=str(0).zfill(4)
        obs.make_2D_projection(filestring, snapno, units, SEEDY=100, PROJ=(0,1))
        print_obspara(filestring, snapno)
	#print k

