import numpy as np
from glob import glob
import blackhole as bh
import extract_observed_prop as obs

path=np.genfromtxt('/projects/b1011/syr904/projects/PULSAR/kickgrid_runs/path_newmodel.dat', dtype='|S')
#path=['/projects/b1011/syr904/cmc/cmc-mpi-06/rundir/3.5e6rv1fb10kick1.0']
#path=['/projects/b1011/syr904/cmc/cmc-mpi-09/rundir/kickgrid/kickgrid_0.11', '/projects/b1011/syr904/cmc/cmc-mpi-09/rundir/kickgrid/kickgrid_0.13']

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


def print_obspara(filestring, snapno, tc):
    #units = obs.read_units(filepath+'/'+'initial')
    f_obs_params = open(filestring+'.snap'+snapno+'.obs_params.dat','w')
    M_total, props = obs.get_obs_props(filestring, snapno, FAC=1.)
    databh=np.genfromtxt(filestring+'.bh.dat')
    for i in range(len(databh[:,1])):
        if databh[:,1][i]==tc:
		N_BH=float(databh[:,2][i]); N_BHsingle=float(databh[:,3][i]); N_BHBH=float(databh[:,5][i]); N_BHnBH=float(databh[:,6][i]); N_BHNS=float(databh[:,7][i]); N_BHWD=float(databh[:,8][i]); N_BHMS=float(databh[:,10][i]); N_BHG=float(databh[:,11][i])
    print>>f_obs_params, '#0.Ltot, 1.M/L, 2.Mave, 3.drc, 4.drhl, 5.dsigmac, 6.dvsigmac_rv, 7.rc, 8.rhl, 9.sigmac, 10.t, 11.vsigmac_rv, 12.M_total, 13.N_BH, 14.N_BHsingle, 15.N_BHBH, 16.N_BHnBH, 17.N_BHNS, 18.N_BHWD, 19.N_BHMS, 20.N_BHG'
    print>>f_obs_params, props['Ltot'], props['M/L'], props['Mave'], props['drc'], props['drhl'], props['dsigmac'], props['dvsigmac_rv'], props['rc'], props['rhl'], props['sigmac'], props['t'], props['vsigmac_rv'], M_total, N_BH, N_BHsingle, N_BHBH, N_BHnBH, N_BHNS, N_BHWD, N_BHMS, N_BHG
    print>>f_obs_params, props['Ltot'], props['M/L'], props['Mave'], props['drc'], props['drhl'], props['dsigmac'], props['dvsigmac_rv'], props['rc'], props['rhl'], props['sigmac'], props['t'], props['vsigmac_rv'], M_total, N_BH, N_BHsingle, N_BHBH, N_BHnBH, N_BHNS, N_BHWD, N_BHMS, N_BHG


for k in range(len(path)):   ### For all the snapshots in a model
    snaps=np.sort(glob(path[k]+'/'+'initial.snap*.dat.gz'))
    snapno_max=len(snaps)-1
    step=int(step_of_loop(snapno_max))
    filestring = path[k]+'/initial'
    units = obs.read_units(filestring)
    for n in range(len(snaps)-1, 0, step):
        t_conv=conv('t',path[k]+'/'+'initial.conv.sh')
        tcode=get_time(snaps[n])
        time=tcode*t_conv
        if time >= 10000.0:
            snapno=str(n).zfill(4)
            obs.make_2D_projection(filestring, snapno, units, SEEDY=100, PROJ=(0,1))
            print_obspara(filestring, snapno, tcode)
        #if n >= 140:
            bh.get_sbp_from_2D_projection(filestring, snapno, BINNO=50, LCUT=15)

    print k

    #if n > 0:
        #snapno=str(0).zfill(4)
	#tcode=get_time(snaps[0])
        #obs.make_2D_projection(filestring, snapno, units, SEEDY=100, PROJ=(0,1))
        #print_obspara(filestring, snapno, tcode)
	#print k


#for k in range(len(path)):   ###For only the last snapshot
#    snaps=np.sort(glob(path[k]+'/'+'initial.snap*.dat.gz'))
#    snapno_max=len(snaps)-1
#    filestring = path[k]+'/initial'
#    units = obs.read_units(filestring)
#    t_conv=conv('t',path[k]+'/'+'initial.conv.sh')
#    tcode=get_time(snaps[snapno_max])
#    time=tcode*t_conv
#    snapno=str(int(snapno_max)).zfill(4)
#    obs.make_2D_projection(filestring, snapno, units, SEEDY=100, PROJ=(0,1))
#    print_obspara(filestring, snapno, tcode)
#    bh.get_sbp_from_2D_projection(filestring, snapno, BINNO=50, LCUT=15)
#    print k
