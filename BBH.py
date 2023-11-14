import numpy as np

savepath = '/fs/lustre/cita/claireshiye/projects/BBH_catalog/'

m0 = []; m1 = []; bbh_type = []; t_mer = []

bbh_col = np.genfromtxt(savepath+'GWcap_BBH_maingrid.dat')
three_col = bbh_col[:,12]
bbh_mer = np.genfromtxt(savepath+'Incluster_BBH_maingrid.dat')
bbh_esc = np.genfromtxt(savepath+'Esc_BBH_maingrid.dat')
t_merger = bbh_esc[:,2]+bbh_esc[:,3]

m0 = m0+list(bbh_col[:,10][three_col==-100])
m1 = m1+list(bbh_col[:,11][three_col==-100])
bbh_type = bbh_type+list(np.full_like(bbh_col[:,10][three_col==-100], 1))
t_mer = t_mer + list(bbh_col[:,2][three_col==-100])

m0 = m0+list(bbh_mer[:,7])
m1 = m1+list(bbh_mer[:,8])
bbh_type = bbh_type+list(np.full_like(bbh_mer[:,7], 2))
t_mer = t_mer + list(bbh_mer[:,2])

m0=m0+list(bbh_esc[:,4][t_merger<=14000.])
m1=m1+list(bbh_esc[:,5][t_merger<=14000.])
bbh_type = bbh_type+list(np.full_like(bbh_esc[:,4][t_merger<=14000.], 3))
t_mer = t_mer + list(t_merger[t_merger<=14000.])

M0 = np.maximum(m0,m1)
M1 = np.minimum(m0,m1)
m0 = np.array(m0); m1 = np.array(m1)
mtot = m0+m1
t_mer = np.array(t_mer)