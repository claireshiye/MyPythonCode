import pandas as pd
import xlrd as xl
import numpy as np
import conversions as conv

df=pd.read_excel('/Users/shiye/Documents/ClusterGroup/MWGC.xlsx')
#print df.columns

d_sun=df['d_sun[kpc]'].values
r_c=df['rc[arcsec]'].values
r_h=df['rh[arcsec]'].values
age=df['age[gyr]'].values

d_sun=np.delete(d_sun, [0])
r_c=np.delete(r_c, [0])
r_h=np.delete(r_h,[0])
age=np.delete(age, [0])

ds=[]; rc=[]; rh=[]; t=[]
for k in range(len(d_sun)):
	x=float(d_sun[k]); y=float(r_c[k]); z=float(r_h[k]); time=float(age[k])
	if x<50.0:
		if time >-100.: 
			ds.append(x); rc.append(y); rh.append(z); t.append(time)

#print len(ds), ds, rc, rh, t

RC=[]; RH=[]
for j in range(len(ds)):
	rcpc=conv.arcsec_to_pc(rc[j], ds[j])
	RC.append(rcpc)
	rhpc=conv.arcsec_to_pc(rh[j], ds[j])
	RH.append(rhpc)

np.savetxt('/Users/shiye/Documents/ClusterGroup/MWGC.dat', np.c_[t, RC, RH], fmt='%f %f %f', header='1.Age(Gyr), 2.Rc(pc), 3. Rh(pc)', comments='#')

#datalmc=np.genfromtxt('/Users/shiye/Documents/ClusterGroup/lmc.dat')
#logt=datalmc[:,2]
#logrc=datalmc[:,0]
#logrh=datalmc[:,1]
dflmc=pd.read_excel('/Users/shiye/Documents/ClusterGroup/lmc_mclaughlin.xlsx')
logt=dflmc['Age'].values
logrc=dflmc['log(rc)'].values
logrh=dflmc['log(rh)'].values


tlmc=[]; rclmc=[]; rhlmc=[]
for i in range(len(logt)):
	tlmc.append(10**float(logt[i])/10.**9); rclmc.append(10**float(logrc[i])); rhlmc.append(10**float(logrh[i]))
np.savetxt('/Users/shiye/Documents/ClusterGroup/LMC.dat', np.c_[tlmc, rclmc, rhlmc], fmt='%f %f %f', header='1.Age(Gyr), 2.Rc(pc), 3. Rh(pc)', comments='#')
