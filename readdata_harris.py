import pandas as pd
import xlrd as xl
import numpy as np
import conversions as conv


def read_mclaughlin_tables():
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

	np.savetxt('/Users/shiye/Documents/ClusterGroup/MWGC.dat', np.c_[t, RC, RH], fmt='%f %f %f', header='	1.Age(Gyr), 2.Rc(pc), 3. Rh(pc)', comments='#')
	
	table8=pd.read_excel('/Users/shiye/Documents/ClusterGroup/table8_mclaughlin.xlsx')
	table11=pd.read_excel('/Users/shiye/Documents/ClusterGroup/table11_mclaughlin.xlsx')
	table12=pd.read_excel('/Users/shiye/Documents/ClusterGroup/table12_mclaughlin.xlsx')
	logt=table8['logage(yr)'].values
	logrc=table11['logrc(pc)'].values
	logrh=table11['logrh(pc)'].values
	logtrh=table12['logtrh(yr)'].values
	print logrc, logrh, logt, logtrh
	
	lrc=[]; lrh=[]; lage=[]; ltrh=[]
	for o in range(0, len(logrc), 3):
		lrc.append(float(logrc[o])); lrh.append(float(logrh[o])); ltrh.append(float(logtrh[o]))
	
	for p in range(0, len(logt), 2):
		lage.append(float(logt[p]))
	
	tlmc=[]; rclmc=[]; rhlmc=[]; trhlmc=[]
	for i in range(len(lage)):
		tlmc.append(10**float(lage[i])/10.**9); rclmc.append(10**float(lrc[i])); rhlmc.append(10**float(lrh[i])); trhlmc.append(10**ltrh[i]/10.**9)
#np.savetxt('/Users/shiye/Documents/ClusterGroup/LMC.dat', np.c_[tlmc, rclmc, rhlmc, trhlmc], fmt='%f %f %f %f', header='1.Age(Gyr), 2.Rc(pc), 3. Rh(pc), 4.trh(Gyr)', comments='#')



def read_harris_catalog():
	df_part3=pd.read_excel('/projects/b1095/syr904/projects/harris_catalog/harris_part3.xlsx')
	df_part1=pd.read_excel('/projects/b1095/syr904/projects/harris_catalog/harris_part1.xlsx')
	#print df_part3.columns

	rc=df_part3['r_c'].values
	rhl=df_part3['r_h'].values
	n_psr=df_part3['Npsr(Freire)'].values
	mass=df_part3['Mass(Msun)(Baumgardt_2018)'].values

	r_sun=df_part1['R_Sun'].values

	M=[]
	for k in range(len(mass)):
		if pd.isna(mass[k]):
			continue
		else:
			x = mass[k].replace(u'\xa0\xb1\xa0', u' ')
			y = x.replace(u'\xd7', u'')
			z = y.replace(u'\xa0', u' ').split()
			vs=float(z[0])
			exps=float(z[-1][-1])
			#print exps
			M.append(vs*10**exps)

	#print M
	#print rc, rhl, r_sun
	#print len(rc), len(rhl), len(r_sun)
	R_c=[]; R_hl=[]; N_psr=[]
	for i in range(len(rc)):
		if pd.isna(rc[i]) or pd.isna(rhl[i]) or pd.isna(r_sun[i]):
			continue
		else:
			R_c.append(conv.arcmin_to_pc(rc[i], r_sun[i]))
			R_hl.append(conv.arcmin_to_pc(rhl[i], r_sun[i]))
			N_psr.append(n_psr[i])


	return R_c, R_hl, N_psr, M


	






