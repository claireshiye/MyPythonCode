import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import LISA_calculations2 as LISA
import LISA_calculations as LISA_original
import scipy.constants as ct
from scipy.interpolate import interp1d
import ecc_calc as LISA2


G = 6.67e-11
Msun = 1.99e30
M_sun = Msun
AU = 1.5e11
d = 16.0
c = 3.0e8

#factor = np.sqrt(32./5.)
factor = 1.0

z = LISA.zAtLuminosityDistance(d)

ax1 = plt.subplot(111) # strain
#ax2 = plt.subplot(212) # t_insp

def get_roche_r(q, a,e):
        R = a*(1-e)*0.49*q**(2./3.) / (0.6*q**(2./3.) + np.log(1+q**(1./3.)))
	return R

def calc_fGW(m1,m2,a,e,z=0):
        return np.sqrt(G*Msun*(m1+m2))/np.pi*(1+e)**1.1954/(a*AU*(1-e**2.))**1.5/(1.+z)

def find_n_peak(m1,m2,a,e):
        Porb = np.sqrt(4*np.pi**2/(G*Msun*(m1+m2))*(a*AU)**3.)
        forb = 1./Porb
        fGW = calc_fGW(m1,m2,a,e,z=0)
        n_peak = fGW/forb
        #n_peak = int(fGW/forb)
        return n_peak

def DECIGO_S(f):
	fp = 7.36 #Hz
	Sh = 7.05*10**(-48)*(1.+(f/fp)**2) + 4.8*10**(-51)* f**(-4)/(1.+(f/fp)**2) + 5.33*10**(-52)*f**(-4)
	return np.sqrt(Sh*f)

f_dec = np.logspace(-3,1,1000)
h_dec = DECIGO_S(f_dec)

data = np.genfromtxt('binaries.dat')

data2 = np.genfromtxt('characteristic_noise_strain.dat')
fn = data2[:,0]
hn_lisa = data2[:,1]

data3 = np.genfromtxt('LIGO_sensitivity_curve.dat')

#3173.05782808 14.0 14.0 25.244361 25.34543 1.1271514 0.95188347 273577.0 648641.0 23.0 5133.58511343 3.0

f_GW_array = []
strain_array = []
count = 0
noise_interp = interp1d(fn, hn_lisa, bounds_error=False, fill_value=1e60)
T_OBS = 4


d_arr = [10,800,16000,1000000]
#d_label = [r'$10\,\rm{kpc}$',r'$800\,\rm{kpc}$',r'$16\,\rm{Mpc}$',r'$1\,\rm{Gpc}$']
d_label = [r'$M_{\rm{BH}}=100\,M_{\odot}$',r'$M_{\rm{BH}}=1000\,M_{\odot}$',r'$M_{\rm{BH}}=10^4\,M_{\odot}$']
mBH_arr = [100,1000,10000]
r_tc_fac = [3.03, 3.37, 3.69]
color_array = ['black','steelblue','forestgreen','gold']


for i in range(len(mBH_arr)):
	for k in range(0,2):
		color = color_array[i]
		label = d_label[i]

		#0.0 14.0 0.001 25.0 1.0 0.7 2 0.1
		f_GW_array = []
		strain_array = []
		strain_array_0 = []

		f_GW_array_RL = []
		strain_array_RL = []

		#d = data[i,8]
		d = 10./1000. #Mpc
		z = LISA.zAtLuminosityDistance(d)
		m1 = 0.6 # Mstar
		#m2 = 1000.0 # Mbh
		m2 = mBH_arr[i]
		#a0 = 0.02
		Radius = 0.012 #stellar radius in solar radii
		r_TDE = (m2/m1)**(1./3.)*Radius

		if k ==0:
			e0 = .99
			a0 = 3.*r_TDE/(1-e0)/215.  # initial SMA in AU
			q = m1/m2
			fGW = calc_fGW(m1,m2,a0,e0,z=z)
			LINESTYLE = '-'
		if k == 1:
			e0 = 0.01
			a0 = 1.0
			q = m1/m2
			#print 'CHECK', np.sqrt(G*Msun*(m1+m2))/np.pi/(a0*AU)**1.5
                        fGW = calc_fGW(m1,m2,a0,e0,z=z)
			LINESTYLE = '--'
		print 'd =', d, m1, m2, 'rTDE =', r_TDE, 'Radius = ', Radius, 'a0 =', a0, e0, 'fGW =',fGW
		if fGW > 5:
			continue

		t, a, e = LISA.t_inspiral_2(a0,e0,m1,m2,t_flag=0,array=1,LIGO=2)

		t_insp = LISA.t_inspiral_2(a0,e0,m1,m2)/1.e6 # In Myr

		delta = 1
		
		LIGO_limit = 5.*2.*G/c**2.*(m1+m2)*Msun
		# Inspiral limit is 5 times sum of schwarschild radii

		r_TD = (m2/m1)**(1./3.)*Radius/215.0 # TD radius in AU

		print i,m1, m2, Radius, r_TD, t_insp
		RL_flag = 0

		for j in range(0,len(t),delta):
			a_temp = a[j]
			e_temp = e[j]

			#if a_temp*AU*(1-e_temp) > LIGO_limit:
			fGW = calc_fGW(m1,m2,a_temp,e_temp,z=z)

			r_peri = a_temp*(1-e_temp)
			RL = get_roche_r(q, a_temp,0)
			#print Radius, r_peri, r_TD
			#print a_temp, e_temp, r_peri, RL, r_TD
			if r_peri > r_TDE/215.:
				n_peak = find_n_peak(m1,m2,a_temp, e_temp)
				#strain = LISA.hn(m1,m2,a_temp, e_temp, d, n_peak)
				#############
				m1_geo = m1*M_sun*ct.G/ct.c**2.
				m2_geo = m2*M_sun*ct.G/ct.c**2.
				a_geo = a[j]*AU
				f_orb_geo = np.sqrt((m1_geo+m2_geo)/(4.*np.pi**2.*a_geo**3.))

				strain = LISA2.hcn_func(m1_geo, m2_geo, z, n_peak, n_peak*f_orb_geo, e[j])
				### Find stationary strain and choose minumum
				Porb = (4.*np.pi**2./(G*Msun*(m1+m2))*(a_temp*AU)**3.)**0.5
				f_orb = 1./Porb
				Tobs = T_OBS*3.15e7
				strain_0 = strain*np.sqrt(LISA_original.fdotn_2(m1,m2,Porb,e_temp,n_peak,z)*Tobs/fGW)
				if RL > Radius/215.:     # cut off when r_peri is less than r_TD
					strain_array.append(factor*np.min([strain,strain_0]))
					f_GW_array.append(fGW)
					a_RL = a_temp
					e_RL = e_temp
				else: #star is filling RL				
					strain_array_RL.append(factor*np.min([strain,strain_0]))
					f_GW_array_RL.append(fGW)
				#################

			else:
				a_temp = a[j-1]
				e_temp = e[j-1]
				#print 'reached LIGO limit!'
				break

		#######

		m1 = m1*Msun
		m2 = m2*Msun
		SNR = LISA2.snr(m1, m2, a0*1.5e11, e0, d, 50, noise_interp, T_OBS)
		SNR_RL = LISA2.snr(m1, m2, a_RL*1.5e11, e_RL, d_arr[i]/1000., 50, noise_interp, T_OBS)
		SNR_stationary = LISA.SNR_calculator(m1,m2,a_RL,e_RL,d_arr[i]/1000.)
		forb = calc_fGW(m1/Msun,m2/Msun,a_RL,e_RL,z=0)
		Mc = pow(m2*m1,3./5.)/pow(m2+m1,1./5.)
		h_o = G/pow(c,2.) * Mc/(d*1e6*3.086e16) * pow(G/pow(c,3.)*np.pi*forb*Mc,2./3.)
		h_f = np.interp(forb, fn, hn_lisa)
		SNR2 =  h_o*np.sqrt(T_OBS*3.15e7*forb)/h_f

		print 'd =', d, forb, 'SNR = ',SNR
		if k == 0:
			ax1.plot(f_GW_array,strain_array, color=color,lw=2,alpha=1,ls=LINESTYLE,label=label)
			ax1.plot(f_GW_array_RL,strain_array_RL, color=color,lw=2,ls=LINESTYLE,alpha=1)
		if k ==1:
			ax1.plot(f_GW_array,strain_array, color=color,lw=2,alpha=1,ls=LINESTYLE)
                        ax1.plot(f_GW_array_RL,strain_array_RL, color=color,lw=2,ls=LINESTYLE,alpha=1)


ax1.plot(data2[:,0], data2[:,1],color='black',lw=2)
#ax1.plot(f_dec, h_dec,color='darkgray',lw=2,ls='-')
#ax1.plot(data3[:,0], data3[:,5],color='gray',lw=3)
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlim(1e-5,3)
ax1.set_ylim(5e-23,1e-15)
#ax1.set_ylim(5e-26,1e-15)
ax1.set_xlabel(r'$f_{\rm{GW}}\,(\rm{Hz})$',fontsize=20)
ax1.set_ylabel(r'$\rm{Characteristic\,Strain}$',fontsize=20)
#plt.setp(ax1.get_xticklabels(), visible=False)
plt.grid(lw=0.5,color='gray',ls='-',zorder=0)

ax1.legend(loc=3,frameon=True, prop={'size': 16},scatterpoints=1)

ax1.xaxis.set_tick_params(width=1,length=6,which='major')
ax1.xaxis.set_tick_params(width=1,length=3.5,which='minor')
ax1.yaxis.set_tick_params(width=1,length=6, which='major')
ax1.yaxis.set_tick_params(width=1,length=3.5, which='minor')

for tick in ax1.xaxis.get_major_ticks():
	tick.label.set_fontsize(14)
for tick in ax1.yaxis.get_major_ticks():
	tick.label.set_fontsize(14)

plt.subplots_adjust(wspace=0, hspace=0.02)
plt.tight_layout()
plt.savefig('BHWD.pdf')
#plt.show()


