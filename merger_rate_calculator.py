import numpy as np
from glob import glob
import collections
from collections import Counter
import os,sys
import subprocess
import gzip
import math
import re
import random
import conversions
import unit_convert as uc
import dynamics as dyn
import scripts
import history_cmc

##Assign weights to the models
def get_weights():
    # Define bins: 108 total bins (3*3*3*4)
    mass_array = [0,5e4,1e5,2e5,7e5]
    rcrh_array = [0,0.4,0.8,1.3]
    rgc_array = [0,5,14,100]
    z_array = [0,0.0005,0.005,0.03]
    #rgc_array = [0,100]
    #z_array = [0,0.03]

    ##Read Data
    data_harris = np.genfromtxt('/projects/b1095/kylekremer/python_code/CMC_Grid_March2019/harris.dat')
    R_sun_obs = data_harris[:,2]; M_bolo = data_harris[:,4]-0.107
    rc_obs = data_harris[:,7]*(R_sun_obs*1000)/3437.75; rh_obs = data_harris[:,8]*(R_sun_obs*1000)/3437.75
    rgc_obs = data_harris[:,3]; z_obs = conversions.metallicity(data_harris[:,5],'fe/htoz')
    mass_obs = 2.0*10**(0.4*(4.74-M_bolo)); rcrh_obs = rc_obs/rh_obs
    #print rc_obs, rh_obs, rcrh_obs

    data_model = np.genfromtxt('/projects/b1095/syr904/projects/SGRB/newruns/finaldata/clusterproperty_12Gyr_maingrid.dat')
    data_path=np.genfromtxt('/projects/b1095/syr904/projects/SGRB/newruns/finaldata/path_allfinished_newruns_maingrid.dat', dtype=str)
    paths=data_path[:,0]
    mass_mod=data_model[:,2]; status_mod=data_model[:,12]; rc_mod=data_model[:,5]; rh_mod=data_model[:,6]
    rcrh_mod=rc_mod/rh_mod
    #print len(rcrh_mod)

    rgc_mod=[]; z_mod=[]
    for o in range(len(paths)):
        s=paths[o].split('/')
        n_star=float(s[-2]); zredshift=float(s[-3][1:]); rg=int(s[-4][2:]); rv=float(s[-5][2:])
        rgc_mod.append(rg); z_mod.append(zredshift)


    ##If use Kyle's /projects/b1095/kylekremer/python_code/CMC_Grid_March2019/cluster_params_12gyr.dat file
    #mass_mod = data_model[:,6]; rc_mod = data_model[:,7]; rh_mod = data_model[:,8]
    #rcrh_mod = rc_mod/rh_mod
    #rgc_mod = data_model[:,4]; z_mod = data_model[:,5]


    ##Calculating the bins
    obs_bin_count=[]
    mod_bin_count=[]
    for i in range(len(mass_array)-1):
        for j in range(len(rcrh_array)-1):
            for k in range(len(rgc_array)-1):
                for l in range(len(z_array)-1):
                    obscount=0; modcount=0
                    for x in range(len(mass_obs)):
                        if mass_array[i]<=mass_obs[x]<mass_array[i+1] and rcrh_array[j]<=rcrh_obs[x]<rcrh_array[j+1] and rgc_array[k]<=rgc_obs[x]<rgc_array[k+1] and z_array[l]<=z_obs[x]<z_array[l+1]:
                            obscount+=1

                    obs_bin_count.append(obscount)
                
                    for y in range(len(mass_mod)):
                        if mass_array[i]<=mass_mod[y]<mass_array[i+1] and rcrh_array[j]<=rcrh_mod[y]<rcrh_array[j+1] and rgc_array[k]<=rgc_mod[y]<rgc_array[k+1] and z_array[l]<=z_mod[y]<z_array[l+1]:
                            if int(status_mod[y])==1:
                                modcount+=1

                    mod_bin_count.append(modcount)

    #print len(obs_bin_count), len(mod_bin_count)
    #print np.sum(obs_bin_count), np.sum(mod_bin_count)
    #print mod_bin_count

    obs_total=0.
    for p in range(len(obs_bin_count)):
        if mod_bin_count[p]!=0:
            obs_total+=obs_bin_count[p]

    #print obs_total

    ##Calculating the weights
    mod_weights=[]
    for m in range(len(mass_mod)):
        theweight=0.
        for i in range(len(mass_array)-1):
            for j in range(len(rcrh_array)-1):
                for k in range(len(rgc_array)-1):
                    for l in range(len(z_array)-1):
                        if mass_array[i]<=mass_mod[m]<mass_array[i+1] and rcrh_array[j]<=rcrh_mod[m]<rcrh_array[j+1] and rgc_array[k]<=rgc_mod[m]<rgc_array[k+1] and z_array[l]<=z_mod[m]<z_array[l+1] and int(status_mod[m])==1:
                            binno=int(i*27+j*9+k*3+l*1)
                            #print binno
                            theweight=float(obs_bin_count[binno])/float(mod_bin_count[binno])/obs_total
        
        mod_weights.append(theweight)

    #print len(mod_weights)

    nomod_weights=[]
    for q in range(len(status_mod)):
        if status_mod[q]==1: nomod_weights.append(1)
        else: nomod_weights.append(0)

    return mod_weights, nomod_weights
    #print np.count_nonzero(mod_weights), np.sum(mod_weights), mod_weights, 
    #return np.where(np.array(mod_weights)==0)[0]




##Generate age distribution according to metallicity
def make_metallicity_bins():
    datamtl=np.genfromtxt('/projects/b1095/syr904/projects/SGRB/AgeDistribution/Mvir1e14.txt')
    ages=datamtl[:,0]; FeH=datamtl[:,1]
    zmetal=[]

    z1=[]; z2=[]; z3=[]
    a1=[]; a2=[]; a3=[]
    f1=[]; f2=[]; f3=[]
    for i in range(len(ages)):
        zmetal.append(uc.metallicity(FeH[i], 'fe/htoz'))

    for j in range(len(zmetal)):
        if zmetal[j]<0.0005:
            z1.append(zmetal[j]); f1.append(FeH[j]); a1.append(ages[j])
        elif zmetal[j]<0.005:
            z2.append(zmetal[j]); f2.append(FeH[j]); a2.append(ages[j])
        else:
            z3.append(zmetal[j]); f3.append(FeH[j]); a3.append(ages[j])


    np.savetxt('/projects/b1095/syr904/projects/SGRB/AgeDistribution/bin1/Metallicity_low.dat', np.c_[a1, f1, z1], fmt='%f %f %f', header='1.Age 2.Fe/H 3.Z', delimiter='', comments='#')
    np.savetxt('/projects/b1095/syr904/projects/SGRB/AgeDistribution/bin1/Metallicity_middle.dat', np.c_[a2, f2, z2], fmt='%f %f %f', header='1.Age 2.Fe/H 3.Z', delimiter='', comments='#')
    np.savetxt('/projects/b1095/syr904/projects/SGRB/AgeDistribution/bin1/Metallicity_high.dat', np.c_[a3, f3, z3], fmt='%f %f %f', header='1.Age 2.Fe/H 3.Z', delimiter='', comments='#')


    z1=[]; z2=[]; z3=[]
    a1=[]; a2=[]; a3=[]
    f1=[]; f2=[]; f3=[]
    for k in range(len(zmetal)):
        if zmetal[k]<=0.00065:
            z1.append(zmetal[k]); f1.append(FeH[k]); a1.append(ages[k])
        elif zmetal[k]<=0.0065:
            z2.append(zmetal[k]); f2.append(FeH[k]); a2.append(ages[k])
        else:
            z3.append(zmetal[k]); f3.append(FeH[k]); a3.append(ages[k])


    np.savetxt('/projects/b1095/syr904/projects/SGRB/AgeDistribution/bin2/Metallicity_low.dat', np.c_[a1, f1, z1], fmt='%f %f %f', header='1.Age 2.Fe/H 3.Z', delimiter='', comments='#')
    np.savetxt('/projects/b1095/syr904/projects/SGRB/AgeDistribution/bin2/Metallicity_middle.dat', np.c_[a2, f2, z2], fmt='%f %f %f', header='1.Age 2.Fe/H 3.Z', delimiter='', comments='#')
    np.savetxt('/projects/b1095/syr904/projects/SGRB/AgeDistribution/bin2/Metallicity_high.dat', np.c_[a3, f3, z3], fmt='%f %f %f', header='1.Age 2.Fe/H 3.Z', delimiter='', comments='#')




##Randomly draws the mergers and make a new merger file with more mergers
def make_random_draw(filepath):
    ##Read merger data
    esc_dns=np.genfromtxt(filepath+'Esc_DNS_maingrid.dat')
    esc_nsbh=np.genfromtxt(filepath+'Esc_NSBH_maingrid.dat')
    in_dns=np.genfromtxt(filepath+'Incluster_DNS_maingrid.dat')
    in_nsbh=np.genfromtxt(filepath+'Incluster_NSBH_maingrid.dat')
    gw_dns=np.genfromtxt(filepath+'GWcap_DNS_maingrid.dat')
    gw_nsbh=np.genfromtxt(filepath+'GWcap_NSBH_maingrid.dat')

    esc_all=[esc_dns, esc_nsbh]
    in_all=[in_dns, in_nsbh]
    gw_all=[gw_dns, gw_nsbh]


    ##Read metallicity data
    metal_low=np.genfromtxt('/projects/b1095/syr904/projects/SGRB/AgeDistribution/bin2/Metallicity_low.dat')
    metal_middle=np.genfromtxt('/projects/b1095/syr904/projects/SGRB/AgeDistribution/bin2/Metallicity_middle.dat')
    metal_high=np.genfromtxt('/projects/b1095/syr904/projects/SGRB/AgeDistribution/bin2/Metallicity_high.dat')

    dict_metallicity = {'0.0002':metal_low,'0.002':metal_middle,'0.02':metal_high}

    data_mod=np.genfromtxt('/projects/b1095/syr904/projects/SGRB/newruns/finaldata/ns_number_9to14Gyr_maingrid.dat')
    metal_mod=data_mod[:,16]
    #metal_mod=[0.002]

    ##Combine the files
    models=[]; tmergers=[]; types=[]; formflag=[]
    models_pri=[]; tmergers_pri=[]; types_pri=[]
    models_dyn=[]; tmergers_dyn=[]; types_dyn=[]
    for i in range(2):
        types.append([])
        types[i]=list(np.full(len(esc_all[i][:,0]), 1))+list(np.full(len(in_all[i][:,0]), 2))+list(np.full(len(gw_all[i][:,0]), 3))
        #print len(types[i])

        models.append([])
        models[i]=list(esc_all[i][:,0])+list(in_all[i][:,0])+list(gw_all[i][:,0])
        #print len(models[i])

        tmergers.append([])
        tesc=esc_all[i][:,2]; tinsp=esc_all[i][:,3]
        tmergers[i]=list(np.add(tesc, tinsp))+list(in_all[i][:,2])+list(gw_all[i][:,2])
        #print len(tmergers[i])

        formflag.append([])
        formflag[i]=list(esc_all[i][:,11])+list(in_all[i][:,11])+list(np.full(len(gw_all[i][:,0]), 0))
        models_pri.append([]); tmergers_pri.append([]); types_pri.append([])
        models_dyn.append([]); tmergers_dyn.append([]); types_dyn.append([])
        for x in range(len(formflag[i])):
            if formflag[i][x]==1: 
                models_pri[i].append(models[i][x])
                tmergers_pri[i].append(tmergers[i][x])
                types_pri[i].append(types[i][x])

            else:
                models_dyn[i].append(models[i][x])
                tmergers_dyn[i].append(tmergers[i][x])
                types_dyn[i].append(types[i][x])


    ##Make the random draws
    fname=['DNS', 'NSBH']
    for j in range(2):
        f=open(filepath+'random_draw_'+fname[j]+'.dat', 'a+', 0)
        f.write('#1.Model, 2.Tmerger, 3.Types(1-esc, 2-incluster, 3-coll)\n')
        for k in range(len(tmergers[j])):
            modelno=int(models[j][k])
            n=0
            while n<500:
                age = np.random.choice(dict_metallicity[str(metal_mod[modelno])][:,0])*1000.0
                t_new=13721.-age+tmergers[j][k]
                #print age, t_new
                if t_new<13721.:
                    f.write('%d %f %d\n'%(models[j][k], t_new, types[j][k]))
                n+=1

        f.close()


        #f1=open(filepath+'random_draw_'+fname[j]+'_pri_more.dat', 'a+', 0)
        #f2=open(filepath+'random_draw_'+fname[j]+'_dyn_more.dat', 'a+', 0)
        #for l in range(len(tmergers_pri[j])):
        #    modelno=int(models_pri[j][l])
        #    n=0
        #    while n<10000:
        #        age = np.random.choice(dict_metallicity[str(metal_mod[modelno])][:,0])*1000.0
        #        t_new=13721.-age+tmergers_pri[j][l]
        #        #print age, t_new
        #        if t_new<13721.:
        #            f1.write('%d %f %d\n'%(models_pri[j][l], t_new, types_pri[j][l]))
        #        n+=1

        #f1.close()

        #for m in range(len(tmergers_dyn[j])):
        #    modelno=int(models_dyn[j][m])
        #    n=0
        #    while n<10000:
        #        age = np.random.choice(dict_metallicity[str(metal_mod[modelno])][:,0])*1000.0
        #        t_new=13721.-age+tmergers_dyn[j][m]
        #        #print age, t_new
        #        if t_new<13721.:
        #            f2.write('%d %f %d\n'%(models_dyn[j][m], t_new, types_dyn[j][m]))
        #       n+=1

        #f2.close()




def get_merger_rate(filepath):
    ##Read data
    data_dns=np.genfromtxt(filepath+'random_draw_DNS_pri.dat')
    data_nsbh=np.genfromtxt(filepath+'random_draw_NSBH_pri.dat')
    data_all=[data_dns, data_nsbh]


    ##Define time bins
    steps=400
    t_array=np.linspace(800., 13721., steps)
    dt=(13721.-800.)/steps*1.e6   ##in years
    rho0=0.33 ##Mpc^-3


    ##Model weights
    weights, noweights=get_weights()
    n_random=500.
    weights[:] = [x / n_random for x in weights]
    modtot=np.sum(noweights)
    print(modtot)
    noweights[:] = [y / n_random/modtot for y in noweights]
    #print modtot, noweights
    #noweights=[1./n_random]
    


    ##Counting number of mergers in each bin
    ##Also calculate the reshift and comoving distance at each bin
    model_all=[]; tmerger_all=[]; count_all=[]

    for i in range(2):
        model_all.append([])
        model_all[i]=data_all[i][:,0]

        tmerger_all.append([])
        tmerger_all[i]=data_all[i][:,1]


        count_all.append([])
        redshift=[]; Dcom=[]
        for j in range(len(t_array)-1):
            counts=0
            for x in range(len(tmerger_all[i])):
                modelno=int(model_all[i][x])
                if t_array[j]<tmerger_all[i][x]<=t_array[j+1]:
                    counts+=1.*noweights[modelno]

            count_all[i].append(counts)


            ##Calculate the redshifts and the comoving distances
            redshift.append(uc.ttoredshift(t_array[j]/1000.))
            Dcom.append(uc.ComovingDistance(redshift[j]).value)

        #print count_all[i]

    redshift.append(uc.ttoredshift(t_array[-1]/1000.))
    Dcom.append(uc.ComovingDistance(redshift[-1]).value)

    #print redshift, Dcom


    ##Calculate the merger rate according to O'leary+2006 equation (13).
    R_all=[]; R_density=[]
    for i in range(2):
        R_all.append([]); R_density.append([])
        zmid=[]
        for k in range(len(t_array)-1):
            zmid.append((redshift[k+1]+redshift[k])/2.)
            R=(count_all[i][k]/dt)*(4.*np.pi/3.)*rho0*(Dcom[k]**3-Dcom[k+1]**3)/(1+zmid[k])
            Rden=(count_all[i][k]/dt)*rho0*(1000.**3)
            #print R, Rden
            R_all[i].append(R); R_density[i].append(Rden)

        R_all[i]=list(reversed(R_all[i])); R_density[i]=list(reversed(R_density[i]))

    zmid=list(reversed(zmid))

    #np.savetxt(filepath+'rates_dns_pri_rho0.33_500draw.dat', np.c_[R_all[0], R_density[0], zmid], fmt='%f %f %f', header='Rates, Rate_density, Redshift', delimiter='', comments='#')
    #np.savetxt(filepath+'rates_nsbh_pri_rho0.33_500draw.dat', np.c_[R_all[1], R_density[1], zmid], fmt='%f %f %f', header='Rates, Rate_density, Redshift', delimiter='', comments='#')


    #return R_all, R_density, zmid



##Calculate merger rate without weight and age distribution
##Unfinished
def get_merger_rate_simple(filepath):
    ##Read merger data
    esc_dns=np.genfromtxt(filepath+'Esc_DNS.dat')
    esc_nsbh=np.genfromtxt(filepath+'Esc_NSBH.dat')
    in_dns=np.genfromtxt(filepath+'Incluster_DNS.dat')
    in_nsbh=np.genfromtxt(filepath+'Incluster_NSBH.dat')
    gw_dns=np.genfromtxt(filepath+'GWcap_DNS.dat')
    gw_nsbh=np.genfromtxt(filepath+'GWcap_NSBH.dat')

    esc_all=[esc_dns, esc_nsbh]
    in_all=[in_dns, in_nsbh]
    gw_all=[gw_dns, gw_nsbh]


    ##Combine the files
    models=[]; tmergers=[]; types=[]; formflag=[]
    models_pri=[]; tmergers_pri=[]; types_pri=[]
    models_dyn=[]; tmergers_dyn=[]; types_dyn=[]
    for i in range(2):
        types.append([])
        types[i]=list(np.full(len(esc_all[i][:,0]), 1))+list(np.full(len(in_all[i][:,0]), 2))+list(np.full(len(gw_all[i][:,0]), 3))
        #print len(types[i])

        models.append([])
        models[i]=list(esc_all[i][:,0])+list(in_all[i][:,0])+list(gw_all[i][:,0])
        #print len(models[i])

        tmergers.append([])
        tesc=esc_all[i][:,2]; tinsp=esc_all[i][:,3]
        tmergers[i]=list(np.add(tesc, tinsp))+list(in_all[i][:,2])+list(gw_all[i][:,2])
        #print len(tmergers[i])

        formflag.append([])
        formflag[i]=list(esc_all[i][:,11])+list(in_all[i][:,11])+list(np.full(len(gw_all[i][:,0]), 0))
        models_pri.append([]); tmergers_pri.append([]); types_pri.append([])
        models_dyn.append([]); tmergers_dyn.append([]); types_dyn.append([])
        for x in range(len(formflag[i])):
            if formflag[i][x]==1: 
                models_pri[i].append(models[i][x])
                tmergers_pri[i].append(tmergers[i][x])
                types_pri[i].append(types[i][x])

            else:
                models_dyn[i].append(models[i][x])
                tmergers_dyn[i].append(tmergers[i][x])
                types_dyn[i].append(types[i][x])


    ##Make the random draws
    fname=['DNS', 'NSBH']
    for j in range(2):
        f=open(filepath+'random_draw_'+fname[j]+'.dat', 'a+', 0)
        f.write('#1.Model, 2.Tmerger, 3.Types(1-esc, 2-incluster, 3-coll)\n')
        for k in range(len(tmergers[j])):
            modelno=int(models[j][k])
            n=0
            while n<500:
                age = np.random.choice(dict_metallicity[str(metal_mod[modelno])][:,0])*1000.0
                t_new=13721.-age+tmergers[j][k]
                #print age, t_new
                if t_new<13721.:
                    f.write('%d %f %d\n'%(models[j][k], t_new, types[j][k]))
                n+=1

        f.close()


        #f1=open(filepath+'random_draw_'+fname[j]+'_pri_more.dat', 'a+', 0)
        #f2=open(filepath+'random_draw_'+fname[j]+'_dyn_more.dat', 'a+', 0)
        #for l in range(len(tmergers_pri[j])):
        #    modelno=int(models_pri[j][l])
        #    n=0
        #    while n<10000:
        #        age = np.random.choice(dict_metallicity[str(metal_mod[modelno])][:,0])*1000.0
        #        t_new=13721.-age+tmergers_pri[j][l]
        #        #print age, t_new
        #        if t_new<13721.:
        #            f1.write('%d %f %d\n'%(models_pri[j][l], t_new, types_pri[j][l]))
        #        n+=1

        #f1.close()

        #for m in range(len(tmergers_dyn[j])):
        #    modelno=int(models_dyn[j][m])
        #    n=0
        #    while n<10000:
        #        age = np.random.choice(dict_metallicity[str(metal_mod[modelno])][:,0])*1000.0
        #        t_new=13721.-age+tmergers_dyn[j][m]
        #        #print age, t_new
        #        if t_new<13721.:
        #            f2.write('%d %f %d\n'%(models_dyn[j][m], t_new, types_dyn[j][m]))
        #       n+=1

        #f2.close()


    ##Read data
    data_dns=np.genfromtxt(filepath+'random_draw_DNS_pri.dat')
    data_nsbh=np.genfromtxt(filepath+'random_draw_NSBH_pri.dat')
    data_all=[data_dns, data_nsbh]


    ##Define time bins
    steps=400
    t_array=np.linspace(800., 13721., steps)
    dt=(13721.-800.)/steps*1.e6   ##in years
    rho0=0.33 ##Mpc^-3


    ##Model weights
    weights, noweights=get_weights()
    n_random=500.
    weights[:] = [x / n_random for x in weights]
    modtot=np.sum(noweights)
    print(modtot)
    noweights[:] = [y / n_random/modtot for y in noweights]
    #print modtot, noweights
    #noweights=[1./n_random]
    


    ##Counting number of mergers in each bin
    ##Also calculate the reshift and comoving distance at each bin
    model_all=[]; tmerger_all=[]; count_all=[]

    for i in range(2):
        model_all.append([])
        model_all[i]=data_all[i][:,0]

        tmerger_all.append([])
        tmerger_all[i]=data_all[i][:,1]


        count_all.append([])
        redshift=[]; Dcom=[]
        for j in range(len(t_array)-1):
            counts=0
            for x in range(len(tmerger_all[i])):
                modelno=int(model_all[i][x])
                if t_array[j]<tmerger_all[i][x]<=t_array[j+1]:
                    counts+=1.*noweights[modelno]

            count_all[i].append(counts)


            ##Calculate the redshifts and the comoving distances
            redshift.append(uc.ttoredshift(t_array[j]/1000.))
            Dcom.append(uc.ComovingDistance(redshift[j]).value)

        #print count_all[i]

    redshift.append(uc.ttoredshift(t_array[-1]/1000.))
    Dcom.append(uc.ComovingDistance(redshift[-1]).value)

    #print redshift, Dcom


    ##Calculate the merger rate according to O'leary+2006 equation (13).
    R_all=[]; R_density=[]
    for i in range(2):
        R_all.append([]); R_density.append([])
        zmid=[]
        for k in range(len(t_array)-1):
            zmid.append((redshift[k+1]+redshift[k])/2.)
            R=(count_all[i][k]/dt)*(4.*np.pi/3.)*rho0*(Dcom[k]**3-Dcom[k+1]**3)/(1+zmid[k])
            Rden=(count_all[i][k]/dt)*rho0*(1000.**3)
            #print R, Rden
            R_all[i].append(R); R_density[i].append(Rden)

        R_all[i]=list(reversed(R_all[i])); R_density[i]=list(reversed(R_density[i]))

    zmid=list(reversed(zmid))

    #np.savetxt(filepath+'rates_dns_pri_rho0.33_500draw.dat', np.c_[R_all[0], R_density[0], zmid], fmt='%f %f %f', header='Rates, Rate_density, Redshift', delimiter='', comments='#')
    #np.savetxt(filepath+'rates_nsbh_pri_rho0.33_500draw.dat', np.c_[R_all[1], R_density[1], zmid], fmt='%f %f %f', header='Rates, Rate_density, Redshift', delimiter='', comments='#')


    #return R_all, R_density, zmid




