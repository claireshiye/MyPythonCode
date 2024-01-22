import numpy as np 
import os,sys
from collections import Counter
import re
import gzip
import scripts
import scripts1
import scripts2
import scripts3
import dynamics as dyn
import unit_convert as uc
import gw_ecc_calc as gwcalc
import LISA_calculations as lisa
import ns_history as nh
import useful_function as uf

yearsc=31557600.
twopi=6.283185307179586
Gconst=6.674*10**-8 ##cm3*g-1*s-2
clight=3*10**10 ##cm/s
Msun=2*10**33 ##gram
AU=1.496*10**13  ##cm
PC=3.086*10**18  ##cm
Kconst=9.87*10**-48 ##yr/G^2
Lsun=4.02*10**16 ##mJy*kpc^2

##Find the number of XX collisions in the models
##startype is a list
def get_XX_collision(pathlist, start, end, startype, savepath, readflag):
    if readflag == 1:
        sourcedir=np.genfromtxt(pathlist, dtype=str)
        filepaths=sourcedir[:,0]; status=sourcedir[:,1]
    else:
        filepaths=pathlist; status=[1]

    #model=[]; model_status=[]; mm=[]; mcom=[]; ktypem=[]; kcom=[]; timem=[]; idm=[]; rm=[]; colltype=[]
    fcoll=open(savepath+'_coll_all.dat', 'w+')
    fcoll.write('#1.Model 2.Time(Myr) 3.IDcoll 4.Radius(pc) 5.Mcoll 6.M0 7.M1 8.M2 9.M3 10.kcoll 11.k0 12.k1 13.k2 14.k3 15.model_status 16.COLLTYPE\n')
    for i in range(len(filepaths)):
        filestr=filepaths[i]+'initial'

        t_conv=dyn.conv('t', filestr+'.conv.sh')
        l_conv=dyn.conv('l', filestr+'.conv.sh')

        collfile=filestr+'.collision.log'
        collfile2=filestr+'2.collision.log'
        colldata=scripts1.readcollfile(collfile)
        if os.path.isfile(collfile2) and os.path.getsize(collfile2) > 0:
            colldata2=scripts1.readcollfile(collfile2)
            colldata=colldata+colldata2

        for j in range(len(colldata)):
            line=colldata[j].split()
            if line[1]=='single-single':  ##Single-single star collision
                colltype='SS'
                if int(line[11]) in startype or int(line[12]) in startype:
                    model=i; model_status=int(status[i]); timem=t_conv*float(line[0])
                    mm=float(line[4]); m0=float(line[6]); m1=float(line[8]); m2=-100; m3=-100
                    ktypem=int(line[10]); ktype0=int(line[11]); ktype1=int(line[12]); ktype2=-100; ktype3=-100
                    idm=int(line[3]); rm=float(line[9])*l_conv

                    fcoll.write('%d %f %d %f %f %f %f %f %f %d %d %d %d %d %d %s\n'%(model, timem, idm, rm, mm, m0, m1, m2, m3, ktypem, ktype0, ktype1, ktype2, ktype3, model_status, colltype))


            if line[1]=='binary-single':   ##Binary-single collision
                colltype='BS'
                if int(line[2])==2:
                    if int(line[11]) in startype or int(line[12]) in startype:
                        model=i; model_status=int(status[i]); timem=t_conv*float(line[0])
                        mm=float(line[4]); m0=float(line[6]); m1=float(line[8]); m2=-100; m3=-100
                        ktypem=int(line[10]); ktype0=int(line[11]); ktype1=int(line[12]); ktype2=-100; ktype3=-100
                        idm=int(line[3]); rm=float(line[9])*l_conv

                        fcoll.write('%d %f %d %f %f %f %f %f %f %d %d %d %d %d %d %s\n'%(model, timem, idm, rm, mm, m0, m1, m2, m3, ktypem, ktype0, ktype1, ktype2, ktype3, model_status, colltype))

                if int(line[2])==3:
                    if int(line[13]) in startype or int(line[14]) in startype or int(line[15]) in startype:
                        model=i; model_status=int(status[i]); timem=t_conv*float(line[0])
                        mm=float(line[4]); m0=float(line[6]); m1=float(line[8]); m2=float(line[10]); m3=-100
                        ktypem=int(line[12]); ktype0=int(line[13]); ktype1=int(line[14]); ktype2=int(line[15]); ktype3=-100
                        idm=int(line[3]); rm=float(line[11])*l_conv

                        fcoll.write('%d %f %d %f %f %f %f %f %f %d %d %d %d %d %d %s\n'%(model, timem, idm, rm, mm, m0, m1, m2, m3, ktypem, ktype0, ktype1, ktype2, ktype3, model_status, colltype))


            if line[1]=='binary-binary':   ##Binary-binary collision
                colltype='BB'
                if int(line[2])==2:
                    if int(line[11]) in startype or int(line[12]) in startype:
                        model=i; model_status=int(status[i]); timem=t_conv*float(line[0])
                        mm=float(line[4]); m0=float(line[6]); m1=float(line[8]); m2=-100; m3=-100
                        ktypem=int(line[10]); ktype0=int(line[11]); ktype1=int(line[12]); ktype2=-100; ktype3=-100
                        idm=int(line[3]); rm=float(line[9])*l_conv

                        fcoll.write('%d %f %d %f %f %f %f %f %f %d %d %d %d %d %d %s\n'%(model, timem, idm, rm, mm, m0, m1, m2, m3, ktypem, ktype0, ktype1, ktype2, ktype3, model_status, colltype))

                if int(line[2])==3:
                    if int(line[13]) in startype or int(line[14]) in startype or int(line[15]) in startype:
                        model=i; model_status=int(status[i]); timem=t_conv*float(line[0])
                        mm=float(line[4]); m0=float(line[6]); m1=float(line[8]); m2=float(line[10]); m3=-100
                        ktypem=int(line[12]); ktype0=int(line[13]); ktype1=int(line[14]); ktype2=int(line[15]); ktype3=-100
                        idm=int(line[3]); rm=float(line[11])*l_conv

                        fcoll.write('%d %f %d %f %f %f %f %f %f %d %d %d %d %d %d %s\n'%(model, timem, idm, rm, mm, m0, m1, m2, m3, ktypem, ktype0, ktype1, ktype2, ktype3, model_status, colltype))

                if int(line[2])==4:
                    if int(line[15]) in startype or int(line[16]) in startype or int(line[17]) in startype or int(line[18]) in startype:
                        model=i; model_status=int(status[i]); timem=t_conv*float(line[0])
                        mm=float(line[4]); m0=float(line[6]); m1=float(line[8]); m2=float(line[10]); m3=float(line[12])
                        ktypem=int(line[14]); ktype0=int(line[15]); ktype1=int(line[16]); ktype2=int(line[17]); ktype3=int(line[18])
                        idm=int(line[3]); rm=float(line[13])*l_conv

                        fcoll.write('%d %f %d %f %f %f %f %f %f %d %d %d %d %d %d %s\n'%(model, timem, idm, rm, mm, m0, m1, m2, m3, ktypem, ktype0, ktype1, ktype2, ktype3, model_status, colltype))


        print(i)

    fcoll.close()



##Find the number of XX collisions in the models
##klim1 and klim2 are lists
def get_coll_NSfromWDWD(pathlist, start, end, k_finl, klim1, klim2, savepath, readflag):
    if readflag == 1:
        sourcedir=np.genfromtxt(pathlist, dtype=str)
        filepaths=sourcedir[:,0]; status=sourcedir[:,1]
    else:
        filepaths=pathlist; status=[1]

    #fcoll=open(savepath+'_coll.dat', 'w+')
    #fcoll.write('#1.Model 2.Time(Myr) 3.IDcoll 4.Radius(pc) 5.Mcoll 6.M0 7.M1 8.M2 9.M3 10.kcoll 11.k0 12.k1 13.k2 14.k3 15.model_status 16.COLLTYPE\n')
    for i in range(start, end):
        fcoll=open(filepaths[i]+'ns_wdwd_coll.dat', 'w+')
        fcoll.write('#1.Model 2.Time(Myr) 3.IDcoll 4.Radius(pc) 5.Mcoll 6.M0 7.M1 8.M2 9.M3 10.kcoll 11.k0 12.k1 13.k2 14.k3 15.model_status 16.COLLTYPE\n')

        filestr=filepaths[i]+'initial'
        
        t_conv=dyn.conv('t', filestr+'.conv.sh')
        l_conv=dyn.conv('l', filestr+'.conv.sh')

        collfile=filestr+'.collision.log'
        collfile2=filestr+'2.collision.log'
        colldata=scripts1.readcollfile(collfile)
        if os.path.isfile(collfile2) and os.path.getsize(collfile2) > 0:
            colldata2=scripts1.readcollfile(collfile2)
            colldata=colldata+colldata2

        for j in range(len(colldata)):
            line=colldata[j].split()
            if line[1]=='single-single':  ##Single-single star collision
                if int(line[10])==k_finl:
                    print(int(line[10]), int(line[11]), int(line[12]))
                colltype='SS'
                if int(line[10])==k_finl and int(line[11]) in klim1 and int(line[12]) in klim2:
                    model=i; model_status=int(status[i]); timem=t_conv*float(line[0])
                    mm=float(line[4]); m0=float(line[6]); m1=float(line[8]); m2=-100; m3=-100
                    ktypem=int(line[10]); ktype0=int(line[11]); ktype1=int(line[12]); ktype2=-100; ktype3=-100
                    idm=int(line[3]); rm=float(line[9])*l_conv

                    fcoll.write('%d %f %d %f %f %f %f %f %f %d %d %d %d %d %d %s\n'%(model, timem, idm, rm, mm, m0, m1, m2, m3, ktypem, ktype0, ktype1, ktype2, ktype3, model_status, colltype))


            if line[1]=='binary-single':   ##Binary-single collision
                colltype='BS'
                if int(line[2])==2:
                    if int(line[10])==k_finl:
                        print(int(line[10]), int(line[11]), int(line[12]))
                    if int(line[10])==k_finl and int(line[11]) in klim1 and int(line[12]) in klim2:
                        model=i; model_status=int(status[i]); timem=t_conv*float(line[0])
                        mm=float(line[4]); m0=float(line[6]); m1=float(line[8]); m2=-100; m3=-100
                        ktypem=int(line[10]); ktype0=int(line[11]); ktype1=int(line[12]); ktype2=-100; ktype3=-100
                        idm=int(line[3]); rm=float(line[9])*l_conv

                        fcoll.write('%d %f %d %f %f %f %f %f %f %d %d %d %d %d %d %s\n'%(model, timem, idm, rm, mm, m0, m1, m2, m3, ktypem, ktype0, ktype1, ktype2, ktype3, model_status, colltype))


            if line[1]=='binary-binary':   ##Binary-binary collision
                colltype='BB'
                if int(line[2])==2:
                    if int(line[10])==k_finl:
                        print(int(line[10]), int(line[11]), int(line[12]))
                    if int(line[10])==k_finl and int(line[11]) in klim1 and int(line[12]) in klim2:
                        model=i; model_status=int(status[i]); timem=t_conv*float(line[0])
                        mm=float(line[4]); m0=float(line[6]); m1=float(line[8]); m2=-100; m3=-100
                        ktypem=int(line[10]); ktype0=int(line[11]); ktype1=int(line[12]); ktype2=-100; ktype3=-100
                        idm=int(line[3]); rm=float(line[9])*l_conv

                        fcoll.write('%d %f %d %f %f %f %f %f %f %d %d %d %d %d %d %s\n'%(model, timem, idm, rm, mm, m0, m1, m2, m3, ktypem, ktype0, ktype1, ktype2, ktype3, model_status, colltype))

        print(i)

    fcoll.close()


##Find the number of XX mergers in the models
##klim1 and klim2 are lists
def get_merger_NSfromWDWD(pathlist, start, end, k_finl, klim1, klim2, savepath, readflag):
    if readflag == 1:
        sourcedir=np.genfromtxt(pathlist, dtype=str)
        filepaths=sourcedir[:,0]; status=sourcedir[:,1].astype(int)
    else:
        filepaths=pathlist; status=[1]

    #fmerge=open(savepath+'_merger.dat', 'w+')
    #fmerge.write('#1.Model 2.Time(Myr) 3.IDmerg 4.Radius(pc) 5.Mmerg 6.M0 7.M1 8.kmerg 9.k0 10.k1 11.model_status\n')
    for i in range(start,end):
        fmerge=open(filepaths[i]+'ns_wdwd_merger.dat', 'w+')
        fmerge.write('#1.Model 2.Time(Myr) 3.IDmerg 4.Radius(pc) 5.Mmerg 6.M0 7.M1 8.kmerg 9.k0 10.k1 11.model_status\n')

        filestr=filepaths[i]+'initial'

        t_conv=dyn.conv('t', filestr+'.conv.sh')
        l_conv=dyn.conv('l', filestr+'.conv.sh')

        semergefile = filestr+'.semergedisrupt.log'
        semergefile2 = filestr+'2.semergedisrupt.log'
        semergedata=scripts2.readmergefile(semergefile)
        if os.path.isfile(semergefile2) and os.path.getsize(semergefile2)>0:
            semergedata2=scripts2.readmergefile(semergefile2)
            semergedata=semergedata+semergedata2

        for hh in range(len(semergedata)):
            line=semergedata[hh].split()
            #print(line)
            if int(line[1])<3:
                if int(line[-1]) in klim1 and int(line[-2]) in klim2 and int(line[-3])==k_finl:
                    fmerge.write('%d %f %d %f %f %f %f %d %d %d %d\n'%(i, t_conv*float(line[0]), 
                    int(line[2]), float(line[8])*l_conv, float(line[3]), float(line[5]), float(line[7]),int(line[-3]), int(line[-2]), int(line[-1]), status[i]))

        print(i)

    fmerge.close()


def get_coll_NSfromNSMS(pathlist, start, end, savepath, readflag):
    if readflag == 1:
        sourcedir=np.genfromtxt(pathlist, dtype=str)
        filepaths=sourcedir[:,0]; status=sourcedir[:,1]
    else:
        filepaths=pathlist; status=[1]

    for i in range(start, end):
        collfile=filepaths[i]+'initial.collision.log'
        collfile2=filepaths[i]+'initial2.collision.log'
        colldata=scripts1.readcollfile(collfile)
        if os.path.isfile(collfile2) and os.path.getsize(collfile2) > 0:
            colldata2=scripts1.readcollfile(collfile2)
            colldata=colldata+colldata2
            
        l_conv = dyn.conv('l', filepaths[i]+'initial.conv.sh')
        t_conv = dyn.conv('t', filepaths[i]+'initial.conv.sh')
            
        fcoll=open(filepaths[i]+'ns_nsms_coll.dat', 'w+')
        fcoll.write('#1.Model 2.Time(Myr) 3.IDcoll 4.Radius(pc) 5.Mcoll 6.M0 7.M1 8.M2 9.M3 10.kcoll 11.k0 12.k1 13.k2 14.k3 15.COLLTYPE\n')
    
        for j in range(len(colldata)):
            line=colldata[j].split()
            if line[1]=='single-single':  ##Single-single star collision
                colltype='SS'
                ktypem=int(line[10])
                if ktypem==13 and ((int(line[11])==13 and int(line[12]) <=1) or (int(line[12])==13 and int(line[11]) <=1)):
                    model=i; timem=t_conv*float(line[0])
                    mm=float(line[4]); m0=float(line[6]); m1=float(line[8]); m2=-100; m3=-100
                    ktype0=int(line[11]); ktype1=int(line[12]); ktype2=-100; ktype3=-100
                    idm=int(line[3]); rm=float(line[9])*l_conv
                    
                    fcoll.write('%d %f %d %f %f %f %f %f %f %d %d %d %d %d %s\n'%(model, timem, idm, rm, mm, m0, m1, m2, m3, ktypem, ktype0, ktype1, ktype2, ktype3, colltype))
    
    
            if line[1]=='binary-single':   ##Binary-single collision
                colltype='BS'
                if int(line[2])==2:
                    ktypem=int(line[10])
                    if ktypem==13 and ((int(line[11])==13 and int(line[12]) <=1) or (int(line[12])==13 and int(line[11]) <=1)):
                        model=i; timem=t_conv*float(line[0])
                        mm=float(line[4]); m0=float(line[6]); m1=float(line[8]); m2=-100; m3=-100
                        ktype0=int(line[11]); ktype1=int(line[12]); ktype2=-100; ktype3=-100
                        idm=int(line[3]); rm=float(line[9])*l_conv
                        
                        fcoll.write('%d %f %d %f %f %f %f %f %f %d %d %d %d %d %s\n'%(model, timem, idm, rm, mm, m0, m1, m2, m3, ktypem, ktype0, ktype1, ktype2, ktype3, colltype))
    
                if int(line[2])==3:
                    ktypem=int(line[12])
                    if ktypem==13:
                        if int(line[13])==13 and int(line[14])<=1 and int(line[15])<=1:
                            model=i; timem=t_conv*float(line[0])
                            mm=float(line[4]); m0=float(line[6]); m1=float(line[8]); m2=float(line[10]); m3=-100
                            ktype0=int(line[13]); ktype1=int(line[14]); ktype2=int(line[15]); ktype3=-100
                            idm=int(line[3]); rm=float(line[11])*l_conv
                        
                            fcoll.write('%d %f %d %f %f %f %f %f %f %d %d %d %d %d %s\n'%(model, timem, idm, rm, mm, m0, m1, m2, m3, ktypem, ktype0, ktype1, ktype2, ktype3, colltype))
                            
                        if int(line[14])==13 and int(line[13])<=1 and int(line[15])<=1:
                            model=i; timem=t_conv*float(line[0])
                            mm=float(line[4]); m0=float(line[6]); m1=float(line[8]); m2=float(line[10]); m3=-100
                            ktype0=int(line[13]); ktype1=int(line[14]); ktype2=int(line[15]); ktype3=-100
                            idm=int(line[3]); rm=float(line[11])*l_conv
                        
                            fcoll.write('%d %f %d %f %f %f %f %f %f %d %d %d %d %d %s\n'%(model, timem, idm, rm, mm, m0, m1, m2, m3, ktypem, ktype0, ktype1, ktype2, ktype3, colltype))
    
                        if int(line[15])==13 and int(line[13])<=1 and int(line[14])<=1:
                            model=i; timem=t_conv*float(line[0])
                            mm=float(line[4]); m0=float(line[6]); m1=float(line[8]); m2=float(line[10]); m3=-100
                            ktype0=int(line[13]); ktype1=int(line[14]); ktype2=int(line[15]); ktype3=-100
                            idm=int(line[3]); rm=float(line[11])*l_conv
                        
                            fcoll.write('%d %f %d %f %f %f %f %f %f %d %d %d %d %d %s\n'%(model, timem, idm, rm, mm, m0, m1, m2, m3, ktypem, ktype0, ktype1, ktype2, ktype3, colltype))
    
    
            if line[1]=='binary-binary':   ##Binary-binary collision
                colltype='BB'
                if int(line[2])==2:
                    ktypem=int(line[10])
                    if ktypem==13 and ((int(line[11])==13 and int(line[12])<=1) or (int(line[12])==13 and int(line[11])<=1)):
                        model=i; timem=t_conv*float(line[0])
                        mm=float(line[4]); m0=float(line[6]); m1=float(line[8]); m2=-100; m3=-100
                        ktype0=int(line[11]); ktype1=int(line[12]); ktype2=-100; ktype3=-100
                        idm=int(line[3]); rm=float(line[9])*l_conv
                        
                        fcoll.write('%d %f %d %f %f %f %f %f %f %d %d %d %d %d %s\n'%(model, timem, idm, rm, mm, m0, m1, m2, m3, ktypem, ktype0, ktype1, ktype2, ktype3, colltype))
    
                        
                if int(line[2])==3:
                    ktypem=int(line[12])
                    if ktypem==13:
                        if int(line[13])==13 and int(line[14])<=1 and int(line[15])<=1:
                            model=i; timem=t_conv*float(line[0])
                            mm=float(line[4]); m0=float(line[6]); m1=float(line[8]); m2=float(line[10]); m3=-100
                            ktype0=int(line[13]); ktype1=int(line[14]); ktype2=int(line[15]); ktype3=-100
                            idm=int(line[3]); rm=float(line[11])*l_conv
                        
                            fcoll.write('%d %f %d %f %f %f %f %f %f %d %d %d %d %d %s\n'%(model, timem, idm, rm, mm, m0, m1, m2, m3, ktypem, ktype0, ktype1, ktype2, ktype3, colltype))
                            
                        if int(line[14])==13 and int(line[13])<=1 and int(line[15])<=1:
                            model=i; timem=t_conv*float(line[0])
                            mm=float(line[4]); m0=float(line[6]); m1=float(line[8]); m2=float(line[10]); m3=-100
                            ktype0=int(line[13]); ktype1=int(line[14]); ktype2=int(line[15]); ktype3=-100
                            idm=int(line[3]); rm=float(line[11])*l_conv
                        
                            fcoll.write('%d %f %d %f %f %f %f %f %f %d %d %d %d %d %s\n'%(model, timem, idm, rm, mm, m0, m1, m2, m3, ktypem, ktype0, ktype1, ktype2, ktype3, colltype))
    
                        if int(line[15])==13 and int(line[13])<=1 and int(line[14])<=1:
                            model=i; timem=t_conv*float(line[0])
                            mm=float(line[4]); m0=float(line[6]); m1=float(line[8]); m2=float(line[10]); m3=-100
                            ktype0=int(line[13]); ktype1=int(line[14]); ktype2=int(line[15]); ktype3=-100
                            idm=int(line[3]); rm=float(line[11])*l_conv
                        
                            fcoll.write('%d %f %d %f %f %f %f %f %f %d %d %d %d %d %s\n'%(model, timem, idm, rm, mm, m0, m1, m2, m3, ktypem, ktype0, ktype1, ktype2, ktype3, colltype))
    
    
                if int(line[2])==4:
                    ktypem=int(line[14])
                    if ktypem==13:
                        if int(line[15])==13 and int(line[16])<=1 and int(line[17])<=1 and int(line[18])<=1:
                            model=i; model_status=int(status[i]); timem=t_conv*float(line[0])
                            mm=float(line[4]); m0=float(line[6]); m1=float(line[8]); m2=float(line[10]); m3=float(line[12])
                            ktype0=int(line[15]); ktype1=int(line[16]); ktype2=int(line[17]); ktype3=int(line[18])
                            idm=int(line[3]); rm=float(line[13])*l_conv
                            
                            fcoll.write('%d %f %d %f %f %f %f %f %f %d %d %d %d %d %s\n'%(model, timem, idm, rm, mm, m0, m1, m2, m3, ktypem, ktype0, ktype1, ktype2, ktype3, colltype))
                        
                        if int(line[16])==13 and int(line[15])<=1 and int(line[17])<=1 and int(line[18])<=1:
                            model=i; model_status=int(status[i]); timem=t_conv*float(line[0])
                            mm=float(line[4]); m0=float(line[6]); m1=float(line[8]); m2=float(line[10]); m3=float(line[12])
                            ktype0=int(line[15]); ktype1=int(line[16]); ktype2=int(line[17]); ktype3=int(line[18])
                            idm=int(line[3]); rm=float(line[13])*l_conv
                            
                            fcoll.write('%d %f %d %f %f %f %f %f %f %d %d %d %d %d %s\n'%(model, timem, idm, rm, mm, m0, m1, m2, m3, ktypem, ktype0, ktype1, ktype2, ktype3, colltype))
                        
                        if int(line[17])==13 and int(line[16])<=1 and int(line[15])<=1 and int(line[18])<=1:
                            model=i; model_status=int(status[i]); timem=t_conv*float(line[0])
                            mm=float(line[4]); m0=float(line[6]); m1=float(line[8]); m2=float(line[10]); m3=float(line[12])
                            ktype0=int(line[15]); ktype1=int(line[16]); ktype2=int(line[17]); ktype3=int(line[18])
                            idm=int(line[3]); rm=float(line[13])*l_conv
                            
                            fcoll.write('%d %f %d %f %f %f %f %f %f %d %d %d %d %d %s\n'%(model, timem, idm, rm, mm, m0, m1, m2, m3, ktypem, ktype0, ktype1, ktype2, ktype3, colltype))


                        if int(line[18])==13 and int(line[16])<=1 and int(line[17])<=1 and int(line[15])<=1:
                            model=i; model_status=int(status[i]); timem=t_conv*float(line[0])
                            mm=float(line[4]); m0=float(line[6]); m1=float(line[8]); m2=float(line[10]); m3=float(line[12])
                            ktype0=int(line[15]); ktype1=int(line[16]); ktype2=int(line[17]); ktype3=int(line[18])
                            idm=int(line[3]); rm=float(line[13])*l_conv
                            
                            fcoll.write('%d %f %d %f %f %f %f %f %f %d %d %d %d %d %s\n'%(model, timem, idm, rm, mm, m0, m1, m2, m3, ktypem, ktype0, ktype1, ktype2, ktype3, colltype))

        fcoll.close()


##Find the number of XX and XX collisions in the models
##klim1 and klim2 are lists
def get_coll_XXXX(pathlist, start, end, klim1, klim2, savepath, readflag):
    if readflag == 1:
        sourcedir=np.genfromtxt(pathlist, dtype=str)
        filepaths=sourcedir[:,0]; status=sourcedir[:,1]
    else:
        filepaths=pathlist; status=[1]

    fcoll=open(savepath+'_coll.dat', 'w+')
    fcoll.write('#1.Model 2.Time(Myr) 3.IDcoll 4.Radius(pc) 5.Mcoll 6.M0 7.M1 8.M2 9.M3 10.kcoll 11.k0 12.k1 13.k2 14.k3 15.model_status 16.COLLTYPE\n')
    for i in range(len(filepaths)):
        filestr=filepaths[i]+'initial'

        t_conv=dyn.conv('t', filestr+'.conv.sh')
        l_conv=dyn.conv('l', filestr+'.conv.sh')

        collfile=filestr+'.collision.log'
        collfile2=filestr+'2.collision.log'
        colldata=scripts1.readcollfile(collfile)
        if os.path.isfile(collfile2) and os.path.getsize(collfile2) > 0:
            colldata2=scripts1.readcollfile(collfile2)
            colldata=colldata+colldata2

        for j in range(len(colldata)):
            line=colldata[j].split()
            if line[1]=='single-single':  ##Single-single star collision
                colltype='SS'
                if (int(line[11]) in klim1 and int(line[12]) in klim2) or (int(line[12]) in klim1 and int(line[11]) in klim2):
                    model=i; model_status=int(status[i]); timem=t_conv*float(line[0])
                    mm=float(line[4]); m0=float(line[6]); m1=float(line[8]); m2=-100; m3=-100
                    ktypem=int(line[10]); ktype0=int(line[11]); ktype1=int(line[12]); ktype2=-100; ktype3=-100
                    idm=int(line[3]); rm=float(line[9])*l_conv

                    fcoll.write('%d %f %d %f %f %f %f %f %f %d %d %d %d %d %d %s\n'%(model, timem, idm, rm, mm, m0, m1, m2, m3, ktypem, ktype0, ktype1, ktype2, ktype3, model_status, colltype))


            if line[1]=='binary-single':   ##Binary-single collision
                colltype='BS'
                if int(line[2])==2:
                    if (int(line[11]) in klim1 and int(line[12]) in klim2) or (int(line[12]) in klim1 and int(line[11]) in klim2):
                        model=i; model_status=int(status[i]); timem=t_conv*float(line[0])
                        mm=float(line[4]); m0=float(line[6]); m1=float(line[8]); m2=-100; m3=-100
                        ktypem=int(line[10]); ktype0=int(line[11]); ktype1=int(line[12]); ktype2=-100; ktype3=-100
                        idm=int(line[3]); rm=float(line[9])*l_conv

                        fcoll.write('%d %f %d %f %f %f %f %f %f %d %d %d %d %d %d %s\n'%(model, timem, idm, rm, mm, m0, m1, m2, m3, ktypem, ktype0, ktype1, ktype2, ktype3, model_status, colltype))


            if line[1]=='binary-binary':   ##Binary-binary collision
                colltype='BB'
                if int(line[2])==2:
                    if (int(line[11]) in klim1 and int(line[12]) in klim2) or (int(line[12]) in klim1 and int(line[11]) in klim2):
                        model=i; model_status=int(status[i]); timem=t_conv*float(line[0])
                        mm=float(line[4]); m0=float(line[6]); m1=float(line[8]); m2=-100; m3=-100
                        ktypem=int(line[10]); ktype0=int(line[11]); ktype1=int(line[12]); ktype2=-100; ktype3=-100
                        idm=int(line[3]); rm=float(line[9])*l_conv

                        fcoll.write('%d %f %d %f %f %f %f %f %f %d %d %d %d %d %d %s\n'%(model, timem, idm, rm, mm, m0, m1, m2, m3, ktypem, ktype0, ktype1, ktype2, ktype3, model_status, colltype))

        print(i)

    fcoll.close()


##Find the number of XX-XX mergers in the models
##klim1 and klim2 are lists
def get_merger_XXXX(pathlist, start, end, klim1, klim2, savepath, readflag):
    if readflag == 1:
        sourcedir=np.genfromtxt(pathlist, dtype=str)
        filepaths=sourcedir[:,0]; status=sourcedir[:,1]
    else:
        filepaths=pathlist; status=[1]

    fmerge=open(savepath+'_merger.dat', 'w+')
    fmerge.write('#1.Model 2.Time(Myr) 3.IDmerg 4.Radius(pc) 5.Mmerg 6.M0 7.M1 8.kmerg 9.k0 10.k1 11.model_status\n')
    for i in range(len(filepaths)):
        filestr=filepaths[i]+'initial'

        t_conv=dyn.conv('t', filestr+'.conv.sh')
        l_conv=dyn.conv('l', filestr+'.conv.sh')

        semergefile = filestr+'.semergedisrupt.log'
        semergefile2 = filestr+'2.semergedisrupt.log'
        semergedata=scripts2.readmergefile(semergefile)
        if os.path.isfile(semergefile2) and os.path.getsize(semergefile2)>0:
            semergedata2=scripts2.readmergefile(semergefile2)
            semergedata=semergedata+semergedata2

        for hh in range(len(semergedata)):
            line=semergedata[hh].split()
            #print(line)
            if int(line[1])<3:
                if int(line[-1]) in klim1 and int(line[-2]) in klim2:
                    fmerge.write('%d %f %d %f %f %f %f %d %d %d %d\n'%(i, t_conv*float(line[0]), 
                    int(line[2]), float(line[8])*l_conv, float(line[3]), float(line[5]), float(line[7]),int(line[-3]), int(line[-2]), int(line[-1]), status[i]))

        print(i)

    fmerge.close()

            
##Find the number of collisions in the models where the product is XX
##startype is a list
def get_XX_collproduct(pathlist, start, end, startype, savepath, readflag):
    if readflag == 1:
        sourcedir=np.genfromtxt(pathlist, dtype=str)
        filepaths=sourcedir[:,0]; status=sourcedir[:,1]
    else:
        filepaths=pathlist; status=[1]

    #model=[]; model_status=[]; mm=[]; mcom=[]; ktypem=[]; kcom=[]; timem=[]; idm=[]; rm=[]; colltype=[]
    #fcoll=open(savepath+'_collproduct_all.dat', 'w+')
    #fcoll.write('#1.Model 2.Time(Myr) 3.IDcoll 4.Radius(pc) 5.Mcoll 6.M0 7.M1 8.M2 9.M3 10.kcoll 11.k0 12.k1 13.k2 14.k3 15.id0 16.id1 17.id2 18.id3 19.model_status 20.COLLTYPE\n')
    for i in range(start, end):
        fcoll=open(filepaths[i]+'ns_collproduct_all.dat', 'w+')
        fcoll.write('#1.Model 2.Time(Myr) 3.IDcoll 4.Radius(pc) 5.Mcoll 6.M0 7.M1 8.M2 9.M3 10.kcoll 11.k0 12.k1 13.k2 14.k3 15.id0 16.id1 17.id2 18.id3 19.model_status 20.COLLTYPE\n')

        filestr=filepaths[i]+'initial'

        t_conv=dyn.conv('t', filestr+'.conv.sh')
        l_conv=dyn.conv('l', filestr+'.conv.sh')

        collfile=filestr+'.collision.log'
        collfile2=filestr+'2.collision.log'
        colldata=scripts1.readcollfile(collfile)
        if os.path.isfile(collfile2) and os.path.getsize(collfile2) > 0:
            colldata2=scripts1.readcollfile(collfile2)
            colldata=colldata+colldata2

        for j in range(len(colldata)):
            line=colldata[j].split()
            if line[1]=='single-single':  ##Single-single star collision
                colltype='SS'
                if int(line[10]) in startype:
                    model=i; model_status=int(status[i]); timem=t_conv*float(line[0])
                    mm=float(line[4]); m0=float(line[6]); m1=float(line[8]); m2=-100; m3=-100
                    ktypem=int(line[10]); ktype0=int(line[11]); ktype1=int(line[12]); ktype2=-100; ktype3=-100
                    idm=int(line[3]); rm=float(line[9])*l_conv
                    id0=int(line[5]); id1=int(line[7]); id2=-100; id3=-100

                    fcoll.write('%d %f %d %f %f %f %f %f %f %d %d %d %d %d %d %d %d %d %d %s\n'%(model, timem, idm, rm, mm, m0, m1, m2, m3, ktypem, ktype0, ktype1, ktype2, ktype3, id0, id1, id2, id3, model_status, colltype))


            if line[1]=='binary-single':   ##Binary-single collision
                colltype='BS'
                if int(line[2])==2:
                    if int(line[10]) in startype:
                        model=i; model_status=int(status[i]); timem=t_conv*float(line[0])
                        mm=float(line[4]); m0=float(line[6]); m1=float(line[8]); m2=-100; m3=-100
                        ktypem=int(line[10]); ktype0=int(line[11]); ktype1=int(line[12]); ktype2=-100; ktype3=-100
                        idm=int(line[3]); rm=float(line[9])*l_conv
                        id0=int(line[5]); id1=int(line[7]); id2=-100; id3=-100

                        fcoll.write('%d %f %d %f %f %f %f %f %f %d %d %d %d %d %d %d %d %d %d %s\n'%(model, timem, idm, rm, mm, m0, m1, m2, m3, ktypem, ktype0, ktype1, ktype2, ktype3, id0, id1, id2, id3, model_status, colltype))

                if int(line[2])==3:
                    if int(line[12]) in startype:
                        model=i; model_status=int(status[i]); timem=t_conv*float(line[0])
                        mm=float(line[4]); m0=float(line[6]); m1=float(line[8]); m2=float(line[10]); m3=-100
                        ktypem=int(line[12]); ktype0=int(line[13]); ktype1=int(line[14]); ktype2=int(line[15]); ktype3=-100
                        idm=int(line[3]); rm=float(line[11])*l_conv
                        id0=int(line[5]); id1=int(line[7]); id2=int(line[9]); id3=-100

                        fcoll.write('%d %f %d %f %f %f %f %f %f %d %d %d %d %d %d %d %d %d %d %s\n'%(model, timem, idm, rm, mm, m0, m1, m2, m3, ktypem, ktype0, ktype1, ktype2, ktype3, id0, id1, id2, id3, model_status, colltype))


            if line[1]=='binary-binary':   ##Binary-binary collision
                colltype='BB'
                if int(line[2])==2:
                    if int(line[10]) in startype:
                        model=i; model_status=int(status[i]); timem=t_conv*float(line[0])
                        mm=float(line[4]); m0=float(line[6]); m1=float(line[8]); m2=-100; m3=-100
                        ktypem=int(line[10]); ktype0=int(line[11]); ktype1=int(line[12]); ktype2=-100; ktype3=-100
                        idm=int(line[3]); rm=float(line[9])*l_conv
                        id0=int(line[5]); id1=int(line[7]); id2=-100; id3=-100

                        fcoll.write('%d %f %d %f %f %f %f %f %f %d %d %d %d %d %d %d %d %d %d %s\n'%(model, timem, idm, rm, mm, m0, m1, m2, m3, ktypem, ktype0, ktype1, ktype2, ktype3, id0, id1, id2, id3, model_status, colltype))

                if int(line[2])==3:
                    if int(line[12]) in startype:
                        model=i; model_status=int(status[i]); timem=t_conv*float(line[0])
                        mm=float(line[4]); m0=float(line[6]); m1=float(line[8]); m2=float(line[10]); m3=-100
                        ktypem=int(line[12]); ktype0=int(line[13]); ktype1=int(line[14]); ktype2=int(line[15]); ktype3=-100
                        idm=int(line[3]); rm=float(line[11])*l_conv
                        id0=int(line[5]); id1=int(line[7]); id2=int(line[9]); id3=-100

                        fcoll.write('%d %f %d %f %f %f %f %f %f %d %d %d %d %d %d %d %d %d %d %s\n'%(model, timem, idm, rm, mm, m0, m1, m2, m3, ktypem, ktype0, ktype1, ktype2, ktype3, id0, id1, id2, id3, model_status, colltype))

                if int(line[2])==4:
                    if int(line[14]) in startype:
                        model=i; model_status=int(status[i]); timem=t_conv*float(line[0])
                        mm=float(line[4]); m0=float(line[6]); m1=float(line[8]); m2=float(line[10]); m3=float(line[12])
                        ktypem=int(line[14]); ktype0=int(line[15]); ktype1=int(line[16]); ktype2=int(line[17]); ktype3=int(line[18])
                        idm=int(line[3]); rm=float(line[13])*l_conv
                        id0=int(line[5]); id1=int(line[7]); id2=int(line[9]); id3=int(line[11])

                        fcoll.write('%d %f %d %f %f %f %f %f %f %d %d %d %d %d %d %d %d %d %d %s\n'%(model, timem, idm, rm, mm, m0, m1, m2, m3, ktypem, ktype0, ktype1, ktype2, ktype3, id0, id1, id2, id3, model_status, colltype))


        print(i)

    fcoll.close()


##Find the number of XX-XX mergers in the models
##klim1 and klim2 are lists
def get_XX_mergerproduct(pathlist, start, end, startype, savepath, readflag):
    if readflag == 1:
        sourcedir=np.genfromtxt(pathlist, dtype=str)
        filepaths=sourcedir[:,0]; status=sourcedir[:,1]
    else:
        filepaths=pathlist; status=[1]

    #fmerge=open(savepath+'_mergerproduct_all.dat', 'w+')
    #fmerge.write('#1.Model 2.Time(Myr) 3.IDmerg 4.Radius(pc) 5.Mmerg 6.M0 7.M1 8.kmerg 9.k0 10.k1 11.id0 12.id1 13.model_status\n')
    for i in range(start, end):
        fmerge=open(filepaths[i]+'ns_mergerproduct_all.dat', 'w+')
        fmerge.write('#1.Model 2.Time(Myr) 3.IDmerg 4.Radius(pc) 5.Mmerg 6.M0 7.M1 8.kmerg 9.k0 10.k1 11.id0 12.id1 13.model_status\n')

        filestr=filepaths[i]+'initial'

        t_conv=dyn.conv('t', filestr+'.conv.sh')
        l_conv=dyn.conv('l', filestr+'.conv.sh')

        semergefile = filestr+'.semergedisrupt.log'
        semergefile2 = filestr+'2.semergedisrupt.log'
        semergedata=scripts2.readmergefile(semergefile)
        if os.path.isfile(semergefile2) and os.path.getsize(semergefile2)>0:
            semergedata2=scripts2.readmergefile(semergefile2)
            semergedata=semergedata+semergedata2

        for hh in range(len(semergedata)):
            line=semergedata[hh].split()
            #print(line)
            if int(line[1])<3:
                if int(line[-3]) in startype:
                    fmerge.write('%d %f %d %f %f %f %f %d %d %d %d %d %d\n'%(i, t_conv*float(line[0]), 
                    int(line[2]), float(line[8])*l_conv, float(line[3]), float(line[5]), float(line[7]),int(line[-3]), int(line[-2]), int(line[-1]), int(line[4]), int(line[6]), int(status[i])))

        print(i)

    fmerge.close()


def coll_merger_hist(modelpath, IDinit_list, Tinit_list, startype):
    '''Input initial IDs, output the coll/mer history (IDs, time of coll/mer)'''
    '''Tinit_list is in code unit'''

    filestr=modelpath+'initial'

    t_conv=dyn.conv('t', filestr+'.conv.sh')

    collfile=filestr+'.collision.log'
    collfile2=filestr+'2.collision.log'
    colldata=scripts1.readcollfile(collfile)
    if os.path.isfile(collfile2) and os.path.getsize(collfile2) > 0:
        colldata2=scripts1.readcollfile(collfile2)
        colldata=colldata+colldata2
     
    timem = []; idm = []; ids = [[],[],[],[]]
    for jj in range(len(colldata)):
        line=colldata[jj].split()
        if line[1]=='single-single':  ##Single-single star collision
                colltype='SS'
                list_star=[int(line[11]), int(line[12])]
                num_star=Counter(list_star)
                if int(line[10])==startype and num_star[startype]==1:
                    #print(line)
                    timem.append(float(line[0]))
                    idm.append(int(line[3]))
                    ids[0].append(int(line[5])); ids[1].append(int(line[7]))
                    ids[2].append(-100); ids[3].append(-100)


        if line[1]=='binary-single':   ##Binary-single collision
            colltype='BS'
            if int(line[2])==2:
                list_star=[int(line[11]), int(line[12])]
                num_star=Counter(list_star)
                if int(line[10])==startype and num_star[startype]==1:
                    #print(line)
                    timem.append(float(line[0]))
                    idm.append(int(line[3]))
                    ids[0].append(int(line[5])); ids[1].append(int(line[7]))
                    ids[2].append(-100); ids[3].append(-100)


            if int(line[2])==3:
                list_star=[int(line[13]), int(line[14]), int(line[15])]
                num_star=Counter(list_star)
                if int(line[12])==startype and num_star[startype]==1:
                    #print(line)
                    timem.append(float(line[0]))
                    idm.append(int(line[3]))
                    ids[0].append(int(line[5])); ids[1].append(int(line[7]))
                    ids[2].append(int(line[9])); ids[3].append(-100)


        if line[1]=='binary-binary':   ##Binary-binary collision
            colltype='BB'
            if int(line[2])==2:
                list_star=[int(line[11]), int(line[12])]
                num_star=Counter(list_star)
                if int(line[10])==startype and num_star[startype]==1:
                    #print(line)
                    timem.append(float(line[0]))
                    idm.append(int(line[3]))
                    ids[0].append(int(line[5])); ids[1].append(int(line[7]))
                    ids[2].append(-100); ids[3].append(-100)


            if int(line[2])==3:
                list_star=[int(line[13]), int(line[14]), int(line[15])]
                num_star=Counter(list_star)
                if int(line[12])==startype and num_star[startype]==1:
                    #print(line)
                    timem.append(float(line[0]))
                    idm.append(int(line[3]))
                    ids[0].append(int(line[5])); ids[1].append(int(line[7]))
                    ids[2].append(int(line[9])); ids[3].append(-100)


            if int(line[2])==4:
                list_star=[int(line[15]), int(line[16]), int(line[17]), int(line[18])]
                num_star=Counter(list_star)
                if int(line[14])==startype and num_star[startype]==1:
                    #print(line)
                    timem.append(float(line[0]))
                    idm.append(int(line[3]))
                    ids[0].append(int(line[5])); ids[1].append(int(line[7]))
                    ids[2].append(int(line[9])); ids[3].append(int(line[11]))


    semergefile = filestr+'.semergedisrupt.log'
    semergefile2 = filestr+'2.semergedisrupt.log'
    semergedata=scripts2.readmergefile(semergefile)
    if os.path.isfile(semergefile2) and os.path.getsize(semergefile2)>0:
        semergedata2=scripts2.readmergefile(semergefile2)
        semergedata=semergedata+semergedata2

    for kk in range(len(semergedata)):
        line=semergedata[kk].split()
        if int(line[1])<3:
            list_star=[int(line[-2]),int(line[-1])]
            num_star=Counter(list_star)
            if int(line[-3])==startype and num_star[startype]==1:
                #print(line)
                timem.append(float(line[0]))
                idm.append(int(line[2]))
                ids[0].append(int(line[4])); ids[1].append(int(line[6]))
                ids[2].append(-100); ids[3].append(-100)

    print(modelpath, 'done extracting data')

    timem = np.array(timem)
    idm = np.array(idm)
    ids[0] = np.array(ids[0]); ids[1] = np.array(ids[1])
    ids[2] = np.array(ids[2]); ids[3] = np.array(ids[3])
    
    tindex_sort = np.argsort(timem)
    timem_sort = timem[tindex_sort]
    idm_sort = idm[tindex_sort]
    id0_sort = ids[0][tindex_sort]
    id1_sort = ids[1][tindex_sort]
    id2_sort = ids[2][tindex_sort]
    id3_sort = ids[3][tindex_sort]

    dict_allid = {}
    for xx in range(len(IDinit_list)):
        if IDinit_list[xx]==-100:
            continue
            
        #print(xx, IDinit_list[xx])

        oldidm = IDinit_list[xx]
        oldtime = Tinit_list[xx]

        timem_tcut = timem_sort[timem_sort>oldtime]
        idm_tcut = idm_sort[timem_sort>oldtime]
        id0_tcut = id0_sort[timem_sort>oldtime]
        id1_tcut = id1_sort[timem_sort>oldtime]
        id2_tcut = id2_sort[timem_sort>oldtime]
        id3_tcut = id3_sort[timem_sort>oldtime]

        #print(timem_tcut, idm_tcut)

        yy=0
        tm_list = []; idm_list = []
        check = 0
        while yy < len(timem_tcut):
            if oldidm in id0_tcut[yy:] or oldidm in id1_tcut[yy:] or oldidm in id2_tcut[yy:] or oldidm in id3_tcut[yy:]:
                check=1
                new_index_list = [np.where(oldidm==id0_tcut[yy:])[0], np.where(oldidm==id1_tcut[yy:])[0],
                                  np.where(oldidm==id2_tcut[yy:])[0], np.where(oldidm==id3_tcut[yy:])[0]]
                #print(new_index_list)
                actual_index_list = []
                for zz in range(len(new_index_list)):
                    if new_index_list[zz].size>0:
                        actual_index_list.append(new_index_list[zz][0])

                new_index = np.min(actual_index_list)
                #print(new_index)

                newtm = timem_tcut[yy:][new_index]
                newidm = idm_tcut[yy:][new_index]
                #print(newidm, newtm)
                idm_list.append(newidm)
                tm_list.append(newtm)

            #if new_index==0:
            #    break

            if not check:
                break
            else:
                oldidm = newidm
                if new_index==0:
                    new_index+=1
                yy=yy+new_index
                check=0

        #print(idm_list, tm_list)
        dict_allid[str(int(IDinit_list[xx]))] = dict(zip(tm_list,idm_list, ))

    return dict_allid