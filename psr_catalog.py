import numpy as np
from glob import glob
import collections
from collections import Counter
import os,sys
import subprocess
import matplotlib
import pandas as pd
import gzip
import math
import re
#import history_maker_full4 as hi4
import history_cmc as hic
import dynamics as dyn
import scripts3
import scripts1
import scripts2
#from scipy import stats

yearsc=31557600.
twopi=6.283185307179586
Gconst=6.674*10**-8 ##cm3*g-1*s-2
clight=3*10**10 ##cm/s
Msun=2*10**33 ##gram
AU=1.496*10**13  ##cm
PC=3.086*10**18  ##cm
Kconst=9.87*10**-48 ##yr/G^2
Lsun=4.02*10**16 ##mJy*kpc^2



def readdata_freire():
    #from astropy.extern.six.moves.urllib import request
    #url = 'http://www.naic.edu/~pfreire/GCpsr.txt'
    #open('/projects/b1095/syr904/projects/PULSAR/GC_psr.txt', 'wb').write(request.urlopen(url).read())
    
    K=5.7*10**19
    p=[]; pdot=[]; binflag=[]; namespin=[] #Ps=[]; Pdots=[]; Bs=[]; Pb=[]; Pdotb=[]; Bb=[]  ##P unit ms, dpdt has a factor of 10^-20
    Ns=0; Nb=0.  ##Calculating the numbers of single and binary pulsars
    period=[]; ecc=[]; mc=[]; names=[]
    
    pall=[]; bfall=[]; nameall=[]

    ntot=0

    with open('/projects/b1095/syr904/projects/PULSAR/data_observed/GC_psr.txt', 'rb') as f:
        for _ in range(4):
            next(f)
        for line in f:
            datanew=line.split()
            #print(data)
            data=[]
            for item in datanew:
                data.append(re.sub(r'\([^)]*\)', '', item.decode('utf-8')))
            #print(datanew)
            #print(data)
            if not data: continue
            if str(data[0][0])=='J' or data[0][0]=='B':
                ntot+=1
                ##Calculate numbers
                if str(data[5])=='i': 
                    Ns+=1; pall.append(float(data[2])); bfall.append(0); nameall.append(data[0])
                elif str(data[5])!='i' and str(data[5])!='*': 
                    Nb+=1; pall.append(float(data[2])); bfall.append(1); nameall.append(data[0])
                else:
                    pall.append(float(data[2])); bfall.append(-100); nameall.append(data[0]) 
             
                ##Extract orbital data
                if str(data[5])!='i' and str(data[5])!='*':
                    #if str(data[5][0])=='<' or str(data[7][0])=='<' or str(data[5][0])=='>' or str(data[7][0])=='>': 
                    #    continue
                    #else: 
                    period.append(data[5]); ecc.append(data[7]); mc.append(data[8])
                    names.append(data[0])


                ##Extract spin period data
                if str(data[3])!='*':
                    if str(data[3][0])=='<':
                        ##May want to include this in the future (there is just one of them)
                        continue        
                    if str(data[5])=='i':
                        namespin.append(data[0]) 
                        binflag.append(0)
                        p.append(float(data[2]))
                        #dpdts=data[3].split('(')
                        #if len(dpdts)>1:
                        #    errs=dpdts[1].split(')')
                        #    if errs[1]!='':
                        #        Pdots.append(float(dpdts[0])*10**-15)
                        #    else: 
                        #        Pdots.append(float(dpdts[0])*10**-20)
                        #else: 
                        #    Pdots.append(float(dpdts[0])*10**-20)
                        dpdts=data[3]
                        if '*' in dpdts:
                            pdot.append(float(dpdts.split('*')[0])*10**-15)
                        else:
                            pdot.append(float(dpdts)*10**-20)
                    else: 
                        namespin.append(data[0]) 
                        binflag.append(1)
                        p.append(float(data[2]))
                        #dpdtb=data[3].split('(')
                        #if len(dpdtb)>1:
                        #    errb=dpdtb[1].split(')')
                        #    if errb[1]!='':
                        #        Pdotb.append(float(dpdtb[0])*10**-15)
                        #    else:
                        #        Pdotb.append(float(dpdtb[0])*10**-20)
                        #else:
                        #    Pdotb.append(float(dpdtb[0])*10**-20)
                        dpdtb=data[3]
                        if '*' in dpdtb:
                            pdot.append(float(dpdtb.split('*')[0])*10**-15)
                        else:
                            pdot.append(float(dpdtb)*10**-20)


    ##print Pdots, Pdotb, Ps, Pb
    #Bs=K*np.sqrt(np.abs(Pdots)*np.array(Ps)*0.001)
    #Bb=K*np.sqrt(np.abs(Pdotb)*np.array(Pb)*0.001)
    #Ps=np.array(Ps); Pb=np.array(Pb); Pdots=np.array(Pdots); Pdotb=np.array(Pdotb)
    #Bs=np.array(Bs); Bb=np.array(Bb)
    ##print Bs, Bb
    #print(pdot)
    print(ntot)
    return p, pdot, binflag, namespin, period, ecc, mc, names, pall, bfall, nameall
