import numpy as np
import os,sys
import subprocess
import gzip
import scripts
import matplotlib.pyplot as plt
import history_cmc

def history_maker(id0, id1, path, snapno, time_i, savepath):
    ID_FLAG = 0
    MTflag = 1    
    id0_str = str(id0)
    f = open(savepath+'/'+id0_str+'_history.dat','w+')
    sys.stdout = f
    units=scripts.read_units(path+'/'+'initial')
    km = units[0]['l_cgs']*1.0e-5
    time_units = units[0]['t_myr']

    star_ids = [id0]
    positions = [1]
    filestring = 'initial'
    path = path
    binary = 1
    h = history_cmc.history_maker(star_ids,positions,filestring,path,binary)
    # First create the time array
    time_array = []
    id0 = id0
    for i in range(len(h[id0]['binint']['binint'])):
        time_array.append(h[id0]['binint']['binint'][i]['interaction']['type']['time']*time_units)
    print(time_array)
    print('t_onset=',time_i)
    time = h[id0]['binint']['binint']
        if snapno <= 300:
                Delta = -1
        if snapno <= 600 and snapno > 300:
                Delta = -4
        if snapno <= 1500 and snapno > 600:
                Delta = -6
        if snapno <= 10000 and snapno > 1500:
                Delta = -8
        if snapno > 10000:
                Delta = -10

    snapno_max = int(snapno)
    units=scripts.read_units(path+'/'+'initial')
    km = units[0]['l_cgs']*1.0e-5
    time_units = units[0]['t_myr']

    flag = snapno
    m0 = []        # 0 is the companion
    m1 = []        # 1 is the BH    
    ID0 = []
    Porb = []
    e = []
    timearr = []
    for j in range(snapno,0,Delta):
        if ID_FLAG == 1:
            print('stopped!')
            break
        N_XRB = 0
        N_BH = 0
        N_BHnBH = 0
        N_BBH = 0

        if j < 10:
            snapno = '000'+str(j)
        if j < 100 and j >= 10:
            snapno = '00'+str(j)
        if j < 1000 and j >= 100:
            snapno = '0'+str(j)
        if j < 10000 and j >= 1000:
            snapno = str(j)
        if j >= 10000:
            snapno = str(j)

        f2 = gzip.open(path+'/'+'initial.snap'+snapno+'.dat.gz','r')
        lines = f2.readlines()
        line = lines[0]
        parsed = line.split('[')
        parsed = parsed[0]
        parsed = parsed.split('=')
        time = float(parsed[1])*time_units      # DEFINE THE TIME OF THE SNAPFILE
        ###############################
        #if time > time_i:    # START RECORDING DATA BEFORE MT BEGINS
        #    break
        ###############################
                # DEFINE END TIME (TIME OF NEXT SNAPSHOT)
                if j-Delta < snapno_max:
                if j-Delta < 10:
                                snapno_f = '000'+str(j-Delta)
                        if j-Delta < 100 and j+1 >= 10:
                                snapno_f = '00'+str(j-Delta)
                        if j-Delta < 1000 and j+1 >= 100:
                                snapno_f = '0'+str(j-Delta)
                        if j-Delta < 10000 and j+1 >= 1000:
                                snapno_f = str(j-Delta)
                        if j-Delta >= 10000:
                                snapno_f = str(j-Delta)

                        f2_f = gzip.open(path+'/'+'initial.snap'+snapno_f+'.dat.gz','r')
                        lines2 = f2_f.readlines()
                        line2 = lines2[0]
                        parsed = line2.split('[')
                        parsed = parsed[0]
                        parsed = parsed.split('=')
                        time_f = float(parsed[1])*time_units      # DEFINE THE END TIME
                else:
                        snapno_max_str = snapno
            f2_f = gzip.open(path+'/'+'initial.snap'+snapno_max_str+'.dat.gz','r')
                        lines2 = f2_f.readlines()
                        line2 = lines2[0]
                        parsed = line2.split('[')
                        parsed = parsed[0]
                        parsed = parsed.split('=')
                        time_f = float(parsed[1])*time_units
                #########################################################################

    ################## determine if you had an interaction ##########    
        t_snap_i = time
        t_snap_f = time_f

    #########################
        TRIPLEflag = 0
        for n in range(len(time_array)-1,-1,-1):
            if time_array[n] >= t_snap_i and time_array[n] < t_snap_f:
                i = n
                timeI = h[id0]['binint']['binint'][i]['interaction']['type']['time']*time_units
                type = h[id0]['binint']['binint'][i]['interaction']['type']['type']
                for m in range(len(h[id0]['binint']['binint'][i]['interaction']['output'])):
                    for k in range(len(h[id0]['binint']['binint'][i]['interaction']['output'][m]['ids'])):
                    try:
                        VAR = float(h[id0]['binint']['binint'][i]['interaction']['output'][m]['ids'][k])
                        if VAR == id0:
                            for l in range(len(h[id0]['binint']['binint'][i]['interaction']['output'][m]['ids'])):
                                if float(h[id0]['binint']['binint'][i]['interaction']['output'][m]['ids'][l]) == id1:
                                    ############  TRIPLE STUFF   ######################
                                    for k in range(len(h[id0]['binint']['binint'][i]['interaction']['output'][m]['m'])):
                                        if k == 2:    # DETERMINE IF A TRIPLE IS FORMED
                                            TRIPLEflag = 1
                                    if TRIPLEflag == 1:
                                        TRIPLEflag = 0
                                        print(timeI, 3)
                                        m_1 = h[id0]['binint']['binint'][i]['interaction']['output'][m]['m'][0]
                                        m_2 = h[id0]['binint']['binint'][i]['interaction']['output'][m]['m'][1]
                                        m_S = h[id0]['binint']['binint'][i]['interaction']['output'][m]['m'][2]
                                        print(m_1, m_2)
                                        aI_i = h[id0]['binint']['binint'][i]['interaction']['output'][m]['a'][0] 
                                        ao = h[id0]['binint']['binint'][i]['interaction']['output'][m]['a'][1] 
                                        #print "aI_i=",aI_i,"ao=",ao,
                                        eo = h[id0]['binint']['binint'][i]['interaction']['output'][m]['e'][1]
                                        eI = h[id0]['binint']['binint'][i]['interaction']['output'][m]['e'][0]
                                        aI_f = pow(m_S*(m_1+m_2)/(m_1*m_2)/ao + 1./aI_i,-1.)    
                                        print(aI_f, eI)
                                        ##### Stability criteria
                                        Stability = 3.3/(1-eo)*pow(2./3.*(1+m_S/(m_1+m_2))*(1+eo)/(1-eo)**0.5,2./5.)    
                                        if ao/aI_i > Stability:
                                            print(MTflag, 1)
                                        else:
                                            print(MTflag, 0)                    


                                    ##########################
                                    print timeI,
                                                                        if type == 'BS':
                                                                                print 1,
                                                                      
                                    if type == 'BB':
                                        print 2,

                                    for k in range(0,len(h[id0]['binint']['binint'][i]['interaction']['output'][m]['m'])):
                                                                                if k <= 1:
                                                                                        print h[id0]['binint']['binint'][i]['interaction']['output'][m]['m'][k],
                                    print h[id0]['binint']['binint'][i]['interaction']['output'][m]['a'][0], h[id0]['binint']['binint'][i]['interaction']['output'][m]['e'][0], MTflag, 0
                                    print timeI, 0,
                                    for k in range(len(h[id0]['binint']['binint'][i]['interaction']['output'][m]['m'])):
                                                                                if k <= 1:
                                            print h[id0]['binint']['binint'][i]['interaction']['output'][m]['m'][k],

                    except Exception:
                        pass
                for m in range(len(h[id0]['binint']['binint'][i]['interaction']['input'])):
                                    for k in range(len(h[id0]['binint']['binint'][i]['interaction']['input'][m]['ids'])):
                                        try:
                                                VAR = float(h[id0]['binint']['binint'][i]['interaction']['input'][m]['ids'][k])
                                                if VAR == id0:
                                                        for l in range(len(h[id0]['binint']['binint'][i]['interaction']['input'][m]['ids'])):
                                                                if float(h[id0]['binint']['binint'][i]['interaction']['input'][m]['ids'][l]) == id1:
                                                                        print h[id0]['binint']['binint'][i]['interaction']['input'][m]['a'], h[id0]['binint']['binint'][i]['interaction']['input'][m]['e'], MTflag, 0
                    except Exception:
                                                pass 
    ###################
    #        if time_array[n] >= t_snap_i and time_array[n] < t_snap_f:
    #            timeI = h[id0]['binint']['binint'][i]['interaction']['type']['time']*time_units
    #            type = h[id0]['binint']['binint'][i]['interaction']['type']['type']
    #            print timeI, type
    #            print '     ','OUTPUT:'
    #            for j in range(len(h[id0]['binint']['binint'][i]['interaction']['output'])):
    #                print '     ',h[id0]['binint']['binint'][i]['interaction']['output'][j]['type'],'(',
    #                for k in range(len(h[id0]['binint']['binint'][i]['interaction']['output'][j]['ids'])):
    #                print h[id0]['binint']['binint'][i]['interaction']['output'][j]['ids'][k],
    #                print ')',
    #                for k in range(len(h[id0]['binint']['binint'][i]['interaction']['output'][j]['m'])):
    #                print 'mass:',h[id0]['binint']['binint'][i]['interaction']['output'][j]['m'][k],
    #                print 'a:',h[id0]['binint']['binint'][i]['interaction']['output'][j]['a'],
    #                print 'e:',h[id0]['binint']['binint'][i]['interaction']['output'][j]['e']
    #            print '     ','INPUT:'
    #            for j in range(len(h[id0]['binint']['binint'][i]['interaction']['input'])):
    #                print '     ',h[id0]['binint']['binint'][i]['interaction']['input'][j]['type'],'(',
    #                for k in range(len(h[id0]['binint']['binint'][i]['interaction']['input'][j]['ids'])):
    #                print h[id0]['binint']['binint'][i]['interaction']['input'][j]['ids'][k],
    #                print ')',
    #                for k in range(len(h[id0]['binint']['binint'][i]['interaction']['input'][j]['m'])):
    #                print 'mass:',h[id0]['binint']['binint'][i]['interaction']['input'][j]['m'][k],
    #                print 'a:',h[id0]['binint']['binint'][i]['interaction']['input'][j]['a'],
    #                print 'e:',h[id0]['binint']['binint'][i]['interaction']['input'][j]['e']
    ##############################
        for i in range(2,len(lines)):
            ID_FLAG = 1
            data = lines[i]
            data = data.split(' ')
            MTflag = 0
            if float(data[7]) == 1.:
                M0 = float(data[8])
                M1 = float(data[9])
                q0 = M0/M1
                q1 = M1/M0
                a = float(data[12])
                ecc = float(data[13])
                RL_Egg_0 = a*0.49*pow(q0,2/3.)/(0.6*pow(q0,2/3.) + np.log(1+pow(q0,1/3.)))
                RL_Egg_1 = a*0.49*pow(q1,2/3.)/(0.6*pow(q1,2/3.) + np.log(1+pow(q1,1/3.)))
                R0 = float(data[21])*0.00465 # Radius of star0 converted from Rsun to AU
                R1 = float(data[22])*0.00465 # " "
                if float(data[10]) == id0 and float(data[11]) == id1:   # asks if still member of same binary
                    if R0 > RL_Egg_0 or R1 > RL_Egg_1: 
                        MTflag = 1
                    m0.append(float(data[8]))
                    m1.append(float(data[9]))
                    ID0.append(float(data[17]))
                    period = pow(4*np.pi*np.pi*pow(float(data[12])*1.5e11,3.)/(6.67e-11*1.99e30*(float(data[8]) + float(data[9]))),1./2.)        
                    Porb.append(period)
                    e.append(float(data[13]))
                    timearr.append(time)
                    ID_FLAG = 0
                    print time, 0, float(data[8]), float(data[9]), data[12], data[13], MTflag, snapno 
                    break
                if float(data[11]) == id0 and float(data[10]) == id1:
                    if R0 > RL_Egg_0 or R1 > RL_Egg_1: 
                        MTflag = 1
                    m0.append(float(data[9]))
                    m1.append(float(data[8]))
                    ID0.append(float(data[18]))
                    period = pow(4*np.pi*np.pi*pow(float(data[12])*1.5e11,3.)/(6.67e-11*1.99e30*(float(data[8]) + float(data[9]))),1./2.)    
                    Porb.append(period)
                    e.append(float(data[13]))
                    timearr.append(time)
                    ID_FLAG = 0
                    print time, 0, float(data[8]), float(data[9]), data[12], data[13], MTflag, snapno
                    break
    m_ZAMS = 0
    return m0[-1], m1[-1], m_ZAMS, ID0[-1], Porb[-1], e[-1], timearr[-1]
