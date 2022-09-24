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
#import ecc_calc as gwcalc
#import LISA_calculations as lisa
#import ns_history as nh
#import useful_function as uf


#sys.path.insert(1, '/projects/b1095/syr904/MyCodes/cmctoolkit')
#import cmctoolkit as cmct


##Output data for NS history
def output_psr(time_code, data_string):
    return {'t_form': time_code, 'id0':int(data_string[3]),'id1':int(data_string[4]), 
            'M0':float(data_string[5]), 'M1':float(data_string[6]),
            'B0':float(data_string[7]), 'B1':float(data_string[8]),
            'P0':float(data_string[9]), 'P1':float(data_string[10]),
            'K0':float(data_string[11]), 'K1':float(data_string[12]),
            'sma':float(data_string[13]), 'ecc':float(data_string[14]),
            'radrol0':float(data_string[15]), 'radrol1':float(data_string[16]),
            'dmdt0':float(data_string[17]), 'dmdt1':float(data_string[18]),
            'r':float(data_string[19]),
            'vr':float(data_string[20]), 'vt':float(data_string[21]),
            'bacc0':float(data_string[22]), 'bacc1':float(data_string[23]),
            'tacc0':float(data_string[24]), 'tacc1':float(data_string[25])}



##Find the formation info of NS and MSP
def output_formation_NS_MSP(data_str, check_ns, check_msp, theid):
    NS_form = {}; MSP_form = {}
    if int(data_str[2])!=1:
        if check_ns == 0 and int(data_str[11])==13 and int(data_str[3])==theid:
            check_ns = 1
            t_nsform = float(data_str[1])
            NS_form =  output_psr(t_nsform, data_str)
           
        Pspin=float(data_str[9])  ##in sec
        B=float(data_str[7])
        deathcut=(Pspin**2)*(0.17*10**12)
        if check_msp == 0 and deathcut<B and Pspin<=0.03 and int(data_str[3])==theid:
            check_msp = 1
            t_mspform = float(data_str[1])
            MSP_form = output_psr(t_mspform, data_str)
    
    else:
        if int(data_str[11]) == 13:
            if check_ns == 0 and int(data_str[3])==theid:
                check_ns = 1
                t_nsform = float(data_str[1])
                NS_form =  output_psr(t_nsform, data_str)

            Pspin=float(data_str[9])  ##in sec
            B=float(data_str[7])
            deathcut=(Pspin**2)*(0.17*10**12)
            if check_msp == 0 and deathcut<B and Pspin<=0.03 and int(data_str[3])==theid:
                check_msp = 1
                t_mspform = float(data_str[1])
                MSP_form = output_psr(t_mspform, data_str)
                            
        if int(data_str[12]) == 13:
            if check_ns == 0 and int(data_str[4])==theid:
                check_ns = 1
                t_nsform = float(data_str[1])
                NS_form =  output_psr(t_nsform, data_str)

            Pspin=float(data_str[10])  ##in sec
            B=float(data_str[8])
            deathcut=(Pspin**2)*(0.17*10**12)
            if check_msp == 0 and deathcut<B and Pspin<=0.03 and int(data_str[4])==theid:
                check_msp = 1
                t_mspform = float(data_str[1])
                MSP_form = output_psr(t_mspform, data_str)

    return check_ns, check_msp, NS_form, MSP_form


##Find all collision involving the NS
def output_coll_NS(datacoll, theid):
    output_coll = {}
    for xx in range(len(datacoll)):
        linecoll = datacoll[xx].split()
        if int(linecoll[3])==theid:
            colltime=float(linecoll[0])
            colltype = linecoll[1]; collnum=int(linecoll[2])
            collmm=float(linecoll[4])
            if collnum==2:
                collids = [int(linecoll[5]), int(linecoll[7])]
                collmktype=int(linecoll[10])
                startypes = [int(linecoll[11]), int(linecoll[12])]
                masses = [float(linecoll[6]), float(linecoll[8])]

            elif collnum==3:
                collids=[int(linecoll[5]), int(linecoll[7]), int(linecoll[9])]
                collmktype=int(linecoll[12])
                startypes = [int(linecoll[13]), int(linecoll[14]), int(linecoll[15])]
                masses = [float(linecoll[6]), float(linecoll[8]), float(linecoll[10])]

            elif collnum==4:
                collids=[int(linecoll[5]), int(linecoll[7]), int(linecoll[9]), int(linecoll[11])]
                collmktype=int(linecoll[14])
                startypes = [int(linecoll[15]), int(linecoll[16]), int(linecoll[17]), int(linecoll[18])]
                masses = [float(linecoll[6]), float(linecoll[8]), float(linecoll[10]), float(linecoll[12])]
    
            output_coll = {'t_coll':colltime, 'type_coll':colltype, 'num_coll':collnum,
            'mid':int(linecoll[3]), 'mm':collmm, 'mk':collmktype,
            'ids':collids, 'masses':masses, 'startypes':startypes}

            break


    return output_coll


##Find all binary merger/disruption of the NS
def output_semer_NS(datasemer, theid):
    output_semer={}
    for xx in range(len(datasemer)):
        linesemer=datasemer[xx].split()
        check_semer = 0
        if int(linesemer[1])<3:
            if int(linesemer[2])==theid:
                check_semer=1
                t_semer=float(linesemer[0]); type_semer=int(linesemer[1])
                mr=float(linesemer[3]); kr=int(linesemer[9])
                semerids = [int(linesemer[4]), int(linesemer[6])]
                semermass = [float(linesemer[5]), float(linesemer[7])]
                rsemer=float(linesemer[8])
                semerktype=[int(linesemer[10]), int(linesemer[11])]

        else:
            if int(linesemer[2])==theid:
                check_semer=1
                t_semer=float(linesemer[0]); type_semer=int(linesemer[1])
                mr=float(linesemer[3]); kr=int(linesemer[7])
                semerids=int(linesemer[4]); semermass=float(linesemer[5])
                rsemer=float(linesemer[6])
                semerktype=int(linesemer[8])

            if int(linesemer[4])==theid:
                check_semer=1
                t_semer=float(linesemer[0]); type_semer=int(linesemer[1])
                mr=float(linesemer[5]); kr=int(linesemer[8])
                semerids=int(linesemer[2]); semermass=float(linesemer[3])
                rsemer=float(linesemer[6])
                semerktype=int(linesemer[7])

        if check_semer==1:
            output_semer[str(t_semer)]={'type':type_semer, 'rid':theid, 'rm':mr, 'rk':kr,
                                        'ids':semerids, 'masses':semermass, 'startypes':semerktype,
                                        'r_pos':rsemer}

    return output_semer


def output_tc_gc_NS(sourcepath, theid):
    output_tc_gc={}
    property_init, property_finl, property_des = find_tc_properties_final(sourcepath)
    for xx in range(len(prop_init['id0'])):
        if int(prop_init['id0'][xx])==theid:
            period = uc.au_to_period(prop_finl['sma'][xx], prop_finl['m0'][xx], prop_finl['m1'][xx])
                
            if prop_init['type'][xx] == 'SS_COLL_TC_P':
                output_tc_gc[str(prop_init['time'][xx])] = {'type':'TC', 'ids':[prop_init['id0'][xx], prop_init['id1'][xx]], 'kinis':[prop_init['k0'][xx],prop_init['k1'][xx]], 'minis':[prop_init['m0'][xx],prop_init['m1'][xx]], 'rperi':prop_init['rperi'][xx], 'v_inf':prop_init['vinf'][xx], 'kfnl':[prop_finl['k0'][xx],prop_finl['k1'][xx]], 'mfnls':[prop_finl['m0'][xx],prop_finl['m1'][xx]], 'sma':prop_finl['sma'][xx], 'ecc':prop_finl['ecc'][xx], 'period':period}

            if prop_init['type'][xx] == 'SS_COLL_Giant':
                output_tc_gc[str(prop_init['time'][xx])] = {'type':'GC', 'ids':[prop_init['id0'][xx], prop_init['id1'][xx]], 'kinis':[prop_init['k0'][xx],prop_init['k1'][xx]], 'minis':[prop_init['m0'][xx],prop_init['m1'][xx]], 'rperi':prop_init['rperi'][xx], 'v_inf':prop_init['vinf'][xx], 'kfnl':[prop_finl['k0'][xx],prop_finl['k1'][xx]], 'mfnls':[prop_finl['m0'][xx],prop_finl['m1'][xx]], 'sma':prop_finl['sma'][xx], 'ecc':prop_finl['ecc'][xx], 'period':period}

        elif int(prop_init['id1'][xx])==theid:
            period = uc.au_to_period(prop_finl['sma'][xx], prop_finl['m0'][xx], prop_finl['m1'][xx])
                
            if prop_init['type'][xx] == 'SS_COLL_TC_P':
                output_tc_gc[str(prop_init['time'][xx])] = {'type':'TC', 'ids':[prop_init['id1'][xx], prop_init['id0'][xx]], 'kinis':[prop_init['k1'][xx],prop_init['k0'][xx]], 'minis':[prop_init['m1'][xx],prop_init['m0'][xx]], 'rperi':prop_init['rperi'][xx], 'v_inf':prop_init['vinf'][xx], 'kfnl':[prop_finl['k1'][xx],prop_finl['k0'][xx]], 'mfnls':[prop_finl['m1'][xx],prop_finl['m0'][xx]], 'sma':prop_finl['sma'][xx], 'ecc':prop_finl['ecc'][xx], 'period':period}

            if prop_init['type'][xx] == 'SS_COLL_Giant':
                output_tc_gc[str(prop_init['time'][xx])] = {'type':'GC', 'ids':[prop_init['id1'][xx], prop_init['id0'][xx]], 'kinis':[prop_init['k1'][xx],prop_init['k0'][xx]], 'minis':[prop_init['m1'][xx],prop_init['m0'][xx]], 'rperi':prop_init['rperi'][xx], 'v_inf':prop_init['vinf'][xx], 'kfnl':[prop_finl['k1'][xx],prop_finl['k0'][xx]], 'mfnls':[prop_finl['m1'][xx],prop_finl['m0'][xx]], 'sma':prop_finl['sma'][xx], 'ecc':prop_finl['ecc'][xx], 'period':period}

    return output_tc_gc


##Find binary encounters throughout the cluster evolution
def output_binint_NS(sourcepath, theid):
    output_binint = {}

    binintfile = sourcepath+'initial.binint.log'
    fbinint=open(binintfile,'r')
    positions=scripts3.find_positions(fbinint)
    for xx in range(len(positions)-1):
        binint=scripts3.read_segment(fbinint,positions[xx])
        bininput = binint['input']
        binoutput = binint['output']
        for yy in range(len(bininput)):
            if int(bininput[yy]['no'])==1 and int(bininput[yy]['ids'][0])==theid:
                t_binint = float(binint['type']['time'])
                for zz in range(len(binoutput)):
                    if int(binoutput[zz]['no'])==1 and ':' in binoutput[zz]['ids'][0]:
                        binoutput_ids = binoutput[zz]['ids'][0].split(':')
                    else:
                        binoutput_ids = [binoutput[zz]['ids'][0]]

                    if int(binoutput[zz]['no'])==2:
                        if ':' in binoutput[zz]['ids'][0]:
                            binoutput_id0 = binoutput[zz]['ids'][0].split(':')
                        else:
                            binoutput_id0 = [binoutput[zz]['ids'][0]]
                        if ':' in binoutput[zz]['ids'][1]:
                            binoutput_id1 = binoutput[zz]['ids'][1].split(':')
                        else:
                            binoutput_id1 = [binoutput[zz]['ids'][1]]

                    if int(binoutput[zz]['no'])==3:
                        if ':' in binoutput[zz]['ids'][0]:
                            binoutput_id0 = binoutput[zz]['ids'][0].split(':')
                        else:
                            binoutput_id0 = [binoutput[zz]['ids'][0]]
                        if ':' in binoutput[zz]['ids'][1]:
                            binoutput_id1 = binoutput[zz]['ids'][1].split(':')
                        else:
                            binoutput_id1 = [binoutput[zz]['ids'][1]]
                        if ':' in binoutput[zz]['ids'][2]:
                            binoutput_id2 = binoutput[zz]['ids'][2].split(':')
                        else:
                            binoutput_id2 = [binoutput[zz]['ids'][2]]

                    if int(binoutput[zz]['no'])==1 and str(theid) in binoutput_ids:
                        output_binint[str(t_binint)] = {'typeout':int(binoutput[zz]['no']), 'kouts':[int(binoutput[zz]['startype'][0]),-100], 'mouts':[float(binoutput[zz]['m'][0]),-100], 'idouts':[binoutput[zz]['ids'][0],str(-100)], 'smaeccout':[-100,-100], 'typein':int(bininput[yy]['no']), 'kins':[int(bininput[yy]['startype'][0]),-100], 'mins':[float(bininput[yy]['m'][0]),-100], 'idins':[int(bininput[yy]['ids'][0]),-100], 'smaeccin':[-100,-100]}
                        
                    elif int(binoutput[zz]['no'])==2 and (str(theid) in binoutput_id0 or str(theid) in binoutput_id1):
                        output_binint[str(t_binint)] = {'typeout':int(binoutput[zz]['no']), 'kouts':binoutput[zz]['startype'], 'mouts':binoutput[zz]['m'], 'idouts':binoutput[zz]['ids'], 'smaeccout':[float(binoutput[zz]['a'][0]),float(binoutput[zz]['e'][0])], 'typein':int(bininput[yy]['no']), 'kins':[int(bininput[yy]['startype'][0]),-100], 'mins':[float(bininput[yy]['m'][0]),-100], 'idins':[int(bininput[yy]['ids'][0]),-100], 'smaeccin':[-100,-100]}

                    elif int(binoutput[zz]['no'])==3 and (str(theid) in binoutput_id0 or str(theid) in binoutput_id1):
                        output_binint[str(t_binint)] = {'typeout':int(binoutput[zz]['no']), 'kouts':binoutput[zz]['startype'][:2], 'mouts':binoutput[zz]['m'][:2], 'idouts':binoutput[zz]['ids'][:2], 'smaeccout':[float(binoutput[zz]['a'][0]),float(binoutput[zz]['e'][0])], 'typein':int(bininput[yy]['no']), 'kins':[int(bininput[yy]['startype'][0]),-100], 'mins':[float(bininput[yy]['m'][0]),-100], 'idins':[int(bininput[yy]['ids'][0]),-100], 'smaeccin':[-100,-100]}

                    elif int(binoutput[zz]['no'])==3 and str(theid) in binoutput_id2:
                        output_binint[str(t_binint)] = {'typeout':int(binoutput[zz]['no']), 'kouts':[binoutput[zz]['startype'][-1],-100], 'mouts':[binoutput[zz]['m'][-1],-100], 'idouts':[binoutput[zz]['ids'][-1],-100], 'smaeccout':[float(binoutput[zz]['a'][1]),float(binoutput[zz]['e'][1])], 'typein':int(bininput[yy]['no']), 'kins':[int(bininput[yy]['startype'][0]),-100], 'mins':[float(bininput[yy]['m'][0]),-100], 'idins':[int(bininput[yy]['ids'][0]),-100], 'smaeccin':[-100,-100]}

            elif int(bininput[yy]['no'])==2 and (int(bininput[yy]['ids'][0])==theid or int(bininput[yy]['ids'][1])==theid):
                t_binint = float(binint['type']['time'])
                for zz in range(len(binoutput)):
                    if int(binoutput[zz]['no'])==1 and ':' in binoutput[zz]['ids'][0]:
                        binoutput_ids = binoutput[zz]['ids'][0].split(':')
                    else:
                        binoutput_ids = [binoutput[zz]['ids'][0]]

                    if int(binoutput[zz]['no'])==2:
                        if ':' in binoutput[zz]['ids'][0]:
                            binoutput_id0 = binoutput[zz]['ids'][0].split(':')
                        else:
                            binoutput_id0 = [binoutput[zz]['ids'][0]]
                        if ':' in binoutput[zz]['ids'][1]:
                            binoutput_id1 = binoutput[zz]['ids'][1].split(':')
                        else:
                            binoutput_id1 = [binoutput[zz]['ids'][1]]

                    if int(binoutput[zz]['no'])==3:
                        if ':' in binoutput[zz]['ids'][0]:
                            binoutput_id0 = binoutput[zz]['ids'][0].split(':')
                        else:
                            binoutput_id0 = [binoutput[zz]['ids'][0]]
                        if ':' in binoutput[zz]['ids'][1]:
                            binoutput_id1 = binoutput[zz]['ids'][1].split(':')
                        else:
                            binoutput_id1 = [binoutput[zz]['ids'][1]]
                        if ':' in binoutput[zz]['ids'][2]:
                            binoutput_id2 = binoutput[zz]['ids'][2].split(':')
                        else:
                            binoutput_id2 = [binoutput[zz]['ids'][2]]

                    if int(binoutput[zz]['no'])==1 and str(theid) in binoutput_ids:
                        output_binint[str(t_binint)] = {'typeout':int(binoutput[zz]['no']), 'kouts':[int(binoutput[zz]['startype'][0]),-100], 'mouts':[float(binoutput[zz]['m'][0]),-100], 'idouts':[int(binoutput[zz]['ids'][0]),-100], 'smaeccout':[-100,-100], 'typein':int(bininput[yy]['no']), 'kins':bininput[yy]['startype'], 'mins':bininput[yy]['m'], 'idins':bininput[yy]['ids'], 'smaeccin':[float(bininput[yy]['a']),float(bininput[yy]['e'])]}

                    if int(binoutput[zz]['no'])==2 and (str(theid) in binoutput_id0 or str(theid) in binoutput_id1):
                        output_binint[str(t_binint)] = {'typeout':int(binoutput[zz]['no']), 'kouts':binoutput[zz]['startype'], 'mouts':binoutput[zz]['m'], 'idouts':binoutput[zz]['ids'], 'smaeccout':[float(binoutput[zz]['a'][0]),float(binoutput[zz]['e'][0])], 'typein':int(bininput[yy]['no']), 'kins':bininput[yy]['startype'], 'mins':bininput[yy]['m'], 'idins':bininput[yy]['ids'], 'smaeccin':[float(bininput[yy]['a']),float(bininput[yy]['e'])]}

                    elif int(binoutput[zz]['no'])==3 and (str(theid) in binoutput_id0 or str(theid) in binoutput_id1):
                        output_binint[str(t_binint)] = {'typeout':int(binoutput[zz]['no']), 'kouts':binoutput[zz]['startype'][:2], 'mouts':binoutput[zz]['m'][:2], 'idouts':binoutput[zz]['ids'][:2], 'smaeccout':[float(binoutput[zz]['a'][0]),float(binoutput[zz]['e'][0])], 'typein':int(bininput[yy]['no']), 'kins':bininput[yy]['startype'], 'mins':bininput[yy]['m'], 'idins':bininput[yy]['ids'], 'smaeccin':[float(bininput[yy]['a']),float(bininput[yy]['e'])]}

                    elif int(binoutput[zz]['no'])==3 and str(theid) in binoutput_id2:
                        output_binint[str(t_binint)] = {'typeout':int(binoutput[zz]['no']), 'kouts':[binoutput[zz]['startype'][-1],-100], 'mouts':[binoutput[zz]['m'][-1],-100], 'idouts':[binoutput[zz]['ids'][-1],-100], 'smaeccout':[float(binoutput[zz]['a'][1]),float(binoutput[zz]['e'][1])], 'typein':int(bininput[yy]['no']), 'kins':bininput[yy]['startype'], 'mins':bininput[yy]['m'], 'idins':bininput[yy]['ids'], 'smaeccin':[float(bininput[yy]['a']),float(bininput[yy]['e'])]}

    fbinint.close()

    return output_binint


#def output_pribinary_info(sourcepath,theid,compid):
#    pribinary = 0
#    snap0 = cmct.Snapshot(fname=sourcepath+'initial.snapshots.h5', 
#                snapshot_name='/0(t=0)', conv=sourcepath+'initial.conv.sh')
#    binflag = np.array(snap0.data['binflag'])
#    id0_bin = np.array(snap0.data['id0'])[binflag == 1]; id1_bin = np.array(snap0.data['id1'])[binflag == 1]
#
#    for xx in range(len(id0_bin)):
#        if (theid ==  int(id0_bin[xx]) and compid == int(id1_bin[xx])) or (theid ==  int(id1_bin[xx]) and compid == int(id0_bin[xx])):
#            pribinary = 1
#
#    return pribinary


def output_pribinary_info_oldcmc(sourcepath,theid,compid):
    pribinary = 0
    output_pribinary = {}
    filestr = sourcepath+'initial'
    firstsnap = filestr+'.snap0000.dat.gz'
    

    bid0 = []; bid1 = []
    with gzip.open(firstsnap, 'r') as fsnap:
        next(fsnap); next(fsnap)
        for line in fsnap:
            datasnap = line.split()
            if int(datasnap[7])==1:
                if (int(datasnap[10])==theid and int(datasnap[11])==compid) or (int(datasnap[11])==theid and int(datasnap[10])==compid):
                    pribinary=1
                    output_pribinary={'priflag':pribinary, 'id0':int(datasnap[10]), 'id1':int(datasnap[11]), 'm0[MSUN]':float(datasnap[8]), 'm1[MSUN]':float(datasnap[9]), 'sma[AU]':float(datasnap[12]), 'ecc':float(datasnap[13])}

                    break
            else:
                if int(datasnap[0])==theid:
                    output_pribinary={'priflag':pribinary, 'id0':int(datasnap[0]), 'id1':-100, 'm0[MSUN]':float(datasnap[1]), 'm1[MSUN]':-100, 'sma[AU]':-100, 'ecc':-100}


    return output_pribinary


def convert_str_to_four_digits(value):
    """
    Take the four digit name of snapshot file, add 3, and return four digit name
    of the next snapshot file
    """
    new_value = int(value)+5
    new_len = len(str(new_value))
    num_zeroes = 4-new_len
    return num_zeroes*'0'+str(new_value)



def check_snapshots_oldcmc(sourcepath, theid, NS_tform):
    """
    Check every 3 snapshot files up to when the NS is formed for the system 
    (timestep < NS_timestep) with IDs theid and compid. 
    Output: return all of the information regarding that system as a list of dictionaries
    """
    filestr = sourcepath+'initial'
    snap = filestr+'.snap0000.dat.gz'
    snap_num = '0000'
    t_current = 0
    print(NS_tform)
    output_snap = {}


    while t_current < NS_tform:
        t_current = dyn.get_time(snap)
        with gzip.open(snap, 'r') as fsnap:
            next(fsnap); next(fsnap)
            for line in fsnap:
                datasnap = line.split()
                if int(datasnap[7])==1:
                    if int(datasnap[10])==theid or int(datasnap[11])==theid:
                        output_snap[str(t_current)]={'id0':int(datasnap[10]), 'id1':int(datasnap[11]), 'm0[MSUN]':float(datasnap[8]), 'm1[MSUN]':float(datasnap[9]), 'sma[AU]':float(datasnap[12]), 'ecc':float(datasnap[13])}
                        break
                else:
                    if int(datasnap[0])==theid:
                        output_snap[str(t_current)]={'id0':int(datasnap[0]), 'id1':-100, 'm0[MSUN]':float(datasnap[1]), 'm1[MSUN]':-100, 'sma[AU]':-100, 'ecc':-100}

        snap_num = convert_str_to_four_digits(snap_num)
        snap = filestr+'.snap%s.dat.gz' %(snap_num)
        print(snap, t_current)

    return output_snap



##Find the full history of a NS or MSP
def get_full_history_NS_id(NSid, Compid, modelpath):
    psrfile = modelpath+'initial.morepulsars.dat'
    binintfile = modelpath+'initial.binint.log'
    collfile = modelpath+'initial.collision.log'
    semergefile = modelpath+'initial.semergedisrupt.log'
    
    t_conv = dyn.conv('t', modelpath+'initial.conv.sh')

    s=modelpath.split('/')
    n_star=s[-2]
    z=s[-3][1:]
    rg=s[-4][2:]
    rv=s[-5][2:]

    #datasnapkey = np.genfromtxt(modelpath+'snap_keys.txt')
    #snapnos = datasnapkey[:,0]; snaptimes = datasnapkey[:,1]

    ##Formation times of the NS and the MSP
    check_nsform = 0; check_mspform = 0
    data_nsform = {}; data_mspform = {}
    with open(psrfile, 'r') as fpsr:
        next(fpsr)
        for line in fpsr:
            datapsr = line.split()
            check_nsform, check_mspform, NS_form_info, MSP_form_info = output_formation_NS_MSP(datapsr, check_nsform, check_mspform, NSid)
            if check_nsform == 1 and not data_nsform:
                data_nsform = NS_form_info
            if check_mspform == 1 and not data_mspform:
                data_mspform = MSP_form_info
            if check_nsform ==1 and check_mspform ==1:
                break
    print('formation time check done')
    fform = open('/projects/b1095/syr904/projects/PULSAR_Catalog/newruns/history_912gyr_maingrid_nondissolved/Formation_rv'+rv+'_rg'+rg+'_z'+z+'_n'+n_star+'_id'+str(NSid)+'.txt', 'w+')
    print(data_nsform, data_mspform, file = fform)
    fform.close()


    ##Is it a primordial binary?
    priinfo=output_pribinary_info_oldcmc(modelpath,NSid,Compid)
    print('primordial binary check doen')
    fpri = open('/projects/b1095/syr904/projects/PULSAR_Catalog/newruns/history_912gyr_maingrid_nondissolved/Pribinary_rv'+rv+'_rg'+rg+'_z'+z+'_n'+n_star+'_id'+str(NSid)+'.txt', 'w+')
    print(priinfo, file = fpri)
    fpri.close()


    ##Collision history
    colldata=scripts1.readcollfile(collfile)
    coll_info = output_coll_NS(colldata, NSid)
    print('collision check done')
    fcoll = open('/projects/b1095/syr904/projects/PULSAR_Catalog/newruns/history_912gyr_maingrid_nondissolved/Coll_rv'+rv+'_rg'+rg+'_z'+z+'_n'+n_star+'_id'+str(NSid)+'.txt', 'w+')
    print(coll_info, file = fcoll)
    fcoll.close()
    
    ##Merger history
    semerdata = scripts2.readmergefile(semergefile)
    semer_info = output_semer_NS(semerdata, NSid)
    print('merger check done')
    fsemer = open('/projects/b1095/syr904/projects/PULSAR_Catalog/newruns/history_912gyr_maingrid_nondissolved/Semerge_rv'+rv+'_rg'+rg+'_z'+z+'_n'+n_star+'_id'+str(NSid)+'.txt', 'w+')
    print(semer_info, file = fsemer)
    fsemer.close()


    ##Tidal capture and giant collision binary formation history
    #tc_gc_info = output_tc_gc_NS(modelpath, NSid)
    #print('tidal capture and giant collision check done')


    ##Binary interaction history
    binint_info = output_binint_NS(modelpath, NSid)
    print('binint check done')
    fbinint = open('/projects/b1095/syr904/projects/PULSAR_Catalog/newruns/history_912gyr_maingrid_nondissolved/Binint_rv'+rv+'_rg'+rg+'_z'+z+'_n'+n_star+'_id'+str(NSid)+'.txt', 'w+')
    print(binint_info, file = fbinint)
    fbinint.close()


    ##Snapshot history
    snap_info = check_snapshots_oldcmc(modelpath, NSid, float(NS_form_info['t_form']))
    fsnaps = open('/projects/b1095/syr904/projects/PULSAR_Catalog/newruns/history_912gyr_maingrid_nondissolved/Snap_rv'+rv+'_rg'+rg+'_z'+z+'_n'+n_star+'_id'+str(NSid)+'.txt', 'w+')
    print(snap_info, file = fsnaps)
    fsnaps.close()



##Find the full collision history of a NS
#output_coll = {'t_coll':colltime, 'type_coll':colltype, 'num_coll':collnum, 'mid':int(linecoll[3]), 'mm':collmm, 'mk':collmktype,'ids':collids, 'masses':masses, 'startypes':startypes}
#def get_full_coll_NS_id(NSid, modelpath):
#    t_conv = dyn.conv('t', modelpath+'initial.conv.sh')
#    collfile = modelpath+'initial.collision.log'
#    colldata=scripts1.readcollfile(collfile)
#    coll_info = output_coll_NS(colldata, NSid)
#    if not coll_info:
#        return 
#    else:
#        for uu in range(len(coll_info['num_coll'])):
#            if 13 not in coll_info['startypes']:
#                
#            if coll_info['startypes'][uu]==13:




                
                            


    
