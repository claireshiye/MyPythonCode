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
    for xx in range(len(colldata)):
        linecoll = colldata[xx].split()
        if int(linecoll[3])==theid:
            colltime=float(linecoll[0])
            colltype = linecoll[1]; collnum=int(linecoll[2])
            collmm=float(linecoll[4])
            if int(linecoll[2])==2:
                collids = [int(linecoll[5]), int(linecoll[7])]
                collmktype=int(linecoll[10])
                startypes = [int(linecoll[11]), int(linecoll[12])]
                masses = [float(linecoll[6]), float(linecoll[8])]

            elif int(linecoll[2])==3:
                collids=[int(linecoll[5]), int(linecoll[7]), int(linecoll[9])]
                collmktype=int(linecoll[12])
                startypes = [int(linecoll[13]), int(linecoll[14]), int(linecoll[15])]
                masses = [float(linecoll[6]), float(linecoll[8]), float(linecoll[10])]

            elif int(linecoll[2])==4:
                collids=[int(linecoll[5]), int(linecoll[7]), int(linecoll[9]), int(linecoll[11])]
                collmktype=int(linecoll[14])
                startypes = [int(linecoll[15]), int(linecoll[16]), int(linecoll[17]), int(linecoll[18])]
                masses = [float(linecoll[6]), float(linecoll[8]), float(linecoll[10]), float(linecoll[12])]
    
            break

    output_coll = {'t_coll':colltime, 'type_coll':colltype, 'num_coll':collnum,
                   'mm':collmm, 'mk':collmktype,
                   'ids':collids, 'masses':masses, 'startypes':startypes}

    return output_coll


##Find all binary merger/disruption of the NS
def output_semer_NS(datasemer, theid):
    output_semer={}
    for xx in range(len(semerdata)):
        linesemer=semerdata[xx].split()
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
            output_semer[str(t_semer)]={'type':type_semer, 'rm':mr, 'rk':kr,
                                        'ids':semerids, 'masses':semermass, 'startypes':semerktype,
                                        'r_pos':rsemer}

    return output_semer



##Find the full history of a NS or MSP
def get_full_history_NS_id(NSid, Compid, modelpath, snflag, tcflag):
    psrfile = modelpath+'initial.morepulsars.dat'
    binintfile = modelpath+'initial.binint.log'
    collfile = modelpath+'initial.collision.log'
    semergefile = modelpath+'initial.semergedisrupt.log'
    snaps=dyn.get_snapshots(modelpath+'initial')
    snap0 = snaps[0]

    t_conv = ns.conv('t', modelpath+'initial.conv.sh')

    dict_ns = {time_conv:t_conv}

    check_nsform = 0; check_mspform = 0
    with.open(psrfile, 'r') as fpsr:
        next(fpsr)
            for line in fpsr:
                datapsr = line.split()
                check_nsform, check_mspform, NS_form_info, MSP_form_info = output_formation_NS_MSP(datapsr, check_nsform, check_mspform)
                if check_nsform == 1 and not NS_form_info:
                    data_nsform = NS_form_info
                if check_mspform == 1 and not MSP_form_info:
                    data_mspform = MSP_form_info
                if check_nsform ==1 and check_msp ==1:
                    break


    colldata=scripts1.readcollfile(collfile)
    coll_info = output_coll_NS(colldata, NSid)


    semerdata = scripts2.readmergefile(semerge)
    semer_info = output_semer_NS(semerdata, NSid)

    
    property_init, property_finl, property_des = find_tc_properties_final(modelpath)
    for xx in range(len(prop_init['id0'])):
        if int(prop_init['id0'][xx])==NSid:

             period = uc.au_to_period(prop_finl['sma'][xx], prop_finl['m0'][xx], prop_finl['m0'][xx])
                
            if prop_finl['k1'][xx]==13 and prop_init['type'][xx] == 'SS_COLL_TC_P':
                period = uc.au_to_period(prop_finl['sma'][xx], prop_finl['m0'][xx], prop_finl['m0'][xx])
                ftc.write('%d %f %d %d %d %d %f %f %f %f %f %f %d %d %f %f %f %f %f\n'%(ii, prop_init['time'][xx]*t_conv, prop_init['id1'][xx], prop_init['id0'][xx], 
                      prop_init['k1'][xx], prop_init['k0'][xx], 
                      prop_init['m1'][xx], prop_init['m0'][xx],
                      prop_init['r1'][xx], prop_init['r0'][xx],
                      prop_init['rperi'][xx], prop_init['vinf'][xx],
                      prop_finl['k1'][xx], prop_finl['k0'][xx], 
                      prop_finl['m1'][xx], prop_finl['m0'][xx],
                      prop_finl['sma'][xx], prop_finl['ecc'][xx],
                      period))

            if prop_finl['k0'][xx]==13 and prop_init['type'][xx] == 'SS_COLL_Giant':
                period = uc.au_to_period(prop_finl['sma'][xx], prop_finl['m0'][xx], prop_finl['m0'][xx])
                fcoll.write('%d %f %d %d %d %d %f %f %f %f %f %f %d %d %f %f %f %f %f\n'%(ii, prop_init['time'][xx]*t_conv, prop_init['id0'][xx], prop_init['id1'][xx], 
                      prop_init['k0'][xx], prop_init['k1'][xx], 
                      prop_init['m0'][xx], prop_init['m1'][xx],
                      prop_init['r0'][xx], prop_init['r1'][xx],
                      prop_init['rperi'][xx], prop_init['vinf'][xx],
                      prop_finl['k0'][xx], prop_finl['k1'][xx], 
                      prop_finl['m0'][xx], prop_finl['m1'][xx],
                      prop_finl['sma'][xx], prop_finl['ecc'][xx],
                      period))
            if prop_finl['k1'][xx]==13 and prop_init['type'][xx] == 'SS_COLL_Giant':
                period = uc.au_to_period(prop_finl['sma'][xx], prop_finl['m0'][xx], prop_finl['m0'][xx])
                fcoll.write('%d %f %d %d %d %d %f %f %f %f %f %f %d %d %f %f %f %f %f\n'%(ii, prop_init['time'][xx]*t_conv, prop_init['id1'][xx], prop_init['id0'][xx], 
                      prop_init['k1'][xx], prop_init['k0'][xx], 
                      prop_init['m1'][xx], prop_init['m0'][xx],
                      prop_init['r1'][xx], prop_init['r0'][xx],
                      prop_init['rperi'][xx], prop_init['vinf'][xx],
                      prop_finl['k1'][xx], prop_finl['k0'][xx], 
                      prop_finl['m1'][xx], prop_finl['m0'][xx],
                      prop_finl['sma'][xx], prop_finl['ecc'][xx],
                      period))



                
                            


    
