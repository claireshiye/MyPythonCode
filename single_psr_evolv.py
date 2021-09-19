import numpy as np

twopi=2.*np.pi
yearsc=3.1557*10**7
bconst=-3000.   ##Myr
K=twopi**2*2.5*10**-49  ##Unit is year/Gauss^2
Bbot=5.*10**7


def single_psr_evolv_onestep(B0, P0, Time):  ##Time in Myr, P0 in sec, B0 in Gauss
    Pyear0=P0/yearsc
    #print t_step

    Pdot=K*B0*B0/Pyear0
    B=B0*np.exp(Time/bconst)+Bbot
    Pyear=Pyear0+Pdot*1e6

    P=Pyear*yearsc

    return B, P


def single_psr_evolv_itgstep(Bini, Pini, Time):  ##Time in Myr, P0 in sec, B0 in Gauss
    B = []; P = []
    for ii in range(len(Time)):
        last_step = Time[ii]%1.
        tot_step = int(Time[ii]-last_step) 
        t_step=list(np.ones(tot_step))
        t_step.append(last_step)

        Pyear0=Pini[ii]/yearsc
        B0 = Bini[ii]
        #print t_step

        for x in range(len(t_step)):
            Pdot=K*B0*B0/Pyear0
            Bt=B0*np.exp((t_step[x])/bconst)+Bbot
            Pyear=Pyear0+Pdot*1e6

            B0=Bt; Pyear0=Pyear
  
        P.append(Pyear*yearsc)
        B.append(Bt)

    return B, P