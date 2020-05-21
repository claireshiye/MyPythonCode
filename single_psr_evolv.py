import numpy as np

twopi=2.*np.pi
yearsc=3.1557*10**7
bconst=-3000.   ##Myr
K=twopi**2*2.5*10**-49  ##Unit is year/Gauss^2
Bbot=5.*10**7


def single_psr_evolv(B0, P0, Time):  ##Time in Myr, P0 in sec, B0 in Gauss
    if Time<=1:
        t_step=[Time]
    elif Time>1:
        t_step=list(np.arange(0., Time, 1.))
        t_step.append(Time)
        del t_step[0]
    P=P0/yearsc
    #print t_step
    for i in range(len(t_step)):
        B=B0*np.exp((t_step[i])/bconst)+Bbot
        Pdot=K*B*B/P
        P=P+Pdot*1e6

    P=P*yearsc

    return B, P