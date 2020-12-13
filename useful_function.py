import numpy as np


##Add a column at the end of file
def add_column(filepath, savepath, new_column):
    dataset = np.genfromtxt(filepath)
    headers = '1.Model 2.Time(Myr) 3.Status 4.r(pc) 5.B(G) 6.P(sec) 7.dmdt0(Msun/yr) 8.dmdt1(Msun/yr) 9.rolrad0 10.rolrad1 11.m0(Msun) 12.m1(Msun) 13.ID0 14.ID1 15.k0 16.k1 17.a(AU) 18.ecc 19.Formation 20.TCflag'
    nc = np.reshape(new_column, (dataset.shape[0],1))   # adjust dimension of the new array
    result = np.append(dataset, nc, 1)         # append as last column
    np.savetxt(savepath, result, delimiter=" ", fmt="%s", header=headers, comments='#')
    

