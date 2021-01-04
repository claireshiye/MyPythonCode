import numpy as np


##Add a column at the end of file
def add_column(filepath, savepath, new_column):
    dataset = np.genfromtxt(filepath)
    headers = '1.Model 2.Time(Myr) 3.Status 4.r(pc) 5.B(G) 6.P(sec) 7.dmdt0(Msun/yr) 8.dmdt1(Msun/yr) 9.rolrad0 10.rolrad1 11.m0(Msun) 12.m1(Msun) 13.ID0 14.ID1 15.k0 16.k1 17.a(AU) 18.ecc 19.Formation 20.TCflag'
    nc = np.reshape(new_column, (dataset.shape[0],1))   # adjust dimension of the new array
    result = np.append(dataset, nc, 1)         # append as last column
    np.savetxt(savepath, result, delimiter=" ", fmt="%s", header=headers, comments='#')
    
##Read ASCII tables from websites
def read_ascii_table(table_name_list, column_needed):
    all_cols = []
    tot_n = 0

    for ii in range(len(table_name_list)):
        col_nums = column_needed[ii]

        for i in range(len(col_nums)):
            all_cols.append([])

        with open(table_name_list[ii], 'r') as f:
            lines = f.readlines()
            #print(lines)
            for x in range(6, len(lines)):
                line = lines[x].strip('\n')
                data = line.split()
                data[1] = data[0]+data[1]
                data.pop(0)

                for y in range(len(col_nums)):
                    index_num = int(y+tot_n)
                    #print(all_cols, index_num)
                    all_cols[index_num].append(data[int(col_nums[y])])

        tot_n += len(col_nums)
    
    print(all_cols)
    print(tot_n)

    fout = open('/Users/shiye/Documents/Research/pulsar/field_redbacks.dat', 'a+')
    fout.write('#1.ID 2.P(ms) 3.Pdot_obs(10^-20) 4.Pb(days) 5.Mns(Msun) 6.Mns_minus(Msun) 7.Mns_plus(Msun) 8.Mns>=? 9.Mc(Msun) 10.Mc_minus(Msun) 11.Mc_plus(Msun) 12.Mc>=?\n')
    for jj in range(len(all_cols[0])):
        for j in range(tot_n):
            if j>0:
                num1 = all_cols[j][jj].split('(')[0]
                #print(num1)

                if '>or=' in num1:
                    num2 = num1.split('=')[1]
                    num34 = 'yes'
                else:
                    num2 =  num1
                    num34 = 'no'
    
                if '{' in num2:
                    num3 = num2.split('{')
                    #print(num3)
                    num31 = num3[1].split('}')[0]; num32 = num3[2].split('}')[0]; num33 = num3[3].split('}')[0]
                else:
                    num31 = num2; num32 = -100; num33 = -100
                        
                if j<=3:
                    if num31!='cdots':
                        fout.write('%g '%(float(num31)))
                    else:
                        fout.write('%s '%(num31))
                elif j>3:
                    if num31!='cdots':
                        print(num31, num32, num33, num34)
                        fout.write('%g %g %g %s '%(float(num31), float(num32), float(num33), num34))
                    else:
                        fout.write('%s %g %g %s '%(num31, num32, num33, num34))

            else:
                num1 = all_cols[j][jj]
                fout.write('%s '%(num1))

        fout.write('\n')

    fout.close()
                

              
         
