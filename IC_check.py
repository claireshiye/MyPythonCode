import numpy as np
import re


path_old = '/projects/b1091/CMC_Grid_March2019/rundir/rv0.5/rg8/z0.0002/8e5/'
path_new = '/projects/b1095/syr904/cmc/CMC-COSMIC/wdmerger_update/CMC-COSMIC_modified/CMC/runs/n8-rv0.5-rg8-z0.0002_ICold_wdmass_nstde_test/'

ini_old = path_old+'initial.cmc.parsed'
ini_new = path_new+'KingProfile.ini'

allflags_new = []
allflagvalues_new = []
with open(ini_new, 'r') as fini_new:
    for line in fini_new:
        if line[0]==';' or line[0]=='' or line[0]=='[':
            continue
        else:
            data = line.split()
            #print(data)
            if not data:
                #print(data)
                continue
            allflags_new.append(data[0])
            allflagvalues_new.append(data[2])
#print(allflags_new, allflagvalues_new)


allflags_old = []
allflagvalues_old = []
with open(ini_old, 'r') as fini_old:
    for line in fini_old:
        if line[0]=='#':
            continue
        else:
            data = line.split()
            allflags_old.append(data[0])
            allflagvalues_old.append(data[1])
#print(allflags_old, allflagvalues_old)


diff_flags = []
diff_values_new = []
diff_values_old = []
for xx in range(len(allflags_new)):
    check=0
    for yy in range(len(allflags_old)):
        if allflags_new[xx]==allflags_old[yy]:
            check=1
            if ',' in allflagvalues_new[xx] or re.search('[a-zA-Z]', allflagvalues_new[xx]):
                diff_flags.append(allflags_new[xx])
                diff_values_new.append(allflagvalues_new[xx])
                diff_values_old.append(allflagvalues_old[yy])
            elif float(allflagvalues_new[xx])!=float(allflagvalues_old[yy]):
                diff_flags.append(allflags_new[xx])
                diff_values_new.append(allflagvalues_new[xx])
                diff_values_old.append(allflagvalues_old[yy])

            break
        else:
            if allflags_old[yy][:3]=='BSE':
                convert_flag = allflags_old[yy][4:].lower()
                if convert_flag == allflags_new[xx]:
                    check=1
                    if ',' in allflagvalues_new[xx] or re.search('[a-zA-Z]', allflagvalues_new[xx]):
                        diff_flags.append(allflags_new[xx])
                        diff_values_new.append(allflagvalues_new[xx])
                        diff_values_old.append(allflagvalues_old[yy])
                    elif float(allflagvalues_new[xx])!=float(allflagvalues_old[yy]):
                        diff_flags.append(allflags_new[xx])
                        diff_values_new.append(allflagvalues_new[xx])
                        diff_values_old.append(allflagvalues_old[yy])

    if check==0:
        diff_flags.append(allflags_new[xx])
        diff_values_new.append(allflagvalues_new[xx])
        diff_values_old.append(-100)

print(diff_flags)
print(diff_values_new)
print(diff_values_old)


