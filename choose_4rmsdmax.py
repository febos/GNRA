
import os

motif_csv = "HL_85603.2.csv"
pdb_folder = "PDB1"
pairwise_file = "GNRA_pairwise_RMSD4.csv"

motifs = []
rmsds  = {}
prev_rmsds = {}

with open("HL_85603.2.csv") as inp:
    for line in inp:

        nts = [x.split('|') for x in line.strip().strip('"').split('","')]

        pdb = nts[0][0]

        nts = ["#{}/{}:_{}".format(x[1],x[2],x[4]) if len(x) < 6 else
               "#{}/{}:_{}{}".format(x[1],x[2],x[4],x[-1])
               for x in nts][1:-1]

        motifs.append((pdb,' '.join(nts)))
        rmsds[(pdb,' '.join(nts))] = []

       
with open(pairwise_file,'w') as outp:

    for i in range(35,36):
        m1 = motifs[i]
        print(i+1, m1)
        for j in range(len(motifs)):
            print(j+1,end=' ')
            m2 = motifs[j]
            
            command = 'python ARTEM/artem.py r={} rres="{}" q={} qres="{}"'.format(os.path.join(pdb_folder,
                                                                                                m1[0]+'.cif1'),
                                                                                   m1[1],
                                                                                   os.path.join(pdb_folder,
                                                                                                m2[0]+'.cif1'),
                                                                                   m2[1],)
            command += " rformat=cif qformat=cif sizemin=4 > artem.tmp"

            if (m1[0],m1[1],m2[0],m2[1]) in prev_rmsds:
                rmsd = prev_rmsds[(m1[0],m1[1],m2[0],m2[1])]
            else:
                os.system(command)
                try:
                    
                    with open("artem.tmp") as inp:
                        rmsd = inp.readlines()[1].split()[2]
                except:
                    rmsd = "10"

            outp.write('\t'.join([*m1,*m2,rmsd])+'\n')
            print('\t'.join([*m1,*m2,rmsd]))
            rmsds[m1].append(float(rmsd))
            rmsds[m2].append(float(rmsd))
        print()

best_motif, best_median = '', 10

N = len(motifs)
for m in motifs:
    print(m)
    print(sorted(rmsds[m]))
    median = sorted(rmsds[m])[N//2]
    if median < best_median:
        best_motif = m
        best_median = median

print("CENTROID:")
print(best_motif)
print("Median RMSD:",round(best_median,3))
    
