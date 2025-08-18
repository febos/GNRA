
import matplotlib.pyplot as plt

highres_threshold = 4.0

resols = []

bad = set()
nmr = set()
good = set()

with open("../pdb_resol.tsv") as inp:
    for line in inp:
        ls = line.strip().split()
        if ls[1] != '\\N':
            resol = float(ls[1])
            if resol <= highres_threshold:
                good.add(ls[0])
            else:
                bad.add(ls[0])
            resols.append(resol)
        else:
            nmr.add(ls[0])




GOOD = 0
BAD  = 0
NMR  = 0

with open("annotated_hits_highres.tsv",'w') as outp:
    with open("../annotated_hits.tsv") as inp:
        for line in inp:
            pdb = line[:4]
            if pdb in good:
                GOOD += 1
                outp.write(line)
            elif pdb in bad:
                BAD += 1
            elif pdb in nmr:
                NMR += 1
        
                
    
print(len(good),len(bad),len(nmr))
print(GOOD, BAD, NMR)
plt.hist(resols,bins=100)
plt.xlim([0,10])
plt.show()
