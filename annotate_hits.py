

from urslib2 import RSS, Rfam


hits = {}
with open("processed_hits.tsv") as inp:
    for line in inp:
        if line.strip():
            linesplit = line.strip().split()
            if linesplit[0] not in hits:
                hits[linesplit[0]] = []
            hits[linesplit[0]].append(linesplit[1:])


with open("annotated_hits.tsv","w") as outp:
    for k,pdb in enumerate(sorted(hits.keys())):
        print(pdb, k+1)
        pdbmodel = "PDB1/{}.cif1".format(pdb)
        outmodel = "PDB1dssr/{}.out1".format(pdb)
        model = RSS.SecStruct(pdbmodel,outmodel)

        for hit in hits[pdb]:

            size, rmsd, rmsdsize, resrmsd, matchmask, rs, bs, crystalmask = *hit[:5],hit[5:11],hit[11:17],hit[17] 

            connects = []
            rels = []
            for n1,n2 in zip(rs[:-1],rs[1:]):
                connect = 0
                rel = 'NA'
                if n1 != '-' and n2 != '-':
                    c1, m1, i1 = model.dssrnucls[n1[2:]]
                    c2, m2, i2 = model.dssrnucls[n2[2:]]
                    if c1==c2 and m1==m2:
                        connect = i2-i1
                    rel = model.NuclRelation(n1[2:], n2[2:])
                connects.append(connect)
                rels.append(rel)

            elems = []
            for n1 in rs:
                elem = 'NA'
                if n1 != '-':
                    elem = model.NuclSS(n1[2:])
                elems.append(elem)

            loopseq,looplen = '-',0

            if set(elems[1:-1]) == {'HC',} and set(rels[1:-1]) == {'SM',}:
                c,m,i = model.dssrnucls[rs[1][2:]]
                t = model.chains[c][m][i]['THREAD']
                loopseq,looplen = model.threads[t-1]['SEQ'], model.threads[t-1]['LEN']
                
#rfam = Rfam.GetRfamInfo()
#rfam['1ffk']['0']

#anti/syn,sugar-pucker,
#additional interactions, U-turns,
#Rfam-families, NRlist-classes, ...
