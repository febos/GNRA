

from urslib2 import RSS, Rfam


rfam = Rfam.GetRfamInfo()

nrclass = {}

with open("nrlist_3.370_all.csv") as inp:
    for line in inp:
        if line.strip():
            ls = line.strip().strip('"').split('","')
            for efu in ls[2].split(','):
                for rna in efu.split('+'):
                    pdb,mod,ch = rna.split('|')
                    nrclass[(pdb,ch)] = ls[0]


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

        uturns = [0 for i in range(len(model.u_turns))]

        for hit in hits[pdb]:

            size, rmsd, rmsdsize, resrmsd, matchmask, rs, bs, crystalmask = *hit[:5],hit[5:11],hit[11:17],hit[17] 

            # consecutive pairwise connectivity + same/local/long-range 
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

            # stems/loops
            elems = []
            for n1 in rs:
                elem = 'NA'
                if n1 != '-':
                    elem = model.NuclSS(n1[2:])
                elems.append(elem)

            # sequence+length if hairpin
            loopseq,looplen = '-',0
            if set(elems[1:-1]) == {'HC',} and set(rels[1:-1]) == {'SM',}:
                c,m,i = model.dssrnucls[rs[1][2:]]
                t = model.chains[c][m][i]['THREAD']
                loopseq,looplen = model.threads[t-1]['SEQ'], model.threads[t-1]['LEN']

            ## conformations + puckers
            confs = []
            for n1 in rs:
                conf = ['NA','NA']
                if n1 != '-' and n1[2:] in model.summary:
                    desc = model.summary[n1[2:]]
                    if "~C2'-endo," in desc:
                        conf[1] = "~C2'-endo"
                    elif "~C3'-endo," in desc:
                        conf[1] = "~C3'-endo"
                    if "anti," in desc:
                        conf[0] = "anti"
                    elif "syn," in desc:
                        conf[0] = "syn"
                confs.append(','.join(conf))

            # GA-like bp
            ga = '-'
            resG = rs[1][2:]
            resA = rs[-2][2:]
            if resG and resA:
                for bp in model.bpairs:
                    pair = bp['NUCL1'][0]+bp['NUCL2'][0]
                    if pair == resG+resA:
                        ga = bp['BOND']+',' +bp['CLASS'][1]
                        break
                    elif pair == resA+resG:
                        ga = '-'.join(bp['BOND'].split('-')[::-1]) + ',' + bp['CLASS'][1][0] +\
                             bp['CLASS'][1][2]+bp['CLASS'][1][1]
                        break

            # additional bps/stacks
            inners = {n1[2:]:k for k,n1 in enumerate(rs)}
            inters = [['-','-'] for x in rs]
            for bp in model.bpairs:
                if bp['NUCL1'][0] in inners and bp['NUCL2'][0] not in inners:
                    inters[inners[bp['NUCL1'][0]]][0] = "basepair"
                elif bp['NUCL2'][0] in inners and bp['NUCL1'][0] not in inners:
                    inters[inners[bp['NUCL2'][0]]][0] = "basepair"

            for stack in model.stacks:
                for n1,n2 in zip(stack[1][:-1],stack[1][1:]):
                    if n1 in inners and n2 not in inners:
                        inters[inners[n1]][1] = "stack"
                    elif n2 in inners and n1 not in inners:
                        inters[inners[n2]][1] = "stack"

            inters = [','.join(x) for x in inters]

            #U-turns
            uts = 0
            for i,ut in enumerate(model.u_turns):
                if ut['NUCL1'] == rs[1][2:] or ut['NUCL2'] == rs[-2][2:]:
                    uturns[i] += 1
                    uts = 1
                    break

            #Rfam + RNA3DHub
            fams = set()
            classes = set()
            for n1 in rs:
                if n1 != '-':
                    c1, m1, i1 = model.dssrnucls[n1[2:]]
                    if pdb.lower() in rfam and c1 in rfam[pdb.lower()]:
                        fams.add(rfam[pdb.lower()][c1][0])
                    if (pdb,c1) in nrclass:
                        classes.add(nrclass[(pdb,c1)])
            fams    = ','.join(sorted(fams)) if fams else '-'
            classes = ','.join(sorted(classes)) if classes else '-'
            #print(classes,fams)

            outp.write('\t'.join([str(x) for x in [pdb, size, rmsd, rmsdsize, resrmsd, matchmask, *rs, *bs, crystalmask,
                                                   *connects, *rels, *elems, loopseq, looplen,
                                                   *confs, ga, *inters, uts, fams, classes]])+'\n')
            
            
            

        #for i in range(len(uturns)):
        #    if uturns[i] != 1:
        #        print(uturns[i],model.u_turns[i])
                    
                    
                
#rfam = Rfam.GetRfamInfo()
#rfam['1ffk']['0']

#Rfam-families, NRlist-classes, ...
