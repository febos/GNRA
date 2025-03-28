

bag = {}

with open("annotated_hits.tsv") as inp:
    cnt = 0
    for line in inp:
        cnt += 1
        if line.strip():
            hit = line.strip().split()
            
            pdb, size, rmsd, rmsdsize, resrmsd, matchmask,\
            rs, bs, crystalmask, conns, rels, elems,\
            loopseq, loopsize, confs, gabp, flbp, inters,\
            ut, rfam, nrclass = *hit[:6],hit[6:12],hit[12:18],hit[18],hit[19:24],hit[24:29],\
                                hit[29:35],hit[35],hit[36],hit[37:43],hit[43],hit[44],hit[45:51],\
                                hit[51],hit[52],hit[53]
            rmsd = float(rmsd)
            desc = ':'.join([str(x) for x in [size, matchmask, *bs, crystalmask,
                                              *conns, *rels, *elems, loopseq,
                                              loopsize, confs, gabp, flbp, inters,
                                              ut, rfam, nrclass]])

            if desc not in bag or bag[desc][0] > rmsd:
                bag[desc] = [rmsd,cnt,line]   
            

with open("nr_hits.tsv","w") as outp:
    for cnt,line in sorted([(v[1],v[2]) for k,v in bag.items()]):
        outp.write(line)
    
bag = {}
counts = {}

with open("nr_hits.tsv") as inp:
    cnt = 0
    for line in inp:
        cnt += 1
        if line.strip():
            hit = line.strip().split()
            
            pdb, size, rmsd, rmsdsize, resrmsd, matchmask,\
            rs, bs, crystalmask, conns, rels, elems,\
            loopseq, loopsize, confs, gabp, flbp, inters,\
            ut, rfam, nrclass = *hit[:6],hit[6:12],hit[12:18],hit[18],hit[19:24],hit[24:29],\
                                hit[29:35],hit[35],hit[36],hit[37:43],hit[43],hit[44],hit[45:51],\
                                hit[51],hit[52],hit[53]
            rmsd = float(rmsd)
            
            desc = ':'.join([str(x) for x in [size, matchmask, *bs, crystalmask,
                                              *conns, *rels, *elems, loopseq,
                                              loopsize, confs, gabp, flbp, inters,
                                              ut]])

            if desc not in counts:
                counts[desc] = [set(), set()]
            counts[desc][0].add(rfam)
            counts[desc][1].add(nrclass)

            if desc not in bag or bag[desc][0] > rmsd:
                bag[desc] = [rmsd,cnt,line]   
            

with open("unique_hits.tsv","w") as outp:
    for cnt,line, desc in sorted([(v[1],v[2],k) for k,v in bag.items()]):
        outp.write(line.strip()+'\t'+str(len(counts[desc][0]))+'\t'+str(len(counts[desc][1]))+'\t'+\
                   ','.join(sorted(counts[desc][0]))+'\t'+','.join(sorted(counts[desc][1]))+'\n')


