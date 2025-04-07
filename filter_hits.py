

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
                bag[desc] = [rmsd,cnt,line, pdb]   
            

with open("nr_hits.tsv","w") as outp:
    for cnt,line,pdb in sorted([(v[1],v[2],v[3]) for k,v in bag.items()]):
        outp.write(line.replace(pdb,"pdb_0000"+pdb.lower()))
    
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

            topo = []
            for i in range(len(bs)-1):
                if bs[i]!='-':
                    topo.append('n{}'.format(i+1))
                    if bs[i+1]!='-':
                        if conns[i]=='1':
                            topo.append('-')
                        else:
                            topo.append('|')
                    else:
                        topo.append('_')
                else:
                    topo.append('___')
            if bs[-1] != '-':
                topo.append('n{}'.format(len(bs)))
            else:
                topo.append('__')
            topo = ''.join(topo)


            link = "https://rna.bgsu.edu/rna3dhub/display3D/unitid/{}"\
                   .format(','.join([pdb.upper()+'|'+x.replace('.','|')[:-1] for x in rs if x != '-']))
            
            desc = ':'.join([str(x) for x in [size, matchmask, *bs, crystalmask,
                                              *conns, *rels, *elems, loopseq,
                                              loopsize, confs, gabp, flbp, inters,
                                              ut]])

            if desc not in counts:
                counts[desc] = [set(), set()]
            counts[desc][0].add(rfam)
            counts[desc][1].add(nrclass)

            if desc not in bag or bag[desc][0] > rmsd:
                bag[desc] = [rmsd,cnt,line,topo,link,pdb]   
            

with open("unique_hits.tsv","w") as outp:
    for cnt,line, desc, topo,link,pdb in sorted([(v[1],v[2],k,v[3],v[4],v[5]) for k,v in bag.items()]):
        outp.write(line.strip()+'\t'+str(len(counts[desc][0]))+'\t'+str(len(counts[desc][1]))+'\t'+\
                   topo+'\t'+\
                   ','.join(sorted(counts[desc][0]))+'\t'+','.join(sorted(counts[desc][1]))+'\t'+link+'\n')

    
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

            topo = []
            for i in range(len(bs)-1):
                if bs[i]!='-':
                    topo.append('n{}'.format(i+1))
                    if bs[i+1]!='-':
                        if conns[i]=='1':
                            topo.append('-')
                        else:
                            topo.append('|')
                    else:
                        topo.append('_')
                else:
                    topo.append('___')
            if bs[-1] != '-':
                topo.append('n{}'.format(len(bs)))
            else:
                topo.append('__')
            topo = ''.join(topo)


            link = "https://rna.bgsu.edu/rna3dhub/display3D/unitid/{}"\
                   .format(','.join([pdb.upper()+'|'+x.replace('.','|')[:-1] for x in rs if x != '-']))
            
            desc = topo

            if desc not in counts:
                counts[desc] = [set(), set()]
            counts[desc][0].add(rfam)
            counts[desc][1].add(nrclass)

            if desc not in bag or bag[desc][0] > rmsd:
                bag[desc] = [rmsd,cnt,line,topo,link,pdb]   
            

with open("unique_topologies.tsv","w") as outp:
    for cnt,line, desc, topo,link,pdb in sorted([(v[1],v[2],k,v[3],v[4],v[5]) for k,v in bag.items()]):
        outp.write(line.strip()+'\t'+str(len(counts[desc][0]))+'\t'+str(len(counts[desc][1]))+'\t'+\
                   topo+'\t'+\
                   ','.join(sorted(counts[desc][0]))+'\t'+','.join(sorted(counts[desc][1]))+'\t'+link+'\n')
