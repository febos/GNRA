
import os


def AtomAtomDistance(coords1, coords2):
    return pow(sum((c1 - c2)**2
                   for c1, c2 in zip(coords1, coords2)), 0.5)

def DSSRidstr(atom):
    '''1.A.G.22.B - DSSR-like id string
       [model.chain.res.resnum.inscode]'''
    return '.'.join([str(x) for x in [atom["pdbx_PDB_model_num"],
                                      atom["auth_asym_id"],
                                      atom["auth_comp_id"],
                                      atom["auth_seq_id"],
                                      atom["pdbx_PDB_ins_code"]]])


def ParseAtomCIF(line, title):

    linesplit = line.strip().split()
    atom = {title[i]:linesplit[i] for i in range(len(title))}

    for frag in ("atom", "comp", "asym", "seq"):

        auth  =  "auth_{}_id".format(frag)
        label = "label_{}_id".format(frag)
        
        if auth not in atom and label in atom:
            atom[auth] = atom[label]
        elif label not in atom and auth in atom:
            atom[label] = atom[auth]

    # Convert all integers from strings
    for int_token in ("id", "auth_seq_id"):
        atom[int_token] = int(atom[int_token]) if int_token in atom else float("nan")

    # By default, the model is set to 1
    atom["pdbx_PDB_model_num"] = int(atom["pdbx_PDB_model_num"])\
                                 if "pdbx_PDB_model_num" in atom and\
                                 atom["pdbx_PDB_model_num"].isdigit() else 1

    # Convert all floats from strings
    for float_token in ("Cartn_x", "Cartn_y", "Cartn_z","occupancy","B_iso_or_equiv"):
        atom[float_token] = float(atom[float_token]) if float_token in atom else float("nan")
    
    if      "auth_atom_id"  not in atom: atom["auth_atom_id"]       = ''
    if     "label_atom_id"  not in atom: atom["label_atom_id"]      = ''
    if      "label_alt_id"  not in atom: atom["label_alt_id"]       = ''
    if      "auth_comp_id"  not in atom: atom["auth_comp_id"]       = ''
    if     "label_comp_id"  not in atom: atom["label_comp_id"]      = ''
    if      "auth_asym_id"  not in atom: atom["auth_asym_id"]       = ''
    if "pdbx_PDB_ins_code"  not in atom: atom["pdbx_PDB_ins_code"]  = ''
    if "pdbx_formal_charge" not in atom: atom["pdbx_formal_charge"] = ''

    if atom["pdbx_PDB_ins_code"] == '?':
        atom["pdbx_PDB_ins_code"] = ''

    if atom["pdbx_formal_charge"] == '?':
        atom["pdbx_formal_charge"] = ''

    if atom["label_alt_id"] == '.':
        atom["label_alt_id"] = ''

    # Strip double quotes for cases like this: "O3'"
    atom["auth_atom_id"] = atom["auth_atom_id"].strip('"') 
    atom["label_atom_id"] = atom["label_atom_id"].strip('"')

    # Convert SPDBV asterisks into prime symbols
    atom["auth_atom_id"]  = atom["auth_atom_id"].replace('*',"'")
    atom["label_atom_id"] = atom["label_atom_id"].replace('*',"'")
    # Convert SPDBV phosphate oxygen format O1P -> OP1, O2P -> OP2
    atom["auth_atom_id"]  = atom["auth_atom_id"].replace('O1P','OP1').replace('O2P','OP2')
    atom["label_atom_id"] = atom["label_atom_id"].replace('O1P','OP1').replace('O2P','OP2')

    return atom


def ParseCIF(filepath):

    residuelst, residuedct = [], {}

    with open(filepath) as file:
        for line in file:
            if line.startswith("data_"):
                model = 0
                if "sym" in line:
                    model = 1
                title  = []
                titlei = None
            if line.startswith("_atom_site."):
                title.append(line.strip().split('.')[-1])

            elif line.startswith('ATOM') or line.startswith('HETATM'):

                atom = ParseAtomCIF(line, title)

                atom["pdbx_PDB_model_num"] += model

                if not titlei:

                    title2 = [k for k in atom if k not in title]
                    finaltitle  = title + title2
                    titlei = {x:i for i,x in enumerate(finaltitle)}

                atom['DSSR'] = DSSRidstr(atom)

                # new residue found
                if atom['DSSR'] not in residuedct:
                    residuelst.append(atom['DSSR'])
                    residuedct[atom['DSSR']] = {'ATOMS': []}

                residue = residuedct[atom['DSSR']]
                # new atom found (otherwise we have a new altloc for an already found atom)
                if atom['auth_atom_id'] not in residue:
                    residue['ATOMS'].append(atom['auth_atom_id'])
                # consider only the last altloc version of the atom and ignore all the previous ones
                residue[atom['auth_atom_id']] = [atom[t] for t in ["Cartn_x", "Cartn_y", "Cartn_z"]]

    Entry = {'IDLIST': residuelst,
             'IDDICT': residuedct}

    return Entry


def CrystalContacts(pdb,hits, thresh = 4.0, farthresh = 20.0):

    idstrs = {'-':'0'}

    for hit in hits:
        for idstr in hit[0]:
            idstrs[idstr] = '0'

    entry = ParseCIF("PDB1cc/{}.cif".format(pdb))

    for idstr in entry['IDLIST']:
        if idstr in idstrs:
            found = False
            for idstr2 in entry['IDLIST']:
                toofar = False
                if idstr2.startswith('1.'):
                    continue
                for atom1 in entry['IDDICT'][idstr]['ATOMS']:
                    for atom2 in entry['IDDICT'][idstr2]['ATOMS']:
                        D = AtomAtomDistance(entry['IDDICT'][idstr][atom1],
                                             entry['IDDICT'][idstr2][atom2])
                        if D <= thresh:
                            found = True
                            break
                        elif D > farthresh:
                            toofar = True
                            break
                    if found or toofar:
                        break
                if found:
                    break
            if found:
                idstrs[idstr] = '1'

    return [''.join([idstrs[idstr] for idstr in hit[0]]) for hit in hits]



hits56 = {}
hits34 = {}

refmask = ['1.1A.C.2374.',
           '1.1A.G.2375.',
           '1.1A.A.2376.',
           '1.1A.A.2377.',
           '1.1A.A.2378.',
           '1.1A.G.2379.']

for file,pool in (("CGAAAG.artem", hits56),
                  ("GAAA.artem", hits34),):

    with open(file) as inp:
        for line in inp:
            if line.startswith("ID"):
                continue
            ls = line.strip().split()
            if not ls:
                continue
            _, size, rmsd, rmsdsize, resrmsd, prim, scnd, ref, qry = ls
            pdb = os.path.basename(ref).replace('.cif','')
            if pdb not in hits56:
                hits56[pdb] = []
            if pdb not in hits34:
                hits34[pdb] = []

            d = {}
            for pair in scnd.split(','):
                v,w = pair.split('=')
                d[w] = v

            hit = [d[idstr] if idstr in d else '-' for idstr in refmask]
            mask = ''.join(['0' if x=='-' else '1' for x in hit])
            bases = [x.split('.')[2] if x!='-' else '-' for x in hit]
            pool[pdb].append([hit,mask,bases,size,rmsd,rmsdsize,resrmsd])

hits = {}

with open("processed_hits.tsv",'w') as outp:

    for pdb in sorted(set(hits56.keys()) | set(hits34.keys())):
        print(pdb)
        hits[pdb] = []

        for hit in hits56[pdb][::-1]:
            hits[pdb].append(hit)

        #add 3/4-hits only if they don't match any of the 5/6-hits
        for hit in hits34[pdb][::-1]:
            for hit2 in hits56[pdb]:
                if all(x==y or x=='-' for x,y in zip(hit[0],hit2[0])):
                    break
            else:
                hits[pdb].append(hit)

        crystalcontacts = CrystalContacts(pdb,hits[pdb])
        
        for lst, crcts in zip(hits[pdb],crystalcontacts):
            outp.write('\t'.join([str(x) for x in [pdb,*lst[3:],"M"+lst[1],*lst[0],*lst[2],"M"+crcts]])+'\n')
