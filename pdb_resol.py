
import os, glob

from urslib2 import mmCIF


if __name__ == "__main__":

    with open("pdb_resol.tsv",'w') as outp:
        for file in sorted(glob.glob("PDB1/*.cif1")):
            pdb = os.path.basename(file)[:4]
            print(file)
            model = mmCIF.Model(file)
            outp.write('\t'.join([pdb,str(model.headers['RESOL']),model.headers['EXPDTA'][0]])+'\n')
            
