
import glob, os
from urslib2 import DSSR

pdbs = set()
with open("processed_hits.tsv") as inp:
    for line in inp:
        if line.strip():
            pdb = line.strip().split()[0]
            pdbs.add(pdb)

modelfolder = "PDB1"
dssrfolder = "PDB1dssr"

path_to_models = os.path.join(modelfolder,"*.cif1")
os.makedirs(dssrfolder,   exist_ok = True)

models = sorted(glob.glob(path_to_models))

for model in models:

    pdb = os.path.basename(model).split('.')[0]
    if pdb not in pdbs:
        continue

    outmodel = model.replace(modelfolder,dssrfolder).replace('.cif1','.out1')
    
    if os.path.exists(outmodel) and os.path.getsize(outmodel) > 0:
        print(model,"(skipped)")
        continue
    print(model)
    DSSR.run(model, dssrfolder)
    
