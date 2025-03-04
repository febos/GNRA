
# Steps to reproduce the results (Linux-based)

##Prepare ARTEM and urslib2

 - "git clone https://github.com/david-bogdan-r/ARTEM"
 
 - "git clone https://github.com/febos/urslib2"
 
 - "mv urslib2 urslib2_repo"
 
 - "ln -s urslib2_repo/urslib2 urslib2"
 
 - set the path to DSSR in urslib2/config.py 

##Prepare the search database and the reference instance 

 - "pip install requests"
 
 - "python3 pdb_download.py" 
 
 - install [UCSF ChimeraX](https://www.cgl.ucsf.edu/chimerax/)
 
 - "python3 pdb_process.py" (keep only the 1st models + add interacting symmetry mates for X-ray structures)
 
 - "python3 choose_GNRA_reference.py"
 
##Obtain, annotate, and filter the GNRA matches

 - "python3 ARTEM/artem.py r=PDB1cc q=8VTW_CGAAAG.cif rres=# rseed=#1 rformat=cif sizemin=5 rmsdsizemax=0.25 rnosub=1 silent=1 > CGAAAG.artem"
 
 - "python3 ARTEM/artem.py r=PDB1cc q=8VTW_GAAA.cif rres=# rseed=#1 rformat=cif sizemin=3 rmsdsizemax=0.25 qrst="/1A:_2375" rnosub=1 silent=1 > GAAA.artem"
 
 - "python3 postprocess_artem.py"
 
 - "python3 run_DSSR.py"
 
 - "python3 annotate_hits.py"
 
 - "python3 filter_hits.py"


