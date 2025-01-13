
## PDB search

Prepare ARTEM and urslib2

 - "git clone https://github.com/david-bogdan-r/ARTEM"
 
 - "git clone https://github.com/febos/urslib2"
 
 - "mv urslib2 urslib2_repo"
 
 - "ln -s urslib2_repo/urslib2 urslib2"
 
 - set the path to DSSR in urslib2/config.py [TODO]()

Prepare a search database

 - "pip install requests"
 
 - "python3 pdb_download.py" 
 
======================TODO==============================
 

 - [process the pdb entries == 1st models + crystal contacts]()
 
Run ARTEM

 - "python3 ARTEM/artem.py r=./PDB/*.cif q=refmotifs/1ffk_0_kt7.cif rres=#1 qres=":_77_82 :_92_94 :_97_100" sizemin=12 rmsdmax=2.0 rnosub=1 qrst=":_79_81 :_93_94 :_97_98" silent=True > 1ffk_hits.tsv"
 
 - "python3 ARTEM/artem.py r=./PDB/*.cif q=refmotifs/3d2g_A_kj.cif rres=#1 qres=":_8_10 :_36_38 :_42 :_44_46 :_69_72" sizemin=12 rmsdmax=2.0 rnosub=1 qrst=":_9 :_37_38 :_42 :_44 :_71_72" silent=True > 3d2g_hits.tsv"

Merge identical hits between the motifs

 - "python3 unique_hits.py"
 
Prepare urslib2
 


Annotate hits

 - see post-processing.ipynb
