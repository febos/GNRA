# GNRA
Exploring GNRA tetraloop-like motifs in nucleic acid 3D structures

## Reference

[J.M. Bujnicki, E.F. Baulin (2025) Exploring GNRA tetraloop-like motifs in nucleic acid 3D structures. bioRxiv. DOI: 10.1101/2025.07.03.663028](https://doi.org/10.1101/2025.07.03.663028)

## Check out [our other developments](https://github.com/febos/wiki)

## Content 

### Figures

Manuscript figures

### Variants

The representative matches of the twelve recurrent backbone topology GNRA-tetraloop motif variants in PDB format

### files

- **8VTW_CGAAAG.cif** - The six-residue reference in mmCIF format
- **8VTW_GAAA.cif** - The four-residue reference in mmCIF format
- **CGAAAG.artem** - Raw ARTEM output for the six-residue reference search
- **FigureRMSDthresh.py** - Python script to plot RMSD distributions
- **GAAA.artem** - Raw ARTEM output for the four-residue reference search
- **GNRA_pairwise_RMSD.csv** - Pairwise RMSD values among six-residue motifs
- **GNRA_pairwise_RMSD4.csv** - Pairwise RMSD values among four-residue motifs
- **HL_85603.2.csv** - The GNRA motif class of the [RNA 3D Motif Atlas](https://rna.bgsu.edu/rna3dhub/motifs)
- **README.md** - This README file
- **REPRODUCE.md** - The steps to reproduce the results
- **TableS1.xlsx** - Table S1. List of 23,283 non-redundant matches of the GNRA motif.
- **TableS2.xlsx** - Table S2. List of 13,729 unique instances of the GNRA tetraloop motif matches
- **TableS3.xlsx** - Table S3. List of 59 unique instances of the GNRA tetraloop motif strand topologies
- **annotate_hits.py** - Python script to annotate ARTEM matches with structural features
- **annotated_hits.tsv** - The list of annotated matches
- **choose_4rmsdmax.py** - Python script to select the RMSD thresholds
- **choose_GNRA_reference.py** - Python script to select the GNRA tetraloop motif reference
- **filter_hits.py** - Python script to filter redundant matches
- **nr_hits.tsv** - The list of non-redundant matches
- **nrlist_3.370_all.csv** - [The BGSU representative set of RNA structures](https://rna.bgsu.edu/rna3dhub/nrlist)
- **pdb_download.py** - Python script to download nucleic acid-containing PDB entries
- **pdb_process.py** - Python script to process the PDB entries
- **postprocess_artem.py** - Python script to parse ARTEM matches
- **processed_hits.tsv** - The list of parsed ARTEM matches
- **run_DSSR.py** - Python script to run DSSR annotations
- **stat.ipynb** - Jupyter Notebook with summary statistics of the dataset
- **unique_hits.tsv** - The list of structurally unique matches
- **unique_topologies.tsv** - The list of unique backbone topology variants of the GNRA tetraloop motif
