# ensem_water_rmsf

##Summary:

This script computes the root mean square fluctuation (RMSF) for each user-specified HETATM (solvent) in an NMR-style ensemble model and writes to either a flat .csv file, or a .pml file in order to add those atoms and their positions as "pseudoatoms" to an existing pymol file for the purposes of displaying them by cartoon putty with a color scale based on B-factor. 

A max RMSF can be provided for .pml output mode so the user can make comparable visualizations across different ensembles.

##Example usage:

phenix.python ensem_water_rmsf.py ensem.pdb O csv > outfile.csv
phenix.python ensem_water_rmsf.py ensem.pdb O pml max_rmsf > outfile.pml
##Conditions to run:

###Software requirements:

phenix==1.9-1692
###Additional requirements:

A multi-MODEL PDB file with HETATMs
