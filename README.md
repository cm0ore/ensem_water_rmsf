# ensem_water_rmsf

##Summary:

This script computes the root mean square fluctuation (RMSF) for each user-specified HETATM (solvent) in an NMR-style ensemble model and writes to:
        -a flat .csv file
        -a .pml file in order to add those solvent atoms and their positions as "pseudoatoms" to an existing pymol session for the purposes of displaying the solvent atoms by cartoon putty with a color scale based on B-factor. 

Notes:

 -This script considers any solvent atoms from different model states that fall within 3A of an existing solvent atom to be the same atom as the atom in the previous state. That is, the identity of each solvent molecule is determined by making a 3A^3 box around each modeled position and collecting the statistics (mean position and rmsf) for solvent molecules modeled inside those 3A^3 boxes. 
 
 -A max RMSF (or B-factor) can be provided for .pml output mode so the user can make comparable visualizations across different ensembles.

##Example usage:

phenix.python ensem_water_rmsf.py ensem.pdb O 
phenix.python ensem_water_rmsf.py ensem.pdb O 

##Conditions to run:

###Software requirements:

phenix==1.9-1692
###Additional requirements:

A multi-MODEL PDB file with HETATMs
