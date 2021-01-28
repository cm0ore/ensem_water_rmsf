
from __future__ import division
import iotbx.pdb
import sys
import os
import math
import iotbx.pdb.amino_acid_codes

aa_resnames = iotbx.pdb.amino_acid_codes.one_letter_given_three_letter

file_name = sys.argv[1]
pdb_obj = iotbx.pdb.hierarchy.input(file_name=file_name)

targ_atom_names_str = sys.argv[2]
targ_atom_names = []
for targ_atom_name in targ_atom_names_str.split(","):
  targ_atom_names.append(targ_atom_name)
print >> sys.stderr, "Calculating RMSF using these atoms:", " ".join(targ_atom_names)

assert sys.argv[3] in ['csv', 'pml']
output = sys.argv[3]

added_atoms = {} # (model id, atom xyz) --> (model id, atom xyz) FOR BOOKKEEPING!
xyzs = {} # (model id, atom xyz) --> list of xyzs in each 3*3*3A box
mean_xyz = {} # (model id, atom xyz) --> mean xyz
atom_count = {} # atom xyz --> # of instances (multi-MODEL or alt confs)
for model in pdb_obj.hierarchy.models():
  for chain in model.chains():
    for rg in chain.residue_groups():
      if rg.atom_groups()[0].resname in aa_resnames: continue
      for ag in rg.atom_groups():
        if ag.resname == 'HOH':
	  for atom in ag.atoms():
	    atom_name = (int(model.id), atom.xyz)
	    if int(model.id) == 1:
	      mean_xyz[atom_name] = atom.xyz 
	      xyzs[atom_name] = [atom.xyz]
	      atom_count[atom_name] = 1
	      added_atoms[atom_name] = atom_name
	    
	    for element in mean_xyz:
	      coords = element[1]
	      distance = tuple(abs(i-j) for i, j in zip(coords, atom.xyz))
	      if distance[0] > 3 or distance[1] > 3 or distance[2] > 3:   #size of each box is 3^3 A^3
		continue
		
	      elif element != atom_name:
		mean_xyz[element] = tuple(i + j for i,j in zip(mean_xyz[element], atom.xyz)) #update the summed coordinates of that box
		added_atoms[atom_name] = atom_name     #this is just registering that we've added the atom to an existing box in the mean_xyz dic
		atom_count[element] += 1               #for mean and rmsf calculations later
		xyzs[element].append(atom.xyz)        #add the new atom's coordinates to the list of coordinates for that box (for the rmsf calculation later)
	 
            if atom_name not in added_atoms:           #append any loners to mean_xyz{} by giving them a new box of their own
	      mean_xyz[atom_name] = atom.xyz
	      xyzs[atom_name] = [atom.xyz]
	      atom_count[atom_name] = 1
		
	    

for atom_name in mean_xyz:
  xyz_sum = mean_xyz[atom_name] #summed coordinates
  mean_xyz[atom_name] = tuple(i / atom_count[atom_name] for i in xyz_sum)

summed_distances = {}
rmsfs = {} # HOH atom coordinates --> dist
for element in mean_xyz:
  for atom in xyzs[element]:                  #xyzs has a list of all the coordinates of Os in each box
    xyz_pairs = zip(atom, mean_xyz[element])
    sq_dist = sum((i - j)**2 for i, j in xyz_pairs)
    if element not in summed_distances:
	summed_distances[element] = sq_dist
    else:
      summed_distances[element] += sq_dist    #summed squared distances of all O atom from the mean for each box

for element in summed_distances: 	      
  summed_distances[element] /= atom_count[element] 	#take an avg of the summed distances
  summed_distances[element] = math.sqrt(summed_distances[element])     #take the root to get the rmsf

max_rmsf = max(summed_distances.values())
if len(sys.argv) > 4:
  max_rmsf = float(sys.argv[4])


file_basename = os.path.basename(file_name).split('.pdb')[0]
if output == 'csv':
  print 'dataset,model_id,atom_count,mean_atomic_coordinates,rmsf'

for element in sorted(summed_distances):
	
  rmsf = summed_distances[element]
  if output == 'csv':
    print '%s,%s,%s,%s,%.4f' % \
      (file_basename, element[0], atom_count[element], mean_xyz[element], rmsf)

  elif output == 'pml':
    bfactor = rmsf ** 2
    maxbfactor = max_rmsf ** 2
    meanx = mean_xyz[element][0]
    meany = mean_xyz[element][1]
    meanz = mean_xyz[element][2]
    #scaled_rmsf = (rmsf / max_rmsf) * 100
    scaled_bfactor = (bfactor / maxbfactor) * 100
    print 'pseudoatom water, chain=ZZ, state=-1, atoms_in_site= %d, b=%.3f, pos=[%.4f,%.4f,%.4f]' % \
      (atom_count[element], scaled_bfactor, meanx, meany, meanz)



