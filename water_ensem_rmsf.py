from __future__ import division
import iotbx.pdb
import sys
import os
import math
import iotbx.pdb.amino_acid_codes

aa_resnames = iotbx.pdb.amino_acid_codes.one_letter_given_three_letter

file_name = sys.argv[1]
pdb_obj = iotbx.pdb.hierarchy.input(file_name=file_name)
#n_models = len(pdb_obj.hierarchy.models()) -- deprecated in favor of atom_count[atom_name]

targ_atom_names_str = sys.argv[2]
targ_atom_names = []
for targ_atom_name in targ_atom_names_str.split(","):
  targ_atom_names.append(targ_atom_name)
print >> sys.stderr, "Calculating RMSF using these atoms:", " ".join(targ_atom_names)

assert sys.argv[3] in ['csv', 'pml']
output = sys.argv[3]

#if output == 'pml':
#  max_rmsf = float(sys.argv[3])

mean_xyz = {} # atom xyz --> xyz
atom_count = {} # atom xyz --> # of instances (multi-MODEL or alt confs)
for model in pdb_obj.hierarchy.models():
  if int(model.id) > 20:
    continue
  for chain in model.chains():
    for rg in chain.residue_groups():
      if rg.atom_groups()[0].resname in aa_resnames: continue
      for ag in rg.atom_groups():
        if ag.resname == 'HOH':
	  for atom in ag.atoms():
	    atom_name = (model.id, atom.xyz)
	    mean_xyz[atom_name] = atom.xyz #seems redundant but actually useful
	    atom_count[atom_name] = 1
            #print('name', atom_name)
	    for element in mean_xyz:
	      if element[0] == model.id:
		continue
	      coords = element[1]
	      distance = tuple(abs(i-j) for i, j in zip(coords, atom.xyz))
	      if distance[0] > 1 or distance[1] > 1 or distance[2] > 1:
		continue
	      #elif distance[0] == 0 and distance[1] == 0 and distance[2] == 0:
		#continue
	      else:
		mean_xyz[element] = tuple(i + j for i,j in zip(mean_xyz[element], atom.xyz))
		atom_count[element] += 1
	        mean_xyz[atom_name] = tuple(i + j for i,j in zip(element[1], atom.xyz))
		atom_count[atom_name] += 1


for atom_name in mean_xyz:
  xyz_sum = mean_xyz[atom_name] #summed coordinates
  mean_xyz[atom_name] = tuple(i / atom_count[atom_name] for i in xyz_sum)


rmsfs = {} # HOH atom coordinates --> dist
for model in pdb_obj.hierarchy.models():
  if int(model.id) > 20:
    continue
  for chain in model.chains():
    for rg in chain.residue_groups():
      if rg.atom_groups()[0].resname in aa_resnames: continue
      for ag in rg.atom_groups():
        if ag.resname == 'HOH':
          for atom in ag.atoms():
            atom_name = (model.id, atom.xyz)

	    xyz_pairs = zip(atom.xyz, mean_xyz[atom_name])
            sq_dist = sum((i - j)**2 for i, j in xyz_pairs)
            #print 'dist: %s' % math.sqrt(sq_dist)
            if atom_name not in rmsfs:
              rmsfs[atom_name] = sq_dist 
            else:
              rmsfs[atom_name] = rmsfs[atom_name] + sq_dist

max_rmsf = max(rmsfs.values())
if len(sys.argv) > 4:
  max_rmsf = float(sys.argv[4])

for atom_name in rmsfs:
  sq_dist_sum = rmsfs[atom_name]
  rmsf = math.sqrt(sq_dist_sum / atom_count[atom_name])
  rmsfs[atom_name] = rmsf

file_basename = os.path.basename(file_name).split('.pdb')[0]
if output == 'csv':
  print 'dataset,model_id,atomic_coordinates,rmsf'
for atom_name in sorted(rmsfs):
  rmsf = rmsfs[atom_name]
  if rmsf < 0.0001:
    continue 
  if output == 'csv':
    print '%s,%s,%s,%.4f' % \
      (file_basename, atom_name[0], atom_name[1], rmsf)
  elif output == 'pml':
    bfactor = rmsf ** 2
    maxbfactor = max_rmsf ** 2
    #scaled_rmsf = (rmsf / max_rmsf) * 100
    scaled_bfactor = (bfactor / maxbfactor) * 100
    print 'pseudoatom water, chain=ZZ, state=%s,  b=%.3f, pos=[%.4f,%.4f,%.4f]' % \
      (atom_name[0], scaled_bfactor, atom_name[1][0], atom_name[1][1], atom_name[1][2])



