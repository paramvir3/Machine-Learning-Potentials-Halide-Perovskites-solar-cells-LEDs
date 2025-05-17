from ase.io import read, write
from ase.io.espresso import read_espresso_out
from ase.visualize import view

in_filename = './movie.xyz'

out_filename = './out.pdb'

all_atoms = read(in_filename, format='extxyz', index=':')

write(out_filename, all_atoms, append=False, format='proteindatabank')

#all_atoms = read('./mol.extxyz', format='extxyz', index=':')

#view(all_atoms)
