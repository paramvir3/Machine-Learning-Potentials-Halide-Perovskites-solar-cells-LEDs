from ase.io import read, write
from ase.build.tools import sort
atoms_cubic = read('POSCAR_cubic').repeat((12,12,12))
atoms_hex   = read('POSCAR_delta').repeat((6,18,4))
atoms_gamma   = read('POSCAR_gamma').repeat((8,9,6))

atoms_hex = sort(atoms_hex)
write('model.extxyz', atoms_hex, format='extxyz')
