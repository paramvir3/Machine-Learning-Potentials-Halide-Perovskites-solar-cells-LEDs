from ase.io import read, write
from ase.build.tools import sort
atoms_cubic = read('POSCAR_cubic').repeat((6,6,6))
atoms_hex   = read('POSCAR_delta').repeat((3,12,3))
atoms_gamma   = read('POSCAR_gamma').repeat((6,6,3))

atoms_hex = sort(atoms_hex)
write('model.extxyz', atoms_hex, format='extxyz')
