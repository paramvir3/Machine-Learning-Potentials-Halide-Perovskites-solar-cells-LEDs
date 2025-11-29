from ase.io import read, write
from ase.build.tools import sort
atoms_cubic = read('POSCAR_cubic').repeat((6,6,12))
atoms_hex   = read('POSCAR_delta').repeat((3,12,3))
atoms_gamma   = read('POSCAR_gamma').repeat((6,6,3))

atoms_cubic = sort(atoms_cubic)
write('model.extxyz', atoms_cubic, format='extxyz')
