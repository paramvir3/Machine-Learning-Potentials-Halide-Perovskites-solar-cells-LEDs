from ase.io import read, write

atoms = read('POSCAR')

write('model.xyz', atoms, format='extxyz')
