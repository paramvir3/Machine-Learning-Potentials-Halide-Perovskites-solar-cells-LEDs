from ase.io import read, write
from ase.build.tools import sort
#atoms =  read('hexagonal_FAPI.cif')
atoms =  read('hexa.vasp')

atoms = sort(atoms)
write('POSCAR', atoms, format='vasp')
