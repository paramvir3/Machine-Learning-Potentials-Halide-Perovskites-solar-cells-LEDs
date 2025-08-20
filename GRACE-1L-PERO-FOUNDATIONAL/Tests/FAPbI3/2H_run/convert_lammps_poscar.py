from ase.io import read, write
from ase.io.espresso import read_espresso_out
from ase.visualize import view
from ase.build import bulk
from ase.build.tools import sort

all_atoms = read('data.relax', format='lammps-data')

supercell = all_atoms.repeat((3,3,2))

all_atoms=sort(supercell)

write('POSCAR.vasp', all_atoms, format='vasp')
