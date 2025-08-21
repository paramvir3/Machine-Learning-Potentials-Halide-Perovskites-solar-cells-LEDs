from ase.io import read, write
from ase.io.espresso import read_espresso_out
from ase.visualize import view
from ase.build import bulk
from ase.build.tools import sort

all_atoms = read('spiro_ss.extxyz')

supercell = all_atoms.repeat((1,1,1))

all_atoms=sort(supercell)

write('structure.data', all_atoms, format='lammps-data')
