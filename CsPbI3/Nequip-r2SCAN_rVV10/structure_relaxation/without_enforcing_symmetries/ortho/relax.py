from ase import Atoms
import sys
from ase.io import read, write
import warnings, os
from typing import Type
import numpy as np
from ase.build import bulk
from ase.build.tools import sort
from ase.io import read, write
from ase.io.vasp import read_vasp
from ase.calculators.calculator import Calculator
import ase
from ase.constraints import FixSymmetry
from ase.filters import UnitCellFilter, ExpCellFilter, StrainFilter,FrechetCellFilter
from ase.optimize import LBFGS, FIRE, MDMin, GPMin
from ase.spacegroup import get_spacegroup
from ase.spacegroup.symmetrize import check_symmetry
from nequip.ase import NequIPCalculator


def opt_with_symmetry(
    atoms_in: Atoms,
    calculator: Calculator,
    fix_symmetry: bool = False,
    hydrostatic_strain: bool = False,
) -> Atoms:
    atoms = atoms_in.copy()
    atoms.calc = calculator
    if fix_symmetry:
        atoms.set_constraint([FixSymmetry(atoms)])
    ecf = FrechetCellFilter(atoms, hydrostatic_strain=hydrostatic_strain)
    #log = 'relax_%d.log' %(x) 
    opt = LBFGS(ecf, logfile='log')
    opt.run(fmax=0.005, steps=10000) 
    cell_diff = (atoms.cell.cellpar() / atoms_in.cell.cellpar() - 1.0) * 100
    #ene_diff = (atoms.get_potential_energy() - atoms_in.get_potential_energy())
    print(f"Optimization finished with steps = {opt.nsteps}")
    print("Optimized Cell         :", atoms.cell.cellpar())
    print("Optimized Cell diff (%):", cell_diff)
    print("Scaled positions       :\n", atoms.get_scaled_positions())
    print(f"Epot after opt: {atoms.get_potential_energy()} eV")
    #print("Energy difference:", ene_diff)
    return atoms

unit_atoms=read("ortho.cif")

crystal_atoms = unit_atoms.repeat((1,2,1))

crystal_atoms = sort(crystal_atoms)

write('POSCAR_initial',crystal_atoms,format='vasp')

print("Initial symmetry at precision 1e-3")

check_symmetry(crystal_atoms, 1.0e-3, verbose=True)

calculator = NequIPCalculator.from_deployed_model(model_path="best.pth", device='cpu')
atoms_new = opt_with_symmetry(crystal_atoms, calculator, fix_symmetry=False, hydrostatic_strain=False)

# We print out the initial symmetry groups at two different precision levels
print("Final symmetry at precision 1e-3")
check_symmetry(atoms_new, 1.0e-3, verbose=True)

new_poscar = 'POSCAR_relaxed'
atoms_new = sort(atoms_new)

write(new_poscar,atoms_new,format='vasp')
