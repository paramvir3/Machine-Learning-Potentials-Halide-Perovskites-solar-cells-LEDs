import os
import datetime
import warnings
from typing import Literal, Any
from collections.abc import Callable
import traceback
from copy import deepcopy
from importlib.metadata import version
import json

import numpy as np
import pandas as pd

from tqdm import tqdm

from ase.constraints import FixSymmetry
from ase.filters import ExpCellFilter, FrechetCellFilter
from ase.optimize import FIRE, LBFGS
from ase.optimize.optimize import Optimizer
from ase.spacegroup import get_spacegroup
from ase import Atoms
from ase.io import read, write

from pynep.calculate import NEP

from k_srme import aseatoms2str,  two_stage_relax, ID, STRUCTURES, NO_TILT_MASK
from k_srme.conductivity import init_phono3py, get_fc2_and_freqs, get_fc3, calculate_conductivity

calc = NEP('nep.txt')

# Relaxation parameters
ase_optimizer : Literal["FIRE", "LBFGS","BFGS"] = "LBFGS"
ase_filter: Literal["frechet", "exp"] = "frechet"
if_two_stage_relax = True # Use two-stage relaxation enforcing symmetries
max_steps = 300
force_max = 1e-4  # Run until the forces are smaller than this in eV/A

# Symmetry parameters
# symmetry precision for enforcing relaxation and conductivity calculation
symprec = 1e-5 
# Enforce symmetry with during relaxation if broken
enforce_relax_symm = True 

prog_bar = True
save_forces = True # Save force sets to file

filter_cls: Callable[[Atoms], Atoms] = {
    "frechet": FrechetCellFilter,
    "exp": ExpCellFilter,
}[ase_filter]
optim_cls: Callable[..., Optimizer] = {"FIRE": FIRE, "LBFGS": LBFGS}[ase_optimizer]

run_params={
    "k_srme_version": version("k_srme"),
    "ase_optimizer": ase_optimizer,
    "ase_filter": ase_filter,
    "if_two_stage_relax": if_two_stage_relax,
    "max_steps": max_steps,
    "force_max": force_max,
    "symprec": symprec,
    "enforce_relax_symm": enforce_relax_symm,
}

atoms = read("POSCAR")

atoms.calc = calc

atoms, relax_dict = two_stage_relax(
    atoms,
    fmax_stage1=force_max,
    fmax_stage2=force_max,
    steps_stage1=max_steps,
    steps_stage2=max_steps,
    Optimizer = optim_cls,
    Filter = filter_cls,
    allow_tilt = True,
    log=True,#f'log.log',
    enforce_symmetry = enforce_relax_symm,
)

write('POSCAR_relaxed', atoms, format='vasp')

atoms.calc = None

atoms.info["fc2_supercell"] = [5,13,3]
atoms.info["fc3_supercell"] = [3,8,2]
atoms.info["q_mesh"] = [5,13,3]

ph3 = init_phono3py(atoms,
                    log=True,
                    symprec=symprec)

ph3, fc2_set, freqs = get_fc2_and_freqs(
                            ph3,
                            calculator = calc,
                            log=True,
                            pbar_kwargs={"leave":False,"disable":not prog_bar} 
                        )

ph3, fc3_set = get_fc3(
                    ph3,
                    calculator = calc,
                    log=True,
                    pbar_kwargs={"leave":False,"disable":not prog_bar} 
                )


# Save force constants to npz file
npz_filename = "force_sets.npz"
np.savez(npz_filename, fc2_set=fc2_set, fc3_set=fc3_set)
print(f"Force constants saved to {npz_filename}")
