
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


from k_srme import aseatoms2str,  two_stage_relax, ID, STRUCTURES, NO_TILT_MASK
from k_srme.conductivity import init_phono3py, get_fc2_and_freqs, get_fc3, calculate_conductivity, calculate_fc3_set, load_force_sets



atoms = read("POSCAR_relaxed")
atoms.info["fc2_supercell"] = [5,13,3]
atoms.info["fc3_supercell"] = [3,8,2]
atoms.info["q_mesh"] = [5,13,3]

init_phono3py

npz_filename = "force_sets.npz"

data = np.load(npz_filename)
fc2_set = data['fc2_set']
fc3_set = data['fc3_set']


ph3 = init_phono3py(atoms,
                    log=True,
                    symprec=1e-5)


ph3 = load_force_sets(ph3,fc2_set,fc3_set)

temp_range = np.arange(100, 601, 25)

ph3, kappa_dict = calculate_conductivity(ph3,temperatures=temp_range,log=True)

###with open("kappa_dict.json", "w") as f:
###    json.dump(kappa_dict, f, indent=2)
###
###filename = "kappa_dict.dat"
###
###with open(filename, "w") as f:
###    for key, value in kappa_dict.items():
###        f.write(f"{key}\t{value}\n") 
