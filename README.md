# <span style="font-size:larger;">ML-PERO-0</span>

[![arXiv](https://img.shields.io/badge/arXiv-1234.56789-b31b1b.svg)](https://arxiv.org/abs/2404.05644)
[![arXiv](https://img.shields.io/badge/arXiv-1234.56789-b31b1b.svg)](https://arxiv.org/abs/2405.11599) 

## General Information
- GRACE-1Layer-FOUNDATIONAL model can be used for:
  
  1. OD, 1D, 2D, and 3D halide perovskites
  2. Any other components of perovskite solar cells including ETLs like TiO2, SnO2, SAMs etc, and including SPIRO and other organic HTMs etc, all kinds of passivating molecules, pseudo-halides mixtures
  3. All kinds of phase transitions can be performed for example titled phase transitions, mixed cation-anion perovskite crystallization from homogeneous mixtures of ions, non-perovskite to perovskite phase transitions, nucleation from solutions etc.
  4. Also perovskite LEDs and their components
     
  Important Note: Testing is under progress for many aspects since last year, finetuning is recommended for targeted problems
     
- Use experimental structures https://github.com/paramvir3/Crystal-Structures-Halide-perovskite of halide perovskites

## Codes Used
* ASE:
* ACE-Julia:
* GRACE:
* k_SRME:https://github.com/MPA2suite/k_SRME.git
* LAMMPS:
* NEP-GPUMD:
* Nequip:
* PET-MAD:
* PLUMED:https://github.com/paramvir3/plumed2
* Quantum espresso:
* VASP:
* 7Net:

## Features
- CsPbBr<sub>3</sub> [1] -- NEQUIP message passing machine learning interatomic potentials

![melt_crystal](https://github.com/ahlawat-paramvir/MLIP-Perovskites/assets/10708344/803ad827-2fea-4ed7-8696-f46d1f5ee1fe)

- CsPbI<sub>3</sub> [2] -- Neuroevolution Potential (NEP) and GPUMD
1. Benchmarking against experiments:

* &Delta; H: relative enthalpy (Kj/mol)
<img src="https://github.com/user-attachments/assets/be4c6bf8-8d44-46df-a1b3-7c2e845c0d32" alt="relative energies" width="200" >


* C<sub>P</sub>: Heat Capacity of tilted phase transitions γ --> β --> α

<img src="https://github.com/user-attachments/assets/a5a565b1-fe0c-436b-b473-e6eed0ce9562" alt="latice_parameters" width="215" >

<img src="https://github.com/user-attachments/assets/1e31df37-bc42-4569-9f81-666180eb6392" alt="heat capacity" width="250" >


* T<sub>m</sub>: melting point (K), T<sub>m</sub> <sub>experiment</sub> = 750K

<img src="https://github.com/user-attachments/assets/b1f918da-f912-4a2a-958d-05749fb9266f" alt="Melting point" width="200" >


* k<sub>T</sub>: lattice thermal conductivity

<img src="https://github.com/user-attachments/assets/b5b976e7-2a8f-44af-b2b9-2c47da025a36" alt="Lattice thermal conductivity" width="200" >


2. Simulations for designing and improving experiments for solar cells and LEDs [11-13]:

* Smaller scale MD simulations

![new_gif](https://github.com/ahlawat-paramvir/MLIP-Perovskites/assets/10708344/1f028241-0ac0-4797-ba8a-91ec38bfbfea)

* Size Matters: device scale million atoms MD simulations revealing Zig-Zag Ruddlesden-Popper (RP) grain boundaries 

![SM6-ezgif com-video-to-gif-converter](https://github.com/user-attachments/assets/712eff9b-64b5-4c7f-9a9f-5f194a9c48c6)




## Acknowledgements
- Swiss National Science Foundation through post-doc mobility Fellowship No. P500PN_206693

## References
1. Lattice matched heterogeneous nucleation eliminate defective buried interface in halide perovskites: https://doi.org/10.1021/acs.chemmater.4c03034 , DFT dataset: https://doi.org/10.5281/zenodo.10975237
2. Size dependent solid-solid crystallization of halide perovskites. https://doi.org/10.48550/arXiv.2404.05644
3. Crystallization of FAPbI3: Polytypes and stacking faults." The Journal of Chemical Physics 159.15 (2023): https://doi.org/10.1063/5.0165285
4. Organic Spacers in 2D Perovskites: General Trends and Structure‐Property Relationships from Computational Studies. Helvetica Chimica Acta 104.4 (2021): e2000232: https://doi.org/10.1002/hlca.202000232
5. A combined molecular dynamics and experimental study of two-step process enabling low-temperature formation of phase-pure α-FAPbI3.Sci. Adv.7,eabe3326(2021) https://www.science.org/doi/10.1126/sciadv.abe3326
6. https://doi.org/10.1002/adma.202001069, https://doi.org/10.1021/acsnano.8b00267
7. https://doi.org/10.1126/science.aax3878
8. https://doi.org/10.1073/pnas.1711744114
9. 10.1524/zpch.1992.175.part_1.063

11. Vertically stacked monolithic perovskite colour photodetectors. Nature 642, 592–598 (2025). https://doi.org/10.1038/s41586-025-09062-3
12. Intragrain 3D perovskite heterostructure for high-performance pure-red perovskite LEDs. Nature 641, 352–357 (2025). https://doi.org/10.1038/s41586-025-08867-6
13. Ruddlesden–Popper Defects Act as a Free Surface: Role in Formation and Photophysical Properties of CsPbI3: https://doi.org/10.1002/adma.202501788
14. Nanoscale heterophase regulation enables sunlight-like full-spectrum white electroluminescence. Nat Commun 16, 3621 (2025). https://doi.org/10.1038/s41467-025-58743-0

## Contact
paramvir.chem@gmail.com

