import numpy as np
from ase import units
from ase.calculators.calculator import Calculator, all_changes
from ase.md.langevin import Langevin
from ase.build import bulk
from ase.io import read
import os
from ase.optimize import FIRE
from pynep.calculate import NEP

# --- CONFIGURATION ---
CUBIC_FILE = 'POSCAR_cubic'
DELTA_FILE = 'POSCAR_delta'

calc_nep = NEP('./nep.txt')

def relax_structure(atoms, real_calc):
    print("    > 🧘 Relaxing atomic positions (FIRE minimization)...")
    atoms.calc = real_calc
    opt = FIRE(atoms, logfile=None)
    opt.run(fmax=0.001, steps=10000)
    return atoms

# --- 1. ROBUST EINSTEIN CRYSTAL ---
class EinsteinCrystal(Calculator):
    """
    Reference Hamiltonian: Harmonic springs tethering atoms to fixed lattice sites.
    """
    implemented_properties = ['energy', 'forces']

    def __init__(self, atoms, k_springs={'Cs': 8.0, 'Pb': 20.0, 'I': 12.0}):
        super().__init__()
        self.k_springs = k_springs
        self.r0 = atoms.get_positions().copy() 
        self.symbols = atoms.get_chemical_symbols()
        
        # Pre-compute k array (N, 1)
        self.k_array = np.array([self.k_springs.get(s, 10.0) for s in self.symbols]).reshape(-1, 1)

    def calculate(self, atoms=None, properties=['energy'], system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        
        positions = atoms.get_positions()
        cell = atoms.get_cell()
        
        if np.linalg.det(cell) < 1e-8:
             inv_cell = np.eye(3)
        else:
             inv_cell = np.linalg.inv(cell)
        
        diff_cartesian = positions - self.r0
        
        # General Minimum Image Convention
        if any(atoms.get_pbc()):
            diff_frac = np.dot(diff_cartesian, inv_cell)
            diff_frac -= np.round(diff_frac)
            diff = np.dot(diff_frac, cell)
        else:
            diff = diff_cartesian

        # Energy: 0.5 * k * dr^2
        energy = 0.5 * np.sum(self.k_array * diff**2)
        forces = -self.k_array * diff
        
        self.results['energy'] = energy
        self.results['forces'] = forces

# --- 2. HYBRID CALCULATOR ---
class HybridCalculator(Calculator):
    implemented_properties = ['energy', 'forces', 'du_dl']

    def __init__(self, calc_ref, calc_real, lamb=0.0):
        super().__init__()
        self.calc_ref = calc_ref       
        self.calc_real = calc_real 
        self.lamb = lamb               
        self.du_dl = 0.0

    def set_lambda(self, lamb):
        self.lamb = lamb

    def calculate(self, atoms=None, properties=['energy', 'forces'], system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        
        # Reference (Springs)
        atoms.calc = self.calc_ref
        e_ref = atoms.get_potential_energy()
        f_ref = atoms.get_forces()
        
        # Real (GRACE)
        atoms.calc = self.calc_real
        e_real = atoms.get_potential_energy()
        f_real = atoms.get_forces()
        
        # Mix
        self.results['energy'] = (1.0 - self.lamb) * e_ref + self.lamb * e_real
        self.results['forces'] = (1.0 - self.lamb) * f_ref + self.lamb * f_real
        
        # Gradient dU/dlambda
        self.results['du_dl'] = e_real - e_ref
        self.du_dl = e_real - e_ref # store for access
        
        atoms.calc = self

# --- 3. ANALYTICAL FREE ENERGY ---
def get_analytical_free_energy(atoms, k_springs, T):
    """
    Absolute free energy of the classical Einstein Crystal reference.
    """
    kb = units.kB
    a_harmonic = 0.0
    symbols = atoms.get_chemical_symbols()
    
    for s in symbols:
        k_val = k_springs.get(s, 10.0)
        # Classical harmonic oscillator term relative to Ideal Gas
        # A = 1.5 * kT * ln( k / (2 * pi * kT) )
        term = 1.5 * kb * T * np.log(k_val / (2.0 * np.pi * kb * T))
        a_harmonic += term

    return a_harmonic

# --- 4. UTILS ---
def constrain_com(atoms, reference_positions):
    """ Removes Center of Mass drift """
    p = atoms.get_momenta()
    masses = atoms.get_masses()
    total_mass = np.sum(masses)
    com_p = np.sum(p, axis=0)
    atoms.set_momenta(p - (com_p / total_mass) * masses[:, np.newaxis])
    
    # Re-center
    current_com = atoms.get_center_of_mass()
    ref_com = np.average(reference_positions, axis=0, weights=masses)
    drift = current_com - ref_com
    if np.linalg.norm(drift) > 1e-5:
        atoms.translate(-drift)

def calculate_optimal_k(atoms, real_calc, temp_k, n_steps=2000):
    """ Optimize Spring Constants via MSD """
    print(f"    > Optimizing Spring Constants (MSD) @ {temp_k} K...")
    atoms_opt = atoms.copy()
    atoms_opt.calc = real_calc
    r0 = atoms_opt.get_positions().copy()
    
    dyn = Langevin(atoms_opt, 0.5 * units.fs, temperature_K=temp_k, friction=0.02)
    dyn.run(n_steps // 4) # Equilibrate
    
    species_msd = {s: [] for s in set(atoms_opt.get_chemical_symbols())}
    
    for _ in range(n_steps):
        dyn.run(1)
        constrain_com(atoms_opt, r0)
        
        diff = atoms_opt.get_positions() - r0
        # Simple MIC check (assuming nearly orthogonal for MSD is okay, or use full)
        cell = atoms_opt.get_cell().diagonal()
        diff -= np.round(diff / cell) * cell
        
        sq_disp = np.sum(diff**2, axis=1)
        symbols = atoms_opt.get_chemical_symbols()
        for i, sym in enumerate(symbols):
            species_msd[sym].append(sq_disp[i])

    k_optimized = {}
    kb = units.kB
    for sym in species_msd:
        avg_msd = np.mean(species_msd[sym])
        avg_msd = max(avg_msd, 1e-5) # Avoid div/0
        k_val = (3.0 * kb * temp_k) / avg_msd
        k_optimized[sym] = k_val
        
    return k_optimized

# --- 5. MAIN SIMULATION ---
def run_ti(atoms, real_calc, temp_k, n_steps=5000, optimize_k=True):
    
    # 1. Optimize Spring Constants
    if optimize_k:
        k_springs = calculate_optimal_k(atoms, real_calc, temp_k, n_steps=2000)
    else:
        k_springs = {'Cs': 8.0, 'Pb': 20.0, 'I': 12.0}

    # 2. Analytical Part
    f_analytical = get_analytical_free_energy(atoms, k_springs, temp_k)
    
    # 3. Setup Hybrid
    ref_calc = EinsteinCrystal(atoms, k_springs=k_springs)
    hybrid = HybridCalculator(ref_calc, real_calc)
    atoms.calc = hybrid

    # 4. Gaussian Quadrature (7 points)
    deg = 7 
    lambdas, weights = np.polynomial.legendre.leggauss(deg)
    lambdas = 0.5 * (lambdas + 1)
    weights = 0.5 * weights
    
    gradients = []
    r0_ref = atoms.get_positions().copy()
    
    print(f"    > Running TI ({n_steps} steps/window)...")
    
    for i, (lam, w) in enumerate(zip(lambdas, weights)):
        hybrid.set_lambda(lam)
        
        dyn = Langevin(atoms, 0.5 * units.fs, temperature_K=temp_k, friction=0.02)
        dyn.run(n_steps // 5) # Equilibrate
        
        win_grads = []
        for _ in range(n_steps):
            dyn.run(1)
            constrain_com(atoms, r0_ref)
            win_grads.append(hybrid.du_dl)
        
        avg_grad = np.mean(win_grads)
        gradients.append(avg_grad)
        
    # 5. Integrate
    f_integral = np.sum(np.array(gradients) * weights)
    f_total = f_analytical + f_integral
    
    return f_total

# --- 6. EXECUTION LOOP ---
if __name__ == "__main__":
    
    # Define Temperature Range (Start, Stop+Step, Step)
    # 250, 300, 350, ..., 700
    TEMPS = range(300,700,50) 
    STEPS = 10000 # Production steps per window
    
    print(f"=== STARTING FREE ENERGY SCAN: {TEMPS[0]}K -> {TEMPS[-1]}K ===")
    
    # Check input files
    if not (os.path.exists(CUBIC_FILE) and os.path.exists(DELTA_FILE)):
        print("Error: POSCAR files not found. Please check paths.")
        exit()

    results = []

    for T in TEMPS:
        print(f"\n{'='*40}\n TEMPERATURE: {T} K\n{'='*40}")
        
        # Load fresh structures for each temperature
        # Ideally, you should RELAX (NPT) these at T before running TI
        # But here we use the provided POSCARs
      #  atoms_cubic = read(CUBIC_FILE).repeat((2,3,2))
      #  atoms_hex   = read(DELTA_FILE).repeat((2,6,1))
        atoms_cubic = read(CUBIC_FILE).repeat((4,4,3))
        atoms_hex   = read(DELTA_FILE).repeat((2,6,1))
        
        atoms_cubic = relax_structure(atoms_cubic, calc_nep)
        atoms_hex = relax_structure(atoms_hex, calc_nep)
        
        # 1. Cubic Phase
        print(f"--- Phase: Cubic (Alpha) ---")
        F_cub = run_ti(atoms_cubic, calc_nep, T, STEPS, optimize_k=True)
        F_cub_norm = F_cub / len(atoms_cubic)
        print(f"    -> G_cubic = {F_cub_norm:.5f} eV/atom")

        # 2. Hexagonal Phase
        print(f"--- Phase: Hexagonal (Delta) ---")
        F_hex = run_ti(atoms_hex, calc_nep, T, STEPS, optimize_k=True)
        F_hex_norm = F_hex / len(atoms_hex)
        print(f"    -> G_hex   = {F_hex_norm:.5f} eV/atom")
        
        # 3. Compare
        dG = F_hex_norm - F_cub_norm
        stable = "HEX (Delta)" if dG < 0 else "CUBIC (Alpha)"
        
        results.append({
            "T": T,
            "G_cub": F_cub_norm,
            "G_hex": F_hex_norm,
            "dG": dG,
            "Stable": stable
        })

    # --- FINAL SUMMARY TABLE ---
    print("\n" + "="*70)
    print(f"{'Temp (K)':<10} | {'G_cubic (eV/at)':<15} | {'G_hex (eV/at)':<15} | {'Delta G':<10} | {'Stable Phase'}")
    print("-" * 70)
    for r in results:
        print(f"{r['T']:<10} | {r['G_cub']:.5f}           | {r['G_hex']:.5f}           | {r['dG']:>6.4f}     | {r['Stable']}")
