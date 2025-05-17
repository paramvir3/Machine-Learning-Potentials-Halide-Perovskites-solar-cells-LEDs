import re
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pl

# Define a more robust function to extract the last total potential energy
def extract_last_total_potential_energy(file_path):
    """
    Reads a file and extracts the last reported total potential energy.

    Parameters:
        file_path (str): Path to the output file.

    Returns:
        float or None: The last total potential energy value found, or None if not found.
    """
    energies = []
    with open(file_path, 'r') as file:
        for line in file:
            if "total_potential" in line and "eV" in line:
                match = re.search(r"total_potential\s*=\s*([-+]?\d*\.\d+|\d+)\s*eV", line)
                if match:
                    energies.append(float(match.group(1)))
    return energies[-1] if energies else None

# Assume these functions have already extracted the energy values from different paths
delta_phase = extract_last_total_potential_energy("delta/output")
ortho_phase = extract_last_total_potential_energy("ortho/output")
beta_phase = extract_last_total_potential_energy("beta/output")
cubic_phase = extract_last_total_potential_energy("cubic/output")
delta_FAPI = extract_last_total_potential_energy("hexagonal_FAPI/output")

# Define the phases and their corresponding energy values
phases = [r'$\delta$-Phase', r'$\gamma$-Phase', r'$\beta$-Phase', r'$\alpha$-Phase', r'$\delta$-FAPI']
energy_values = [delta_phase, ortho_phase, beta_phase, cubic_phase, delta_FAPI]

# Conversion factor (if necessary, assuming it converts from eV to another unit, like kJ/mol)
ev_to_kjmol = 96.485332 / 4

# Calculate relative values with respect to delta_phase, applying the conversion factor consistently
# Filter out None values to avoid errors in calculation
valid_energies = [(name, energy) for name, energy in zip(phases, energy_values) if energy is not None]

if not valid_energies or valid_energies[0][1] is None:
    print("Error: Could not extract a valid reference energy for the delta phase.")
else:
    delta_reference_energy = valid_energies[0][1]
    relative_energy_values = [(energy * ev_to_kjmol) - (delta_reference_energy * ev_to_kjmol) for _, energy in valid_energies]
    valid_phases = [name for name, _ in valid_energies]

    # Create a figure with one subplot
    fig, ax2 = pl.subplots(1, 1, figsize=(6, 6))

    # Define font size for all labels and titles
    font_size = 16

    # Plot 2: Relative energy values to delta phase
    ax2.bar(valid_phases, relative_energy_values, color=['skyblue', 'salmon', 'lightgreen', 'orange', 'blue'][:len(valid_phases)])
    ax2.set_xlabel('CsPbI$_3$ polymorphs', fontsize=font_size)
    ax2.set_ylabel('Relative Energy (KJ/mol)', fontsize=font_size)
    ax2.set_title('Relative Energy Values with respect to $\delta$-CsPbI$_3$', fontsize=font_size)
    ax2.tick_params(axis='x', rotation=45, labelsize=font_size - 2)
    ax2.tick_params(axis='y', labelsize=font_size - 2)
    ax2.grid(axis='y', linestyle='--', alpha=0.7)

    # Adjust layout and display
    pl.tight_layout()
    pl.show()

    fig.savefig("rE.pdf", dpi=300, bbox_inches="tight", transparent=True)
    fig.savefig("rE.png", dpi=300, bbox_inches="tight", transparent=False)
