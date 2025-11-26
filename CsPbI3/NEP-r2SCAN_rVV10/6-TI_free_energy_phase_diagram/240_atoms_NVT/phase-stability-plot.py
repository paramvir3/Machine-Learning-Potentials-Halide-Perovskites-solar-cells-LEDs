import sys
import os
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams

# Use 'Agg' backend for non-interactive plotting
matplotlib.use('Agg')

# -----------------------------------------------------------------------------
# 1. PRX PLOT STYLE SETTINGS
# -----------------------------------------------------------------------------
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Tahoma'] 

matplotlib.rcParams.update({
    'font.size': 12,             
    'text.usetex': True,         
    'axes.labelsize': 12,        
    'xtick.labelsize': 12,       
    'ytick.labelsize': 12,       
    'legend.fontsize': 12,       
    'axes.linewidth': 0.8,       
    'xtick.direction': 'in',     
    'ytick.direction': 'in',     
    'xtick.major.size': 4,       
    'ytick.major.size': 4,       
    'xtick.minor.size': 2,       
    'ytick.minor.size': 2,       
    'xtick.top': True,           
    'ytick.right': True,         
    'legend.frameon': False,     
    'savefig.dpi': 300,          
    'figure.figsize': (5, 5),    
    'grid.linestyle': '--',      
    'grid.alpha': 0.6,           
    'lines.linewidth': 2,        
    'xtick.major.width': 1.2,    
    'ytick.major.width': 1.2,    
    'xtick.minor.visible': True, 
    'ytick.minor.visible': True, 
    'xtick.minor.width': 1,      
    'ytick.minor.width': 1,      
    'axes.spines.top': True,     
    'axes.spines.right': True,   
    'axes.spines.bottom': True,  
    'axes.spines.left': True,    
    'axes.edgecolor': '0.1',     
})

# LaTeX preamble
matplotlib.rcParams['text.latex.preamble'] = '\n'.join([
    r'\usepackage{siunitx}',     
    r'\sisetup{detect-all}',     
    r'\usepackage{amsmath, amssymb}', 
])

# -----------------------------------------------------------------------------
# 2. COLOR PALETTE
# -----------------------------------------------------------------------------
# Defining specific colors for the phases
COLOR_CUBIC = '#5DA5DA'  # Blue
COLOR_HEX = '#F15854'    # Red
COLOR_LINE = '#4D4D4D'   # Dark Gray for the main line

# -----------------------------------------------------------------------------
# 3. DATA PARSING
# -----------------------------------------------------------------------------

def parse_delta_g(filename):
    """
    Parses the summary table to extract Temperature and Delta G.
    Returns numpy arrays for easier plotting manipulation.
    """
    temp_list = []
    delta_g_list = []
    
    try:
        with open(filename, 'r') as f:
            lines = f.readlines()
            
        start_reading = False
        
        for line in lines:
            if "Temp (K)" in line and "Delta G" in line:
                start_reading = True
                continue
            
            if start_reading and "---" in line:
                continue
            
            if start_reading:
                parts = line.split('|')
                if len(parts) >= 4:
                    try:
                        t = float(parts[0].strip())
                        dg = float(parts[3].strip())
                        temp_list.append(t)
                        delta_g_list.append(dg)
                    except ValueError:
                        continue
                        
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
        sys.exit(1)

    return np.array(temp_list), np.array(delta_g_list)

# -----------------------------------------------------------------------------
# 4. PLOTTING FUNCTION
# -----------------------------------------------------------------------------

def plot_phase_stability(temp, delta_g, output_filename="phase_stability_filled.png"):
    fig, ax = plt.subplots()

    # Labels
    ax.set_xlabel(r'Temperature [K]', fontsize=22)
    ax.set_ylabel(r'$\Delta G$ [eV/atom]', fontsize=22)
    ax.set_title(r'Phase Stability Map', fontsize=14)

    # 1. Plot the main Delta G line
    ax.plot(temp, delta_g, color=COLOR_LINE, linewidth=2, label=r'$\Delta G = G_{\text{hex}} - G_{\text{cubic}}$')

    # 2. Add a horizontal line at 0
    ax.axhline(0, color='black', linestyle='-', linewidth=0.8, alpha=0.8)

    # 3. Fill logic
    # Fill RED where Hexagonal is stable (Delta G < 0)
    ax.fill_between(temp, delta_g, 0, 
                    where=(delta_g <= 0), 
                    interpolate=True, 
                    color=COLOR_HEX, 
                    alpha=0.2, 
                    label=r'Hexagonal ($\delta$) Stable')

    # Fill BLUE where Cubic is stable (Delta G > 0)
    ax.fill_between(temp, delta_g, 0, 
                    where=(delta_g >= 0), 
                    interpolate=True, 
                    color=COLOR_CUBIC, 
                    alpha=0.2, 
                    label=r'Cubic ($\alpha$) Stable')

    # 4. Annotate the Transition Temperature (approximate zero crossing)
    # Find index where sign changes
    idx = np.where(np.diff(np.sign(delta_g)))[0]
    if len(idx) > 0:
        # Simple linear interpolation for the zero crossing
        x1, x2 = temp[idx[0]], temp[idx[0]+1]
        y1, y2 = delta_g[idx[0]], delta_g[idx[0]+1]
        t_c = x1 - y1 * (x2 - x1) / (y2 - y1)
        
        # Plot a vertical dashed line at the transition
        ax.axvline(t_c, color='black', linestyle='--', linewidth=1, alpha=0.6)
        
        # Add text label for Tc
        ax.text(t_c, max(delta_g)*0.1, r'$T_c \approx \SI{' + f"{t_c:.0f}" + r'}{K}$', 
                fontsize=12, ha='right', rotation=90, verticalalignment='bottom')

    # Legend
    ax.legend(loc='upper left', frameon=False)
    
    # Grid and Layout
    ax.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()

    # Save
    print(f"Saving plot to {output_filename} and .pdf version...")
    plt.savefig(output_filename, format='png', dpi=600)
    plt.savefig(output_filename.replace('.png', '.pdf'), format='pdf', dpi=600)

# -----------------------------------------------------------------------------
# 5. MAIN
# -----------------------------------------------------------------------------

def main():
    input_file = "output"
    
    if not os.path.exists(input_file):
        print(f"File '{input_file}' not found.")
        return

    print(f"Reading data from {input_file}...")
    temp, delta_g = parse_delta_g(input_file)
    
    if len(temp) == 0:
        print("No data found.")
        return

    plot_phase_stability(temp, delta_g)
    print("Done!")

if __name__ == "__main__":
    main()
