import sys
import os
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from matplotlib import rcParams

# Use 'Agg' backend for non-interactive plotting
matplotlib.use('Agg')

# -----------------------------------------------------------------------------
# 1. PRX / PUBLICATION PLOT STYLE SETTINGS
# -----------------------------------------------------------------------------
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Tahoma'] 

matplotlib.rcParams.update({
    'font.size': 12,              
    'text.usetex': True,          
    'axes.labelsize': 14,         
    'xtick.labelsize': 12,        
    'ytick.labelsize': 12,        
    'legend.fontsize': 11,        
    'axes.linewidth': 1.0,        
    'xtick.direction': 'in',      
    'ytick.direction': 'in',      
    'xtick.major.size': 5,        
    'ytick.major.size': 5,        
    'xtick.minor.size': 2.5,      
    'ytick.minor.size': 2.5,      
    'xtick.top': True,            
    'ytick.right': True,          
    'legend.frameon': False,      
    'savefig.dpi': 300,           
    'figure.figsize': (4, 3.5),   
    'grid.linestyle': '--',       
    'grid.alpha': 0.5,            
    'lines.linewidth': 2,         
    'axes.edgecolor': '0.1',      
})

# LaTeX preamble
matplotlib.rcParams['text.latex.preamble'] = '\n'.join([
    r'\usepackage{siunitx}',      
    r'\sisetup{detect-all}',      
    r'\usepackage{amsmath, amssymb}', 
])

# -----------------------------------------------------------------------------
# 2. CONFIGURATION & COLORS
# -----------------------------------------------------------------------------
FILE_CUBIC = "rs_free_energy_cubic.csv"
FILE_DELTA = "rs_free_energy_delta.csv"

# Experimental Value
T_EXPERIMENTAL = 600  # Kelvin

# Colors
COLOR_CUBIC = '#5DA5DA'  # Blue
COLOR_DELTA = '#F15854'  # Red/Orange
COLOR_BLACK = '#4D4D4D'
COLOR_EXP   = '#60BD68'  # Green for experimental line
COLOR_CALC  = '#4D4D4D'  # Dark Grey/Black for calculated line

# [CRITICAL] Atom counts for normalization
N_ATOMS_CUBIC = 1
N_ATOMS_DELTA = 1

# -----------------------------------------------------------------------------
# 3. DATA PROCESSING LOGIC
# -----------------------------------------------------------------------------
def load_data(filename):
    try:
        if not os.path.exists(filename):
            print(f"Error: File '{filename}' not found.")
            return None, None

        df = pd.read_csv(filename, sep=r'\s+|,', engine='python')
        df.columns = df.columns.str.strip()
        
        t_col = [c for c in df.columns if 'temp' in c.lower() or 'T' in c]
        g_col = [c for c in df.columns if 'free' in c.lower() or 'energy' in c.lower() or 'G' in c]
        
        if t_col and g_col:
            return df[t_col[0]].values, df[g_col[0]].values
        else:
            return df.iloc[:, 0].values, df.iloc[:, 1].values
            
    except Exception as e:
        print(f"Error reading {filename}: {e}")
        return None, None

def find_crossing_point(T, G1, G2):
    diff = G2 - G1
    f_diff = interp1d(T, diff, kind='cubic')
    t_fine = np.linspace(T.min(), T.max(), 5000)
    diff_fine = f_diff(t_fine)
    zero_crossings = np.where(np.diff(np.sign(diff_fine)))[0]
    
    if len(zero_crossings) > 0:
        idx = zero_crossings[0]
        t1, t2 = t_fine[idx], t_fine[idx+1]
        d1, d2 = diff_fine[idx], diff_fine[idx+1]
        t_cross = t1 + (0 - d1) * (t2 - t1) / (d2 - d1)
        return t_cross
    return None

def main():
    # 1. Load Data
    print("--- Loading Phase Data ---")
    T_cub, G_cub = load_data(FILE_CUBIC)
    T_del, G_del = load_data(FILE_DELTA)
    
    if T_cub is None or T_del is None: 
        return

    # 2. Normalize
    if np.mean(G_cub) < -20 and N_ATOMS_CUBIC > 1:
        G_cub = G_cub / N_ATOMS_CUBIC
        
    if np.mean(G_del) < -20 and N_ATOMS_DELTA > 1:
        G_del = G_del / N_ATOMS_DELTA

    # Align Grids
    if not np.allclose(T_cub, T_del):
        f_del = interp1d(T_del, G_del, kind='linear', fill_value="extrapolate")
        G_del = f_del(T_cub)
        T_del = T_cub 

    # 3. Calc Delta G and Tc
    Delta_G = G_del - G_cub
    Tc = find_crossing_point(T_cub, G_cub, G_del)

    # -----------------------------------------------------------------------------
    # 4. PLOTTING
    # -----------------------------------------------------------------------------
    print("--- Generating Plot ---")
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
    plt.subplots_adjust(hspace=0.05)

    # --- Plot 1: Absolute G ---
    ax1.plot(T_cub, G_cub, 'o-', color=COLOR_CUBIC, label=r'Cubic ($\alpha$)', 
             markevery=100, markersize=0.1, markerfacecolor='white', markeredgewidth=0.1)
    ax1.plot(T_del, G_del, 's-', color=COLOR_DELTA, label=r'Hexagonal ($\delta$)', 
             markevery=100, markersize=0.1, markerfacecolor='white', markeredgewidth=0.1)
    
    # Vertical Line: Experimental
    ax1.axvline(x=T_EXPERIMENTAL, color=COLOR_EXP, linestyle='--', linewidth=1.5, alpha=0.8)

    ax1.set_ylabel(r'$G$ [\si{\electronvolt}/atom]')
    ax1.grid(True)
    ax1.legend(loc='best')

    # --- Plot 2: Delta G ---
    ax2.plot(T_cub, Delta_G, '-', color=COLOR_BLACK, lw=2, label=r'$\Delta G$')
    
    # Fill areas
    ax2.fill_between(T_cub, 0, Delta_G, where=(Delta_G >= 0), 
                     color=COLOR_CUBIC, alpha=0.15, label=r'Cubic Stable')
    ax2.fill_between(T_cub, 0, Delta_G, where=(Delta_G <= 0), 
                     color=COLOR_DELTA, alpha=0.15, label=r'Hexagonal Stable')
    
    # Zero line
    ax2.axhline(0, color='black', ls=':', lw=1)

    # Calculate reasonable Y position for text labels (bottom 20% of the plot range)
    y_limits = np.abs(Delta_G).max()
    text_y_pos = -y_limits * 0.8 

    # --- 1. EXPERIMENTAL LINE (Green) ---
    ax2.axvline(x=T_EXPERIMENTAL, color=COLOR_EXP, linestyle='--', linewidth=1.5, alpha=0.8)
    ax2.text(T_EXPERIMENTAL - 55, text_y_pos, r'$T_{\mathrm{exp}}$', 
             color=COLOR_EXP, fontsize=16, fontweight='bold', ha='left')

    # --- 2. CALCULATED LINE (Black/Grey) ---
    if Tc:
        print(f"\n TRANSITION TEMPERATURE FOUND: {Tc:.2f} K")
        
        # Draw Line on both plots
        ax1.axvline(x=Tc, color=COLOR_CALC, linestyle='-.', linewidth=1.5, alpha=0.8)
        ax2.axvline(x=Tc, color=COLOR_CALC, linestyle='-.', linewidth=1.5, alpha=0.8)
        
        # Add Dot at crossing
        ax2.plot(Tc, 0, 'o', color='white', markeredgecolor='black', markeredgewidth=1.5, markersize=6, zorder=10)
        
        # Add Label (Shifted left to distinguish from Experimental if they are close)
        ax2.text(Tc + 65, text_y_pos, r'$T_{\mathrm{calc}}$', 
                 color='k', fontsize=16, fontweight='bold', ha='right')
                 
    ax2.set_ylabel(r'$\Delta G$ [\si{\electronvolt}/atom]')
    ax2.set_xlabel(r'Temperature [\si{\kelvin}]')
    ax2.grid(True)
    
    # Save
    output_filename = "phase_diagram_plot.png"
    plt.savefig(output_filename, format='png', dpi=600, bbox_inches='tight')
    plt.savefig(output_filename.replace('.png', '.pdf'), format='pdf', dpi=600, bbox_inches='tight')
    print(f"📊 Plot saved to '{output_filename}' and PDF.")

if __name__ == "__main__":
    main()
