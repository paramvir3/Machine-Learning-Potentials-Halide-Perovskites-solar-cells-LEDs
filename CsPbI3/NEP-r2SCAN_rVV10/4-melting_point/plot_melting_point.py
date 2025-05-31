import string, re, struct, sys, math, os
import time
import types
#import datemax_time # Commented out as it's not a standard library and might cause import errors
#from pylab import * # Avoid using 'import *' in general practice
from sys import argv
from shutil import move
from os import remove, close
from subprocess import PIPE, Popen
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg') # Use 'Agg' backend for non-interactive plotting, suitable for saving figures
import matplotlib.pyplot as plt # Keeping plt as per the user's original script
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib import rc
from matplotlib import rcParams
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from scipy.optimize import curve_fit
# from ase.io import read, write # Commented out as not used in the provided snippet
from typing import Tuple
import scipy
import matplotlib.cm as cm # Import the colormap module


# --- PRX Plot Style Settings (Combined and updated with user's preferences) ---
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Tahoma'] # User's preferred sans-serif font

matplotlib.rcParams.update({
    'font.size': 12,             # User's overall font size
    'text.usetex': True,         # Enable LaTeX rendering for text
    'axes.labelsize': 12,        # Font size for axis labels (matches overall font size)
    'xtick.labelsize': 12,       # Font size for x-tick labels (matches overall font size)
    'ytick.labelsize': 12,       # Font size for y-tick labels (matches overall font size)
    'legend.fontsize': 12,       # Slightly smaller legend font size for clarity
    'axes.linewidth': 0.8,       # Line width of the axes frame
    'xtick.direction': 'in',     # Inward ticks for x-axis
    'ytick.direction': 'in',     # Inward ticks for y-axis
    'xtick.major.size': 4,       # Major tick size for x-axis
    'ytick.major.size': 4,       # Major tick size for y-axis
    'xtick.minor.size': 2,       # Minor tick size for x-axis
    'ytick.minor.size': 2,       # Minor tick size for y-axis
    'xtick.top': True,           # Show ticks on top x-axis
    'ytick.right': True,         # Show ticks on right y-axis
    'legend.frameon': False,     # No frame around the legend (as per PRX style)
    'savefig.dpi': 300,          # Resolution for saved figures
    'figure.figsize': (5, 5),    # Standard single-column width (3.37 inches) for PRX style
    'grid.linestyle': '--',      # User's grid linestyle
    'grid.alpha': 0.6,           # User's grid alpha
    'lines.linewidth': 2,        # User's line width
    'xtick.major.width': 1.2,    # User's major tick width
    'ytick.major.width': 1.2,    # User's major tick width
    'xtick.minor.visible': True, # User's minor tick visibility
    'ytick.minor.visible': True, # User's minor tick visibility
    'xtick.minor.width': 1,      # User's minor tick width
    'ytick.minor.width': 1,      # User's minor tick width
    'axes.spines.top': True,     # User's spine visibility
    'axes.spines.right': True,   # User's spine visibility
    'axes.spines.bottom': True,  # User's spine visibility
    'axes.spines.left': True,    # User's spine visibility
    'axes.edgecolor': '0.1',     # User's axes edge color
})

# LaTeX preamble for specific packages
matplotlib.rcParams['text.latex.preamble'] = '\n'.join([
    r'\usepackage{siunitx}',     # For upright micro symbols and unit formatting
    r'\sisetup{detect-all}',     # Forces siunitx to use your font settings
    r'\usepackage{amsmath, amssymb}', # For mathematical symbols and environments
])
# --- End PRX Plot Style Settings ---


# ----------------------------------------
#               COLORS
# ---------------------------------------
# Using the _Colors class as it's designed for indexing and provides a consistent palette.
class _Colors(object):
    """Helper class with different colors for plotting"""
    red = '#F15854'
    blue = '#5DA5DA'
    orange = '#FAA43A'
    green = '#60BD68'
    pink = '#F17CB0'
    brown = '#B2912F'
    purple = '#B276B2'
    yellow = '#DECF3F'
    gray = '#4D4D4D'
    cyan = '#00FFFF'
    rebecca_purple = '#663399'
    chartreuse = '#7FFF00'
    dark_red = '#8B0000'

    def __getitem__(self, i):
        color_list = [
            self.red,
            self.orange,
            self.green,
            self.blue,
            self.pink,
            self.brown,
            self.purple,
            self.yellow,
            self.gray,
            self.cyan,
            self.rebecca_purple,
            self.chartreuse,
            self.dark_red
        ]
        return color_list[i % len(color_list)]

from matplotlib.colors import LinearSegmentedColormap
Colors = _Colors() # Instantiate the Colors class
colormap_colors = ['#000000', Colors[0]] # Using Colors[0] (red) for colormap
colormap = LinearSegmentedColormap.from_list("Custom", colormap_colors, N=256)


# Function to read potential energy from thermo.out file
def read_potential_energy(filename):
    """
    Reads the potential energy (3rd column) from a thermo.out file.

    Args:
        filename (str): The name of the thermo.out file.

    Returns:
        list: A list of potential energy values, or an empty list if an error occurs.
    """
    potential_energies = []
    try:
        with open(filename, 'r') as f:
            for line in f:
                # Split the line into columns
                columns = line.split()
                # Check if the line has at least 3 columns
                if len(columns) >= 3:
                    try:
                        # Convert the 3rd column to float and append to the list
                        potential_energies.append(float(columns[2]))
                    except ValueError:
                        print(f"Error: Could not convert 3rd column to float in file {filename}. Skipping line: {line.strip()}")
    except FileNotFoundError:
        print(f"Error: File not found: {filename}")
        return []
    except Exception as e:
        print(f"An error occurred while reading the file {filename}: {e}")
        return []
    return potential_energies

# Function to plot potential energy vs time
def plot_potential_energy(data_dict, output_filename="potential_energy_plot.png"):
    """
    Plots the potential energy vs time for multiple temperatures on the same graph,
    using the 'tab20' colormap to represent increasing temperatures.

    Args:
        data_dict (dict): A dictionary where keys are temperatures (e.g., "750K")
                          and values are lists of potential energy values.
        output_filename (str, optional): The name of the file to save the plot.
                                         Defaults to "potential_energy_plot.png".
    """

    # Figure and axes creation, now using rcParams for figure size
    fig, ax = plt.subplots()

    # Set labels using LaTeX for proper formatting
    # Changed y-axis label to reflect the new scale
    ax.set_xlabel(r'Time [ns]', fontsize=22)
    ax.set_ylabel(r'Potential energy [eV $\times 10^4$]', fontsize=22) # Updated label for clarity
    ax.set_title(r'Experiment $\sim$ 753K')

    # Sort the temperatures to ensure the colormap aligns with increasing temperature
    sorted_temps = sorted(data_dict.keys(), key=lambda t: int(t[:-1])) # Sort numerically, removing the 'K'

    # Choose the 'tab20' colormap
    cmap = cm.get_cmap('tab20')

    # Normalize the temperature values to the range of the colormap
    norm = plt.Normalize(vmin=float(int(sorted_temps[0][:-1])), vmax=float(int(sorted_temps[-1][:-1])))

    # Store plot handles and labels for reverse legend
    lines = []
    labels = []

    # Plot data for each temperature
    for temp in sorted_temps:
        energies = data_dict[temp]
        if energies:  # Only plot if there is data.
            # Divide energies by 10000 as requested
            energies_scaled = [e / 10000 for e in energies]
            time_ns = [t * 0.5 / 1000 for t in range(len(energies_scaled))]  # Convert time to ns
            color = cmap(norm(int(temp[:-1])) + 0.3)
            line, = ax.plot(time_ns[::100], energies_scaled[::100], label=temp, color=color)
            lines.append(line)
            labels.append(temp)
        else:
            print(f"Skipping {temp} due to missing or invalid data.")

#     # --- Add the extra legend entry for T_experiment = 750K without drawing a line ---
#     # Create an invisible line to serve as a proxy for the legend entry
#     exp_line, = ax.plot([], [],
#                         linestyle='--', color=Colors.dark_red, linewidth=1.5,
#                         label=r'T$_{\text{experiment}}$ = 750K') # Using LaTeX for subscript
#     lines.append(exp_line)
#     labels.append(r'Experiment $\sim$ 753K')
    # --- End of extra legend entry ---

    ax.set_xlim(xmin=2, xmax=50)
    # Updated y-axis limits to reflect the division by 10000
    ax.set_ylim(ymax=-16.2240, ymin=-16.2800)

    # Display legend in reverse order, using bbox_to_anchor for external placement
    # The new experimental line will be at the bottom of the reversed legend
    ax.legend(reversed(lines), reversed(labels), loc='center left', bbox_to_anchor=(1.0, 0.5))
    ax.grid(True)
    plt.tight_layout(rect=[0, 0, 1.02, 1.02]) # Adjust tight_layout to make space for the legend

    plt.savefig('Melting_point.pdf', format='pdf', dpi=600)
    plt.savefig('Melting_point.png', format='png', dpi=1200)
    # plt.show() # Commented out for non-interactive 'Agg' backend

def main():
    """
    Main function to orchestrate the reading and plotting of potential energy data.
    """
    # Define the directories to search for thermo.out files
    directories = ["740K", "750K", "760K", "770K", "775K", "780K", "790K", "800K", "810K", "820K"]

    # Dictionary to store potential energy data for each temperature
    potential_energy_data = {}

    # Loop through each directory, find thermo.out, and read the data
    for temp_dir in directories:
        thermo_file = os.path.join(temp_dir, "thermo.out")
        if os.path.exists(thermo_file):  # Check if the file exists
            energies = read_potential_energy(thermo_file)
            potential_energy_data[temp_dir] = energies
        else:
            print(f"Warning: thermo.out file not found in directory {temp_dir}. Skipping.")

    # Plot the potential energy data
    if potential_energy_data:
        plot_potential_energy(potential_energy_data)
    else:
        print("No data to plot. Exiting.")

if __name__ == "__main__":
    main()

