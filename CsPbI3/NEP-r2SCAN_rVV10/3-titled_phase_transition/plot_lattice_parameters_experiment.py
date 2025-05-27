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


# --- PRX Plot Style Settings (Combined and updated with user's preferences) ---
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Tahoma'] # User's preferred sans-serif font

matplotlib.rcParams.update({
    'font.size': 20,             # User's overall font size
    'text.usetex': True,         # Enable LaTeX rendering for text
    'axes.labelsize': 12,        # Font size for axis labels (matches overall font size)
    'xtick.labelsize': 14,       # Font size for x-tick labels (matches overall font size)
    'ytick.labelsize': 14,       # Font size for y-tick labels (matches overall font size)
    'legend.fontsize': 14,       # Legend font size (as per user's last request)
    'axes.linewidth': 0.8,       # Line width of the axes frame
    'xtick.direction': 'in',     # Inward ticks for x-axis
    'ytick.direction': 'in',     # Inward ticks for y-axis
    'xtick.major.size': 6,       # Major tick size for x-axis
    'ytick.major.size': 6,       # Major tick size for y-axis
    'xtick.minor.size': 2,       # Minor tick size for x-axis
    'ytick.minor.size': 2,       # Minor tick size for y-axis
    'xtick.top': True,           # Show ticks on top x-axis
    'ytick.right': True,         # Show ticks on right y-axis
    'legend.frameon': False,     # No frame around the legend (as per PRX style)
    'savefig.dpi': 300,          # Resolution for saved figures
    'figure.figsize': (6, 5),    # Standard figure size (5x5 inches as per user's last request)
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
    'axes.prop_cycle': plt.cycler( # Removed the default color cycle here to allow explicit colors below
        'color', ['#1f77b4', '#2ca02c', '#d62728'] # Keeping these if they are intended for the lines
    ),
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


def plot_thermo_out(filename="thermo.out"):
    """
    Reads data from a 'thermo.out' file, parses it, and plots column 1 against
    columns 10, 14, and 18 on the same chart, with enhanced visual styling
    featuring larger labels and ticks while maintaining a smaller figure size.

    Args:
        filename (str, optional): The name of the 'thermo.out' file.
            Defaults to "thermo.out".

    Returns:
        None. The function generates a plot and displays it.
    """
    try:
        with open(filename, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found. Please ensure the file exists and the path is correct.")
        return
    except Exception as e:
        print(f"An error occurred while reading the file: {e}")
        return

    data = []
    for line in lines:
        try:
            row_data = [float(x) for x in line.split()]
            data.append(row_data)
        except ValueError:
            print(f"Skipping line due to invalid data: {line.strip()}")
            continue
        except Exception as e:
            print(f"An unexpected error occurred while processing data: {e}")
            return

    if not data:
        print("Error: No valid data found in the file.")
        return

    data_array = np.array(data)
    num_columns = data_array.shape[1]
    if num_columns < 18:
        print(f"Error: The file '{filename}' does not contain enough columns. It has {num_columns} columns, but the script requires at least 18.")
        return

    x = data_array[:, 0]
    y1 = data_array[:, 9] / 16 / 1.415
    y2 = data_array[:, 13] / 16 / 1.415
    y3 = data_array[:, 17] / 12 / 2

    # Load experimental data
    try:
        exp_a_data = np.loadtxt('experimental_a_all.dat')
        x_exp_a = exp_a_data[:, 0]
        y_exp_a = exp_a_data[:, 1]
    except FileNotFoundError:
        print("Warning: 'experimental_a_all.dat' not found. Skipping experimental 'a' data.")
        x_exp_a, y_exp_a = np.array([]), np.array([])

    try:
        exp_b_data = np.loadtxt('experimental_b_all.dat')
        x_exp_b = exp_b_data[:, 0]
        y_exp_b = exp_b_data[:, 1]
    except FileNotFoundError:
        print("Warning: 'experimental_b_all.dat' not found. Skipping experimental 'b' data.")
        x_exp_b, y_exp_b = np.array([]), np.array([])

    try:
        exp_c_data = np.loadtxt('experimental_c_all.dat')
        x_exp_c = exp_c_data[:, 0]
        y_exp_c = exp_c_data[:, 1]
    except FileNotFoundError:
        print("Warning: 'experimental_c_all.dat' not found. Skipping experimental 'c' data.")
        x_exp_c, y_exp_c = np.array([]), np.array([])


    # Create figure and axes using the global rcParams for styling
    fig, ax = plt.subplots()

    # Plot the data
    # Using explicit colors as they were defined in the original script
    ax.plot(x[::2], y2[::2], label=r'a/$\sqrt{2}$', linestyle='-', color='#1f77b4', markersize=6)

    ax.scatter(x_exp_a, y_exp_a, label=r'a/$\sqrt{2}$ (Exp)',
               marker='s', s=50, facecolors='none', edgecolors='#1f77b4', linewidth=1.2) # Added linewidth for consistency

    ax.plot(x[::2], y1[::2], label=r'b/$\sqrt{2}$', linestyle='--', color='#2ca02c', markersize=6)

    ax.scatter(x_exp_b, y_exp_b, label=r'b/$\sqrt{2}$ (Exp)',
               marker='D', s=50, facecolors='none', edgecolors='#2ca02c', linewidth=1.2) # Added linewidth for consistency

    ax.plot(x[::2], y3[::2], label=r'c/2', linestyle=':', color='#d62728', markersize=6)

    ax.scatter(x_exp_c, y_exp_c, label=r'c/2 (Exp)',
               marker='h', s=50, facecolors='none', edgecolors='#d62728', linewidth=1.2) # Added linewidth for consistency

    # Set labels using LaTeX for proper formatting
    ax.set_xlabel(r'Temperature [K]') # Fontsize is set globally by rcParams
    ax.set_ylabel(r'Normalized lattice constants [\AA]') # Fontsize is set globally by rcParams

    # Set limits
    ax.set_xlim(left=100)
    ax.set_ylim(bottom=6.05)

    # Add legend
    ax.legend(loc='best', borderaxespad=0.) # Legend frame and fontsize are now set globally via rcParams

    # Grid is enabled via rcParams, but explicit call here is fine.
    ax.grid(True)
    plt.tight_layout() # Adjust layout to prevent labels from overlapping

    # Save the figure
    plt.savefig('latice_parameters.pdf', format='pdf', dpi=600)
    plt.savefig('latice_parameters.png', format='png', dpi=1200)
    # plt.show() # Commented out for non-interactive 'Agg' backend

if __name__ == "__main__":
    plot_thermo_out()

