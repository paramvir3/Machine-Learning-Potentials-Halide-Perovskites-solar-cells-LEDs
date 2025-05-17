import string, re, struct, sys, math, os
import time
import types
#import datemax_time
#from pylab import *
from sys import argv
from shutil import move
from os import remove, close
from subprocess import PIPE, Popen
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pl
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib import rc
from matplotlib import rcParams
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from scipy.optimize import curve_fit
from matplotlib import rc
from matplotlib import rcParams
from ase.io import read, write
from typing import Tuple
import scipy


rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Tahoma']

matplotlib.rcParams.update({'font.size': 12})
matplotlib.rcParams['text.latex.preamble'] = '\n'.join([
       r'\usepackage{siunitx}',   # i need upright \micro symbols, but you need...
       r'\sisetup{detect-all}',   # ...this to force siunitx to actually use your fonts
       r'\usepackage{amsmath, amssymb}',
])


# ----------------------------------------
#              COLORS
# ---------------------------------------

c_orange=np.array([172./255., 90./255., 22./255.,1.0])
c_red_nature=np.array([206./255., 30./255., 65./255.,1.0])
c_green_nature=np.array([96./255., 172./255., 63./255.,1.0])
c_blue_nature =np.array([54./255., 79./255., 156./255.,1.0])
c_purple_nature=np.array([245./255., 128./255., 32./255.,1.0])#np.array([192./255., 98./255., 166./255.,1.0])
c_black1=np.array([50./255., 50./255., 50./255.,1.0])
c_black2=np.array([100./255., 100./255.,100./255.,1.0])
c_black3=np.array([150./255., 150./255., 150./255.,1.0])
c_black4=np.array([200./255., 200./255., 200./255.,1.0])
c_cyna='#00aeff' #42cef4'
c_dark_green='#004d00'
c_pink='#ff1aff'
c_purple='#990099'
c_dark_orange='#db5f00'
c_orange='orange'


c_black=np.array([40./255., 41./255., 35./255.,1.0])
c_black=np.array([116./255., 112./255., 93./255.,1.0])

c_red=np.array([249./255., 36./255.,114./255.,1.0])
c_purple=np.array([172./255., 128./255., 255./255.,1.0])
c_orange=np.array([253./255., 150./255., 33./255.,1.0])


c_red=np.array([206./255., 30./255., 65./255.,1.0])

#c_green=np.array([62./255., 208./255., 102./255.,1.0])
c_green=np.array([96./255., 172./255., 63./255.,1.0])
c_blue=np.array([26./255., 97./255., 191./255.,1.0])
c_pink=np.array([144./255., 110./255., 209./255.,1.0])#'#AE5BB3'  #F379FB'
c_orange=np.array([245./255., 128./255., 32./255.,1.0])
c_cyan=np.array([30./255., 165./255., 180./255.,1.0])
c_black=np.array([0./255., 0./255., 0./255.,1.0])
c_green=np.array([96./255., 172./255., 63./255.,1.0])
c_blue=np.array([26./255., 97./255., 191./255.,1.0])

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
Colors = _Colors()
colormap_colors = ['#000000', Colors[0]]
colormap = LinearSegmentedColormap.from_list("Custom", colormap_colors, N=256)

import re

def extract_epot_value(file_path):
    """
    Reads a file and extracts the energy value after 'Epot after opt:' if present.

    Parameters:
        file_path (str): Path to the file.

    Returns:
        float: The extracted energy value, or None if not found.
    """
    # Define the pattern to match: 'Epot after opt:' followed by a number
    pattern = r"Epot after opt:\s*([-+]?\d*\.\d+|\d+)"
    with open(file_path, 'r') as file:
        for line in file:
            # Search for the pattern in each line
            match = re.search(pattern, line)
            if match:
                # Convert the matched number to a float and return it
                return float(match.group(1))
    
    # Return None if the pattern is not found
    return None

# Assume these functions have already extracted the energy values from different paths
delta_phase = extract_epot_value("delta/output")
ortho_phase = extract_epot_value("ortho/output")
beta_phase = extract_epot_value("beta/output")
cubic_phase = extract_epot_value("cubic/output") 
delta_FAPI = extract_epot_value("hexagonal_FAPI/output") 

# Define the phases and their corresponding energy values

phases = [r'$\delta$-Phase', r'$\gamma$-Phase', r'$\beta$-Phase', r'$\alpha$-Phase', r'$\delta$-FAPI']

energy_values = [delta_phase, ortho_phase, beta_phase, cubic_phase, delta_FAPI]

# Conversion factor (if necessary, assuming it converts from eV to another unit, like kJ/mol)
ev_to_kjmol = 96.485332 / 8

# Calculate relative values with respect to delta_phase, applying the conversion factor consistently
relative_energy_values = [(value * ev_to_kjmol ) - (delta_phase * ev_to_kjmol) for value in energy_values]

# Create a figure with two subplots
fig, (ax2) = pl.subplots(1, 1, figsize=(6, 6))


# Define font size for all labels and titles
font_size = 16

# Plot 2: Relative energy values to delta phase
ax2.bar(phases, relative_energy_values, color=['skyblue', 'salmon', 'lightgreen', 'orange', 'blue'])
ax2.set_xlabel('CsPbI$_3$ polymorphs', fontsize=font_size)
ax2.set_ylabel('Relative Energy (KJ/mol)', fontsize=font_size)
ax2.set_title('Relative Energy Values with respect to $\delta$-CsPbI$_3$', fontsize=font_size)
ax2.tick_params(axis='x', rotation=45, labelsize=font_size - 2)
ax2.tick_params(axis='y', labelsize=font_size - 2)
ax2.grid(axis='y', linestyle='--', alpha=0.7)

# Adjust layout and display
pl.tight_layout()
pl.show()

fig.savefig("rE"+".pdf", dpi=300, bbox_inches="tight",transparent=True)
fig.savefig("rE"+".png", dpi=300, bbox_inches="tight",transparent=False)

