import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter

def calculate_heat_capacity_smooth_derivative(filename="thermo.out", window_length=50, polyorder=2):
    """
    Calculates the heat capacity (Cp) using a smoothed derivative method (Savitzky-Golay filter)
    and plots it against temperature.

    Args:
        filename (str, optional): The name of the thermo.output file. Defaults to "thermo.out".
        window_length (int, optional): The length of the filter window (must be odd). Defaults to 101.
        polyorder (int, optional): The order of the polynomial used to fit the samples. Defaults to 2.

    Returns:
        None. Displays a plot of the heat capacity vs. temperature.
    """
    try:
        # Load data from file
        data = np.loadtxt(filename)
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
        return
    except Exception as e:
        print(f"Error loading file: {e}")
        return

    # Extract temperature (1st column) and potential energy (3rd column)
    temperature = data[:, 0]
    potential_energy = data[:, 2]

    # Ensure there are enough data points for the Savitzky-Golay filter
    if len(temperature) < window_length:
        print(f"Error: Not enough data points for the given window length ({window_length}).")
        return

    # Calculate the derivative of potential energy with respect to temperature using Savitzky-Golay filter
    try:
        cp = savgol_filter(potential_energy, window_length, polyorder, deriv=1, delta=np.diff(temperature).mean())
    except ValueError as e:
        print(f"Error during Savitzky-Golay filtering: {e}")
        return

   # Enhanced Plot Styling with Smaller Figure and Larger Labels/Ticks
    plt.style.use('seaborn-v0_8-whitegrid')

    plt.rcParams.update({
        'font.family': 'Arial',
        'font.size': 10,  # Base font size
        'axes.labelsize': 20,  # Larger axis labels
        'axes.titlesize': 14,
        'xtick.labelsize': 14,  # Larger tick labels
        'ytick.labelsize': 14,  # Larger tick labels
        'legend.fontsize': 18,
        'figure.figsize': (7, 5),  # Smaller figure size
        'axes.linewidth': 1,
        'grid.linestyle': '--',
        'grid.alpha': 0.6,
        'lines.linewidth': 2,
        'legend.frameon': True,
        'legend.edgecolor': '0.1',
        'legend.shadow': False,
        'legend.loc': 'best',
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
        'xtick.major.size': 8,
        'ytick.major.size': 8,
        'xtick.major.width': 1.2,
        'ytick.major.width': 1.2,
        'xtick.minor.visible': True,
        'ytick.minor.visible': True,
        'xtick.minor.size': 4,
        'ytick.minor.width': 1,
        'ytick.minor.width': 1,
        'axes.spines.top': True,
        'axes.spines.right': True,
        'axes.spines.bottom': True,
        'axes.spines.left': True,
        'axes.edgecolor': '0.1',
        'axes.prop_cycle': plt.cycler(
            'color', ['#1f77b4', '#2ca02c', '#d62728']
        ),
    })

    fig, ax = plt.subplots()

    ax.plot(temperature, cp, label="dU/dT", color='tomato', linestyle='-',  markersize=5)
   

    ax.set_xlabel('Temperature [K]',  fontsize=20)
    ax.set_ylabel('Heat capacity', fontsize=20)
    # ax.set_title('Temperature Dependence of Lattice Parameters', fontweight='bold', fontsize=14)

    ax.set_xlim([min(temperature), max(temperature)])
    #ax.set_ylim([np.min(cp) - 1, np.max(cp) + 1])

    ax.legend()
    ax.grid(True)
    plt.tight_layout()


    plt.savefig('Cp_plot.pdf', format='pdf', dpi=600)
    plt.savefig('Cp_plot.png', format='png', dpi=1200)
    plt.show()


if __name__ == "__main__":
    # Call the function with appropriate arguments
    calculate_heat_capacity_smooth_derivative(window_length=101, polyorder=2)

