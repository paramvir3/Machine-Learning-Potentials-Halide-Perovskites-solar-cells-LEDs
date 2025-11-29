import yaml
import numpy as np
import matplotlib.pyplot as plt
from pandas import read_csv, DataFrame
import sys

# Robust import for cumtrapz (handles SciPy v1.14+ deprecation)
try:
    from scipy.integrate import cumulative_trapezoid as cumtrapz
except ImportError:
    from scipy.integrate import cumtrapz

from ase.units import kB

def analyze_and_save_rs(spring_yaml="ti_spring.yaml", rs_csv="ti_rs.csv", output_csv="rs_free_energy_cubic.csv"):
    try:
        # 1. Load Reference State (G0, T0)
        # --------------------------------
        with open(spring_yaml, "r") as f:
            y = yaml.safe_load(f)

        # Handle keys that might vary slightly by GPUMD version
        T0 = y.get("T") or y.get("Temperature")
        G0 = y.get("G") or y.get("Free Energy (eV/atom)")

        if T0 is None or G0 is None:
            print(f"Error: Could not find 'T' or 'G' in {spring_yaml}")
            return

        print(f"Reference State: T0 = {T0} K, G0 = {G0} eV/atom")

        # 2. Load and Process RS Trajectory
        # ---------------------------------
        print(f"Reading {rs_csv}...")
        rs = read_csv(rs_csv)
        
        # --- FIX START: ROBUST SPLITTING LOGIC ---
        n = int(len(rs)/2)
        
        # Create initial slices
        forward = rs.iloc[:n].copy()
        backward = rs.iloc[n:].copy() 
        
        # Check lengths and trim the excess
        # This handles cases where len(rs) is odd (e.g. 101 rows)
        if len(backward) > len(forward):
            backward = backward.iloc[:len(forward)]
        elif len(forward) > len(backward):
            forward = forward.iloc[:len(backward)]
            
        # Now safely reverse the backward path
        backward = backward.iloc[::-1].copy()
        backward.reset_index(inplace=True, drop=True)
        # --- FIX END ---

        # Determine enthalpy column name
        if "enthalpy" in rs.columns:
            h_col = "enthalpy"
        elif "Potential Energy" in rs.columns:
            h_col = "Potential Energy"
        else:
            print(f"Error: Could not find enthalpy column. Columns: {rs.columns}")
            return

        # Extract values (Now guaranteed to be same shape)
        dl = forward["dlambda"].values
        l = forward["lambda"].values
        H1 = forward[h_col].values
        H2 = backward[h_col].values
        
        # Debug print to verify shapes
        print(f"Data shapes -> Lambda: {l.shape}, H1: {H1.shape}, H2: {H2.shape}")

        # Calculate Temperature array based on lambda scaling
        T = T0 / l

        # 3. Calculate Free Energy (Your Logic)
        # -------------------------------------
        # Calculate work (w) using cumulative trapezoidal rule
        w1 = cumtrapz(H1, l, initial=0)
        w2 = cumtrapz(H2, l, initial=0)
        w = (w1 + w2) * 0.5

        # Apply Reversible Scaling Free Energy Formula
        G = (G0 + 1.5 * kB * T0 * np.log(l) + w) / l

        # 4. Save Final Values to CSV
        # ---------------------------
        df_results = DataFrame({
            "Temperature_K": T,
            "Free_Energy_eV_atom": G,
            "Lambda": l,
            "Work_w": w
        })

        df_results.to_csv(output_csv, index=False)
        print(f"✅ Success! Final G and T values saved to: {output_csv}")

        # 5. Plotting
        # ------------------------------------
        plt.figure(figsize=(8, 6))
        plt.plot(T, G, label="RS Free Energy", linewidth=2)
        plt.xlabel("Temperature (K)")
        plt.ylabel("G (eV/atom)")
        plt.title(f"Free Energy vs Temperature (Ref: {T0}K)")
        plt.legend()
        plt.grid(True, alpha=0.3)
        plot_filename = "G_vs_T_plot.png"
        plt.savefig(plot_filename)
        print(f"Plot saved to: {plot_filename}")

    except FileNotFoundError as e:
        print(f"File Error: {e}")
    except ValueError as e:
        print(f"Math/Shape Error: {e}")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    analyze_and_save_rs()
