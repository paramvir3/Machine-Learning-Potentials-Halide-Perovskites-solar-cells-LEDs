import yaml
import numpy as np
import matplotlib.pyplot as plt
from pandas import read_csv
from scipy.integrate import cumtrapz
from ase.units import kB

# you need to run a ti_spring simulation at t_min
with open("ti_spring.yaml", "r") as f:
        y =  yaml.safe_load(f)

        T0 = y["T"]
        G0 = y["G"]

        rs = read_csv("ti_rs.csv")
        n = int(len(rs)/2)
        forward = rs[:n]
        backward = rs[n:][::-1]
        backward.reset_index(inplace=True)
        dl = forward["dlambda"]
        l = forward["lambda"]
        H1 = forward["enthalpy"]
        H2 = backward["enthalpy"]
        T = T0/l

        w = (cumtrapz(H1,l,initial=0) + cumtrapz(H2,l,initial=0))*0.5

        G = (G0 + 1.5*kB*T0*np.log(l) + w)/l
        plt.plot(T, G, label="RS")
        plt.legend()
        plt.xlabel("T (K)")
        plt.ylabel("G (eV/atom)")
