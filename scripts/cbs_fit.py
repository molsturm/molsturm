from math import exp, sqrt
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# The function to fit against is in all cases
# See
# Frank Jensen
# Estimating the Hartree-Fock limit from finite basis set calculations

def EnergyFunc(x, Einf, A, B):
    return Einf + A* np.exp(-B* np.sqrt(x))

def energy_fit(orders,energies):
    # Make fitting easier by shifting the energies beforehand:
    eshift = int(min(energies))
    energies = [ e - eshift for e in energies ]

    print("Fitting curve E(order) = E(CBS) + A * exp(-B * sqrt(order))")
    popt, pcov = curve_fit(EnergyFunc,orders,energies)

    # Undo energy shift:
    popt[0] += eshift

    print("Estimated parameters:  value (stddev)")
    print("     E(CBS) == {0:16.14G} ({1: >9.4G})".format(popt[0],sqrt(pcov[0][0])))
    print("          A == {0:16.14G} ({1: >9.4G})".format(popt[1],sqrt(pcov[1][1])))
    print("          B == {0:16.14G} ({1: >9.4G})".format(popt[2],sqrt(pcov[2][2])))
    return popt

def cbs_fit_gaussian(orders, energies):
    popt = energy_fit(orders, energies)
    plt.figure()
    plt.plot(orders, energies, 'ko', label="cc-pVnZ")
    x = np.linspace(min(orders)-0.1,8,200)
    plt.plot(x, EnergyFunc(x, *popt), 'r-', label="Fitted convergence")
    plt.legend()
    plt.show()
