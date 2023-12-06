import sys
import numpy as np
import matplotlib.pyplot as plt


def plt_coverage(path):
    with open(f"{path}/specnum_output.txt", "r") as infile:
        header = infile.readline().split()
    data = np.loadtxt(f"{path}/specnum_output.txt", skiprows=1)
    for i in range(5, len(header)):
        plt.plot(data[:, 2], data[:, i], label=header[i])
    plt.xlabel(header[2])
    plt.ylabel("Number of species")
    plt.legend()
    plt.show()


def plt_production(path, n_surf_species, molecule):
    with open(f"{path}/specnum_output.txt", "r") as infile:
        header = infile.readline().split()
    data = np.loadtxt(f"{path}/specnum_output.txt", skiprows=1)
    if molecule is None:
        for i in range(5 + n_surf_species, len(header)):
            if data[-1, i] > 0:
                plt.plot(data[:, 2], data[:, i], label=header[i])
    else:
        if molecule not in header:
            sys.exit(f"ERROR: {molecule} not found")
        plt.plot(data[:, 2], data[:, header.index(molecule)], label=header[header.index(molecule)])
    plt.xlabel("KMC time (s)")
    plt.ylabel("Number of species")
    plt.legend()
    plt.show()


def plt_tof(path, area, molecule):
    with open(f"{path}/specnum_output.txt", "r") as infile:
        header = infile.readline().split()
    if molecule not in header:
        sys.exit(f"ERROR: {molecule} not found")
    data = np.loadtxt(f"{path}/specnum_output.txt", skiprows=1)
    kmc_time = data[-1, 2]
    tof = np.diff(data[:, header.index(molecule)] / area, data[:, 2], 0.1)
    plt.plot(data[:, 2], tof, label=header[header.index(molecule)])
    #tof = dxdt(x=data[:, header.index(molecule)] / area, t=data[:, 2], kind="finite_difference", k=10)
    plt.plot(data[:, 2], tof, label=header[header.index(molecule)])
    plt.xlabel("KMC time (s)")
    plt.ylabel("TOF (molec·s-1·Å-2)")
    plt.legend()
    plt.show()
