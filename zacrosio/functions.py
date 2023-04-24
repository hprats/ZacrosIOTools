import ast
import sys
import numpy as np
import matplotlib.pyplot as plt
from math import sqrt, exp
from scipy.constants import pi, N_A, k, h, physical_constants

k_eV = physical_constants["Boltzmann constant in eV/K"][0]
atomic_mass = physical_constants["atomic mass constant"][0]


def get_q_vib(T, vib_list):
    """Calculates the vibrational partition function (including ZPE).

    Arguments:
        T (float): The temperature in K
        vib_list (str): List of all vibrational modes in meV
    """
    q_vib = 1.0
    vib_list = ast.literal_eval(vib_list)
    for v in vib_list:
        q_vib = q_vib * exp(- v / (1000 * 2 * k_eV * T)) / (1 - exp(- v / (1000 * k_eV * T)))
    return q_vib


def get_q_rot(T, inertia_list, sym_number):
    """Calculates the rotational partition function.

    Arguments: T (float): The temperature in K
    inertia_moments (str): List of inertia moments (1 for linear, 3 for non-linear) in amu*Å2
    sym_number (int): Symmetry number of the molecule
    """
    inertia_list = ast.literal_eval(inertia_list)
    if len(inertia_list) == 1:  # linear
        i = inertia_list[0] * atomic_mass / 1.0e20  # from amu*Å2 to kg*m2
        q_rot_gas = 8 * pi ** 2 * i * k * T / (sym_number * h ** 2)
    elif len(inertia_list) == 3:  # non-linear
        i_a = inertia_list[0] * atomic_mass / 1.0e20
        i_b = inertia_list[1] * atomic_mass / 1.0e20
        i_c = inertia_list[2] * atomic_mass / 1.0e20
        q_rot_gas = (sqrt(pi * i_a * i_b * i_c) / sym_number) * (8 * pi ** 2 * k * T / h ** 2) ** (3 / 2)
    else:
        sys.exit(f"Invalid inertia_list")
    return q_rot_gas


def calc_ads(A_site, molec_mass, T, vib_list_is, vib_list_ts, vib_list_fs, inertia_list, sym_number, degeneracy):
    """Calculates the forward and reverse pre-exponential factors for a reversible activated adsorption."""
    A_site = A_site * 1.0e-20  # Å^2 to m^2
    m = molec_mass * 1.0e-3 / N_A  # g/mol to kg/molec
    q_vib_gas = get_q_vib(T=T, vib_list=vib_list_is)
    q_rot_gas = get_q_rot(T=T, inertia_list=inertia_list, sym_number=sym_number)
    q_trans_2d_gas = A_site * 2 * pi * m * k * T / h ** 2
    q_el_gas = degeneracy
    q_vib_ads = get_q_vib(T=T, vib_list=vib_list_fs)
    if not vib_list_ts:  # non-activated
        pe_fwd = A_site / sqrt(2 * pi * m * k * T)
        pe_rev = (q_el_gas * q_vib_gas * q_rot_gas * q_trans_2d_gas / q_vib_ads) * (k * T / h)
    else:  # activated
        q_vib_ts = get_q_vib(T=T, vib_list=vib_list_ts)
        pe_fwd = (q_vib_ts / (q_el_gas * q_vib_gas * q_rot_gas * q_trans_2d_gas)) * (A_site / sqrt(2 * pi * m * k * T))
        pe_rev = (q_vib_ts / q_vib_ads) * (k * T / h)
    return pe_fwd, pe_rev


def calc_surf_proc(T, vib_list_is, vib_list_ts, vib_list_fs):
    """Calculates the forward and reverse pre-exponential factors for a reversible surface process."""
    q_vib_initial = get_q_vib(T=T, vib_list=vib_list_is)
    q_vib_ts = get_q_vib(T=T, vib_list=vib_list_ts)
    q_vib_final = get_q_vib(T=T, vib_list=vib_list_fs)
    pe_fwd = (q_vib_ts / q_vib_initial) * (k * T / h)
    pe_rev = (q_vib_ts / q_vib_final) * (k * T / h)
    return pe_fwd, pe_rev


def plot_coverage(path):
    f = open(f"{path}/specnum_output.txt")
    header = f.readline().split()
    f.close()
    data = np.loadtxt("{path}/specnum_output.txt", skiprows=1)
    for i in range(5, len(header)):
        plt.plot(xdata=data[:, 2], ydata=data[:, i], label=header[i])
    plt.xlabel(header[2])
    plt.ylabel("Number of species")
    plt.legend()
    plt.show()
