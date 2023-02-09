import ast
import sys
from math import sqrt, exp
from scipy.constants import pi, N_A, k, h, physical_constants

k_eV = physical_constants["Boltzmann constant in eV/K"][0]


def get_q_vib(T, vib_list, include_zpe=False):
    """Calculates the vibrational partition function (including ZPE).

    Arguments:
        T (float): The temperature in K
        vib_list (str): List of all vibrational modes in meV
        include_zpe (bool): Include ZPE contribution in the partition function
    """
    q_vib = 1.0
    vib_list = ast.literal_eval(vib_list)
    for v in vib_list:
        if include_zpe:
            q_vib = q_vib * exp(v / (1000 * 2 * k_eV * T)) / (1 - exp(v / (1000 * k_eV * T)))
        else:
            q_vib = q_vib * 1 / (1 - exp(v / (1000 * k_eV * T)))
    return q_vib


def get_q_rot(T, rot_temperatures, sym_number):
    """Calculates the rotational partition function.

    Arguments: T (float): The temperature in K
    rot_temperatures (str): List of rotational temperatures (1 for linear,
    3 for non-linear). They can be calculated from the inertia moments as rot_t = h**2 / (8 * pi**2 * i * k).
    sym_number (int): Symmetry number of the molecule
    """
    rot_temperatures = ast.literal_eval(rot_temperatures)
    if len(rot_temperatures) == 1:  # linear
        rot_t = rot_temperatures[0]
        q_rot_gas = T / (sym_number * rot_t)
    elif len(rot_temperatures) == 3:  # non-linear
        rot_t_a, rot_t_b, rot_t_c = rot_temperatures[0], rot_temperatures[1], rot_temperatures[2]
        q_rot_gas = (sqrt(pi) / sym_number) * sqrt(T**3 / (rot_t_a * rot_t_b * rot_t_c))
    else:
        sys.exit(f"Invalid inertia_list")
    return q_rot_gas


def calc_non_act_ads(A_site, molec_mass, T, vib_list_ads, vib_list_gas, rot_temperatures, sym_number):
    """Calculates the forward and reverse pre-exponential factors for a reversible non-activated adsorption."""
    m = molec_mass / 1000 / N_A
    q_vib_ads = get_q_vib(T=T, vib_list=vib_list_ads)
    q_vib_gas = get_q_vib(T=T, vib_list=vib_list_gas)
    q_rot_gas = get_q_rot(T=T, rot_temperatures=rot_temperatures, sym_number=sym_number)
    q_trans_2d_gas = A_site * 2 * pi * m * k * T / h ** 2
    pe_fwd = A_site / sqrt(2 * pi * m * k * T)
    pe_rev = (q_vib_gas * q_rot_gas * q_trans_2d_gas / q_vib_ads) * (k * T / h)
    return pe_fwd, pe_rev


def calc_act_ads(A_site, molec_mass, T, vib_list_ads, vib_list_gas, vib_list_ts, rot_temperatures, sym_number):
    """Calculates the forward and reverse pre-exponential factors for a reversible activated adsorption."""
    m = molec_mass / 1000 / N_A
    q_vib_ads = get_q_vib(T=T, vib_list=vib_list_ads)
    q_vib_gas = get_q_vib(T=T, vib_list=vib_list_gas)
    q_vib_ts = get_q_vib(T=T, vib_list=vib_list_ts)
    q_rot_gas = get_q_rot(T=T, rot_temperatures=rot_temperatures, sym_number=sym_number)
    q_trans_2d_gas = A_site * 2 * pi * m * k * T / h ** 2
    pe_fwd = (q_vib_ts / (q_vib_gas * q_rot_gas * q_trans_2d_gas)) * (A_site / sqrt(2 * pi * m * k * T))
    pe_rev = (q_vib_ts / q_vib_ads) * (k * T / h)
    return pe_fwd, pe_rev


def calc_surf_react(T, vib_list_initial, vib_list_ts, vib_list_final):
    """Calculates the forward and reverse pre-exponential factors for a reversible surface reaction."""
    q_vib_initial = get_q_vib(T=T, vib_list=vib_list_initial)
    q_vib_ts = get_q_vib(T=T, vib_list=vib_list_ts)
    q_vib_final = get_q_vib(T=T, vib_list=vib_list_final)
    pe_fwd = (q_vib_ts / q_vib_initial) * (k * T / h)
    pe_rev = (q_vib_ts / q_vib_final) * (k * T / h)
    return pe_fwd, pe_rev
