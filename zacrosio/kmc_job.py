import os
import ast
import sys
from scipy.constants import pi, N_A, k, h, physical_constants
from math import sqrt, exp

k_eV = physical_constants["Boltzmann constant in eV/K"][0]


def get_q_vib(T, vib_list):
    q_vib = 1
    for v in vib_list:
        q_vib = q_vib * exp(-h * v / (2 * k * T)) / (1 - exp(-h * v / (k * T)))
    return q_vib


def get_q_rot(T, inertia_list, sym_number):
    if len(inertia_list) == 1:  # linear
        i = inertia_list[0]
        q_rot_gas = 8 * pi ** 2 * i * k * T / (sym_number * h ** 2)
    elif len(inertia_list) == 3:  # non-linear
        i_a, i_b, i_c = inertia_list[0], inertia_list[1], inertia_list[2]
        q_rot_gas = (sqrt(pi * i_a * i_b * i_c) / sym_number) * (8 * pi ** 2 * k * T / h ** 2) ** (3 / 2)
    else:
        sys.exit(f"Invalid inertia_list")
    return q_rot_gas


def calc_non_act_ads(A_site, molec_mass, T, vib_list_ads, vib_list_gas, inertia_list, sym_number):
    m = molec_mass / 1000 / N_A
    q_vib_ads = get_q_vib(T=T, vib_list=vib_list_ads)
    q_vib_gas = get_q_vib(T=T, vib_list=vib_list_gas)
    q_rot_gas = get_q_rot(T=T, inertia_list=inertia_list, sym_number=sym_number)
    q_trans_2d_gas = A_site * 2 * pi * m * k * T / h ** 2
    pe_fwd = A_site / sqrt(2 * pi * m * k * T)
    pe_rev = (q_vib_gas * q_rot_gas * q_trans_2d_gas / q_vib_ads) * (k * T / h)
    return pe_fwd, pe_rev


def calc_act_ads(A_site, molec_mass, T, vib_list_ads, vib_list_gas, vib_list_ts, inertia_list, sym_number):
    m = molec_mass / 1000 / N_A
    q_vib_ads = get_q_vib(T=T, vib_list=vib_list_ads)
    q_vib_gas = get_q_vib(T=T, vib_list=vib_list_gas)
    q_vib_ts = get_q_vib(T=T, vib_list=vib_list_ts)
    q_rot_gas = get_q_rot(T=T, inertia_list=inertia_list, sym_number=sym_number)
    q_trans_2d_gas = A_site * 2 * pi * m * k * T / h ** 2
    pe_fwd = (q_vib_ts / (q_vib_gas * q_rot_gas * q_trans_2d_gas)) * (A_site / sqrt(2 * pi * m * k * T))
    pe_rev = (q_vib_ts / q_vib_ads) * (k * T / h)
    return pe_fwd, pe_rev


def calc_surf_react(T, vib_list_initial, vib_list_ts, vib_list_final):
    q_vib_initial = get_q_vib(T=T, vib_list=vib_list_initial)
    q_vib_ts = get_q_vib(T=T, vib_list=vib_list_ts)
    q_vib_final = get_q_vib(T=T, vib_list=vib_list_final)
    pe_fwd = (q_vib_ts / q_vib_initial) * (k * T / h)
    pe_rev = (q_vib_ts / q_vib_final) * (k * T / h)
    return pe_fwd, pe_rev


class NewKMCJob:
    """A class that represents a new KMC job with ZACROS."""
    def __init__(self, path, simulation_tags, df_mechanism, df_energetics, lattice_path):
        self.path = path
        self.name = path.split('/')[-1]
        self.simulation_tags = simulation_tags
        self.df_mechanism = df_mechanism
        self.df_energetics = df_energetics
        self.lattice_path = lattice_path

    def create_job_dir(self, T=300):
        """Creates a new directory and writes there the ZACROS input files"""
        if not os.path.exists(self.path):
            os.mkdir(self.path)
            self.write_simulation(T=T)
            self.write_mechanism(T=T)
            self.write_energetics()
            self.write_lattice()
        else:
            print(f'{self.path} already exists (nothing done)')

    def write_simulation(self, T=300):
        """Writes the simulation_input.dat file"""
        self.simulation_tags['temperature'] = T
        with open(f"{self.path}/simulation_input.dat", 'w') as infile:
            for tag in self.simulation_tags:
                infile.write(f"{tag} = {self.simulation_tags[tag]}\n")
            infile.write(f"no_restart\n")
            infile.write(f"finish\n")

    def write_mechanism(self, T=300):
        """Writes the mechanism_input.dat file"""
        self.write_header(file_name="mechanism_input.dat")
        with open(f"{self.path}/mechanism_input.dat", 'a') as infile:
            infile.write('\nmechanism\n\n')
            infile.write('############################################################################s\n\n')
            for step in self.df_mechanism.index:
                infile.write(f"reversible_step {step}\n\n")
                step_type = self.df_mechanism.loc[step, 'type']
                if 'adsorption' in step_type:
                    infile.write(f"  gas_reacs_prods {self.df_mechanism.loc[step, 'gas_reacs_prods']}\n")
                sites = self.df_mechanism.loc[step, 'sites']
                infile.write(f"  sites {sites}\n")
                if sites > 1:
                    infile.write(f"  neighboring {self.df_mechanism.loc[step, 'neighboring']}\n")
                infile.write(f"  initial\n")
                initial_state_list = ast.literal_eval(self.df_mechanism.loc[step, 'initial'])
                for element in initial_state_list:
                    infile.write(f"    {element}\n")
                infile.write(f"  final\n")
                final_state_list = ast.literal_eval(self.df_mechanism.loc[step, 'final'])
                for element in final_state_list:
                    infile.write(f"    {element}\n")
                pre_expon, pe_ratio = self.get_pre_expon(step=step, T=T)
                infile.write(f"\n  pre_expon {pre_expon:.3e}\n")
                infile.write(f"  activ_eng {self.df_mechanism.loc[step, 'activ_eng']:.2f}\n")
                infile.write(f"  pe_ratio {pe_ratio:.3e}\n\n")
                infile.write(f"end_reversible_step\n\n")
                infile.write('############################################################################s\n\n')
                infile.write(f"end_mechanism\n")

    def write_energetics(self):
        """Writes the energetics_input.dat file"""
        self.write_header(file_name="energetics_input.dat")
        with open(f"{self.path}/energetics_input.dat", 'a') as infile:
            infile.write('\nenergetics\n\n')
            infile.write('############################################################################s\n\n')
            for cluster in self.df_energetics.index:
                infile.write(f"cluster {cluster}\n\n")
                sites = self.df_energetics.loc[cluster, 'sites']
                infile.write(f"  sites {sites}\n")
                if sites > 1:
                    infile.write(f"  neighboring {self.df_energetics.loc[cluster, 'neighboring']}\n")
                infile.write(f"  lattice_state\n")
                lattice_state_list = ast.literal_eval(self.df_energetics.loc[cluster, 'lattice_state'])
                for element in lattice_state_list:
                    infile.write(f"    {element}\n")
                infile.write(f"\n  cluster_eng {self.df_energetics.loc[cluster, 'cluster_eng']:.2f}\n\n")
                infile.write(f"end_cluster\n\n")
                infile.write('############################################################################s\n\n')
                infile.write(f"end_energetics\n")

    def write_lattice(self):
        """Writes the lattice_input.dat file"""
        os.system(f"cp {self.lattice_path} {self.path}/lattice_input.dat")

    def write_header(self, file_name):
        with open(f"{self.path}/{file_name}", 'w') as infile:
            infile.write('############################################################################\n')
            infile.write('# Zacros Input File generated with the ZacrosIOTools                       #\n')
            infile.write('# https://github.com/hprats/ZacrosIOTools.git                              #\n')
            infile.write('#                                                                          #\n')
            infile.write('# Hector Prats, PhD                                                        #\n')
            infile.write('############################################################################\n')

    def get_pre_expon(self, step, T):
        """Calculates the forward pre-exponential and the pre-exponential ratio, required for the mechanism_input.dat
        file """
        step_type = self.df_mechanism.loc[step, 'type']
        if step_type == 'non_activated_adsorption':
            pe_fwd, pe_rev = calc_non_act_ads(A_site=self.df_mechanism.loc[step, 'A_site'],
                                              molec_mass=self.df_mechanism.loc[step, 'molec_mass'],
                                              T=T,
                                              vib_list_ads=self.df_mechanism.loc[step, 'vib_list_ads'],
                                              vib_list_gas=self.df_mechanism.loc[step, 'vib_list_gas'],
                                              inertia_list=self.df_mechanism.loc[step, 'inertia_list'],
                                              sym_number=self.df_mechanism.loc[step, 'sym_number'])
        elif step_type == 'activated_adsorption':
            pe_fwd, pe_rev = calc_act_ads(A_site=self.df_mechanism.loc[step, 'A_site'],
                                          molec_mass=self.df_mechanism.loc[step, 'molec_mass'],
                                          T=T,
                                          vib_list_ads=self.df_mechanism.loc[step, 'vib_list_ads'],
                                          vib_list_gas=self.df_mechanism.loc[step, 'vib_list_gas'],
                                          vib_list_ts=self.df_mechanism.loc[step, 'vib_list_ts'],
                                          inertia_list=self.df_mechanism.loc[step, 'inertia_list'],
                                          sym_number=self.df_mechanism.loc[step, 'sym_number'])
        elif step_type == 'surface_reaction':
            pe_fwd, pe_rev = calc_surf_react(T=T,
                                             vib_list_initial=self.df_mechanism.loc[step, 'vib_list_initial'],
                                             vib_list_ts=self.df_mechanism.loc[step, 'vib_list_ts'],
                                             vib_list_final=self.df_mechanism.loc[step, 'vib_list_final'])
        else:
            sys.exit(f"Invalid step type: {step_type}")
        pe_ratio = pe_fwd / pe_rev
        return pe_fwd, pe_ratio

