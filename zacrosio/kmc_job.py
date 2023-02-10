import os
import pandas as pd
from random import randint
from zacrosio.functions import *


class NewKMCJob:
    """A class that represents a new KMC job with ZACROS.

        Attributes:
            path (str): The path of the job including the job name. Will be used as the name of the folder.
            simulation_tags (dict): A dictionary including all tags for the simulation_input.dat except random_seed and
            temperature.
            df_mechanism: A Pandas dataframe including the information for the mechanism_input.dat.
            df_energetics: A Pandas dataframe including the information for the energetics_input.dat.
            lattice_path (str): The path of the lattice_input.dat (already created).

        Examples:
        >>> import pandas as pd
        >>> from zacrosio.kmc_job import NewKMCJob
        >>> simulation_tags = {'pressure': 1.0, 'n_gas_species': 5}
        >>> my_job = NewKMCJob(
        >>>    path='/home/my_job',
        >>>    simulation_tags=simulation_tags,
        >>>    df_mechanism=pd.read_csv("/home/input_files/mechanism.csv", index_col=0),
        >>>    df_energetics=pd.read_csv("/home/input_files/energetics.csv", index_col=0),
        >>>    lattice_path='/home/input_files/lattice_input.dat')
        >>> my_job.create_job_dir(T=400, p=2)
        """

    def __init__(self, path, simulation_tags, df_mechanism, df_energetics, lattice_path):
        self.path = path
        self.name = path.split('/')[-1]
        self.simulation_tags = simulation_tags
        self.df_mechanism = df_mechanism
        self.df_energetics = df_energetics
        self.lattice_path = lattice_path

    def create_job_dir(self, T=300, p=1):
        """Creates a new directory and writes there the ZACROS input files"""
        if not os.path.exists(self.path):
            os.mkdir(self.path)
            self.write_simulation(T=T, p=p)
            self.write_mechanism(T=T)
            self.write_energetics()
            self.write_lattice()
        else:
            print(f'{self.path} already exists (nothing done)')

    def write_simulation(self, T=300, p=1):
        """Writes the simulation_input.dat file"""
        gas_specs_names = [x.replace('_gas', '') for x in self.df_energetics.index if '_gas' in x]
        surf_specs_names = [x.replace('_point', '') for x in self.df_energetics.index if '_point' in x]
        surf_specs_names = [x+'*'*int(self.df_energetics.loc[f'{x}_point', 'sites']) for x in surf_specs_names]
        surf_specs_dent = [x.count('*') for x in surf_specs_names]
        self.write_header(file_name="simulation_input.dat")
        with open(f"{self.path}/simulation_input.dat", 'a') as infile:
            infile.write('random_seed\t'.expandtabs(26) + str(randint(100000, 999999)) + '\n')
            infile.write('temperature\t'.expandtabs(26) + str(float(T)) + '\n')
            infile.write('pressure\t'.expandtabs(26) + str(float(p)) + '\n')
            infile.write('n_gas_species\t'.expandtabs(26) + str(len(gas_specs_names)) + '\n')
            infile.write('gas_specs_names\t'.expandtabs(26) + " ".join(str(x) for x in gas_specs_names) + '\n')
            for tag in ['gas_energies', 'gas_molar_fracs', 'gas_molec_weights']:
                tag_list = [self.df_energetics.loc[f"{molecule}_gas", tag] for molecule in gas_specs_names]
                infile.write(f'{tag}\t'.expandtabs(26) + " ".join(str(x) for x in tag_list) + '\n')
            infile.write('n_surf_species\t'.expandtabs(26) + str(len(surf_specs_names)) + '\n')
            infile.write('surf_specs_names\t'.expandtabs(26) + " ".join(str(x) for x in surf_specs_names) + '\n')
            infile.write('surf_specs_dent\t'.expandtabs(26) + " ".join(str(x) for x in surf_specs_dent) + '\n')
            for tag in self.simulation_tags:
                infile.write((tag + '\t').expandtabs(26) + str(self.simulation_tags[tag]) + '\n')
            infile.write(f"no_restart\n")
            infile.write(f"finish\n")

    def write_mechanism(self, T=300):
        """Writes the mechanism_input.dat file"""
        self.write_header(file_name="mechanism_input.dat")
        with open(f"{self.path}/mechanism_input.dat", 'a') as infile:
            infile.write('mechanism\n\n')
            infile.write('############################################################################s\n\n')
            for step in self.df_mechanism.index:
                infile.write(f"reversible_step {step}\n\n")
                step_type = self.df_mechanism.loc[step, 'type']
                if 'adsorption' in step_type:
                    infile.write(f"  gas_reacs_prods {self.df_mechanism.loc[step, 'gas_reacs_prods']}\n")
                infile.write(f"  sites {int(self.df_mechanism.loc[step, 'sites'])}\n")
                if not pd.isnull(self.df_mechanism.loc[step, 'neighboring']):
                    infile.write(f"  neighboring {self.df_mechanism.loc[step, 'neighboring']}\n")
                infile.write(f"  initial\n")
                initial_state_list = ast.literal_eval(self.df_mechanism.loc[step, 'initial'])
                for element in initial_state_list:
                    infile.write(f"    {element}\n")
                infile.write(f"  final\n")
                final_state_list = ast.literal_eval(self.df_mechanism.loc[step, 'final'])
                for element in final_state_list:
                    infile.write(f"    {element}\n")
                infile.write(f"  site_types {self.df_mechanism.loc[step, 'site_types']}\n")
                pre_expon, pe_ratio = self.get_pre_expon(step=step, T=T)
                infile.write(f"  pre_expon {pre_expon:.3e}\n")
                infile.write(f"  pe_ratio {pe_ratio:.3e}\n")
                infile.write(f"  activ_eng {self.df_mechanism.loc[step, 'activ_eng']:.2f}\n")
                for keyword in ['prox_factor', 'angles']:  # optional keywords
                    if not pd.isnull(self.df_mechanism.loc[step, keyword]):
                        infile.write(f"  {keyword} {self.df_mechanism.loc[step, keyword]}\n")
                # todo: add possibility of including no_mirror_images (maybe column with additional flags)
                infile.write(f"\nend_reversible_step\n\n")
                infile.write('############################################################################s\n\n')
            infile.write(f"end_mechanism\n")

    def write_energetics(self):
        """Writes the energetics_input.dat file"""
        self.write_header(file_name="energetics_input.dat")
        with open(f"{self.path}/energetics_input.dat", 'a') as infile:
            infile.write('energetics\n\n')
            infile.write('############################################################################s\n\n')
            for cluster in [x for x in self.df_energetics.index if '_gas' not in x]:
                infile.write(f"cluster {cluster}\n\n")
                infile.write(f"  sites {int(self.df_energetics.loc[cluster, 'sites'])}\n")
                if not pd.isnull(self.df_energetics.loc[cluster, 'neighboring']):
                    infile.write(f"  neighboring {self.df_energetics.loc[cluster, 'neighboring']}\n")
                infile.write(f"  lattice_state\n")
                lattice_state_list = ast.literal_eval(self.df_energetics.loc[cluster, 'lattice_state'])
                for element in lattice_state_list:
                    infile.write(f"    {element}\n")
                infile.write(f"  site_types {self.df_energetics.loc[cluster, 'site_types']}\n")
                if not pd.isnull(self.df_energetics.loc[cluster, 'graph_multiplicity']):
                    infile.write(f"  graph_multiplicity {int(self.df_energetics.loc[cluster, 'graph_multiplicity'])}\n")
                if not pd.isnull(self.df_energetics.loc[cluster, 'angles']):
                    infile.write(f"  angles {self.df_energetics.loc[cluster, 'angles']}\n")
                # todo: add column optional keywords (e.g. no_mirror_images
                infile.write(f"  cluster_eng {self.df_energetics.loc[cluster, 'cluster_eng']:.2f}\n\n")
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
            infile.write('############################################################################\n\n')

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
                                              sym_number=int(self.df_mechanism.loc[step, 'sym_number']),
                                              degeneracy=int(self.df_mechanism.loc[step, 'degeneracy']))
        elif step_type == 'activated_adsorption':
            pe_fwd, pe_rev = calc_act_ads(A_site=self.df_mechanism.loc[step, 'A_site'],
                                          molec_mass=self.df_mechanism.loc[step, 'molec_mass'],
                                          T=T,
                                          vib_list_ads=self.df_mechanism.loc[step, 'vib_list_ads'],
                                          vib_list_gas=self.df_mechanism.loc[step, 'vib_list_gas'],
                                          vib_list_ts=self.df_mechanism.loc[step, 'vib_list_ts'],
                                          inertia_list=self.df_mechanism.loc[step, 'inertia_list'],
                                          sym_number=int(self.df_mechanism.loc[step, 'sym_number']),
                                          degeneracy=int(self.df_mechanism.loc[step, 'degeneracy']))
        elif step_type == 'surface_process':
            pe_fwd, pe_rev = calc_surf_proc(T=T,
                                            vib_list_initial=self.df_mechanism.loc[step, 'vib_list_initial'],
                                            vib_list_ts=self.df_mechanism.loc[step, 'vib_list_ts'],
                                            vib_list_final=self.df_mechanism.loc[step, 'vib_list_final'])
        else:
            sys.exit(f"Invalid step type: {step_type}")
        if not pd.isnull(self.df_mechanism.loc[step, 'scaling_factor']):
            pe_fwd = pe_fwd * self.df_mechanism.loc[step, 'scaling_factor']
            pe_rev = pe_rev * self.df_mechanism.loc[step, 'scaling_factor']
        pe_ratio = pe_fwd / pe_rev
        return pe_fwd, pe_ratio
