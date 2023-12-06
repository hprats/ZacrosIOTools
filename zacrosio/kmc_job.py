import os
import pandas as pd
from random import randint
from zacrosio.input_files.energetics_input import EnergeticModel
from zacrosio.input_files.mechanism_input import ReactionModel
from zacrosio.input_files.lattice_input import LatticeModel
from zacrosio.input_functions import write_header
from zacrosio.read_functions import read_file, check_finished, get_data_from_general_output
from zacrosio.plot_functions import plt_production, plt_tof


class NewKMCJob:
    """A class that represents a new KMC job with ZACROS.

        Attributes: path (str): The path of the job including the job name. Will be used as the name of the folder.
        simulation_tags (dict): A dictionary including keywords relatred to the frequency of sampling and stopping
        criteria, for the simulation_input.dat.
        df_mechanism: A Pandas dataframe including the information for the mechanism_input.dat.
        df_energetics: A Pandas dataframe including the information for the energetics_input.dat. lattice_path (str):
        The path of the lattice_input.dat (already created).

        Example:
        >>> import pandas as pd
        >>> from zacrosio.kmc_job import NewKMCJob
        >>>
        >>> my_job = NewKMCJob(
        >>>     df_gas=pd.read_csv("gas_data.csv", index_col=0),
        >>>     df_mechanism=pd.read_csv("mechanism_data.csv", index_col=0),
        >>>     df_energetics=pd.read_csv("energetics_data.csv", index_col=0),
        >>>     lattice_path="lattice_input.dat")
        >>>
        >>> simulation_tags = {'snapshots': 'on time 5.e-1', 'process_statistics': 'on time 1.e-2',
        >>>                    'species_numbers': 'on time 5.e-3', 'event_report': 'off', 'max_steps': 1000000000,
        >>>                    'max_time': 2.0, 'wall_time': 86400}
        >>>
        >>> for temp in range(400, 1000, 100):
        >>>     my_job.create_job_dir(
        >>>         path=f'./co_oxidation_{temp}K',
        >>>         T=temp,  # in K
        >>>         simulation_tags=simulation_tags,
        >>>         dict_pressure={'CO': 1.2, 'O2': 0.01},  # in bar
        >>>         repeat_cell=[10, 10],
        >>>         dict_manual_scaling={'O_diffusion': 1e-3},
        >>>         list_auto_scaling=['CO_adsorption'])
        """

    def __init__(self, df_gas, df_mechanism, df_energetics, lattice_path):
        self.path = None
        self.df_gas = df_gas
        self.reaction_model = ReactionModel(df=df_mechanism)
        self.energetic_model = EnergeticModel(df=df_energetics)
        self.lattice_model = LatticeModel(path=lattice_path)

    def create_job_dir(self, path, T, simulation_tags, dict_pressure, repeat_cell=None, dict_manual_scaling=None,
                       list_auto_scaling=None):
        """Creates a new directory and writes there the ZACROS input files"""
        if list_auto_scaling is None:
            list_auto_scaling = []
        if dict_manual_scaling is None:
            dict_manual_scaling = {}
        self.path = path
        if not os.path.exists(self.path):
            os.mkdir(self.path)
            self.write_simulation(T=T, simulation_tags=simulation_tags, dict_pressure=dict_pressure,
                                  list_auto_scaling=list_auto_scaling)
            self.reaction_model.write(path=self.path, T=T, df_gas=self.df_gas, dict_manual_scaling=dict_manual_scaling,
                                      list_auto_scaling=list_auto_scaling)
            self.energetic_model.write(path=self.path)
            if repeat_cell is not None:
                self.lattice_model.update_size(size=repeat_cell)
            self.lattice_model.write(path=self.path)
        else:
            print(f'{self.path} already exists (nothing done)')

    def write_simulation(self, T, simulation_tags, dict_pressure, list_auto_scaling):
        """Writes the simulation_input.dat file"""
        gas_specs_names = [x for x in self.df_gas.index]
        surf_specs_names = [x.replace('_point', '') for x in self.energetic_model.df.index if '_point' in x]
        surf_specs_names = [x + '*' * int(self.energetic_model.df.loc[f'{x}_point', 'sites']) for x in surf_specs_names]
        surf_specs_dent = [x.count('*') for x in surf_specs_names]
        write_header(f"{self.path}/simulation_input.dat")
        with open(f"{self.path}/simulation_input.dat", 'a') as infile:
            infile.write('random_seed\t'.expandtabs(26) + str(randint(100000, 999999)) + '\n')
            infile.write('temperature\t'.expandtabs(26) + str(float(T)) + '\n')
            p_tot = sum(dict_pressure.values())
            infile.write('pressure\t'.expandtabs(26) + str(p_tot) + '\n')
            infile.write('n_gas_species\t'.expandtabs(26) + str(len(gas_specs_names)) + '\n')
            infile.write('gas_specs_names\t'.expandtabs(26) + " ".join(str(x) for x in gas_specs_names) + '\n')
            tags_dict = ['gas_energy', 'gas_molec_weight']
            tags_zacros = ['gas_energies', 'gas_molec_weights']
            for tag1, tag2 in zip(tags_dict, tags_zacros):
                tag_list = [self.df_gas.loc[x, tag1] for x in gas_specs_names]
                infile.write(f'{tag2}\t'.expandtabs(26) + " ".join(str(x) for x in tag_list) + '\n')
            gas_molar_frac_list = [dict_pressure[x] / p_tot for x in gas_specs_names]
            infile.write(f'gas_molar_fracs\t'.expandtabs(26) + " ".join(str(x) for x in gas_molar_frac_list) + '\n')
            infile.write('n_surf_species\t'.expandtabs(26) + str(len(surf_specs_names)) + '\n')
            infile.write('surf_specs_names\t'.expandtabs(26) + " ".join(str(x) for x in surf_specs_names) + '\n')
            infile.write('surf_specs_dent\t'.expandtabs(26) + " ".join(str(x) for x in surf_specs_dent) + '\n')
            for tag in simulation_tags:
                infile.write((tag + '\t').expandtabs(26) + str(simulation_tags[tag]) + '\n')
            if len(list_auto_scaling) > 0:
                infile.write(f"enable_stiffness_scaling\n")
                infile.write(f"check_every 5000\n")
                infile.write(f"min_separation 200.0\n")
                infile.write(f"max_separation 400.0\n")
                infile.write(f"scaling_factor 3.0\n")
            infile.write(f"finish\n")


class KMCJob:
    """A class that represents a finished KMC job with ZACROS."""

    def __init__(self, path):
        self.path = path
        self.finished = check_finished(path)

        if self.finished:
            # Read files
            self.lattice_input = read_file(f"{path}/lattice_input.dat")
            # Get basic info
            n_surf_species, n_sites, area = get_data_from_general_output(path)
            self.n_surf_species = n_surf_species
            self.n_sites = n_sites
            self.area = area

    def plot_production(self, molecule=None):
        plt_production(self.path, self.n_surf_species, molecule)

    def plot_tof(self, molecule):
        plt_tof(self.path, self.area, molecule)
