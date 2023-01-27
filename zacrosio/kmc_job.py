import os
import ast


class NewKMCJob:
    """A class that represents a new KMC job with ZACROS.
    """

    def __init__(self, path, simulation_tags, df_mechanism, df_energetics, lattice_path):
        self.path = path
        self.name = path.split('/')[-1]
        self.simulation_tags = simulation_tags
        self.df_mechanism = df_mechanism
        self.df_energetics = df_energetics
        self.lattice_path = lattice_path

    def create_job_dir(self, temperature=300):
        """Creates a new directory and writes there the ZACROS input files"""
        if not os.path.exists(self.path):
            os.mkdir(self.path)
            self.write_simulation(temperature=temperature)
            self.write_mechanism(temperature=temperature)
            self.write_energetics()
            self.write_lattice()
        else:
            print(f'{self.path} already exists (nothing done)')

    def write_header(self, file_name):
        """Description"""
        with open(f"{self.path}/{file_name}", 'w') as infile:
            infile.write('############################################################################\n')
            infile.write('# Zacros Input File generated with the ZacrosIOTools                       #\n')
            infile.write('# https://github.com/hprats/ZacrosIOTools.git                              #\n')
            infile.write('#                                                                          #\n')
            infile.write('# Hector Prats, PhD                                                        #\n')
            infile.write('############################################################################\n')


    def write_simulation(self, temperature=300):
        """Description"""
        self.simulation_tags['temperature'] = temperature
        with open(f"{self.path}/simulation_input.dat", 'w') as infile:
            for tag in self.simulation_tags:
                infile.write(f"{tag} = {self.simulation_tags[tag]}\n")
            infile.write(f"no_restart\n")
            infile.write(f"finish\n")

    def write_mechanism(self, temperature=300):
        """Description"""
        self.write_header(file_name="mechanism_input.dat")
        with open(f"{self.path}/mechanism_input.dat", 'a') as infile:
            infile.write('energetics\n\n')
            infile.write('############################################################################s\n\n')
            for cluster in self.df_energetics.index:
                infile.write(f"cluster {cluster}\n\n")
                sites = self.df_energetics.loc[cluster, 'sites']
                infile.write(f"\tsites {sites}\n")
                if sites > 1:
                    infile.write(f"\tneighboring {self.df_energetics.loc[cluster, 'neighboring']}\n")
                infile.write(f"\tlattice_state\n")
                lattice_state_list = ast.literal_eval(self.df_energetics.loc[cluster, 'lattice_state'])
                for element in lattice_state_list:
                    infile.write(f"\t\t{element}\n\n")
                infile.write(f"\tcluster_eng {self.df_energetics.loc[cluster, 'cluster_eng']}\n\n")
                infile.write(f"end_cluster\n\n")
                infile.write('############################################################################s\n\n')

    def write_energetics(self):
        """Description"""
        self.write_header(file_name="energetics_input.dat")
        with open(f"{self.path}/energetics_input.dat", 'a') as infile:
            infile.write('energetics\n\n')
            infile.write('############################################################################s\n\n')
            for cluster in self.df_energetics.index:
                infile.write(f"cluster {cluster}\n\n")
                sites = self.df_energetics.loc[cluster, 'sites']
                infile.write(f"\tsites {sites}\n")
                if sites > 1:
                    infile.write(f"\tneighboring {self.df_energetics.loc[cluster, 'neighboring']}\n")
                infile.write(f"\tlattice_state\n")
                lattice_state_list = ast.literal_eval(self.df_energetics.loc[cluster, 'lattice_state'])
                for element in lattice_state_list:
                    infile.write(f"\t\t{element}\n\n")
                infile.write(f"\tcluster_eng {self.df_energetics.loc[cluster, 'cluster_eng']}\n\n")
                infile.write(f"end_cluster\n\n")
                infile.write('############################################################################s\n\n')

    def write_lattice(self):
        """Description"""
        os.system(f"cp {self.lattice_path} {self.path}/lattice_input.dat")
