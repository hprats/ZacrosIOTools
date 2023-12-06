import ast
import pandas as pd
from zacrosio.input_functions import write_header


class EnergeticModel:
    """A class that represents a KMC energetic model."""

    def __init__(self, df):
        self.df = df

    def write(self, path):
        """Writes the energetics_input.dat file"""
        write_header(f"{path}/energetics_input.dat")
        with open(f"{path}/energetics_input.dat", 'a') as infile:
            infile.write('energetics\n\n')
            infile.write('############################################################################\n\n')
            for cluster in [x for x in self.df.index if '_gas' not in x]:
                infile.write(f"cluster {cluster}\n\n")
                infile.write(f"  sites {int(self.df.loc[cluster, 'sites'])}\n")
                if not pd.isnull(self.df.loc[cluster, 'neighboring']):
                    infile.write(f"  neighboring {self.df.loc[cluster, 'neighboring']}\n")
                infile.write(f"  lattice_state\n")
                lattice_state_list = ast.literal_eval(self.df.loc[cluster, 'lattice_state'])
                for element in lattice_state_list:
                    infile.write(f"    {element}\n")
                infile.write(f"  site_types {self.df.loc[cluster, 'site_types']}\n")
                if not pd.isnull(self.df.loc[cluster, 'graph_multiplicity']):
                    infile.write(f"  graph_multiplicity {int(self.df.loc[cluster, 'graph_multiplicity'])}\n")
                if not pd.isnull(self.df.loc[cluster, 'angles']):
                    infile.write(f"  angles {self.df.loc[cluster, 'angles']}\n")
                # todo: add column optional keywords (e.g. no_mirror_images)
                infile.write(f"  cluster_eng {self.df.loc[cluster, 'cluster_eng']:.2f}\n\n")
                infile.write(f"end_cluster\n\n")
                infile.write('############################################################################\n\n')
            infile.write(f"end_energetics\n")
