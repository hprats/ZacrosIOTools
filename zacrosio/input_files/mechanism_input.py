import ast
import pandas as pd
from zacrosio.input_functions import write_header, calc_ads, calc_surf_proc


class ReactionModel:
    """A class that represents a KMC reaction model."""

    def __init__(self, df):
        self.df = df

    def write(self, path, T, df_gas, dict_manual_scaling, list_auto_scaling):
        """Writes the mechanism_input.dat file"""
        write_header(f"{path}/mechanism_input.dat")
        with open(f"{path}/mechanism_input.dat", 'a') as infile:
            infile.write('mechanism\n\n')
            infile.write('############################################################################\n\n')
            for step in self.df.index:
                infile.write(f"reversible_step {step}\n\n")
                if not pd.isna(self.df.loc[step, 'molecule']):
                    infile.write(f"  gas_reacs_prods {self.df.loc[step, 'molecule']} -1\n")
                infile.write(f"  sites {int(self.df.loc[step, 'sites'])}\n")
                if not pd.isnull(self.df.loc[step, 'neighboring']):
                    infile.write(f"  neighboring {self.df.loc[step, 'neighboring']}\n")
                infile.write(f"  initial\n")
                initial_state_list = ast.literal_eval(self.df.loc[step, 'initial'])
                for element in initial_state_list:
                    infile.write(f"    {element}\n")
                infile.write(f"  final\n")
                final_state_list = ast.literal_eval(self.df.loc[step, 'final'])
                for element in final_state_list:
                    infile.write(f"    {element}\n")
                infile.write(f"  site_types {self.df.loc[step, 'site_types']}\n")
                pre_expon, pe_ratio = self.get_pre_expon(step=step, T=T, df_gas=df_gas,
                                                         dict_manual_scaling=dict_manual_scaling)
                if step in dict_manual_scaling:
                    infile.write(f"  pre_expon {pre_expon:.3e}   # scaled 1e-{dict_manual_scaling[step]}\n")
                else:
                    infile.write(f"  pre_expon {pre_expon:.3e}\n")
                infile.write(f"  pe_ratio {pe_ratio:.3e}\n")
                infile.write(f"  activ_eng {self.df.loc[step, 'activ_eng']:.2f}\n")
                for keyword in ['prox_factor', 'angles']:  # optional keywords
                    if not pd.isnull(self.df.loc[step, keyword]):
                        infile.write(f"  {keyword} {self.df.loc[step, keyword]}\n")
                if step in list_auto_scaling:
                    infile.write(f"  stiffness_scalable \n")
                # todo: add possibility of including no_mirror_images (maybe column with additional flags)
                infile.write(f"\nend_reversible_step\n\n")
                infile.write('############################################################################\n\n')
            infile.write(f"end_mechanism\n")

    def get_pre_expon(self, step, T, df_gas, dict_manual_scaling):
        """Calculates the forward pre-exponential and the pre-exponential ratio, required for the mechanism_input.dat
        file """
        if not pd.isna(self.df.loc[step, 'molecule']):  # adsorption
            molecule = self.df.loc[step, 'molecule']
            molec_mass = df_gas.loc[molecule, 'gas_molec_weight']
            pe_fwd, pe_rev = calc_ads(A_site=self.df.loc[step, 'A_site'],
                                      molec_mass=molec_mass,
                                      T=T,
                                      vib_list_is=self.df.loc[step, 'vib_list_is'],
                                      vib_list_ts=self.df.loc[step, 'vib_list_ts'],
                                      vib_list_fs=self.df.loc[step, 'vib_list_fs'],
                                      inertia_list=df_gas.loc[molecule, 'inertia_list'],
                                      sym_number=int(df_gas.loc[molecule, 'sym_number']),
                                      degeneracy=int(df_gas.loc[molecule, 'degeneracy']))
        else:  # surface process
            pe_fwd, pe_rev = calc_surf_proc(T=T,
                                            vib_list_is=self.df.loc[step, 'vib_list_is'],
                                            vib_list_ts=self.df.loc[step, 'vib_list_ts'],
                                            vib_list_fs=self.df.loc[step, 'vib_list_fs'])

        if step in dict_manual_scaling:
            pe_fwd = pe_fwd * 10 ** (-dict_manual_scaling[step])
            pe_rev = pe_rev * 10 ** (-dict_manual_scaling[step])
        pe_ratio = pe_fwd / pe_rev
        return pe_fwd, pe_ratio
