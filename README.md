ZacrosIOTools
===========

This project is a collective of tools for the preparation of input files and
post-processing of output files from Zacros.

### Installation

Clone the repo to your preferred location:

    git clone https://github.com/hprats/ZacrosIOTools.git

Add to your .bashrc file a PYTHONPATH to your repository, i.e.

    export PYTHONPATH=path-to-your-repo/ZacrosIOTools:$PYTHONPATH

### Dependencies

ZacrosIOTools relies on a few Python packages, please be sure to have available:

    - pandas 1.5 or higher 
    - scipy 1.10 or higher
    
### Example

With ZacrosIOTools, the input files for a ZACROS job can be prepared as follows:

    from zacrosio.kmc_job import NewKMCJob

    df_mechanism = pd.read_csv("mechanism_data.csv", index_col=0)
    df_energetics = pd.read_csv("energetics_data.csv", index_col=0)
    lattice_path = "./lattice_input.dat"

    simulation_tags = {'snapshots': 'on time 5.e-1',
                       'process_statistics': 'on time 1.e-2',
                       'species_numbers': 'on time 5.e-3',
                       'event_report': 'off',
                       'max_steps': 1000000000,
                       'max_time': 2.0,
                       'wall_time': 86400}

    job = NewKMCJob(path="./new_job",
                    simulation_tags=simulation_tags,
                    df_mechanism=df_mechanism,
                    df_energetics=df_energetics,
                    lattice_path=lattice_path)

    job.create_job_dir(T=1000, p=2)

In this example, ZacrosIOTools will create a new folder named "new_job" and write there all 4 input files at the desired temperature and pressure. The information required to write each input file is described in the following.

#### 1. simulation_input.dat

This file contains information about the species involved, the operating conditions, as well as parameters that specify the behavior of the program, namely when to take samples, what are the stopping criteria. Here, only the keyords related to the frequency of sampling and stopping criteria are required, and should be passed to ZacroIOTools as a dictionary (e.g. see example above).

#### 2. lattice_input.dat

This file defines the lattice model. Currently, the generation of the lattice_input.dat file by ZacrosIOTools is not implemented. This file has to be generated manually and provide its path, so that ZacrosIOTools can copy this file to the new job directory. 

#### 3. energetics_input.dat

This file defines the cluster expansion Hamiltonian to be used for calculating the energy of a given lattice configuration. This information must be given in the form of a Pandas dataframe, where each row of corresponds to a cluster. The clusters are classified as gas-phase molecules (X_gas), point clusters (X_point) or pairwise lateral interactions (X+Y_pair). 
The entries corresponding to gas_phase molecules must end in '_gas' (e.g. CH4_gas) and the following columns are required:
- gas_energy (float, in eV), e.g. 1.40
- gas_molar_frac (float), e.g. 0.65
- gas_molec_weight (float, in g/mol), e.g. 16.04

The entries corresponding to point clusters and pairwise lateral interactions must end in '_point' and '_pair', respectively (e.g. CH3_point, CH2+H_pair) and the following columns are required:
- sites (int), e.g. 2
- site_types (str), e.g. tM tC
- lattice_state (list), e.g. ['1 CH3** 1', '1 CH3** 2']
- cluster_eng (float, in eV), e.g. -0.42

Optional columns:
- neighboring (str), e.g. 1-2
- graph_multiplicity (int), e.g. 2
- angles (str), e.g. 1-2-3:180

#### 4. mechanism_input.dat

This file defines the reaction mechanism. This information must be given in the form of a Pandas dataframe, where each row of corresponds to a reversible elementary step (e.g. adsorption, desorption, diffusion or surface reaction). Irreversible steps can not be defined in the current implementation.
The entries corresponding to elementary steps can be given any name (e.g. CO_dissociation) and the following columns are required:
- type (str): choose between 'non_activated_adsorption', 'activated_adsorption' or 'surface_process'
- sites (int), e.g. 2
- site_types (str), e.g. tM tC
- initial (list), e.g. ['1 * 1','2 CH3** 1','2 CH3** 2']
- final (list), e.g. ['1 H_tC* 1','2 CH2** 1','2 CH2** 2']
- activ_eng (float, in eV), e.g. 1.02

Additional required columns for 'non_activated_adsorption' and 'activated_adsorption' steps:
- gas_reacs_prods (str), e.g. CO -1
- A_site (float), e.g. 4.28
- vib_list_ads (list of floats, in meV), e.g. [249.1, 82.8, 63.8, 63.8, 8.9, 8.3]
- vib_list_gas (list of floats, in meV), e.g. [249.2]
- inertia_list (list of floats, in amu*Å2), 1/3 elements for linear/non-linear molecules, e.g. [8.9]
- sym_number (int), e.g. 1
- degeneracy (int), e.g. 1

Additional required columns for 'activated_adsorption' and 'surface_process' steps:
- vib_list_ts (list of floats, in meV), e.g. [332.7, 196.2, 70.5, 53.9, 37.7]

Additional required columns for 'surface_process' steps:
- vib_list_initial (list of floats, in meV), e.g. [332.7, 196.2, 70.5, 53.9, 37.7, 10.0]
- vib_list_final (list of floats, in meV), e.g. [332.7, 196.2, 70.5, 53.9, 37.7, 10.0]

Optional columns:
- neighboring (str), e.g. 1-2
- angles (str), e.g. 1-2-3:180
- scaling fator (float), e.g. 0.001
- prox_factor (str), e.g. 0.3 # default is 0.5
