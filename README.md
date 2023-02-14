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
- site_types (string), e.g. tM tC
- lattice_state (list of one element), e.g. ['1 CH3** 1', '1 CH3** 2']
- cluster_eng (float, in eV), e.g. -0.42
Optional columns:
- neighboring (string), e.g. 1-2
- graph_multiplicity (int), e.g. 2
- angles (string), e.g. 1-2-3:180

#### 4. mechanism_input.dat

This file defines the reaction mechanism. This information must be given in the form of a Pandas dataframe, where each row of corresponds to an elementary step (e.g. adsorption, desorption, diffusion or surface reaction). (Todo ...)
