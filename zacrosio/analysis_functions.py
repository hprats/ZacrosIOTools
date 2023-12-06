import numpy as np


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]


def get_tof(path, molecule, area, ignore=0.2):
    with open(f"{path}/specnum_output.txt", "r") as infile:
        header = infile.readline().split()
    data = np.loadtxt(f"{path}/specnum_output.txt", skiprows=1)
    production = data[:, header.index(molecule)]
    time = data[:, 2]
    final_time = data[-1, 2]
    index = np.where(time == find_nearest(time, final_time * ignore))[0][0]
    production = np.delete(production, np.s_[0:index], axis=0)
    time = np.delete(time, np.s_[0:index], axis=0)
    tof = np.polyfit(time, production, 1)[0] / area
    if tof < 1.0e-6:
        tof = 1.0e-6
    return tof


def get_selectivity(path, main, secondary, minimum=0.0, ignore=0.2, return_tof=False, area=None):
    with open(f"{path}/specnum_output.txt", "r") as infile:
        header = infile.readline().split()
    data = np.loadtxt(f"{path}/specnum_output.txt", skiprows=1)
    production_main = data[:, header.index(main)]
    production_secondary = np.zeros(len(data))
    for molecule in secondary:
        production_secondary += data[:, header.index(molecule)]
    time = data[:, 2]
    final_time = data[-1, 2]
    index = np.where(time == find_nearest(time, final_time * ignore))[0][0]
    production_main = np.delete(production_main, np.s_[0:index], axis=0)
    production_secondary = np.delete(production_secondary, np.s_[0:index], axis=0)
    production_total = production_main + production_secondary
    time = np.delete(time, np.s_[0:index], axis=0)
    tof_total = np.polyfit(time, production_total, 1)[0] / area
    if tof_total < 1.0e-6:
        selectivity = float('NaN')
    else:
        mol_main = production_main[-1] - production_main[0]
        mol_total = production_total[-1] - production_total[0]
        selectivity = mol_main / mol_total * 100
    if return_tof:
        if selectivity > minimum:
            tof = get_tof(path=path, molecule=main, area=area, ignore=ignore)
        else:
            tof = float('NaN')
        return tof
    else:
        if selectivity > 100:
            print(selectivity)
    return selectivity




