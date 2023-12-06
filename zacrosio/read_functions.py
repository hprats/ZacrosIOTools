import subprocess
import sys
import os


def check_finished(path):
    return 'Normal termination' in str(subprocess.check_output(f"tail -n4 {path}/general_output.txt", shell=True))


def read_file(file_path):
    with open(file_path, 'r') as infile:
        lines = infile.readlines()
    return lines


def read_last_line(file_path):
    with open(file_path, 'rb') as f:
        try:  # catch OSError in case of a one line file
            f.seek(-2, os.SEEK_END)
            while f.read(1) != b'\n':
                f.seek(-2, os.SEEK_CUR)
        except OSError:
            f.seek(0)
        return f.readline().decode()


def get_data_from_general_output(path):
    lines = read_file(f"{path}/general_output.txt")
    n_surf_species, n_sites, area = 0, 0, 0
    for line in lines:
        if 'Number of surface species' in line:
            n_surf_species = int(line.split()[-1])
        if 'Total number of lattice sites' in line:
            n_sites = int(line.split()[-1])
        if 'Lattice surface area' in line:
            area = float(line.split()[-1])
    if n_surf_species == 0:
        sys.exit("ERROR: n_surf_species = 0")
    elif n_sites == 0:
        sys.exit("ERROR: n_sites = 0")
    elif area == 0:
        sys.exit("ERROR: area = 0")
    return n_surf_species, n_sites, area


def parse_line(dmatch, line):
    for key in dmatch.items():
        if key in line:
            return key, True
    return None, None


def parse_general_output(path):
    dmatch = ['Number of surface species',
              'Total number of lattice sites',
              'Lattice surface area',
              'Site type names and total number of sites of that type']
    data = {}
    with open(f"{path}/general_output.txt", 'r') as file_object:
        line = file_object.readline()
        while len(dmatch) != 0:
            if 'Number of surface species' in line:
                data['n_surf_species'] = int(line.split()[-1])
                dmatch.remove('Number of surface species')
            if 'Total number of lattice sites' in line:
                data['n_sites'] = int(line.split()[-1])
                dmatch.remove('Total number of lattice sites')
            if 'Lattice surface area' in line:
                data['area'] = float(line.split()[-1])
                dmatch.remove('Lattice surface area')
            if 'Site type names and total number of sites of that type' in line:
                line = file_object.readline()
                site_types = {}
                while line.strip():
                    num_sites_of_given_type = int(line.strip().split(' ')[1].replace('(', '').replace(')', ''))
                    site_types[line.strip().split(' ')[0]] = num_sites_of_given_type
                    line = file_object.readline()
                data['site_types'] = site_types
                dmatch.remove('Site type names and total number of sites of that type')
            line = file_object.readline()
        return data
