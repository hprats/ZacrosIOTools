
class LatticeModel:
    """Read an existing lattice_input.dat file and store it as a LatticeModel object"""

    def __init__(self, path):
        with open(path, 'r') as infile:
            lines = infile.readlines()
        self._lines = lines

    def write(self, path):
        """Write the lattice_input.dat file."""
        with open(f"{path}/lattice_input.dat", 'w') as infile:
            for line in self._lines:
                infile.write(line)

    def update_size(self, size):
        """Modify the value of repeat_cell tag."""
        for i, line in enumerate(self._lines):
            if 'repeat_cell' in line:
                self._lines[i] = f'   repeat_cell {size[0]} {size[1]}\n'
