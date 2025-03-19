class MDInputReader:
    def __init__(self, filename):
        self.filename = filename
        self.config = None
        self.defaults = {
            'Simulation': {
                'time_step': 0.00001,
                'n_steps': 20000,
                'ensemble': 'NVE'
            },
            'System': {
                'box_length': 10.0,
                'n_particles': 20,
                'boundary': 'reflecting'
            },
            'Potential': {
                'potential_type': 'Lennard-Jones',
                'epsilon': 30.0,
                'sigma': 2.5,
                'cutoff': 2.5,
                'coulomb_k': 9e9
            },
            'Output': {
                'output_frequency': 200,
                'output_file': 'trajectory.xyz'
            }
        }

    def read(self):
        if self.config is None:
            self.config = self.defaults.copy()
            self._parse_file()
        return self.config

    def _parse_file(self):
        current_section = None
        with open(self.filename, 'r') as file:
            for line in file:
                line = line.strip()
                if line.startswith('#'):
                    continue
                if line.startswith('[') and line.endswith(']'):
                    current_section = line[1:-1]
                    if current_section not in self.config:
                        self.config[current_section] = {}
                elif '=' in line and current_section:
                    key, value = [part.strip() for part in line.split('=', 1)]
                    key = key.strip('"')
                    value = value.split('#')[0].strip()
                    if current_section == 'Atoms':
                        try:
                            self.config[current_section][key] = [float(x) for x in value.split()]
                        except ValueError:
                            self.config[current_section][key] = value.split()
                    else:
                        try:
                            self.config[current_section][key] = float(value)
                        except ValueError:
                            try:
                                self.config[current_section][key] = int(value)
                            except ValueError:
                                self.config[current_section][key] = value

    def get_time_step(self):
        return self.config['Simulation']['time_step']

    def get_n_steps(self):
        return self.config['Simulation']['n_steps']

    def get_ensemble(self):
        return self.config['Simulation']['ensemble']

    def get_box_length(self):
        return self.config['System']['box_length']

    def get_n_particles(self):
        return self.config['System']['n_particles']

    def get_boundary(self):
        return self.config['System']['boundary']

    def get_potential_type(self):
        return self.config['Potential']['potential_type']

    def get_epsilon(self):
        return self.config['Potential']['epsilon']

    def get_sigma(self):
        return self.config['Potential']['sigma']

    def get_cutoff(self):
        return self.config['Potential']['cutoff']

    def get_coulomb_k(self):
        return self.config['Potential']['coulomb_k']

    def get_output_frequency(self):
        return self.config['Output']['output_frequency']

    def get_output_file(self):
        return self.config['Output']['output_file']
