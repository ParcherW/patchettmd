class MDInputReader:
    def __init__(self, filename):
        self.filename = filename
        self.config = self.parse_file()

    def read(self):
        return self.config

    def parse_file(self):
        config = {}
        current_section = None
        current_line = ''

        try:
            with open(self.filename, 'r') as file:
                for line in file:
                    line = line.strip()
                    if line.startswith('['):
                        current_section = line[1:-1]
                        config[current_section] = {}
                    elif line.startswith('#'):
                        continue
                    elif '=' in line:
                        key, value = line.split('=', 1)
                        key = key.strip()
                        value = value.strip()
                        if current_section:
                            try:
                                config[current_section][key] = float(value)
                            except ValueError:
                                config[current_section][key] = value
                        else:
                            raise ValueError(f"Key '{key}' is not part of a section")
                    else:
                        current_line += line + '\n'

                if current_line:
                    current_section = None
                    for line in current_line.split('\n'):
                        line = line.strip()
                        if line.startswith('['):
                            current_section = line[1:-1]
                            config[current_section] = {}
                        elif line.startswith('#'):
                            continue
                        elif '=' in line:
                            key, value = line.split('=', 1)
                            key = key.strip()
                            value = value.strip()
                            if current_section:
                                try:
                                    config[current_section][key] = int(value)
                                except ValueError:
                                    config[current_section][key] = value
                            else:
                                raise ValueError(f"Key '{key}' is not part of a section")
                        else:
                            if current_section:
                                try:
                                    config[current_section][key] = [int(x) for x in line.split()]
                                except ValueError:
                                    config[current_section][key] = [x for x in line.split()]
                            else:
                                raise ValueError(f"Key '{key}' is not part of a section")

        except FileNotFoundError:
            raise ValueError(f"File '{self.filename}' not found")

        # Add defaults
        defaults = {
            'Simulation': {
                'time_step': float(0.00001),
                'n_steps': int(20000),
                'ensemble': 'NVE'
            },
            'System': {
                'box_length': float(10.0),
                'n_particles': int(20),
                'boundary': 'reflecting'
            },
            'Potential': {
                'potential_type': 'Lennard-Jones',
                'epsilon': float(30.0),
                'sigma': float(2.5),
                'cutoff': float(2.5),
                'coulomb_k': float(9e9)
            },
            'Output': {
                'output_frequency': int(200),
                'output_file': 'trajectory.xyz'
            }
        }
        for section, values in defaults.items():
            if section not in config:
                config[section] = values

        return config

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
