class MDInputReader:
    def __init__(self, filename):
        self.filename = filename
        self.config = {
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
        current_section = None
        with open(self.filename, 'r') as file:
            for line in file:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                if line.startswith('[') and line.endswith(']'):
                    current_section = line[1:-1].strip()
                    if current_section not in self.config:
                        self.config[current_section] = {}
                    continue
                key_value_line = line.split('#', 1)[0].strip()
                if '=' not in key_value_line:
                    continue
                key_part, value_part = key_value_line.split('=', 1)
                key = key_part.strip().strip('"\'')
                value = value_part.strip()
                if current_section == 'Atoms':
                    numbers = []
                    for part in value.split():
                        try:
                            num = int(part)
                        except ValueError:
                            try:
                                num = float(part)
                            except ValueError:
                                raise ValueError(f"Invalid number '{part}' in [Atoms] section for key '{key}'")
                        numbers.append(num)
                    self.config[current_section][key] = numbers
                else:
                    converted_value = None
                    try:
                        converted_value = int(value)
                    except ValueError:
                        try:
                            converted_value = float(value)
                        except ValueError:
                            converted_value = value
                    self.config[current_section][key] = converted_value
        return self.config

    # Accessor methods for Simulation section
    def get_time_step(self):
        return self.config['Simulation']['time_step']

    def get_n_steps(self):
        return self.config['Simulation']['n_steps']

    def get_ensemble(self):
        return self.config['Simulation']['ensemble']

    # Accessor methods for System section
    def get_box_length(self):
        return self.config['System']['box_length']

    def get_n_particles(self):
        return self.config['System']['n_particles']

    def get_boundary(self):
        return self.config['System']['boundary']

    # Accessor methods for Potential section
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

    # Accessor methods for Output section
    def get_output_frequency(self):
        return self.config['Output']['output_frequency']

    def get_output_file(self):
        return self.config['Output']['output_file']
