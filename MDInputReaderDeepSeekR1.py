class MDInputReader:
    def __init__(self, filename):
        self.filename = filename
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
        self.config = {}
        for section in self.defaults:
            self.config[section] = self.defaults[section].copy()

    def read(self):
        with open(self.filename, 'r') as f:
            current_section = None
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                if line.startswith('[') and line.endswith(']'):
                    current_section = line[1:-1].strip()
                    if current_section not in self.config:
                        self.config[current_section] = {}
                    continue
                if '#' in line:
                    line_part = line.split('#', 1)[0].strip()
                else:
                    line_part = line.strip()
                if '=' not in line_part:
                    continue
                key_part, value_part = line_part.split('=', 1)
                key = key_part.strip().strip('\'"')
                value_str = value_part.strip()

                if current_section == 'Atoms':
                    numbers = []
                    for part in value_str.split():
                        try:
                            num = int(part)
                        except ValueError:
                            num = float(part)
                        numbers.append(num)
                    self.config[current_section][key] = numbers
                else:
                    if current_section in self.defaults and key in self.defaults[current_section]:
                        default_type = type(self.defaults[current_section][key])
                        try:
                            if default_type == int:
                                value = int(value_str)
                            elif default_type == float:
                                value = float(value_str)
                            else:
                                value = value_str
                        except ValueError:
                            value = value_str
                    else:
                        try:
                            value = int(value_str)
                        except ValueError:
                            try:
                                value = float(value_str)
                            except ValueError:
                                value = value_str
                    if current_section not in self.config:
                        self.config[current_section] = {}
                    self.config[current_section][key] = value
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
