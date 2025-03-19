class MDInputReader:
    def __init__(self, filename):
        self.filename = filename
        self.default_config = {
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
        for section, params in self.default_config.items():
            self.config[section] = params.copy()

    def read(self):
        current_section = None
        with open(self.filename, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                if line.startswith('[') and line.endswith(']'):
                    current_section = line[1:-1].strip()
                    if current_section not in self.config:
                        self.config[current_section] = {}
                    continue
                if '=' in line:
                    key_part, value_part = line.split('=', 1)
                    key = key_part.strip().strip('"').strip("'").strip()
                    value_str = value_part.split('#', 1)[0].strip()
                    if current_section == 'Atoms':
                        parts = value_str.split()
                        parsed = []
                        for part in parts:
                            try:
                                parsed.append(int(part))
                            except ValueError:
                                try:
                                    parsed.append(float(part))
                                except ValueError:
                                    parsed.append(part)
                        value = parsed
                    else:
                        try:
                            value = int(value_str)
                        except ValueError:
                            try:
                                value = float(value_str)
                            except ValueError:
                                value = value_str
                    if current_section in self.config:
                        self.config[current_section][key] = value
                    else:
                        self.config[current_section] = {key: value}
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
