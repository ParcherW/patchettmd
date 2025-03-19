class MDInputReader:
    def __init__(self, filename):
        self.filename = filename
        # Initialize default configuration
        self.config = {
            'Simulation': {
                'time_step': 0.00001,
                'n_steps': 20000,
                'ensemble': 'NVE',
            },
            'System': {
                'box_length': 10.0,
                'n_particles': 20,
                'boundary': 'reflecting',
            },
            'Potential': {
                'potential_type': 'Lennard-Jones',
                'epsilon': 30.0,
                'sigma': 2.5,
                'cutoff': 2.5,
                'coulomb_k': 9e9,
            },
            'Output': {
                'output_frequency': 200,
                'output_file': 'trajectory.xyz',
            },
        }

    def read(self):
        with open(self.filename, 'r') as file:
            lines = file.read().splitlines()

        current_section = None
        for line in lines:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            if line.startswith('[') and line.endswith(']'):
                current_section = line[1:-1].strip()
                if current_section not in self.config:
                    self.config[current_section] = {}
                continue

            if '=' not in line:
                continue

            key_part, value_part = line.split('=', 1)
            key = key_part.strip().strip('"\'')
            value_str = value_part.split('#', 1)[0].strip()

            if current_section == 'Atoms':
                value = []
                for part in value_str.split():
                    try:
                        num = int(part)
                    except ValueError:
                        try:
                            num = float(part)
                        except ValueError:
                            num = part
                    value.append(num)
            else:
                try:
                    value = int(value_str)
                except ValueError:
                    try:
                        value = float(value_str)
                    except ValueError:
                        value = value_str

            self.config[current_section][key] = value

        return self.config

    # Accessor methods for default parameters
    def get_time_step(self):
        return float(self.config['Simulation']['time_step'])

    def get_n_steps(self):
        return int(self.config['Simulation']['n_steps'])

    def get_ensemble(self):
        return str(self.config['Simulation']['ensemble'])

    def get_box_length(self):
        return float(self.config['System']['box_length'])

    def get_n_particles(self):
        return int(self.config['System']['n_particles'])

    def get_boundary(self):
        return str(self.config['System']['boundary'])

    def get_potential_type(self):
        return str(self.config['Potential']['potential_type'])

    def get_epsilon(self):
        return float(self.config['Potential']['epsilon'])

    def get_sigma(self):
        return float(self.config['Potential']['sigma'])

    def get_cutoff(self):
        return float(self.config['Potential']['cutoff'])

    def get_coulomb_k(self):
        return float(self.config['Potential']['coulomb_k'])

    def get_output_frequency(self):
        return int(self.config['Output']['output_frequency'])

    def get_output_file(self):
        return str(self.config['Output']['output_file'])
