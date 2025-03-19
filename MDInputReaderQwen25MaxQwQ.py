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
            lines = file.readlines()

        current_section = None

        for line in lines:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            # Handle section headers
            if line.startswith('[') and line.endswith(']'):
                current_section = line[1:-1].strip()
                if current_section not in self.config:
                    self.config[current_section] = {}
                continue

            # Skip lines without key-value pairs
            if '=' not in line:
                continue

            # Parse key-value pairs
            key_part, value_part = line.split('=', 1)
            key = key_part.strip().strip("'\"")  # Handle quoted keys
            value = value_part.split('#', 1)[0].strip()  # Remove inline comments

            # Special handling for Atoms section
            if current_section == 'Atoms':
                elements = value.split()
                converted = []
                for elem in elements:
                    try:
                        converted.append(int(elem))
                    except ValueError:
                        try:
                            converted.append(float(elem))
                        except ValueError:
                            converted.append(elem)  # Keep as string if conversion fails
                self.config[current_section][key] = converted
            else:
                # Handle type conversion based on defaults
                if (current_section in self.config and
                        key in self.config[current_section]):
                    default = self.config[current_section][key]
                    try:
                        if isinstance(default, int):
                            self.config[current_section][key] = int(value)
                        elif isinstance(default, float):
                            self.config[current_section][key] = float(value)
                        else:
                            self.config[current_section][key] = value.strip()
                    except ValueError:
                        self.config[current_section][key] = value.strip()
                else:
                    self.config[current_section][key] = value.strip()

        return self.config

    # Accessor methods for default parameters
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
