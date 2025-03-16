class MDInputReader:
    def __init__(self, input_source=None):
        # Default values
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

        # Initialize a dictionary to hold the parsed data
        self.data = {section: {} for section in self.defaults}

        # Load input if provided
        if input_source:
            self.load_input(input_source)

    def load_input(self, source):
        if isinstance(source, str):
            if '\n' in source:
                lines = source.split('\n')
            else:
                with open(source, 'r') as file:
                    lines = file.readlines()
        else:
            raise ValueError("Input source must be a string or a valid file path.")

        current_section = None
        for line in lines:
            line = line.strip()

            # Skip empty lines and comments
            if not line or line.startswith('#'):
                continue

            # Handle section headers
            if line.startswith('[') and line.endswith(']'):
                section_name = line[1:-1].strip()
                if section_name in self.defaults:
                    current_section = section_name
                    continue
                else:
                    raise ValueError(f"Unknown section: {section_name}")

            # Parse key-value pairs
            if '=' in line:
                key, value = line.split('=', 1)
                key = key.strip().strip('"').strip("'")
                value = value.split('#')[0].strip()

                try:
                    # Convert to appropriate type
                    if current_section == 'Atoms':
                        value = [self.convert_to_number(v) for v in value.split()]
                    elif '.' in value or 'e' in value.lower():
                        value = float(value)
                    else:
                        value = int(value)
                except ValueError:
                    pass  # Keep as string if conversion fails

                self.data[current_section][key] = value

        # Update defaults with parsed data
        for section, values in self.data.items():
            self.defaults[section].update(values)

    def convert_to_number(self, value):
        try:
            return int(value)
        except ValueError:
            return float(value)

    # Accessor functions
    def get_time_step(self):
        return self.defaults['Simulation']['time_step']

    def get_n_steps(self):
        return self.defaults['Simulation']['n_steps']

    def get_ensemble(self):
        return self.defaults['Simulation']['ensemble']

    def get_box_length(self):
        return self.defaults['System']['box_length']

    def get_n_particles(self):
        return self.defaults['System']['n_particles']

    def get_boundary(self):
        return self.defaults['System']['boundary']

    def get_potential_type(self):
        return self.defaults['Potential']['potential_type']

    def get_epsilon(self):
        return self.defaults['Potential']['epsilon']

    def get_sigma(self):
        return self.defaults['Potential']['sigma']

    def get_cutoff(self):
        return self.defaults['Potential']['cutoff']

    def get_coulomb_k(self):
        return self.defaults['Potential']['coulomb_k']

    def get_output_frequency(self):
        return self.defaults['Output']['output_frequency']

    def get_output_file(self):
        return self.defaults['Output']['output_file']

# Example usage:
# reader = MDInputReader("input_deck.txt")
# print(reader.get_time_step())

