class MDInputReader:
    """
    A reader for a molecular dynamics simulation input deck.
    
    The input deck is structured in sections like:
        [SectionName]
        key = value  # optional inline comment

    Lines starting with '#' are treated as comments.
    Defaults are set to ensure required parameters always exist.
    """

    # Define guaranteed default parameters
    DEFAULTS = {
        "Simulation": {
            "time_step": 0.0001,
            "n_steps": 10000,
            "ensemble": "NVE"
        },
        "System": {
            "box_length": 30.0,
            "n_particles": 40,
            "boundary": "reflecting"
        },
        "Potential": {
            "potential_type": "Lennard-Jones",
            "epsilon": 2.0,
            "sigma": 1.0,
            "cutoff": 2.5,
            "coulomb_k":9e9
        },
        "Output": {
            "output_frequency": 200,
            "output_file": "trajectory.xyz"
        }
    }

    def __init__(self, filename):
        """
        Initialize the MDInputReader and set up data structures.
        """
        self.filename = filename
        self.sections = {section: defaults.copy() for section, defaults in self.DEFAULTS.items()}

    def _convert_value(self, value):
        """
        Convert a string value to int, float, or return as string if conversion fails.
        """
        try:
            return int(value)
        except ValueError:
            try:
                return float(value)
            except ValueError:
                return value.strip()  # Keep as string if not numeric

    def read(self):
        """
        Parse the input deck file, updating sections with user-defined values.
        """
        current_section = None

        try:
            with open(self.filename, 'r') as f:
                for line in f:
                    line = line.strip()

                    # Ignore empty lines or comments
                    if not line or line.startswith('#'):
                        continue

                    # Identify section headers (e.g., [Simulation])
                    if line.startswith('[') and line.endswith(']'):
                        current_section = line[1:-1].strip()
                        if current_section not in self.sections:
                            self.sections[current_section] = {}  # Allow custom sections

                    else:
                        # Process key-value pairs
                        if '=' in line:
                            key, value = line.split('=', 1)
                            key = key.strip()
                            value = value.split('#', 1)[0].strip()  # Remove inline comments

                            # Convert and store the value
                            value = self._convert_value(value)

                            if current_section is None:
                                raise ValueError("Key-value pair found outside any section: " + line)

                            self.sections[current_section][key] = value
                        else:
                            raise ValueError("Line not in key-value format: " + line)

        except FileNotFoundError:
            raise FileNotFoundError(f"Input file '{self.filename}' not found.")

        return self.sections

    def get(self, section, key):
        """
        Retrieve a parameter from a specific section, returning a default if not explicitly set.
        """
        return self.sections.get(section, {}).get(key, None)

    # Individual accessors for convenience
    def get_time_step(self):
        return self.get("Simulation", "time_step")

    def get_n_steps(self):
        return self.get("Simulation", "n_steps")

    def get_ensemble(self):
        return self.get("Simulation", "ensemble")

    def get_box_length(self):
        return self.get("System", "box_length")

    def get_n_particles(self):
        return self.get("System", "n_particles")

    def get_boundary(self):
        return self.get("System", "boundary")

    def get_potential_type(self):
        return self.get("Potential", "potential_type")

    def get_epsilon(self):
        return self.get("Potential", "epsilon")

    def get_sigma(self):
        return self.get("Potential", "sigma")

    def get_cutoff(self):
        return self.get("Potential", "cutoff")

    def get_output_frequency(self):
        return self.get("Output", "output_frequency")

    def get_output_file(self):
        return self.get("Output", "output_file")
