class MDInputReader:
    def __init__(self, filename):
        """Initialize the MDInputReader with a filename."""
        self.filename = filename
        self.config = self._get_default_config()  # Load default values

    def read(self):
        """Reads the input file and populates the configuration dictionary."""
        current_section = None

        try:
            with open(self.filename, 'r') as file:
                for line in file:
                    line = line.strip()

                    # Ignore empty lines and full-line comments
                    if not line or line.startswith("#"):
                        continue

                    # Check if the line defines a new section
                    if line.startswith("[") and line.endswith("]"):
                        current_section = line[1:-1].strip()
                        if current_section not in self.config:
                            self.config[current_section] = {}
                        continue

                    # Process key-value pairs
                    if "=" in line and current_section is not None:
                        key, value = line.split("=", 1)
                        key = key.strip().strip('"')  # Remove quotes around keys if any
                        value = value.split("#", 1)[0].strip()  # Remove inline comments

                        # Special handling for [Atoms] section
                        if current_section == "Atoms":
                            values = value.split()
                            parsed_values = []
                            for v in values:
                                try:
                                    parsed_values.append(int(v) if '.' not in v else float(v))
                                except ValueError:
                                    parsed_values.append(v)
                            self.config[current_section][key] = parsed_values
                        else:
                            # Attempt to convert the value to int or float
                            try:
                                converted_value = int(value) if '.' not in value else float(value)
                            except ValueError:
                                converted_value = value  # Keep as string if conversion fails

                            self.config[current_section][key] = converted_value

        except FileNotFoundError:
            print(f"Warning: File '{self.filename}' not found. Using default configuration.")

        return self.config

    def _get_default_config(self):
        """Returns the default configuration dictionary."""
        return {
            "Simulation": {
                "time_step": 0.00001,
                "n_steps": 20000,
                "ensemble": "NVE",
            },
            "System": {
                "box_length": 10.0,
                "n_particles": 20,
                "boundary": "reflecting",
            },
            "Potential": {
                "potential_type": "Lennard-Jones",
                "epsilon": 30.0,
                "sigma": 2.5,
                "cutoff": 2.5,
                "coulomb_k": 9e9,
            },
            "Output": {
                "output_frequency": 200,
                "output_file": "trajectory.xyz",
            },
        }

    # Accessor methods for default values
    def get_time_step(self):
        return self.config["Simulation"]["time_step"]

    def get_n_steps(self):
        return self.config["Simulation"]["n_steps"]

    def get_ensemble(self):
        return self.config["Simulation"]["ensemble"]

    def get_box_length(self):
        return self.config["System"]["box_length"]

    def get_n_particles(self):
        return self.config["System"]["n_particles"]

    def get_boundary(self):
        return self.config["System"]["boundary"]

    def get_potential_type(self):
        return self.config["Potential"]["potential_type"]

    def get_epsilon(self):
        return self.config["Potential"]["epsilon"]

    def get_sigma(self):
        return self.config["Potential"]["sigma"]

    def get_cutoff(self):
        return self.config["Potential"]["cutoff"]

    def get_coulomb_k(self):
        return self.config["Potential"]["coulomb_k"]

    def get_output_frequency(self):
        return self.config["Output"]["output_frequency"]

    def get_output_file(self):
        return self.config["Output"]["output_file"]
