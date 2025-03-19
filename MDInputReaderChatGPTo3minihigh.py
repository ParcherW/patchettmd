class MDInputReader:
    def __init__(self, filename):
        """
        Initialize the MDInputReader with a filename.
        The configuration dictionary is preloaded with defaults.
        
        Parameters:
            filename (str): Name of the input deck file to be parsed.
        """
        self.filename = filename
        # Nested dictionary to store configuration by section.
        self.config = {}
        self._load_defaults()

    def _load_defaults(self):
        """
        Preload the config dictionary with default sections and key–value pairs.
        These defaults are provided even if the user does not specify them
        in the input deck.
        """
        self.config["Simulation"] = {
            "time_step": 0.00001,   # float: time step in simulation time units
            "n_steps": 20000,       # int: total number of simulation steps
            "ensemble": "NVE"       # string: simulation ensemble
        }
        self.config["System"] = {
            "box_length": 10.0,     # float: simulation box length (cubic box)
            "n_particles": 20,      # int: number of particles in the simulation
            "boundary": "reflecting"  # string: either "periodic" or "reflecting"
        }
        self.config["Potential"] = {
            "potential_type": "Lennard-Jones",  # string: type of potential
            "epsilon": 30.0,    # float: parameter epsilon
            "sigma": 2.5,       # float: parameter sigma
            "cutoff": 2.5,      # float: cutoff distance
            "coulomb_k": 9e9    # float: Coulomb constant (note the quoted key in input decks)
        }
        self.config["Output"] = {
            "output_frequency": 200,    # int: how often to write simulation data
            "output_file": "trajectory.xyz"  # string: output filename
        }
        # The Atoms section is optional and will be created if encountered.

    def read(self):
        """
        Read and parse the input deck file.
        
        The file is expected to use sections defined by square brackets [SectionName],
        with key–value pairs separated by an '='. Full-line comments or inline comments
        start with a hash (#) and are ignored.
        
        Special handling:
          - In the [Atoms] section, values that are whitespace-separated lists
            of numbers are parsed as lists of ints/floats.
          - For all other sections, an attempt is made to convert values to int or float.
            If conversion fails, the value is stored as a string.
          - Quoted keys (e.g., "coulomb_k") have their quotes removed.
          - Whitespace is trimmed from keys and values.
        
        If the file cannot be read, only the defaults remain in the config.
        
        Returns:
            dict: The nested configuration dictionary.
        """
        try:
            f = open(self.filename, 'r')
            lines = f.readlines()
            f.close()
        except Exception:
            # If the file is not found or cannot be read, return defaults only.
            return self.config

        current_section = None
        for line in lines:
            line = line.strip()
            # Skip empty lines and full-line comments.
            if not line or line.startswith("#"):
                continue

            # Identify section headers [SectionName]
            if line.startswith('[') and line.endswith(']'):
                current_section = line[1:-1].strip()
                if current_section not in self.config:
                    self.config[current_section] = {}
                continue

            # Remove inline comments (everything after '#')
            comment_index = line.find("#")
            if comment_index != -1:
                line = line[:comment_index].strip()
            if not line:
                continue

            # Process key-value pair, split only on the first '='.
            if '=' in line:
                parts = line.split('=', 1)
                key = parts[0].strip()
                value_str = parts[1].strip()

                # Remove quotes around the key if present.
                if key.startswith('"') and key.endswith('"'):
                    key = key[1:-1]

                # Special handling for the [Atoms] section:
                if current_section == "Atoms":
                    # Split by whitespace and try to convert each token to int or float.
                    tokens = value_str.split()
                    converted_list = []
                    for token in tokens:
                        try:
                            # Determine whether to use int or float.
                            if '.' in token or 'e' in token.lower():
                                converted_list.append(float(token))
                            else:
                                converted_list.append(int(token))
                        except Exception:
                            # If conversion fails, keep token as string.
                            converted_list.append(token)
                    self.config[current_section][key] = converted_list
                else:
                    # Attempt conversion: try int, then float, else retain string.
                    try:
                        converted = int(value_str)
                    except Exception:
                        try:
                            converted = float(value_str)
                        except Exception:
                            converted = value_str
                    self.config[current_section][key] = converted

        return self.config

    # Accessor Functions for Simulation section
    def get_time_step(self):
        """Return the simulation time step (float)."""
        return self.config.get("Simulation", {}).get("time_step")

    def get_n_steps(self):
        """Return the total number of simulation steps (int)."""
        return self.config.get("Simulation", {}).get("n_steps")

    def get_ensemble(self):
        """Return the simulation ensemble (string)."""
        return self.config.get("Simulation", {}).get("ensemble")

    # Accessor Functions for System section
    def get_box_length(self):
        """Return the simulation box length (float)."""
        return self.config.get("System", {}).get("box_length")

    def get_n_particles(self):
        """Return the number of particles (int)."""
        return self.config.get("System", {}).get("n_particles")

    def get_boundary(self):
        """Return the boundary type (string)."""
        return self.config.get("System", {}).get("boundary")

    # Accessor Functions for Potential section
    def get_potential_type(self):
        """Return the potential type (string)."""
        return self.config.get("Potential", {}).get("potential_type")

    def get_epsilon(self):
        """Return the epsilon parameter (float)."""
        return self.config.get("Potential", {}).get("epsilon")

    def get_sigma(self):
        """Return the sigma parameter (float)."""
        return self.config.get("Potential", {}).get("sigma")

    def get_cutoff(self):
        """Return the cutoff distance (float)."""
        return self.config.get("Potential", {}).get("cutoff")

    def get_coulomb_k(self):
        """Return the Coulomb constant (float)."""
        return self.config.get("Potential", {}).get("coulomb_k")

    # Accessor Functions for Output section
    def get_output_frequency(self):
        """Return the output frequency (int)."""
        return self.config.get("Output", {}).get("output_frequency")

    def get_output_file(self):
        """Return the output file name (string)."""
        return self.config.get("Output", {}).get("output_file")

