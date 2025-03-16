class MDInputReader:
    def __init__(self, input_source):
        """
        Initialize the MDInputReader.
        :param input_source: A string that is either a file path to the input deck or the raw input deck content.
        """
        self.input_source = input_source
        self.config = {}
        self._load_defaults()

    def _load_defaults(self):
        """
        Preload the default configuration.
        """
        self.config["Simulation"] = {
            "time_step": 0.00001,   # float
            "n_steps": 20000,       # int
            "ensemble": "NVE"       # string
        }
        self.config["System"] = {
            "box_length": 10.0,     # float
            "n_particles": 20,      # int
            "boundary": "reflecting"  # string ("periodic" or "reflecting")
        }
        self.config["Potential"] = {
            "potential_type": "Lennard-Jones",  # string
            "epsilon": 30.0,       # float
            "sigma": 2.5,          # float
            "cutoff": 2.5,         # float
            "coulomb_k": 9e9       # float
        }
        self.config["Output"] = {
            "output_frequency": 200,      # int
            "output_file": "trajectory.xyz"  # string
        }
        # Prepare an empty Atoms section for special handling (if provided in the input)
        self.config["Atoms"] = {}

    def read(self):
        """
        Read and parse the input deck. Returns the config nested dictionary.
        The input deck is assumed to be structured into sections [SectionName]
        with key = value pairs. Comments start with '#' and may appear full-line or inline.
        """
        # Attempt to open the file using input_source as a filename;
        # if that fails, treat input_source as raw string input.
        try:
            with open(self.input_source, "r") as f:
                lines = f.readlines()
        except Exception:
            lines = self.input_source.splitlines()

        current_section = None

        for line in lines:
            line = line.strip()
            if not line or line.startswith("#"):
                continue  # Skip empty or comment-only lines

            # Check for a section header (e.g., [SectionName])
            if line.startswith("[") and line.endswith("]"):
                section_name = line[1:-1].strip()
                if section_name not in self.config:
                    self.config[section_name] = {}
                current_section = section_name
                continue

            # Remove inline comments (everything after '#')
            if "#" in line:
                line = line.split("#", 1)[0].strip()
            if not line:
                continue

            # Expecting a key = value pair
            if "=" in line:
                parts = line.split("=", 1)
                key = parts[0].strip()
                value = parts[1].strip()

                # Remove quotes around key if present
                if (key.startswith('"') and key.endswith('"')) or (key.startswith("'") and key.endswith("'")):
                    key = key[1:-1]

                # Special handling for the Atoms section: parse a list of numbers.
                if current_section == "Atoms":
                    tokens = value.split()
                    parsed_tokens = []
                    for token in tokens:
                        # Try converting to int; if that fails, try float.
                        try:
                            num = int(token)
                        except ValueError:
                            try:
                                num = float(token)
                            except ValueError:
                                num = token
                        parsed_tokens.append(num)
                    value = parsed_tokens
                else:
                    # For other sections, if the key exists in defaults, convert to the same type.
                    if current_section in self.config and key in self.config[current_section]:
                        default_val = self.config[current_section][key]
                        if isinstance(default_val, int):
                            try:
                                value = int(value)
                            except ValueError:
                                pass
                        elif isinstance(default_val, float):
                            try:
                                value = float(value)
                            except ValueError:
                                pass
                        # If the default value is a string, leave the value as is.

                # Update the configuration dictionary.
                if current_section is not None:
                    self.config[current_section][key] = value

        return self.config

    # Accessor functions for the Simulation section
    def get_time_step(self):
        """Return the simulation time step (float)."""
        return float(self.config["Simulation"].get("time_step", 0.00001))

    def get_n_steps(self):
        """Return the total number of simulation steps (int)."""
        return int(self.config["Simulation"].get("n_steps", 20000))

    def get_ensemble(self):
        """Return the simulation ensemble (string)."""
        return str(self.config["Simulation"].get("ensemble", "NVE"))

    # Accessor functions for the System section
    def get_box_length(self):
        """Return the simulation box length (float)."""
        return float(self.config["System"].get("box_length", 10.0))

    def get_n_particles(self):
        """Return the number of particles in the simulation (int)."""
        return int(self.config["System"].get("n_particles", 20))

    def get_boundary(self):
        """Return the type of boundary conditions (string)."""
        return str(self.config["System"].get("boundary", "reflecting"))

    # Accessor functions for the Potential section
    def get_potential_type(self):
        """Return the potential type (string)."""
        return str(self.config["Potential"].get("potential_type", "Lennard-Jones"))

    def get_epsilon(self):
        """Return the epsilon parameter for the potential (float)."""
        return float(self.config["Potential"].get("epsilon", 30.0))

    def get_sigma(self):
        """Return the sigma parameter for the potential (float)."""
        return float(self.config["Potential"].get("sigma", 2.5))

    def get_cutoff(self):
        """Return the cutoff parameter for the potential (float)."""
        return float(self.config["Potential"].get("cutoff", 2.5))

    def get_coulomb_k(self):
        """Return the Coulomb constant parameter (float)."""
        return float(self.config["Potential"].get("coulomb_k", 9e9))

    # Accessor functions for the Output section
    def get_output_frequency(self):
        """Return the output frequency (int)."""
        return int(self.config["Output"].get("output_frequency", 200))

    def get_output_file(self):
        """Return the output file name (string)."""
        return str(self.config["Output"].get("output_file", "trajectory.xyz"))
