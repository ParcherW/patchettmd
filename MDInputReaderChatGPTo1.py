class MDInputReader:
    """
    A minimal molecular dynamics input deck reader that:
      - Takes a filename in its constructor.
      - Has a read() method to parse the file and populate nested dictionaries.
      - Includes default sections ([Simulation], [System], [Potential], [Output])
        with guaranteed default key–value pairs.
      - Stores the final configuration in self.config, a nested dictionary:
          {
            "<section>": {
               "<key>": <value>,
               ...
            },
            ...
          }
      - Has dedicated accessor methods for each default parameter (e.g., get_time_step()).
      - Handles comments (lines starting with '#', or '#' inline).
      - Handles quoted keys (e.g., "coulomb_k"), removing the surrounding quotes.
      - Special-cases the [Atoms] section so that values are parsed as lists of numbers.
        Each token in the value is attempted to be converted to int, then float,
        and if that fails, remains a string.
      - Elsewhere, each key–value pair is a single token. Attempts int, then float,
        and finally falls back to string if parsing fails.
      - If the user does not supply any input file or the file is missing keys,
        the defaults are still provided in the config.
    """

    def __init__(self, filename):
        """
        Constructor takes the name of the input deck file. The file is not read here.
        """
        self.filename = filename

        # Preload default configuration
        self.config = {
            "Simulation": {
                "time_step": 0.00001,  # float
                "n_steps": 20000,      # int
                "ensemble": "NVE"      # string
            },
            "System": {
                "box_length": 10.0,   # float
                "n_particles": 20,    # int
                "boundary": "reflecting"  # string
            },
            "Potential": {
                "potential_type": "Lennard-Jones",  # string
                "epsilon": 30.0,                   # float
                "sigma": 2.5,                      # float
                "cutoff": 2.5,                     # float
                "coulomb_k": 9e9                   # float (9.0 x 10^9)
            },
            "Output": {
                "output_frequency": 200,   # int
                "output_file": "trajectory.xyz"  # string
            }
        }

    def read(self):
        """
        Reads the input deck file specified by self.filename.
        Returns the nested dictionary containing the final configuration.

        Sections are denoted by [SectionName].
        Key–value pairs are in the form: key = value
        Comments start with '#' and may appear inline or as an entire line.

        Special handling of the [Atoms] section:
          - Values are stored as lists of numbers (int or float where possible).
        Elsewhere:
          - Values are stored as a single int, float, or string (in that priority).
        """
        current_section = None

        try:
            with open(self.filename, 'r') as file_obj:
                for raw_line in file_obj:
                    line = raw_line.strip()
                    # Skip empty lines
                    if not line:
                        continue
                    # Skip full-line comments
                    if line.startswith('#'):
                        continue

                    # Strip inline comments
                    hash_index = line.find('#')
                    if hash_index != -1:
                        line = line[:hash_index].strip()

                    # Identify section headers: [SectionName]
                    if line.startswith('[') and line.endswith(']'):
                        section_name = line[1:-1].strip()
                        current_section = section_name

                        # Create section in config if it doesn't exist yet
                        if current_section not in self.config:
                            self.config[current_section] = {}
                        continue

                    # Handle key=value pairs
                    if '=' in line:
                        key_part, value_part = line.split('=', 1)
                        key = key_part.strip()
                        value_str = value_part.strip()

                        # Remove surrounding quotes from the key if present
                        if ((key.startswith('"') and key.endswith('"')) or
                                (key.startswith("'") and key.endswith("'"))):
                            key = key[1:-1]

                        # Special handling for the [Atoms] section
                        if current_section == "Atoms":
                            # Split the value into tokens
                            tokens = value_str.split()
                            parsed_list = []
                            for t in tokens:
                                # Try parsing as int, then float, else string
                                try:
                                    parsed_list.append(int(t))
                                except ValueError:
                                    try:
                                        parsed_list.append(float(t))
                                    except ValueError:
                                        parsed_list.append(t)
                            self.config[current_section][key] = parsed_list

                        else:
                            # Attempt int, then float, else string
                            parsed_value = None
                            try:
                                parsed_value = int(value_str)
                            except ValueError:
                                try:
                                    parsed_value = float(value_str)
                                except ValueError:
                                    parsed_value = value_str

                            self.config[current_section][key] = parsed_value

        except FileNotFoundError:
            # If file is not found, the config remains at defaults only
            pass

        return self.config

    # ----------------------------------------------------------------
    # Accessor methods for the default parameters
    # (Each accessor returns the parameter from the config dictionary.)
    # ----------------------------------------------------------------

    # [Simulation] section
    def get_time_step(self):
        return self.config["Simulation"]["time_step"]

    def get_n_steps(self):
        return self.config["Simulation"]["n_steps"]

    def get_ensemble(self):
        return self.config["Simulation"]["ensemble"]

    # [System] section
    def get_box_length(self):
        return self.config["System"]["box_length"]

    def get_n_particles(self):
        return self.config["System"]["n_particles"]

    def get_boundary(self):
        return self.config["System"]["boundary"]

    # [Potential] section
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

    # [Output] section
    def get_output_frequency(self):
        return self.config["Output"]["output_frequency"]

    def get_output_file(self):
        return self.config["Output"]["output_file"]
