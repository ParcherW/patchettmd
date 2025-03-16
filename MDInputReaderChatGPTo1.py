import copy
import re

class MDInputReader:
    """
    A class to read and store parameters from an MD input deck.

    The input deck is assumed to have sections like:
        [SectionName]
        key = value     # possible inline comment
    Comments can appear as full lines (starting with #) or inline.
    Keys may be in quotes ("key_name"), and values may be numeric or strings.

    Special handling for the [Atoms] section:
        If the value is a whitespace-separated list of numbers, it will be
        parsed and stored as a list of floats.
    """

    def __init__(self):
        """
        Initialize the reader with default values for:
          - [Simulation]
          - [System]
          - [Potential]
          - [Output]
        """
        self._defaults = {
            "Simulation": {
                "time_step":      "0.00001",
                "n_steps":        "20000",
                "ensemble":       "NVE",
            },
            "System": {
                "box_length":     "10.0",
                "n_particles":    "20",
                "boundary":       "reflecting",  # "periodic" or "reflecting"
            },
            "Potential": {
                "potential_type": "Lennard-Jones",
                "epsilon":        "30.0",
                "sigma":          "2.5",
                "cutoff":         "2.5",
                "coulomb_k":      "9e9",
            },
            "Output": {
                "output_frequency": "200",
                "output_file":      "trajectory.xyz",
            },
            # You can add other default sections/keys as needed
        }

        # This is where the final parsed data will be stored.
        # Start with a deep copy of the defaults so that
        # the defaults are always available.
        self._data = copy.deepcopy(self._defaults)

    def parse_file(self, file_path):
        """
        Parse an input deck from a file.
        """
        with open(file_path, 'r') as f:
            content = f.read()
        self.parse_string(content)

    def parse_string(self, input_str):
        """
        Parse an input deck from a raw string.
        """
        current_section = None
        lines = input_str.splitlines()

        for line in lines:
            # Strip leading/trailing whitespace
            line = line.strip()
            # Skip empty lines and full-line comments
            if not line or line.startswith('#'):
                continue

            # Identify new section if line starts with '[' and ends with ']'
            if line.startswith('[') and line.endswith(']'):
                current_section = line[1:-1].strip()
                # Ensure we have a place in self._data for this section
                if current_section not in self._data:
                    self._data[current_section] = {}
                continue

            # Handle key-value lines
            # Remove inline comments by splitting at the first '#'
            # (assuming '#' is not quoted or escaped)
            inline_comment_index = line.find('#')
            if inline_comment_index != -1:
                line = line[:inline_comment_index].strip()

            # Split by '=' into key and value
            if '=' in line:
                key, value = line.split('=', 1)
                key = key.strip()
                value = value.strip()

                # Remove surrounding quotes around the key, if present
                # e.g., "coulomb_k" -> coulomb_k
                if (key.startswith('"') and key.endswith('"')) or \
                   (key.startswith("'") and key.endswith("'")):
                    key = key[1:-1].strip()

                # Special handling for [Atoms] section:
                # parse the value as a list of numbers if possible.
                if current_section == "Atoms":
                    # Attempt to parse as a list of floats
                    # e.g., "1.0 2.0 3.0" -> [1.0, 2.0, 3.0]
                    value_parts = value.split()
                    if len(value_parts) > 1:
                        # If there's more than one token, parse them as floats
                        try:
                            float_list = [float(v) for v in value_parts]
                            value = float_list
                        except ValueError:
                            # Fallback: store as string if parse fails
                            pass
                    else:
                        # Single token: still might be numeric
                        try:
                            value_f = float(value)
                            value = [value_f]
                        except ValueError:
                            # Non-numeric single token -> store as string
                            pass
                else:
                    # For other sections, store raw string if not recognized as numeric
                    # but we can do basic checks if needed:
                    # If you'd like numeric interpretation for known keys, do it here
                    pass

                # Store the value in the data dictionary
                if current_section not in self._data:
                    self._data[current_section] = {}
                self._data[current_section][key] = value

    # --------------------------------------------------------------
    # Accessors for default parameters
    # --------------------------------------------------------------

    # [Simulation] accessors
    def get_time_step(self) -> float:
        return float(self._data["Simulation"].get("time_step", self._defaults["Simulation"]["time_step"]))

    def get_n_steps(self) -> int:
        return int(self._data["Simulation"].get("n_steps", self._defaults["Simulation"]["n_steps"]))

    def get_ensemble(self) -> str:
        return str(self._data["Simulation"].get("ensemble", self._defaults["Simulation"]["ensemble"]))

    # [System] accessors
    def get_box_length(self) -> float:
        return float(self._data["System"].get("box_length", self._defaults["System"]["box_length"]))

    def get_n_particles(self) -> int:
        return int(self._data["System"].get("n_particles", self._defaults["System"]["n_particles"]))

    def get_boundary(self) -> str:
        return str(self._data["System"].get("boundary", self._defaults["System"]["boundary"]))

    # [Potential] accessors
    def get_potential_type(self) -> str:
        return str(self._data["Potential"].get("potential_type", self._defaults["Potential"]["potential_type"]))

    def get_epsilon(self) -> float:
        return float(self._data["Potential"].get("epsilon", self._defaults["Potential"]["epsilon"]))

    def get_sigma(self) -> float:
        return float(self._data["Potential"].get("sigma", self._defaults["Potential"]["sigma"]))

    def get_cutoff(self) -> float:
        return float(self._data["Potential"].get("cutoff", self._defaults["Potential"]["cutoff"]))

    def get_coulomb_k(self) -> float:
        return float(self._data["Potential"].get("coulomb_k", self._defaults["Potential"]["coulomb_k"]))

    # [Output] accessors
    def get_output_frequency(self) -> int:
        return int(self._data["Output"].get("output_frequency", self._defaults["Output"]["output_frequency"]))

    def get_output_file(self) -> str:
        return str(self._data["Output"].get("output_file", self._defaults["Output"]["output_file"]))

    # --------------------------------------------------------------
    # Example convenience function for retrieving a value generally.
    # --------------------------------------------------------------
    def get_value(self, section: str, key: str, default=None):
        """
        A generic getter that returns the value stored under [section][key].
        If the section/key is not found, returns `default`.
        """
        if section in self._data and key in self._data[section]:
            return self._data[section][key]
        return default

