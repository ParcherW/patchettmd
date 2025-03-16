class MDInputReader:
    """
    A class to read and parse a molecular dynamics simulation input deck.
    
    The input deck is structured into sections such as [Simulation], [System],
    [Potential], and [Output] with key–value pairs. Comments start with '#'
    and can be inline. The class loads a set of default values and then
    overrides them if the user provides an input deck (either as a raw string or
    via a file). For the [Atoms] section, a whitespace‐separated list of numbers
    is parsed and returned as a list of int/float values.
    """
    
    def __init__(self, input_data=None):
        # Predefined defaults for required sections and keys.
        self.data = {
            "Simulation": {
                "time_step": 0.00001,    # float: time step in time units
                "n_steps": 20000,        # int: total number of simulation steps
                "ensemble": "NVE",       # string: simulation ensemble
            },
            "System": {
                "box_length": 10.0,      # float: simulation box length (cubic)
                "n_particles": 20,       # int: number of particles
                "boundary": "reflecting",# string: "periodic" or "reflecting"
            },
            "Potential": {
                "potential_type": "Lennard-Jones",  # string: potential type
                "epsilon": 30.0,         # float: Lennard-Jones epsilon
                "sigma": 2.5,            # float: Lennard-Jones sigma
                "cutoff": 2.5,           # float: cutoff distance
                "coulomb_k": 9e9,        # float: Coulomb constant
            },
            "Output": {
                "output_frequency": 200, # int: simulation data output frequency
                "output_file": "trajectory.xyz",  # string: output file name
            }
            # Note: The Atoms section is not preloaded with defaults.
        }
        
        # If an input deck is provided, parse it.
        if input_data is not None:
            self._parse(input_data)
    
    def _parse(self, input_text: str):
        """
        Parses the input deck text.
        
        Processes section headers, key–value pairs, handles inline comments,
        quoted keys, and trims whitespace. For the [Atoms] section, converts
        whitespace-separated numbers into a list of numeric types.
        """
        current_section = None
        
        for line in input_text.splitlines():
            # Remove leading/trailing whitespace.
            line = line.strip()
            # Skip empty lines and full-line comments.
            if not line or line.startswith('#'):
                continue
            
            # Check if the line defines a new section.
            if line.startswith('[') and line.endswith(']'):
                section = line[1:-1].strip()
                # Create section if not already present.
                if section not in self.data:
                    self.data[section] = {}
                current_section = section
                continue
            
            # Remove inline comments (everything after a '#' is ignored).
            if '#' in line:
                line = line.split('#', 1)[0].strip()
            
            # Skip if no '=' is present.
            if '=' not in line:
                continue
            
            key, value = line.split('=', 1)
            key = key.strip().strip('"').strip("'")  # Remove whitespace and any quotes around key.
            value = value.strip()
            
            if current_section == "Atoms":
                # For the Atoms section, assume the value is a whitespace-separated list of numbers.
                numbers = value.split()
                parsed_numbers = []
                for num in numbers:
                    try:
                        if '.' in num or 'e' in num.lower():
                            parsed_numbers.append(float(num))
                        else:
                            parsed_numbers.append(int(num))
                    except ValueError:
                        parsed_numbers.append(num)
                self.data[current_section][key] = parsed_numbers
            else:
                # If the key exists in the defaults, convert to the proper type.
                if current_section in self.data and key in self.data[current_section]:
                    default_value = self.data[current_section][key]
                    if isinstance(default_value, int):
                        try:
                            converted_value = int(value)
                        except ValueError:
                            converted_value = value
                    elif isinstance(default_value, float):
                        try:
                            converted_value = float(value)
                        except ValueError:
                            converted_value = value
                    else:
                        converted_value = value
                    self.data[current_section][key] = converted_value
                else:
                    # If the key is not a default, attempt a conversion (fallback to string if conversion fails).
                    try:
                        if '.' in value or 'e' in value.lower():
                            converted_value = float(value)
                        else:
                            converted_value = int(value)
                        self.data[current_section][key] = converted_value
                    except ValueError:
                        self.data[current_section][key] = value
    
    def load_from_file(self, filepath: str):
        """
        Reads the input deck from a file and parses its contents.
        """
        with open(filepath, "r") as f:
            content = f.read()
        self._parse(content)
    
    # Accessor functions for Simulation parameters
    def get_time_step(self) -> float:
        return self.data["Simulation"]["time_step"]

    def get_n_steps(self) -> int:
        return self.data["Simulation"]["n_steps"]

    def get_ensemble(self) -> str:
        return self.data["Simulation"]["ensemble"]

    # Accessor functions for System parameters
    def get_box_length(self) -> float:
        return self.data["System"]["box_length"]

    def get_n_particles(self) -> int:
        return self.data["System"]["n_particles"]

    def get_boundary(self) -> str:
        return self.data["System"]["boundary"]

    # Accessor functions for Potential parameters
    def get_potential_type(self) -> str:
        return self.data["Potential"]["potential_type"]

    def get_epsilon(self) -> float:
        return self.data["Potential"]["epsilon"]

    def get_sigma(self) -> float:
        return self.data["Potential"]["sigma"]

    def get_cutoff(self) -> float:
        return self.data["Potential"]["cutoff"]

    def get_coulomb_k(self) -> float:
        return self.data["Potential"]["coulomb_k"]

    # Accessor functions for Output parameters
    def get_output_frequency(self) -> int:
        return self.data["Output"]["output_frequency"]

    def get_output_file(self) -> str:
        return self.data["Output"]["output_file"]

    # Additional accessor for the Atoms section, if present.
    def get_atoms(self, key: str):
        """
        Returns the list of numbers for a given key in the Atoms section,
        or None if the key is not found.
        """
        if "Atoms" in self.data and key in self.data["Atoms"]:
            return self.data["Atoms"][key]
        return None


# === Example Usage ===
if __name__ == "__main__":
    # Example input deck as a raw string.
    example_input = """
    [Simulation]
    time_step = 0.00002   # Overriding default time_step
    n_steps = 25000
    ensemble = NVT

    [System]
    box_length = 12.0
    n_particles = 50
    boundary = periodic

    [Potential]
    potential_type = Lennard-Jones
    epsilon = 35.0
    sigma = 2.8
    cutoff = 3.0
    "coulomb_k" = 8e9

    [Output]
    output_frequency = 100
    output_file = new_trajectory.xyz

    [Atoms]
    positions = 1.0 2.0 3.0 4 5 6
    """

    # Create an MDInputReader instance with the raw string input.
    reader = MDInputReader(example_input)

    # Access default and overridden parameters.
    print("Time Step:", reader.get_time_step())
    print("Number of Steps:", reader.get_n_steps())
    print("Ensemble:", reader.get_ensemble())
    print("Box Length:", reader.get_box_length())
    print("Number of Particles:", reader.get_n_particles())
    print("Boundary:", reader.get_boundary())
    print("Potential Type:", reader.get_potential_type())
    print("Epsilon:", reader.get_epsilon())
    print("Sigma:", reader.get_sigma())
    print("Cutoff:", reader.get_cutoff())
    print("Coulomb Constant:", reader.get_coulomb_k())
    print("Output Frequency:", reader.get_output_frequency())
    print("Output File:", reader.get_output_file())
    print("Atom Positions:", reader.get_atoms("positions"))
