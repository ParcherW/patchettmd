class MDInputReader:
    def __init__(self, filename):
        self.filename = filename
        # Define default configuration
        self.config = {
            "Simulation": {
                "time_step": 0.00001,
                "n_steps": 20000,
                "ensemble": "NVE"
            },
            "System": {
                "box_length": 10.0,
                "n_particles": 20,
                "boundary": "reflecting"
            },
            "Potential": {
                "potential_type": "Lennard-Jones",
                "epsilon": 30.0,
                "sigma": 2.5,
                "cutoff": 2.5,
                "coulomb_k": 9e9
            },
            "Output": {
                "output_frequency": 200,
                "output_file": "trajectory.xyz"
            }
        }

    def read(self):
        try:
            with open(self.filename, 'r') as file:
                lines = file.readlines()
        except FileNotFoundError:
            raise Exception(f"File {self.filename} not found.")
        
        current_section = None
        for line in lines:
            stripped_line = line.strip()
            if stripped_line.startswith('#') or not stripped_line:
                continue
            
            if stripped_line.startswith('[') and stripped_line.endswith(']'):
                current_section = stripped_line[1:-1]
                if current_section not in self.config:
                    self.config[current_section] = {}
                continue
            
            if current_section is None:
                raise Exception("Key-value pair found outside of a section.")
            
            key, value = stripped_line.split('=')
            key = key.strip().strip('"')
            value = value.split('#')[0].strip()
            
            # Attempt to convert to float or int
            try:
                converted_value = float(value)
                if converted_value.is_integer():
                    converted_value = int(converted_value)
            except ValueError:
                converted_value = value
            
            # Special handling for specific keys that must be integers
            if current_section == "System" and key == "n_particles":
                try:
                    converted_value = int(value)
                except ValueError:
                    raise Exception(f"Invalid integer value for n_particles: {value}")
            
            # Special handling for Atoms section
            if current_section == "Atoms":
                converted_value = [self._convert_to_number(v) for v in value.split()]
            
            self.config[current_section][key] = converted_value
        
        return self.config

    def _convert_to_number(self, value):
        try:
            return float(value)
        except ValueError:
            try:
                return int(value)
            except ValueError:
                return value

    # Accessor functions for defaults
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
