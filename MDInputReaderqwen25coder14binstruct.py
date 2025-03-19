class MDInputReader:
    def __init__(self, filename=None):
        self.filename = filename
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
        if self.filename:
            with open(self.filename, 'r') as file:
                for line in file:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    if line.startswith('[') and line.endswith(']'):
                        current_section = line[1:-1].strip()
                    else:
                        try:
                            key, value = map(str.strip, line.split('=', 1))
                            self._parse_value(current_section, key, value)
                        except ValueError:
                            continue
        return self.config

    def _parse_value(self, section, key, raw_value):
        if section not in self.config:
            self.config[section] = {}
        
        stripped_value = raw_value.split('#')[0].strip()
        if '"' in stripped_value or "'" in stripped_value:
            value = stripped_value.replace('"', '').replace("'", '')
        else:
            try:
                value = int(stripped_value)
            except ValueError:
                try:
                    value = float(stripped_value)
                except ValueError:
                    value = stripped_value
        
        if section == "Atoms":
            values_list = stripped_value.split()
            try:
                values_list = [float(v) for v in values_list]
            except ValueError:
                pass
            self.config[section][key] = values_list
        else:
            self.config[section][key] = value

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

