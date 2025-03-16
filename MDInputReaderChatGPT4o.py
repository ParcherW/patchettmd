import re
from typing import Union, List, Dict

class MDInputReader:
    """
    A class to parse and store Molecular Dynamics (MD) input parameters from an input deck file or string.
    """

    # Default parameters
    DEFAULTS = {
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

    def __init__(self, input_source: Union[str, None] = None, from_file: bool = True):
        """
        Initializes the MDInputReader by parsing an input file or raw string.

        :param input_source: Path to the input deck file (if from_file=True) or raw input string.
        :param from_file: Whether the input_source is a file path (True) or a raw string (False).
        """
        self.parameters = {section: values.copy() for section, values in self.DEFAULTS.items()}
        self.atoms_data = {}  # Special case for the [Atoms] section

        if input_source:
            if from_file:
                with open(input_source, "r") as f:
                    content = f.read()
            else:
                content = input_source
            self._parse_input(content)

    def _parse_input(self, content: str):
        """
        Parses the input deck and stores parameters.
        """
        current_section = None
        for line in content.splitlines():
            line = line.strip()

            # Ignore blank lines or comments
            if not line or line.startswith("#"):
                continue

            # Section header detection: [SectionName]
            section_match = re.match(r"\[(.+?)\]", line)
            if section_match:
                current_section = section_match.group(1)
                if current_section not in self.parameters and current_section != "Atoms":
                    self.parameters[current_section] = {}
                continue

            # Key-value parsing (handling inline comments)
            key_value_match = re.match(r'(["\w]+)\s*=\s*(.+?)(\s*#.*)?$', line)
            if key_value_match and current_section:
                key = key_value_match.group(1).strip().strip('"')  # Remove quotes around keys
                value = key_value_match.group(2).strip()

                # Special handling for the [Atoms] section (list of numbers)
                if current_section == "Atoms":
                    self.atoms_data[key] = self._parse_list_of_numbers(value)
                else:
                    self.parameters[current_section][key] = self._convert_value(value)

    def _convert_value(self, value: str):
        """
        Converts a string value to int, float, or keeps it as a string.
        """
        if value.isdigit():
            return int(value)
        try:
            return float(value)
        except ValueError:
            return value  # Keep as a string if conversion fails

    def _parse_list_of_numbers(self, value: str) -> List[Union[int, float]]:
        """
        Converts a space-separated string of numbers into a list of int/float values.
        """
        numbers = []
        for item in value.split():
            if item.isdigit():
                numbers.append(int(item))
            else:
                try:
                    numbers.append(float(item))
                except ValueError:
                    pass  # Ignore non-numeric values
        return numbers

    def get(self, section: str, key: str, default=None):
        """
        Retrieves a parameter value, falling back to a default if missing.
        """
        return self.parameters.get(section, {}).get(key, default)

    # Accessor methods for default parameters
    def get_time_step(self) -> float:
        return self.get("Simulation", "time_step", 0.00001)

    def get_n_steps(self) -> int:
        return self.get("Simulation", "n_steps", 20000)

    def get_ensemble(self) -> str:
        return self.get("Simulation", "ensemble", "NVE")

    def get_box_length(self) -> float:
        return self.get("System", "box_length", 10.0)

    def get_n_particles(self) -> int:
        return self.get("System", "n_particles", 20)

    def get_boundary(self) -> str:
        return self.get("System", "boundary", "reflecting")

    def get_potential_type(self) -> str:
        return self.get("Potential", "potential_type", "Lennard-Jones")

    def get_epsilon(self) -> float:
        return self.get("Potential", "epsilon", 30.0)

    def get_sigma(self) -> float:
        return self.get("Potential", "sigma", 2.5)

    def get_cutoff(self) -> float:
        return self.get("Potential", "cutoff", 2.5)

    def get_coulomb_k(self) -> float:
        return self.get("Potential", "coulomb_k", 9e9)

    def get_output_frequency(self) -> int:
        return self.get("Output", "output_frequency", 200)

    def get_output_file(self) -> str:
        return self.get("Output", "output_file", "trajectory.xyz")

    def get_atoms_data(self) -> Dict[str, List[Union[int, float]]]:
        """
        Returns the parsed atoms data as a dictionary of lists of numbers.
        """
        return self.atoms_data

