class MDInputReader:
    """
    A simple reader for a molecular dynamics simulation input deck.
    The input deck is assumed to be structured in sections like:
    
        [SectionName]
        key = value
        
    Lines starting with '#' are treated as comments.
    """
    
    def __init__(self, filename):
        self.filename = filename
        self.sections = {}
    
    def _convert_value(self, value):
        """
        Attempt to convert a string value to an int or float.
        If conversion fails, return the string unchanged.
        """
        try:
            # Try integer conversion first.
            return int(value)
        except ValueError:
            try:
                # If not integer, try float conversion.
                return float(value)
            except ValueError:
                # Return the original string if both conversions fail.
                return value
    
    def read(self):
        """
        Read the input deck file and return a dictionary with sections and key-value pairs.
        """
        current_section = None

        try:
            with open(self.filename, 'r') as f:
                for line in f:
                    line = line.strip()

                    # Ignore empty lines or full-line comments
                    if not line or line.startswith('#'):
                        continue

                    # Check for a section header (e.g., [Simulation])
                    if line.startswith('[') and line.endswith(']'):
                        current_section = line[1:-1].strip()
                        self.sections[current_section] = {}

                    else:
                        # Process key-value pairs (expects the format key = value)
                        if '=' in line:
                            key, value = line.split('=', 1)
                            key = key.strip()

                            # Remove inline comments
                            value = value.split('#', 1)[0].strip()

                            # Convert value
                            value = self._convert_value(value)

                            if current_section is None:
                                raise ValueError("Key-value pair found outside any section: " + line)

                            self.sections[current_section][key] = value
                        else:
                            raise ValueError("Line not in key-value format: " + line)

        except FileNotFoundError:
            raise FileNotFoundError(f"Input file '{self.filename}' not found.")

        return self.sections


# Example usage:
if __name__ == "__main__":
    # Replace 'md_input.txt' with the path to your actual input file.
    reader = MDInputReader("md_input.txt")
    config = reader.read()
    
    # Print out the parsed configuration.
    for section, params in config.items():
        print(f"Section: {section}")
        for key, value in params.items():
            print(f"  {key}: {value}")

