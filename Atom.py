import numpy as np

# Define an Atom class to store properties for each particle.
class Atom:
    def __init__(self, position, velocity, atom_type=0, mass=1.0, charge=1.0):
        self.position = np.array(position, dtype=float)  # 3D position
        self.velocity = np.array(velocity, dtype=float)  # 3D velocity
        self.atom_type = atom_type
        self.mass = mass
        self.charge = charge