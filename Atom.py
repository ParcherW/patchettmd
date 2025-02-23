import numpy as np

# Define an Atom class to store properties for each particle.
class Atom:
    def __init__(self, position, velocity, atom_type=0, mass=1.0, charge=1.6e-19):
        self.position = np.array(position, dtype=float)  # 3D position
        self.velocity = np.array(velocity, dtype=float)  # 3D velocity
        self.atom_type = atom_type
        self.mass = mass
        self.charge = charge
        # Bond stretching between two atoms.
class BondStretch:
    def __init__(self, atom1, atom2, r0, k_stretch):
        self.atom1 = atom1  # index of first atom
        self.atom2 = atom2  # index of second atom
        self.r0 = r0        # equilibrium bond length
        self.k = k_stretch  # force constant

# Angle bending among three atoms (with atom2 as the central atom).
class AngleBend:
    def __init__(self, atom1, atom2, atom3, theta0, k_bend):
        self.atom1 = atom1  # first atom
        self.atom2 = atom2  # central atom
        self.atom3 = atom3  # third atom
        self.theta0 = theta0  # equilibrium angle (radians)
        self.k = k_bend       # force constant

# Torsion (dihedral) interaction among four atoms.
class Torsion:
    def __init__(self, atom1, atom2, atom3, atom4, V, n, gamma):
        self.atom1 = atom1  # first atom
        self.atom2 = atom2  # second atom
        self.atom3 = atom3  # third atom
        self.atom4 = atom4  # fourth atom
        self.V = V          # amplitude (energy barrier)
        self.n = n          # periodicity
        self.gamma = gamma  # phase offset (radians)




        