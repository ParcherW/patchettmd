import numpy as np
import matplotlib.pyplot as plt

# Define an Atom class to store properties for each particle.
class Atom:
    def __init__(self, position, velocity, atom_type='A', mass=1.0):
        self.position = np.array(position, dtype=float)  # 3D position
        self.velocity = np.array(velocity, dtype=float)  # 3D velocity
        self.atom_type = atom_type                       # e.g., 'A', 'B', etc.
        self.mass = mass

# Define a Bond class for bonded interactions between atoms.
class Bond:
    def __init__(self, atom_index1, atom_index2, r0, k):
        self.atom1 = atom_index1  # index of the first atom in the bond
        self.atom2 = atom_index2  # index of the second atom
        self.r0 = r0              # equilibrium bond length
        self.k = k                # bond force constant

def compute_lj_params(type1, type2, lj_params):
    """
    Determine the mixed Lennard-Jones parameters using Lorentz-Berthelot rules.
    """
    sigma1 = lj_params[type1]['sigma']
    sigma2 = lj_params[type2]['sigma']
    epsilon1 = lj_params[type1]['epsilon']
    epsilon2 = lj_params[type2]['epsilon']
    sigma_ij = 0.5 * (sigma1 + sigma2)
    epsilon_ij = np.sqrt(epsilon1 * epsilon2)
    return epsilon_ij, sigma_ij

def compute_forces(atoms, L, lj_params, bonds=None):
    """
    Compute forces on each atom due to nonbonded Lennard-Jones interactions 
    and bonded (harmonic) interactions if bonds are provided.
    
    Parameters:
      atoms    : List of Atom objects.
      L        : Box length (cubic box of side L).
      lj_params: Dictionary mapping atom types to their LJ parameters.
      bonds    : Optional list of Bond objects for bonded interactions.
    
    Returns:
      forces   : NumPy array of shape (N, 3) with force vectors.
      potential: Total potential energy of the system.
    """
    N = len(atoms)
    forces = np.zeros((N, 3))
    potential = 0.0

    # Nonbonded Lennard-Jones interactions over unique pairs.
    for i in range(N):
        for j in range(i+1, N):
            # Compute displacement vector with periodic boundary conditions
            r_vec = atoms[i].position - atoms[j].position
            r_vec -= L * np.round(r_vec / L)
            r = np.linalg.norm(r_vec)
            if r == 0:
                continue
            # Get mixed LJ parameters based on atom types.
            epsilon_ij, sigma_ij = compute_lj_params(atoms[i].atom_type,
                                                     atoms[j].atom_type,
                                                     lj_params)
            sr6 = (sigma_ij / r) ** 6
            sr12 = sr6 ** 2
            # Force magnitude is the derivative of the Lennard-Jones potential.
            f_mag = 24 * epsilon_ij / r * (2 * sr12 - sr6)
            force_ij = f_mag * (r_vec / r)
            forces[i] += force_ij
            forces[j] -= force_ij
            potential += 4 * epsilon_ij * (sr12 - sr6)

    # Bonded interactions (if bonds are provided).
    if bonds is not None:
        for bond in bonds:
            i = bond.atom1
            j = bond.atom2
            r_vec = atoms[i].position - atoms[j].position
            # Apply periodic boundary conditions to bonds as well.
            r_vec -= L * np.round(r_vec / L)
            r = np.linalg.norm(r_vec)
            # Harmonic bond potential: U = 0.5 * k * (r - r0)^2.
            delta = r - bond.r0
            potential += 0.5 * bond.k * delta**2
            # Force: F = -dU/dr = -k * (r - r0) in the direction of r_vec.
            f_mag = - bond.k * delta
            force_bond = f_mag * (r_vec / r)
            forces[i] += force_bond
            forces[j] -= force_bond

    return forces, potential

def velocity_verlet(atoms, forces, dt, L, lj_params, bonds=None):
    """
    Advance the simulation by one time step using the velocity Verlet algorithm.
    
    Parameters:
      atoms    : List of Atom objects.
      forces   : Current forces on the atoms.
      dt       : Time step.
      L        : Box length.
      lj_params: Lennard-Jones parameters.
      bonds    : Optional list of Bond objects.
    
    Returns:
      new_forces: Updated forces after the step.
      potential : Potential energy at the new positions.
    """
    N = len(atoms)
    # Update positions.
    for i in range(N):
        # r(t+dt) = r(t) + v(t)*dt + 0.5*(F(t)/m)*dt^2.
        atoms[i].position += atoms[i].velocity * dt + 0.5 * forces[i] / atoms[i].mass * dt**2
        # Enforce periodic boundaries.
        atoms[i].position %= L
    # Compute new forces based on updated positions.
    new_forces, potential = compute_forces(atoms, L, lj_params, bonds)
    # Update velocities.
    for i in range(N):
        atoms[i].velocity += 0.5 * (forces[i] + new_forces[i]) / atoms[i].mass * dt
    return new_forces, potential

def kinetic_energy(atoms):
    """
    Calculate the total kinetic energy of the system.
    """
    ke = 0.0
    for atom in atoms:
        ke += 0.5 * atom.mass * np.dot(atom.velocity, atom.velocity)
    return ke

def run_simulation_custom(custom_atoms=None, bonds=None, L=10.0, dt=0.005, n_steps=10000, lj_params=None):
    """
    Run the 3D molecular dynamics simulation.
    
    Parameters:
      custom_atoms: List of Atom objects. If None, random atoms are generated.
      bonds       : List of Bond objects for bonded interactions.
      L           : Size of the cubic simulation box.
      dt          : Time step.
      n_steps     : Number of integration steps.
      lj_params   : Dictionary mapping atom types to LJ parameters.
                   Default uses type 'A' with epsilon=1.0 and sigma=1.0.
    
    The simulation prints the total energy periodically and plots energy vs. time.
    """
    # Set default LJ parameters if not provided.
    if lj_params is None:
        lj_params = {'A': {'epsilon': 1.0, 'sigma': 1.0}}
    
    # If no custom atoms are provided, create a random set.
    if custom_atoms is None:
        N = 10
        np.random.seed(42)
        custom_atoms = []
        for _ in range(N):
            pos = np.random.rand(3) * L
            vel = np.random.randn(3)
            custom_atoms.append(Atom(position=pos, velocity=vel, atom_type='A', mass=1.0))
    
    forces, potential = compute_forces(custom_atoms, L, lj_params, bonds)
    energies = []
    time_arr = []
    
    for step in range(n_steps):
        forces, potential = velocity_verlet(custom_atoms, forces, dt, L, lj_params, bonds)
        ke = kinetic_energy(custom_atoms)
        total_energy = ke + potential
        energies.append(total_energy)
        time_arr.append(step * dt)
        if step % 1000 == 0:
            print(f"Step {step}, Total Energy: {total_energy:.3f}")
    
    # Plot energy vs. time to check for energy conservation.
    plt.plot(time_arr, energies)
    plt.xlabel('Time')
    plt.ylabel('Total Energy')
    plt.title('Total Energy vs Time')
    plt.show()

if __name__ == '__main__':
    # === Example 1: Running a simulation with randomly generated atoms in 3D ===
    # Uncomment the next line to run with random atoms.
    # run_simulation_custom()

    # === Example 2: Running a simulation with custom atoms and a bonded pair (a diatomic molecule) ===
    # Define two atoms with custom positions, velocities, and type 'A'.
    atom1 = Atom(position=[2.0, 2.0, 2.0], velocity=[0.1, 0.0, 0.0], atom_type='A', mass=1.0)
    atom2 = Atom(position=[2.5, 2.0, 2.0], velocity=[-0.1, 0.0, 0.0], atom_type='A', mass=1.0)
    custom_atoms = [atom1, atom2]
    # Define a bond between these two atoms with an equilibrium length of 0.5 and a stiff force constant.
    bonds = [Bond(atom_index1=0, atom_index2=1, r0=0.5, k=100.0)]
    # LJ parameters for atom type 'A'
    lj_params = {'A': {'epsilon': 1.0, 'sigma': 1.0}}
    run_simulation_custom(custom_atoms=custom_atoms, bonds=bonds, L=10.0, dt=0.001, n_steps=5000, lj_params=lj_params)

