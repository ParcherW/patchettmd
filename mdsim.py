import numpy as np
import matplotlib.pyplot as plt
import MDInputReader

# Define an Atom class to store properties for each particle.
class Atom:
    def __init__(self, position, velocity, atom_type=0, mass=1.0):
        self.position = np.array(position, dtype=float)  # 3D position
        self.velocity = np.array(velocity, dtype=float)  # 3D velocity
        self.atom_type = atom_type
        self.mass = mass

# Define a Bond class for bonded interactions between atoms.
class Bond:
    def __init__(self, atom_index1, atom_index2, r0, k):
        self.atom1 = atom_index1  # index of the first atom in the bond
        self.atom2 = atom_index2  # index of the second atom
        self.r0 = r0              # equilibrium bond length
        self.k = k                # bond force constant

import numpy as np

def create_random_atoms(num_atoms=1, seed=1, atom_type=0, mass=1.0, 
                        bounds=(-10, 10, -10, 10, -10, 10), max_velocity=3):
    """
    Creates a list of randomly initialized Atom objects.

    Args:
        num_atoms (int): Number of atoms to create. Defaults to 1.
        seed (int): Random seed for reproducibility. Defaults to 1.
        atom_type (int): Type of the atom. Defaults to 0.
        mass (float): Mass of the atom. Defaults to 1.0.
        bounds (tuple): A six-tuple (xmin, xmax, ymin, ymax, zmin, zmax) defining 
                        the range within which positions should be generated.
        max_velocity (float): Maximum absolute velocity value. Defaults to 3.

    Returns:
        list: A list of Atom objects with random positions and velocities.
    """
    np.random.seed(seed)  # Set random seed

    xmin, xmax, ymin, ymax, zmin, zmax = bounds  # Unpack bounds
    atoms = []

    for _ in range(num_atoms):
        position = np.array([
            np.random.uniform(xmin, xmax),
            np.random.uniform(ymin, ymax),
            np.random.uniform(zmin, zmax)
        ])

        velocity = np.random.uniform(-max_velocity, max_velocity, size=3)

        atoms.append(Atom(position, velocity, atom_type, mass))

    return atoms


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

def compute_forces(atoms, L, lj_params, bonds=None, boundary="periodic"):
    """
    Compute forces on each atom due to nonbonded Lennard-Jones interactions 
    and bonded (harmonic) interactions if bonds are provided.
    
    For periodic boundaries, the minimum image convention is applied.
    For reflecting boundaries, the direct distance is used.
    """
    N = len(atoms)
    forces = np.zeros((N, 3))
    potential = 0.0
    eps = 1e-12  # small number to avoid division by zero

    # Nonbonded Lennard-Jones interactions.
    for i in range(N):
        for j in range(i+1, N):
            r_vec = atoms[i].position - atoms[j].position
            if boundary == "periodic":
                # Apply the minimum image convention for periodic boundaries.
                r_vec -= L * np.round(r_vec / L)
            # For reflecting, we use the direct difference.
            r = np.linalg.norm(r_vec)
            if r < eps:
                r = eps
            epsilon_ij, sigma_ij = compute_lj_params(atoms[i].atom_type,
                                                     atoms[j].atom_type,
                                                     lj_params)
            sr6 = (sigma_ij / r) ** 6
            sr12 = sr6 ** 2
            f_mag = 24 * epsilon_ij / r * (2 * sr12 - sr6)
            force_ij = f_mag * (r_vec / r)
            forces[i] += force_ij
            forces[j] -= force_ij
            potential += 4 * epsilon_ij * (sr12 - sr6)

    # Bonded interactions.
    if bonds is not None:
        for bond in bonds:
            i = bond.atom1
            j = bond.atom2
            r_vec = atoms[i].position - atoms[j].position
            if boundary == "periodic":
                r_vec -= L * np.round(r_vec / L)
            r = np.linalg.norm(r_vec)
            r_safe = r if r >= eps else eps
            delta = r - bond.r0
            potential += 0.5 * bond.k * delta**2
            f_mag = -bond.k * delta
            force_bond = f_mag * (r_vec / r_safe)
            forces[i] += force_bond
            forces[j] -= force_bond

    return forces, potential

def velocity_verlet(atoms, forces, dt, L, lj_params, bonds=None, boundary="periodic"):
    """
    Advance the simulation by one time step using the velocity Verlet algorithm.
    
    For periodic boundaries, positions are wrapped using modulo L.
    For reflecting boundaries, positions are reflected off the walls.
    """
    N = len(atoms)
    for i in range(N):
        # Update positions.
        atoms[i].position += atoms[i].velocity * dt + 0.5 * forces[i] / atoms[i].mass * dt**2

        if boundary == "periodic":
            # Wrap positions for periodic boundaries.
            atoms[i].position %= L
        elif boundary == "reflecting":
            # Reflecting boundaries: bounce off walls.
            for dim in range(3):
                if atoms[i].position[dim] < 0:
                    atoms[i].position[dim] = -atoms[i].position[dim]
                    atoms[i].velocity[dim] = -atoms[i].velocity[dim]
                elif atoms[i].position[dim] > L:
                    atoms[i].position[dim] = 2 * L - atoms[i].position[dim]
                    atoms[i].velocity[dim] = -atoms[i].velocity[dim]

    # Compute new forces.
    new_forces, potential = compute_forces(atoms, L, lj_params, bonds, boundary=boundary)
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

def write_xyz_frame(f, atoms, step, L):
    """
    Write one frame in extended XYZ format to file f.
    The comment line includes lattice information for OVITO.
    """
    f.write(f"{len(atoms)}\n")
    # Include lattice information for a cubic box: Lattice="L 0 0 0  L 0 0 0  L"
    f.write(f"Step {step} Lattice=\"{L:.5f} 0.0 0.0 0.0 {L:.5f} 0.0 0.0 0.0 {L:.5f}\"\n")
    for atom in atoms:
        x, y, z = atom.position
        f.write(f"{atom.atom_type} {x:.5f} {y:.5f} {z:.5f}\n")

def run_simulation_custom(custom_atoms=None, bonds=None, L=10.0, dt=0.005, n_steps=10000, 
                          lj_params=None, output_interval=100, boundary="periodic"):
    """
    Run the 3D molecular dynamics simulation.
    
    In addition to computing energies, this function writes an extended XYZ trajectory file 
    (with lattice information) that can be opened in OVITO for visualization.
    
    The 'boundary' parameter can be either "periodic" or "reflecting".
    """
    # Set default LJ parameters.
    if lj_params is None:
        lj_params = {0: {'epsilon': 1.0, 'sigma': 1.0}}
    # Create random atoms if none are provided.
    if custom_atoms is None:
        N = 10
        np.random.seed(42)
        custom_atoms = []
        for _ in range(N):
            pos = np.random.rand(3) * L
            vel = np.random.randn(3)
            custom_atoms.append(Atom(position=pos, velocity=vel, atom_type=0, mass=1.0))
    
    forces, potential = compute_forces(custom_atoms, L, lj_params, bonds, boundary=boundary)
    energies = []
    time_arr = []
    # Open file for trajectory output.
    traj_file = open("trajectory.xyz", "w")
    # Main simulation loop.
    for step in range(n_steps):
        forces, potential = velocity_verlet(custom_atoms, forces, dt, L, lj_params, bonds, boundary=boundary)
        ke = kinetic_energy(custom_atoms)
        total_energy = ke + potential
        energies.append(total_energy)
        time_arr.append(step * dt)
        # Write trajectory at specified intervals.
        if step % output_interval == 0:
            write_xyz_frame(traj_file, custom_atoms, step, L)
        if step % 1000 == 0:
            print(f"Step {step}, Total Energy: {total_energy:.3f}")
    traj_file.close()
    print("Trajectory saved as 'trajectory.xyz'.")
    
    # Save the energy vs. time plot.
    plt.figure(figsize=(8, 6))
    plt.plot(time_arr, energies, label='Total Energy')
    plt.xlabel('Time')
    plt.ylabel('Total Energy')
    plt.title('Total Energy vs Time')
    plt.legend()
    plt.savefig('energy.png', dpi=300)
    plt.close()
    print("Energy plot saved as 'energy.png'.")

if __name__ == '__main__':
    inputdeck = MDInputReader.MDInputReader("md_input.txt")
    config = inputdeck.read()
    # Print out the parsed configuration.
    for section, params in config.items():
        print(f"Section: {section}")
        for key, value in params.items():
            print(f"  {key}: {value}")
    print("john " + str(config["Potential"]["sigma"]))
    print(config)
    
    if len(config['Atoms']) != config['System']['n_particles']:
        custom_atoms=create_random_atoms(config['System']['n_particles'])
    else:
      # Construct atoms from the input deck.
      custom_atoms = []
      for atomid, atomdescription in config['Atoms'].items():
        custom_atoms.append(
            Atom(position=atomdescription[0:3],
                 velocity=atomdescription[3:6],
                 atom_type=atomdescription[6],
                 mass=atomdescription[7])
        )
    
    # Define a bond (if any).
    bonds = [Bond(atom_index1=0, atom_index2=1, r0=0.5, k=100.0)]
    
    # LJ parameters for atom type 0 from the config.
    lj_params = {0: {'epsilon': config["Potential"]["epsilon"], 'sigma': config["Potential"]["sigma"]}}

    # Choose boundary condition: "periodic" or "reflecting"
    boundary_choice = config["System"]["boundary"]
    
    # Run the simulation.
    run_simulation_custom(custom_atoms=custom_atoms,
                          bonds=bonds,
                          L=float(config["System"]["box_length"]),
                          dt=config["Simulation"]["time_step"],
                          n_steps=config["Simulation"]["n_steps"],
                          lj_params=lj_params,
                          boundary=boundary_choice)
