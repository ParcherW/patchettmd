import numpy as np
import matplotlib.pyplot as plt
import MDInputReader

# Define an Atom class to store properties for each particle.
class Atom:
    def __init__(self, position, velocity, atom_type='A', mass=1.0):
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
    """
    N = len(atoms)
    forces = np.zeros((N, 3))
    potential = 0.0
    eps = 1e-12  # small number to avoid division by zero

    # Nonbonded Lennard-Jones interactions.
    for i in range(N):
        for j in range(i+1, N):
            r_vec = atoms[i].position - atoms[j].position
            r_vec -= L * np.round(r_vec / L)  # minimum image convention
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

def velocity_verlet(atoms, forces, dt, L, lj_params, bonds=None):
    """
    Advance the simulation by one time step using the velocity Verlet algorithm.
    """
    N = len(atoms)
    for i in range(N):
        atoms[i].position += atoms[i].velocity * dt + 0.5 * forces[i] / atoms[i].mass * dt**2
        atoms[i].position %= L  # apply periodic boundaries
    new_forces, potential = compute_forces(atoms, L, lj_params, bonds)
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

def write_xyz_frame(f, atoms, step):
    """
    Write one frame in XYZ format to file f.
    """
    f.write(f"{len(atoms)}\n")
    f.write(f"Step {step}\n")
    for atom in atoms:
        x, y, z = atom.position
        f.write(f"{atom.atom_type} {x:.5f} {y:.5f} {z:.5f}\n")

def run_simulation_custom(custom_atoms=None, bonds=None, L=10.0, dt=0.005, n_steps=10000, lj_params=None, output_interval=100):
    """
    Run the 3D molecular dynamics simulation.
    In addition to computing energies, this function writes an XYZ trajectory file 
    that can be opened in OVITO for visualization.
    """
    # Set default LJ parameters.
    if lj_params is None:
        lj_params = {'A': {'epsilon': 1.0, 'sigma': 1.0}}
    # Create random atoms if none are provided.
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
    # Open file for trajectory output.
    traj_file = open("trajectory.xyz", "w")
    # Main simulation loop.
    for step in range(n_steps):
        forces, potential = velocity_verlet(custom_atoms, forces, dt, L, lj_params, bonds)
        ke = kinetic_energy(custom_atoms)
        total_energy = ke + potential
        energies.append(total_energy)
        time_arr.append(step * dt)
        # Write to trajectory file at specified intervals.
        if step % output_interval == 0:
            write_xyz_frame(traj_file, custom_atoms, step)
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
    # === Example: Running a simulation with custom atoms and a bonded pair (diatomic molecule) ===
    # Define two atoms with custom positions, velocities, and type 'A'.
    # atom1 = Atom(position=[20.0, 20.0, 20.0], velocity=[0.0, 0.0, 0.0], atom_type='A', mass=1.0)
    # atom2 = Atom(position=[20.5, 20.0, 20.0], velocity=[0.0, 0.0, 0.0], atom_type='A', mass=1.0)
    # atom3 = Atom(position=[19.5, 19.6, 19.0], velocity=[0.0, 0.0, 0.0], atom_type='A', mass=1.0)
    atom1 = Atom(position=[15.0, 15.0, 15.0], velocity=[0.0, 0.0, 0.0], atom_type='A', mass=1.0)
    atom2 = Atom(position=[17.0, 15.0, 15.0], velocity=[0.0, 0.0, 0.0], atom_type='A', mass=1.0)
    atom3 = Atom(position=[14.7, 14.6, 14.5], velocity=[0.0, 0.0, 0.0], atom_type='A', mass=1.0)
    custom_atoms = [atom1, atom2, atom3]
    # custom_atoms = [atom1, atom2, atom3]
    # Define a bond between these two atoms with an equilibrium length of 0.5 and a stiff force constant.
    bonds = [Bond(atom_index1=0, atom_index2=1, r0=0.5, k=100.0)]
    # LJ parameters for atom type 'A'
    #lj_params = {'A': {'epsilon': 1.0, 'sigma': 1.0}}
    lj_params = {'A': {'epsilon':config["Potential"]["epsilon"], 'sigma': config["Potential"]["sigma"]}}

    #run_simulation_custom(custom_atoms=custom_atoms, bonds=bonds, L=10.0, dt=0.00001, n_steps=5000, lj_params=lj_params)
    run_simulation_custom(custom_atoms=custom_atoms, bonds=bonds, L=float(config["System"]["box_length"]), dt=config["Simulation"]["time_step"], n_steps=config["Simulation"]["n_steps"], lj_params=lj_params)

