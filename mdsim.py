import numpy as np
import matplotlib.pyplot as plt
import MDInputReader
from Atom import *

# ------------------------------
# Atom and Bonded Interaction Classes
# ------------------------------

# class Atom:
#     def __init__(self, position, velocity, atom_type=0, mass=1.0, charge=1.0):
#         self.position = np.array(position, dtype=float)  # 3D position
#         self.velocity = np.array(velocity, dtype=float)  # 3D velocity
#         self.atom_type = atom_type
#         self.mass = mass
#         self.charge = charge

# # Bond stretching between two atoms.
# class BondStretch:
#     def __init__(self, atom1, atom2, r0, k_stretch):
#         self.atom1 = atom1  # index of first atom
#         self.atom2 = atom2  # index of second atom
#         self.r0 = r0        # equilibrium bond length
#         self.k = k_stretch  # force constant

# # Angle bending among three atoms (with atom2 as the central atom).
# class AngleBend:
#     def __init__(self, atom1, atom2, atom3, theta0, k_bend):
#         self.atom1 = atom1  # first atom
#         self.atom2 = atom2  # central atom
#         self.atom3 = atom3  # third atom
#         self.theta0 = theta0  # equilibrium angle (radians)
#         self.k = k_bend       # force constant

# # Torsion (dihedral) interaction among four atoms.
# class Torsion:
#     def __init__(self, atom1, atom2, atom3, atom4, V, n, gamma):
#         self.atom1 = atom1  # first atom
#         self.atom2 = atom2  # second atom
#         self.atom3 = atom3  # third atom
#         self.atom4 = atom4  # fourth atom
#         self.V = V          # amplitude (energy barrier)
#         self.n = n          # periodicity
#         self.gamma = gamma  # phase offset (radians)

# ------------------------------
# Utility Functions (Random atoms, etc.)
# ------------------------------

def create_random_atoms(num_atoms=1, seed=1, atom_type=0, mass=1.0, charge=1.0, 
                        bounds=(-10, 10, -10, 10, -10, 10), max_velocity=3):
    np.random.seed(seed)
    xmin, xmax, ymin, ymax, zmin, zmax = bounds
    atoms = []
    for _ in range(num_atoms):
        position = np.array([
            np.random.uniform(xmin, xmax),
            np.random.uniform(ymin, ymax),
            np.random.uniform(zmin, zmax)
        ])
        velocity = np.random.uniform(-max_velocity, max_velocity, size=3)
        atoms.append(Atom(position, velocity, atom_type, mass, charge))
    return atoms

# ------------------------------
# Nonbonded Interactions (LJ & Coulomb)
# ------------------------------

def compute_lj_params(type1, type2, lj_params):
    sigma1 = lj_params[type1]['sigma']
    sigma2 = lj_params[type2]['sigma']
    epsilon1 = lj_params[type1]['epsilon']
    epsilon2 = lj_params[type2]['epsilon']
    sigma_ij = 0.5 * (sigma1 + sigma2)
    epsilon_ij = np.sqrt(epsilon1 * epsilon2)
    return epsilon_ij, sigma_ij

def compute_nonbonded_forces(atoms, L, lj_params, boundary="periodic", coulomb_k=1.0):
    N = len(atoms)
    forces = np.zeros((N, 3))
    potential = 0.0
    eps = 1e-12
    for i in range(N):
        for j in range(i+1, N):
            r_vec = atoms[i].position - atoms[j].position
            if boundary == "periodic":
                r_vec -= L * np.round(r_vec / L)
            r = np.linalg.norm(r_vec)
            if r < eps: r = eps
            # Lennard-Jones
            epsilon_ij, sigma_ij = compute_lj_params(atoms[i].atom_type,
                                                     atoms[j].atom_type,
                                                     lj_params)
            sr6 = (sigma_ij / r)**6
            sr12 = sr6**2
            f_lj_mag = 24 * epsilon_ij / r * (2*sr12 - sr6)
            force_lj = f_lj_mag * (r_vec / r)
            forces[i] += force_lj
            forces[j] -= force_lj
            potential += 4 * epsilon_ij * (sr12 - sr6)
            # Coulomb (if enabled)
            if coulomb_k is not None and coulomb_k != 0:
                potential += coulomb_k * atoms[i].charge * atoms[j].charge / r
                force_coulomb = coulomb_k * atoms[i].charge * atoms[j].charge / (r**2) * (r_vec / r)
                forces[i] += force_coulomb
                forces[j] -= force_coulomb
    return forces, potential

# ------------------------------
# Bonded Interaction Functions
# ------------------------------

def compute_stretch(bond, atoms, L, boundary="periodic"):
    i = bond.atom1
    j = bond.atom2
    r_vec = atoms[i].position - atoms[j].position
    if boundary == "periodic":
        r_vec -= L * np.round(r_vec / L)
    r = np.linalg.norm(r_vec)
    dr = r - bond.r0
    energy = 0.5 * bond.k * dr**2
    f_mag = -bond.k * dr
    f_vec = f_mag * (r_vec / (r if r>1e-12 else 1e-12))
    # Return forces on atoms i and j
    return energy, {i: f_vec, j: -f_vec}

def compute_angle(bend, atoms, L, boundary="periodic"):
    i = bend.atom1
    j = bend.atom2  # central atom
    k = bend.atom3
    r_ij = atoms[i].position - atoms[j].position
    r_kj = atoms[k].position - atoms[j].position
    if boundary == "periodic":
        r_ij -= L * np.round(r_ij / L)
        r_kj -= L * np.round(r_kj / L)
    r_i = np.linalg.norm(r_ij)
    r_k = np.linalg.norm(r_kj)
    cos_theta = np.dot(r_ij, r_kj) / ((r_i * r_k) if r_i*r_k>1e-12 else 1e-12)
    cos_theta = np.clip(cos_theta, -1.0, 1.0)
    theta = np.arccos(cos_theta)
    dtheta = theta - bend.theta0
    energy = 0.5 * bend.k * dtheta**2
    dU_dtheta = bend.k * dtheta
    sin_theta = np.sin(theta)
    if abs(sin_theta) < 1e-6:
        sin_theta = 1e-6
    # Forces on atoms i and k:
    F_i = - dU_dtheta/(r_i*sin_theta) * (r_kj/ r_k - cos_theta*(r_ij/r_i))
    F_k = - dU_dtheta/(r_k*sin_theta) * (r_ij/ r_i - cos_theta*(r_kj/r_k))
    F_j = -(F_i + F_k)
    return energy, {i: F_i, j: F_j, k: F_k}

def compute_torsion(tors, atoms, L, boundary="periodic"):
    # For four atoms i, j, k, l, compute the dihedral angle phi.
    i = tors.atom1; j = tors.atom2; k = tors.atom3; l = tors.atom4
    # Vectors between atoms:
    b1 = atoms[j].position - atoms[i].position
    b2 = atoms[k].position - atoms[j].position
    b3 = atoms[l].position - atoms[k].position
    if boundary == "periodic":
        b1 -= L * np.round(b1 / L)
        b2 -= L * np.round(b2 / L)
        b3 -= L * np.round(b3 / L)
    # Normal vectors to the planes:
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)
    # Normalize normals:
    n1_norm = np.linalg.norm(n1)
    n2_norm = np.linalg.norm(n2)
    if n1_norm < 1e-12: n1_norm = 1e-12
    if n2_norm < 1e-12: n2_norm = 1e-12
    n1_unit = n1 / n1_norm
    n2_unit = n2 / n2_norm
    # Dihedral angle:
    cos_phi = np.dot(n1_unit, n2_unit)
    cos_phi = np.clip(cos_phi, -1.0, 1.0)
    phi = np.arccos(cos_phi)
    # Determine sign of phi using b2.
    if np.dot(np.cross(n1_unit, n2_unit), b2) < 0:
        phi = -phi
    energy = tors.V * (1 + np.cos(tors.n * phi - tors.gamma))
    # For forces, a complete analytical derivation is complex.
    # Here we provide a placeholder that returns zero forces.
    # In a production code, you should replace this with proper force expressions.
    zero_force = np.zeros(3)
    return energy, {i: zero_force, j: zero_force, k: zero_force, l: zero_force}

def compute_bonded_forces(atoms, L, bond_stretches, angle_bends, torsions, boundary="periodic"):
    N = len(atoms)
    forces = np.zeros((N, 3))
    bonded_potential = 0.0
    # Stretching contributions.
    for bond in bond_stretches:
        energy, f_dict = compute_stretch(bond, atoms, L, boundary)
        bonded_potential += energy
        for idx, f in f_dict.items():
            forces[idx] += f
    # Angle bending contributions.
    for angle in angle_bends:
        energy, f_dict = compute_angle(angle, atoms, L, boundary)
        bonded_potential += energy
        for idx, f in f_dict.items():
            forces[idx] += f
    # Torsional contributions.
    for tors in torsions:
        energy, f_dict = compute_torsion(tors, atoms, L, boundary)
        bonded_potential += energy
        for idx, f in f_dict.items():
            forces[idx] += f
    return forces, bonded_potential

# ------------------------------
# Overall Force Calculation & Integration
# ------------------------------

def compute_forces_total(atoms, L, lj_params, bond_stretches, angle_bends, torsions,
                         boundary="periodic", coulomb_k=1.0):
    # Nonbonded contributions.
    forces_nb, potential_nb = compute_nonbonded_forces(atoms, L, lj_params, boundary, coulomb_k)
    # Bonded contributions.
    forces_bonded, potential_bonded = compute_bonded_forces(atoms, L, bond_stretches, angle_bends, torsions, boundary)
    total_forces = forces_nb + forces_bonded
    total_potential = potential_nb + potential_bonded
    return total_forces, total_potential

def velocity_verlet(atoms, forces, dt, L, lj_params, bond_stretches, angle_bends, torsions,
                    boundary="periodic", coulomb_k=1.0):
    N = len(atoms)
    for i in range(N):
        atoms[i].position += atoms[i].velocity * dt + 0.5 * forces[i] / atoms[i].mass * dt**2
        if boundary == "periodic":
            atoms[i].position %= L
        elif boundary == "reflecting":
            for dim in range(3):
                if atoms[i].position[dim] < 0:
                    atoms[i].position[dim] = -atoms[i].position[dim]
                    atoms[i].velocity[dim] = -atoms[i].velocity[dim]
                elif atoms[i].position[dim] > L:
                    atoms[i].position[dim] = 2 * L - atoms[i].position[dim]
                    atoms[i].velocity[dim] = -atoms[i].velocity[dim]
    new_forces, potential = compute_forces_total(atoms, L, lj_params, bond_stretches, angle_bends, torsions,
                                                  boundary=boundary, coulomb_k=coulomb_k)
    for i in range(N):
        atoms[i].velocity += 0.5 * (forces[i] + new_forces[i]) / atoms[i].mass * dt
    return new_forces, potential

def kinetic_energy(atoms):
    ke = 0.0
    for atom in atoms:
        ke += 0.5 * atom.mass * np.dot(atom.velocity, atom.velocity)
    return ke

def write_xyz_frame(f, atoms, step, L):
    f.write(f"{len(atoms)}\n")
    f.write(f"Step {step} Lattice=\"{L:.5f} 0.0 0.0 0.0 {L:.5f} 0.0 0.0 0.0 {L:.5f}\"\n")
    for atom in atoms:
        x, y, z = atom.position
        f.write(f"{atom.atom_type} {x:.5f} {y:.5f} {z:.5f}\n")

# ------------------------------
# Main Simulation Function
# ------------------------------

def run_simulation_custom(custom_atoms=None, bond_stretches=None, angle_bends=None, torsions=None,
                          L=10.0, dt=0.005, n_steps=10000, lj_params=None, coulomb_k=1.0,
                          output_interval=100, boundary="periodic"):
    if lj_params is None:
        lj_params = {0: {'epsilon': 1.0, 'sigma': 1.0}}
    if custom_atoms is None:
        N = 10
        np.random.seed(42)
        custom_atoms = []
        for _ in range(N):
            pos = np.random.rand(3) * L
            vel = np.random.randn(3)
            custom_atoms.append(Atom(position=pos, velocity=vel, atom_type=0, mass=1.0, charge=1.0))
    # If no bonded interactions are provided, use empty lists.
    if bond_stretches is None: bond_stretches = []
    if angle_bends is None: angle_bends = []
    if torsions is None: torsions = []
    
    forces, potential = compute_forces_total(custom_atoms, L, lj_params, bond_stretches, angle_bends, torsions,
                                              boundary=boundary, coulomb_k=coulomb_k)
    energies = []
    time_arr = []
    traj_file = open("trajectory.xyz", "w")
    for step in range(n_steps):
        forces, potential = velocity_verlet(custom_atoms, forces, dt, L, lj_params,
                                            bond_stretches, angle_bends, torsions,
                                            boundary=boundary, coulomb_k=coulomb_k)
        ke = kinetic_energy(custom_atoms)
        total_energy = ke + potential
        energies.append(total_energy)
        time_arr.append(step * dt)
        if step % output_interval == 0:
            write_xyz_frame(traj_file, custom_atoms, step, L)
        if step % 1000 == 0:
            print(f"Step {step}, Total Energy: {total_energy:.3f}")
    traj_file.close()
    print("Trajectory saved as 'trajectory.xyz'.")
    plt.figure(figsize=(8, 6))
    plt.plot(time_arr, energies, label='Total Energy')
    plt.xlabel('Time')
    plt.ylabel('Total Energy')
    plt.title('Total Energy vs Time')
    plt.legend()
    plt.savefig('energy.png', dpi=300)
    plt.close()
    print("Energy plot saved as 'energy.png'.")

# ------------------------------
# Main Block (Input Reading, etc.)
# ------------------------------

if __name__ == '__main__':
    inputdeck = MDInputReader.MDInputReader("md_input.txt")
    config = inputdeck.read()
    for section, params in config.items():
        print(f"Section: {section}")
        for key, value in params.items():
            print(f"  {key}: {value}")
    print("john " + str(config["Potential"]["sigma"]))
    print(config)
    
    bounds = [0, config["System"]["box_length"],
              0, config["System"]["box_length"],
              0, config["System"]["box_length"]]
    if len(config['Atoms']) != config['System']['n_particles']:
        custom_atoms = create_random_atoms(config['System']['n_particles'], bounds=bounds, max_velocity=6)
    else:
        custom_atoms = []
        for atomid, atomdescription in config['Atoms'].items():
            custom_atoms.append(
                Atom(position=atomdescription[0:3],
                     velocity=atomdescription[3:6],
                     atom_type=atomdescription[6],
                     mass=atomdescription[7],
                     charge=atomdescription[8])
            )
    
    # Example: define bonded interactions.
    # For instance, one bond stretch between atoms 0 and 1:
    bond_stretches = [BondStretch(atom1=0, atom2=1, r0=0.5, k_stretch=100.0)]
    # One angle bend among atoms 0, 1, 2 (with atom 1 central):
    angle_bends = [AngleBend(atom1=0, atom2=1, atom3=2, theta0=np.deg2rad(109.5), k_bend=20.0)]
    # One torsion among atoms 0, 1, 2, 3:
    torsions = [Torsion(atom1=0, atom2=1, atom3=2, atom4=3, V=5.0, n=3, gamma=0.0)]
    
    lj_params = {0: {'epsilon': config["Potential"]["epsilon"], 'sigma': config["Potential"]["sigma"]}}
    boundary_choice = config["System"]["boundary"]
    coulomb_k = config.get("Coulomb_k", 1.0)
    
    run_simulation_custom(custom_atoms=custom_atoms,
                          bond_stretches=bond_stretches,
                          angle_bends=angle_bends,
                          torsions=torsions,
                          L=float(config["System"]["box_length"]),
                          dt=config["Simulation"]["time_step"],
                          n_steps=config["Simulation"]["n_steps"],
                          lj_params=lj_params,
                          coulomb_k=coulomb_k,
                          output_interval=config["Output"]["output_frequency"],
                          boundary=boundary_choice)
