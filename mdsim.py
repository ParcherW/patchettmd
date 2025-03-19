import numpy as np
import matplotlib.pyplot as plt
#from MDInputReader import *
#from MDInputReaderChatGPTo3minihigh import *
#from MDInputReaderqwen25coder32binstruct import *
#from MDInputReaderChatGPT4o import *
#from MDInputReaderChatGPTo1 import *
from MDInputReaderqwen25coder14binstruct import *
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
# class :
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
                        bounds=(-10, 10, -10, 10, -10, 10), max_velocity=9):
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


def compute_stretch(bond, atoms, L, boundary="periodic", k_stretch=0.1, gamma_stretch=0.1):
    """
    Compute the bond-stretching energy and forces for a bond between atoms i and j,
    and add a damping force proportional to the relative velocity along the bond.

    Returns:
        energy (float): Stretching energy.
        f_dict (dict): A dictionary with force contributions on atoms i and j.
    """
    i = bond.atom1
    j = bond.atom2
    r_vec = atoms[i].position - atoms[j].position
    if boundary == "periodic":
        r_vec -= L * np.round(r_vec / L)
    r = np.linalg.norm(r_vec)
    dr = r - bond.r0
    #energy = 0.5 * bond.k * dr**2
    energy = k_stretch*(dr)**2
    # Conservative force from the harmonic potential:
    f_pot_mag = -bond.k * dr
    f_pot = f_pot_mag * (r_vec / (r if r > 1e-12 else 1e-12))
    # Damping force: proportional to the projection of the relative velocity along the bond.
    rel_vel = atoms[i].velocity - atoms[j].velocity
    f_damp = - gamma_stretch * np.dot(rel_vel, r_vec/(r if r > 1e-12 else 1e-12)) * (r_vec/(r if r > 1e-12 else 1e-12))
    total_force = f_pot + f_damp
    return energy, {i: total_force, j: -total_force}


# def compute_angle(bend, atoms, L, boundary="periodic", k_bend=0.1, gamma_angle=0.1):
#     """
#     Compute the angle-bending energy and forces for three atoms i, j, k (with j central)
#     and add a damping term that opposes rapid changes in the angle.
    
#     This function uses a simple approximation to estimate an angular velocity.
    
#     Returns:
#         energy (float): Bending energy.
#         f_dict (dict): A dictionary with force contributions on atoms i, j, and k.
#     """
#     i = bend.atom1
#     j = bend.atom2  # central atom
#     k = bend.atom3
#     r_ij = atoms[i].position - atoms[j].position
#     r_kj = atoms[k].position - atoms[j].position
#     if boundary == "periodic":
#         r_ij -= L * np.round(r_ij / L)
#         r_kj -= L * np.round(r_kj / L)
#     r_i = np.linalg.norm(r_ij)
#     r_k = np.linalg.norm(r_kj)
#     cos_theta = np.dot(r_ij, r_kj) / ((r_i * r_k) if r_i*r_k > 1e-12 else 1e-12)
#     cos_theta = np.clip(cos_theta, -1.0, 1.0)
#     theta = np.arccos(cos_theta)
#     dtheta = theta - bend.theta0
#     #energy = 0.5 * bend.k * dtheta**2
#     energy = k_bend*(dtheta)**2
#     dU_dtheta = bend.k * dtheta
#     sin_theta = np.sin(theta)
#     if abs(sin_theta) < 1e-6:
#         sin_theta = 1e-6
#     # Conservative forces (from chain rule for the angle):
#     F_i = - dU_dtheta/(r_i*sin_theta) * (r_kj/r_k - cos_theta*(r_ij/r_i))
#     F_k = - dU_dtheta/(r_k*sin_theta) * (r_ij/r_i - cos_theta*(r_kj/r_k))
#     F_j = -(F_i + F_k)
    
#     # Estimate angular velocity (very approximate):
#     # For a rough estimate, project the tangential parts of the velocities.
#     v_i_tan = atoms[i].velocity - np.dot(atoms[i].velocity, r_ij/r_i) * (r_ij/r_i)
#     v_k_tan = atoms[k].velocity - np.dot(atoms[k].velocity, r_kj/r_k) * (r_kj/r_k)
#     dtheta_dt = np.linalg.norm(v_i_tan - v_k_tan) / ((r_i + r_k)/2.0)
#     # Damping torque (acts to reduce dtheta/dt):
#     damping = - gamma_angle * dtheta_dt
#     # Distribute the damping effect using the same geometric factors as the conservative force.
#     F_i_damp = damping/(r_i*sin_theta) * (r_kj/r_k - cos_theta*(r_ij/r_i))
#     F_k_damp = damping/(r_k*sin_theta) * (r_ij/r_i - cos_theta*(r_kj/r_k))
#     F_j_damp = -(F_i_damp + F_k_damp)
    
#     F_i_total = F_i + F_i_damp
#     F_j_total = F_j + F_j_damp
#     F_k_total = F_k + F_k_damp
#     return energy, {i: F_i_total, j: F_j_total, k: F_k_total}

# def compute_angle(bend, atoms, L, boundary="periodic", k_bend=0.1, gamma_angle=0.1):
#     """
#     Compute the bending (angle) energy and forces for an angle defined by three atoms:
#       - atom i (an end atom)
#       - atom j (the central atom)
#       - atom k (the other end atom)
#     with an equilibrium angle theta0 and force constant k (in bend.theta0 and bend.k).

#     This function computes the conservative force using the chain rule:
#         F_i = - (dU/dtheta)/(ar1*sin(theta)) [r2/ar2 - cos(theta)*(r1/ar1)]
#     (and similarly for atom k, with atom j receiving the negative sum of the end forces).

#     In addition, we estimate an angular velocity (dtheta/dt) based on the difference 
#     in the tangential velocities of atoms i and k and add a damping term:
#         F_damp ~ - gamma_angle * (dtheta/dt) * (chain rule factor).

#     Parameters:
#         bend: An AngleBend object with attributes:
#               - atom1, atom2, atom3 (indices for atoms i, j, k)
#               - theta0 (equilibrium angle, in radians)
#               - k (bending force constant)
#         atoms: list of Atom objects
#         L: box length (used if boundary=="periodic")
#         boundary: "periodic" or "reflecting"
#         gamma_angle: damping coefficient for the angular degree of freedom

#     Returns:
#         energy (float): The bending energy, 0.5*k*(theta - theta0)^2.
#         f_dict (dict): A dictionary mapping atom indices (i, j, k) to force vectors.
#     """
#     # Unpack the indices.
#     i = bend.atom1
#     j = bend.atom2  # central atom
#     k = bend.atom3

#     # Compute bond vectors from central atom to each end.
#     r1 = atoms[i].position - atoms[j].position  # from j to i
#     r2 = atoms[k].position - atoms[j].position  # from j to k

#     # Apply periodic adjustment if needed.
#     if boundary == "periodic":
#         r1 -= L * np.round(r1 / L)
#         r2 -= L * np.round(r2 / L)

#     # Compute bond lengths.
#     ar1 = np.linalg.norm(r1)
#     ar2 = np.linalg.norm(r2)
#     if ar1 < 1e-12: ar1 = 1e-12
#     if ar2 < 1e-12: ar2 = 1e-12

#     # Compute the cosine and sine of the angle.
#     cos_theta = np.dot(r1, r2) / (ar1 * ar2)
#     cos_theta = np.clip(cos_theta, -1.0, 1.0)
#     theta = np.arccos(cos_theta)
#     sin_theta = np.sin(theta)
#     if abs(sin_theta) < 1e-6: sin_theta = 1e-6  # avoid division by zero

#     # Bending energy and its derivative.
#     dtheta = theta - bend.theta0
#     energy = k_bend*(dtheta)**2
#     dU_dtheta = bend.k * dtheta

#     # Estimate angular velocity (dtheta/dt).
#     # Compute tangential (perpendicular) components of velocities for atoms i and k.
#     r1_unit = r1 / ar1
#     r2_unit = r2 / ar2
#     v_i_tan = atoms[i].velocity - np.dot(atoms[i].velocity, r1_unit) * r1_unit
#     v_k_tan = atoms[k].velocity - np.dot(atoms[k].velocity, r2_unit) * r2_unit
#     # A rough estimate: difference in tangential speeds divided by average bond length.
#     dtheta_dt = np.linalg.norm(v_i_tan - v_k_tan) / ((ar1 + ar2) / 2.0 + 1e-12)

#     # Damping contribution: we assume damping acts to oppose the angular rate.
#     damping_term = gamma_angle * dtheta_dt

#     # Total effective derivative in the angle coordinate.
#     effective = dU_dtheta + damping_term

#     # Compute the chain rule factor (gradient of theta with respect to positions).
#     # For the end atom i:
#     grad_theta_i = (1.0 / (ar1 * sin_theta)) * (r2/ar2 - cos_theta * r1/ar1)
#     # For the other end atom k:
#     grad_theta_k = (1.0 / (ar2 * sin_theta)) * (r1/ar1 - cos_theta * r2/ar2)
#     # The central atom j gets the negative sum of the end contributions.
#     grad_theta_j = - (grad_theta_i + grad_theta_k)

#     # Now, the force on each atom is:
#     F_i = - effective * grad_theta_i
#     F_k = - effective * grad_theta_k
#     F_j = - effective * grad_theta_j

#     # Return the bending energy and the forces.
#     return energy, {i: F_i, j: F_j, k: F_k}

def compute_angle(bend, atoms, L, boundary="periodic", k_bend=0.1, gamma_angle=0.1):
    """
    Compute the angle bending energy and forces for three atoms:
      - atom1 (an end atom)
      - atom2 (the central atom)
      - atom3 (the other end atom)

    The bending energy is given by:
        U = e0*(theta - theta0)**2,
    where e0 is the bending modulus (bend.k) and theta0 is the equilibrium angle (bend.theta0).

    The force contributions (returned in a dict) follow the scheme from the provided example.
    """
    import numpy as np
    import math

    # Unpack atom indices (using your naming: atom1, atom2 (central), atom3)
    i = bend.atom1
    j = bend.atom2  # central atom
    k = bend.atom3

    # Get positions from atoms list
    r = [atom.position for atom in atoms]

    # Compute bond vectors.
    # (Always define the bonds as from the central atom (j) to the end atoms.)
    if boundary == "periodic":
        r1 = r[i] - r[j]
        r1 -= L * np.round(r1 / L)
        r2 = r[k] - r[j]
        r2 -= L * np.round(r2 / L)
    else:
        r1 = r[i] - r[j]
        r2 = r[k] - r[j]

    # Compute bond lengths (with a safety check against division by zero).
    ar1 = np.linalg.norm(r1)
    ar2 = np.linalg.norm(r2)
    if ar1 < 1e-12: ar1 = 1e-12
    if ar2 < 1e-12: ar2 = 1e-12

    # Compute the cosine of the current angle.
    dot = np.dot(r1, r2)
    ndot = dot / (ar1 * ar2)
    ndot = np.clip(ndot, -1.0, 1.0)
    th = math.acos(ndot)  # current angle in radians

    # Equilibrium angle and bending modulus.
    th0 = bend.theta0
    e0 = bend.k

    # Compute bending energy: U = e0*(th - th0)**2.
    energy = e0 * (th - th0) ** 2

    # Compute the derivative dU/dtheta.
    # (The provided code uses: dUdth = -2.0*e0*(th - th0))
    dUdth = -2.0 * e0 * (th - th0)

    # A common factor in the force calculations.
    denom = math.sqrt(1.0 - ndot ** 2)
    if denom < 1e-8:
        denom = 1e-8

    # Compute force contributions on the end atoms.
    # For atom i:
    F_i = dUdth * ((r2 / (ar1 * ar2)) - (dot / (2.0 * (ar1 ** 3) * ar2))) / denom
    # For atom k:
    F_k = dUdth * ((r1 / (ar1 * ar2)) - (dot / (2.0 * ar1 * (ar2 ** 3)))) / denom

    # For the central atom j, the force is the negative sum of the contributions.
    F_j = dUdth * ((-(r1 + r2) + (dot * r1) / (ar1 ** 2) + (dot * r2) / (ar2 ** 2)) / (ar1 * ar2)) / denom

    # Return the bending energy and a dictionary mapping atom indices to force vectors.
    return energy, {i: F_i, j: F_j, k: F_k}



def compute_torsion(tors, atoms, L, boundary="periodic", gamma_torsion=0.1):
    """
    Compute the torsional (dihedral) energy for four atoms i, j, k, l
    using the cosine dihedral potential:
    
        U = V [1 + cos(n*phi - gamma)]
    
    and add a damping term that opposes rapid changes in the dihedral angle.
    
    Because the analytical derivation of the forces is quite complex,
    here we include only a simple placeholder damping term.
    
    Returns:
        energy (float): Torsional energy.
        f_dict (dict): A dictionary with force contributions on atoms i, j, k, and l.
                      (For a complete implementation, these forces should be computed from the gradient of phi.)
    """
    i, j, k, l = tors.atom1, tors.atom2, tors.atom3, tors.atom4
    b1 = atoms[j].position - atoms[i].position
    b2 = atoms[k].position - atoms[j].position
    b3 = atoms[l].position - atoms[k].position
    if boundary == "periodic":
        b1 -= L * np.round(b1 / L)
        b2 -= L * np.round(b2 / L)
        b3 -= L * np.round(b3 / L)
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)
    n1_norm = np.linalg.norm(n1)
    n2_norm = np.linalg.norm(n2)
    if n1_norm < 1e-12: n1_norm = 1e-12
    if n2_norm < 1e-12: n2_norm = 1e-12
    n1_unit = n1 / n1_norm
    n2_unit = n2 / n2_norm
    cos_phi = np.dot(n1_unit, n2_unit)
    cos_phi = np.clip(cos_phi, -1.0, 1.0)
    phi = np.arccos(cos_phi)
    if np.dot(np.cross(n1_unit, n2_unit), b2) < 0:
        phi = -phi
    energy = tors.V * (1 + np.cos(tors.n * phi - tors.gamma))
    
    # A simple damping term: estimate dphi/dt (placeholder) as difference in velocities of atoms j and k projected on b2.
    dphi_dt = np.linalg.norm(atoms[j].velocity - atoms[k].velocity) / (np.linalg.norm(b2) if np.linalg.norm(b2) > 1e-12 else 1e-12)
    damping = - gamma_torsion * dphi_dt
    # For a complete implementation, one would compute the gradient of phi with respect to each atom.
    # Here we provide a placeholder that returns zero force plus the damping effect (which is not fully implemented).
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

         # Write bond connectivity: first, the number of bonds...
  #  num_bonds = len(bond_stretches)
  #  f.write(f"{num_bonds}\n")
    # ...then one line per bond listing the indices of the two bonded atoms.
    # (Make sure that your bond indices are consistent with the order of atoms in the file.)
  #  for bond in bond_stretches:
 #       f.write(f"{bond.atom1} {bond.atom2}\n")

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
    inputdeck = MDInputReader("md_input.txt")
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
    bond_stretches = [#BondStretch(atom1=0, atom2=1, r0=0.9611, k_stretch=1480000.0), 
    #                   BondStretch(atom1=1, atom2=2, r0=0.9611, k_stretch=1480000.0),
    #                   #BondStretch(atom1=2, atom2=3, r0=5.0, k_stretch=1000.0),
                       ]
    # # One angle bend among atoms 0, 1, 2 (with atom 1 central):
    angle_bends = [#AngleBend(atom1=0, atom2=1, atom3=2, theta0=np.deg2rad(109.5), k_bend=353000.0),
    #                #AngleBend(atom1=1, atom2=2, atom3=3, theta0=np.deg2rad(109.5), k_bend=2000.0),
                    ]
    # One torsion among atoms 0, 1, 2, 3:
    torsions = [#Torsion(atom1=0, atom2=1, atom3=2, atom4=3, V=5.0, n=3, gamma=0.0)
                ]
    
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
