import sys
import numpy as np
from math import degrees, acos
from multiprocessing import Pool, cpu_count

# Constants for filtering atoms and bonds
CHITOSAN_DONOR_ACCEPTOR_TYPES = {7, 8, 10}  # {N7, O8, O10} in chitosan can act as donors or acceptors
WATER_DONOR_ACCEPTOR_TYPE = 12              # O12 in water can act as donor or acceptor
CHITOSAN_HYDROGEN_TYPE = 6                  # H6 in chitosan acts as hydrogen
WATER_HYDROGEN_TYPE = 11                    # H11 in water acts as hydrogen
DISTANCE_CRITERIA = 2.5                     # Maximum distance for hydrogen bond (in Å)
ANGLE_CRITERIA = (150, 180)                 # Acceptable angle range for hydrogen bond (in degrees)

# Map for atom type to labeled identifier
ATOM_LABELS = {
    6: 'H6',
    7: 'N7',
    8: 'O8',
    10: 'O10',
    11: 'H11',
    12: 'O12'
}

# Load bond information from file
def load_bond_info(filename):
    bonds = []
    with open(filename, 'r') as file:
        for line in file:
            parts = line.split()
            bond_id = int(parts[0])
            bond_type = int(parts[1])
            atom1_id = int(parts[2])
            atom2_id = int(parts[3])
            bonds.append((bond_id, bond_type, atom1_id, atom2_id))
    return bonds

# Load trajectory data for multiple timesteps from file
def load_trajectory(filename):
    timesteps = {}
    current_timestep = None
    atoms = []

    with open(filename, 'r') as file:
        reading_atoms = False
        for line in file:
            if "ITEM: TIMESTEP" in line:
                if current_timestep is not None:
                    timesteps[current_timestep] = atoms
                    atoms = []
                current_timestep = int(next(file).strip())
            elif "ITEM: ATOMS" in line:
                reading_atoms = True
            elif "ITEM:" in line:
                reading_atoms = False
            elif reading_atoms:
                parts = line.split()
                atom = {
                    'id': int(parts[0]),
                    'mol': int(parts[1]),  # Molecule ID to distinguish chitosan from water
                    'type': int(parts[2]),
                    'q': float(parts[3]),
                    'x': float(parts[4]),
                    'y': float(parts[5]),
                    'z': float(parts[6]),
                    'ix': int(parts[7]),
                    'iy': int(parts[8]),
                    'iz': int(parts[9])
                }
                atoms.append(atom)

        if current_timestep is not None:
            timesteps[current_timestep] = atoms

    return timesteps

# Calculate angle between donor, hydrogen, and acceptor atoms
def calculate_angle(donor_pos, hydrogen_pos, acceptor_pos):
    vec_donor_hydrogen = hydrogen_pos - donor_pos
    vec_hydrogen_acceptor = acceptor_pos - hydrogen_pos
    cos_angle = np.dot(vec_donor_hydrogen, vec_hydrogen_acceptor) / (np.linalg.norm(vec_donor_hydrogen) * np.linalg.norm(vec_hydrogen_acceptor))
    angle = degrees(acos(np.clip(cos_angle, -1.0, 1.0)))
    return angle

# Identify hydrogen bonds in a single timestep
def identify_hydrogen_bonds(args):
    bonds, atoms, timestep = args
    atom_dict = {atom['id']: atom for atom in atoms}
    hydrogen_bonds = []

    for bond_id, bond_type, atom1_id, atom2_id in bonds:
        atom1 = atom_dict.get(atom1_id)
        atom2 = atom_dict.get(atom2_id)

        if not atom1 or not atom2:
            continue

        # Identify donor-hydrogen-acceptor relationships based on the specified rules
        if (atom1['type'] in CHITOSAN_DONOR_ACCEPTOR_TYPES and atom2['type'] == CHITOSAN_HYDROGEN_TYPE):
            # Case: Chitosan donor {N7, O8, O10} with H6 as hydrogen
            donor, hydrogen = atom1, atom2
            potential_acceptors = [a for a in atoms if a['type'] == WATER_DONOR_ACCEPTOR_TYPE and a['mol'] != donor['mol']]
        elif (atom1['type'] == CHITOSAN_HYDROGEN_TYPE and atom2['type'] in CHITOSAN_DONOR_ACCEPTOR_TYPES):
            # Case: Chitosan donor {N7, O8, O10} with H6 as hydrogen
            donor, hydrogen = atom2, atom1
            potential_acceptors = [a for a in atoms if a['type'] == WATER_DONOR_ACCEPTOR_TYPE and a['mol'] != donor['mol']]
        elif (atom1['type'] == WATER_DONOR_ACCEPTOR_TYPE and atom2['type'] == WATER_HYDROGEN_TYPE):
            # Case: Water donor (O12) with H11 as hydrogen
            donor, hydrogen = atom1, atom2
            potential_acceptors = [a for a in atoms if a['type'] in CHITOSAN_DONOR_ACCEPTOR_TYPES and a['mol'] != donor['mol']]
        elif (atom1['type'] == WATER_HYDROGEN_TYPE and atom2['type'] == WATER_DONOR_ACCEPTOR_TYPE):
            # Case: Water donor (O12) with H11 as hydrogen
            donor, hydrogen = atom2, atom1
            potential_acceptors = [a for a in atoms if a['type'] in CHITOSAN_DONOR_ACCEPTOR_TYPES and a['mol'] != donor['mol']]
        else:
            # Skip if the relationship does not match specified donor-hydrogen pairs
            continue

        hydrogen_pos = np.array([hydrogen['x'], hydrogen['y'], hydrogen['z']])
        donor_pos = np.array([donor['x'], donor['y'], donor['z']])

        for acceptor in potential_acceptors:
            acceptor_pos = np.array([acceptor['x'], acceptor['y'], acceptor['z']])
            distance = np.linalg.norm(hydrogen_pos - acceptor_pos)
            if distance <= DISTANCE_CRITERIA:
                angle = calculate_angle(donor_pos, hydrogen_pos, acceptor_pos)
                if ANGLE_CRITERIA[0] <= angle <= ANGLE_CRITERIA[1]:
                    hydrogen_bonds.append({
                        'donor_id': donor['id'],
                        'donor_type': ATOM_LABELS.get(donor['type'], 'Unknown'),
                        'hydrogen_id': hydrogen['id'],
                        'hydrogen_type': ATOM_LABELS.get(hydrogen['type'], 'Unknown'),
                        'acceptor_id': acceptor['id'],
                        'acceptor_type': ATOM_LABELS.get(acceptor['type'], 'Unknown'),
                        'distance': distance,
                        'angle': angle,
                        'timestep': timestep
                    })

    return hydrogen_bonds

# Main code execution
trajectory_file = sys.argv[1]
bond_info_file = sys.argv[2]
output_file = trajectory_file.rsplit('.', 1)[0] + '-hydrogen-bonds.dat'

bonds = load_bond_info(bond_info_file)
timesteps = load_trajectory(trajectory_file)

# User-defined process count
max_processes = int(sys.argv[3])
max_processes = min(max_processes, cpu_count())

# Prepare arguments for parallel processing
args = [(bonds, atoms, timestep) for timestep, atoms in timesteps.items()]

# Use multiprocessing Pool with user-defined process count
with Pool(processes=max_processes) as pool:
    results = pool.map(identify_hydrogen_bonds, args)

# Collect all hydrogen bonds from results
all_hydrogen_bonds = [hb for result in results for hb in result]

# Define a format string for the output
output_format = "{:<15} {:<15} {:<15} {:<15} {:<15}\n"

# Write all identified hydrogen bonds across timesteps to the output file
with open(output_file, 'w') as f:
    for timestep in sorted(set(hb['timestep'] for hb in all_hydrogen_bonds)):
        f.write(f"Timestep {timestep}:\n")
        f.write(output_format.format("Donor[ID]", "Hydrogen[ID]", "Acceptor[ID]", "Distance (Å)", "Angle(°)"))
        for hb in [bond for bond in all_hydrogen_bonds if bond['timestep'] == timestep]:
            f.write(output_format.format(
                f"{hb['donor_type']}[{hb['donor_id']}]",
                f"{hb['hydrogen_type']}[{hb['hydrogen_id']}]",
                f"{hb['acceptor_type']}[{hb['acceptor_id']}]",
                f"{hb['distance']:.2f}",
                f"{hb['angle']:.2f}"
            ))
        f.write("\n")

print(f"Hydrogen bond identification complete. Results saved to {output_file}.")

