import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# Load and process the trajectory file
if len(sys.argv) != 3:
    print("Usage: python " + sys.argv[0] + " <input_file.lmptrj> <Nwater>")
    exit()

# Parse input arguments
trajectory_file_path = sys.argv[1]
try:
    Nwater = int(sys.argv[2])
    if Nwater <= 0:
        raise ValueError("Nwater must be positive.")
except ValueError:
    print("Error: Nwater must be a positive number.")
    exit()

if not os.path.exists(trajectory_file_path):
    print(f"Error: The file '{trajectory_file_path}' was not found.")
    exit()

# Directory for saving analysis files
output_dir = f'water-density'                       
os.makedirs(output_dir, exist_ok=True)

def parse_lammps_trajectory(trajectory_file_path):
    """Parse the LAMMPS trajectory file."""
    with open(trajectory_file_path, 'r') as file:
        lines = file.readlines()

    data = {}
    i = 0

    while i < len(lines):
        if "ITEM: TIMESTEP" in lines[i]:
            timestep = int(lines[i + 1].strip())
            data[timestep] = {"atoms": []}
            i += 2

        elif "ITEM: NUMBER OF ATOMS" in lines[i]:
            num_atoms = int(lines[i + 1].strip())
            data[timestep]["num_atoms"] = num_atoms
            i += 2

        elif "ITEM: BOX BOUNDS" in lines[i]:
            bounds = []
            for j in range(3):
                bounds.append([float(x) for x in lines[i + 1 + j].strip().split()])
            data[timestep]["box_bounds"] = bounds
            i += 4

        elif "ITEM: ATOMS" in lines[i]:
            atom_lines = lines[i + 1 : i + 1 + num_atoms]
            atom_data = []
            for line in atom_lines:
                atom_data.append([float(x) if '.' in x or 'e' in x else int(x) for x in line.split()])
            data[timestep]["atoms"] = np.array(atom_data)
            i += 1 + num_atoms

        else:
            i += 1

    return data

def compute_density(data, timestep, shell_thickness=5.0, atom_type=12, mol_threshold=25):
    """Compute number density of water molecules in concentric shells."""
    timestep_data = data[timestep]
    atoms = timestep_data["atoms"]
    bounds = timestep_data["box_bounds"]
    
    # Calculate the center of the simulation box
    center = [(b[0] + b[1]) / 2 for b in bounds]
    
    # Filter water molecules (mol > 25 and type = 12)
    water_atoms = atoms[(atoms[:, 2] == atom_type) & (atoms[:, 1] > mol_threshold)]
    
    # Calculate distances from the center
    distances = np.linalg.norm(water_atoms[:, 4:7] - center, axis=1)
    
    # Determine the maximum distance to cover the entire box
    max_distance = np.sqrt(3) * (bounds[0][1] - bounds[0][0]) / 2
    
    # Define shell edges and compute densities
    shell_edges = np.arange(0, max_distance + shell_thickness, shell_thickness)
    densities = []
    for i in range(len(shell_edges) - 1):
        shell_min = shell_edges[i]
        shell_max = shell_edges[i + 1]
        in_shell = (distances >= shell_min) & (distances < shell_max)
        # volume = (4 / 3) * np.pi * (shell_max**3 - shell_min**3)
        densities.append(len(distances[in_shell]) / Nwater)
    
    midpoints = (shell_edges[:-1] + shell_edges[1:]) / 2
    return midpoints, densities

def compute_density_all_timesteps(data, shell_thickness=5.0, atom_type=12, mol_threshold=25):
    """Compute densities for all timesteps and store results."""
    results = {}
    for timestep in data:
        midpoints, densities = compute_density(data, timestep, shell_thickness, atom_type, mol_threshold)
        results[timestep] = {'midpoints': midpoints, 'densities': densities}
    return results

def plot_density(midpoints, densities, timestep):
    """Plot the density vs distance."""
    plt.figure(figsize=(8, 6))
    plt.plot(midpoints, densities, marker='o', linestyle='-')
    plt.xlabel('Distance from Center (Å)')
    plt.ylabel('Number Density (molecules/Å³)')
    plt.title(f'Water Molecule Density vs Distance (Timestep {timestep})')
    plt.grid(True)
    plt.show()

def save_density_results_aligned(results, trajectory_file_path, shell_thickness):
    """Save the density results in aligned format."""
    with open(trajectory_file_path, 'w') as file:
        for timestep, data in results.items():
            file.write(f"Timestep: {timestep}\n")
            file.write(f"Shell Thickness: {shell_thickness:.1f} Å\n")
            file.write(f"{'Distance (Å)':<15}{'Number of Water Molecules':<26}{'Density (molecules/Å³)':<25}\n")
            for distance, density in zip(data['midpoints'], data['densities']):
                num_water = int(density * Nwater)
                file.write(f"{distance:<15.2f}{num_water:<26}{density:<25.6f}\n")
            file.write("\n")
        print(f"Density results saved to {output_file_path}")

trajectory_data = parse_lammps_trajectory(trajectory_file_path)

# Compute density for all timesteps
all_density_results = compute_density_all_timesteps(trajectory_data)

# Save results to a neatly aligned text file
output_file_path = os.path.join(output_dir, trajectory_file_path.rsplit('.', 1)[0] + '-water-density-vs-distance.txt')     
save_density_results_aligned(all_density_results, output_file_path, shell_thickness=5.0)

# Plot for a specific timestep
specific_timestep = list(trajectory_data.keys())[0]  # Example: first timestep
plot_density(
    all_density_results[specific_timestep]['midpoints'],
    all_density_results[specific_timestep]['densities'],
    specific_timestep
)
