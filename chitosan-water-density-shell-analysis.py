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

def compute_density(data, timestep, shell_thickness=5.0, is_chitosan=False):
    """Compute atom counts and normalized number density for water or chitosan molecules in concentric shells."""
    timestep_data = data[timestep]
    atoms = timestep_data["atoms"]  # Assuming atom columns are ordered appropriately
    bounds = timestep_data["box_bounds"]
    
    # Calculate the center of the simulation box
    center = [(b[0] + b[1]) / 2 for b in bounds]
    
    if is_chitosan:
        # Filter chitosan atoms (mol <= 25 and atom type in range 1 to 10)
        chitosan_atoms = atoms[(atoms[:, 1] <= 25) & (np.isin(atoms[:, 2], range(1, 11)))]
        filtered_atoms = chitosan_atoms
        total_atoms = len(chitosan_atoms)
    else:
        # Filter water atoms (mol > 25 and atom type == 12)
        water_atoms = atoms[(atoms[:, 1] > 25) & (atoms[:, 2] == 12)]
        filtered_atoms = water_atoms
        total_atoms = len(water_atoms)

    if total_atoms == 0:
        raise ValueError("No atoms match the criteria. Check input data or criteria.")

    # Calculate distances from the center
    distances = np.linalg.norm(filtered_atoms[:, 4:7] - center, axis=1)
    
    # Determine the maximum distance to cover the entire box
    max_distance = np.sqrt(3) * (bounds[0][1] - bounds[0][0]) / 2
    
    # Define shell edges and compute densities
    shell_edges = np.arange(0, max_distance + shell_thickness, shell_thickness)
    counts = []
    densities = []
    for i in range(len(shell_edges) - 1):
        shell_min = shell_edges[i]
        shell_max = shell_edges[i + 1]
        in_shell = (distances >= shell_min) & (distances < shell_max)
        count = len(distances[in_shell])
        counts.append(count)
        densities.append(count / total_atoms)  # Normalize by total atoms
    
    midpoints = (shell_edges[:-1] + shell_edges[1:]) / 2
    return midpoints, counts, densities

def compute_density_all_timesteps(data, shell_thickness=5.0):
    """Compute atom counts and densities for both water and chitosan for all timesteps and store results."""
    results = {}
    for timestep in data:
        midpoints_water, counts_water, densities_water = compute_density(data, timestep, shell_thickness, is_chitosan=False)
        midpoints_chitosan, counts_chitosan, densities_chitosan = compute_density(data, timestep, shell_thickness, is_chitosan=True)
        results[timestep] = {
            'midpoints': midpoints_water, 
            'counts_water': counts_water,
            'densities_water': densities_water,
            'counts_chitosan': counts_chitosan,
            'densities_chitosan': densities_chitosan
        }
    return results

def save_density_results(results, output_dir):
    """Save the density results to a file."""
    output_file_path = os.path.join(output_dir, trajectory_file_path.rsplit('.', 1)[0] + '-chitosan-water-density-vs-distance.txt')
    with open(output_file_path, 'w') as file:
        for timestep, data in results.items():
            file.write(f"Timestep: {timestep}\n")
            file.write(f"{'Distance (Å)':<15}{'Water Count':<15}{'Water Density':<15}{'Chitosan Count':<15}{'Chitosan Density':<15}\n")
            for dist, count_water, density_water, count_chitosan, density_chitosan in zip(
                data['midpoints'], data['counts_water'], data['densities_water'], 
                data['counts_chitosan'], data['densities_chitosan']):
                file.write(f"{dist:<15.2f}{count_water:<15}{density_water:<15.6f}{count_chitosan:<15}{density_chitosan:<15.6f}\n")
            file.write("\n")
    print(f"Density results saved to {output_file_path}")

def plot_density_combined(midpoints, densities_water, densities_chitosan, timestep):
    """Plot the densities of water and chitosan."""
    plt.figure(figsize=(8, 6))
    plt.plot(midpoints, densities_water, marker='o', linestyle='--', label='Water')
    plt.plot(midpoints, densities_chitosan, marker='s', linestyle='--', label='Chitosan')
    plt.xlabel('Distance from Center (Å)')
    plt.ylabel('Normalized Number Density')
    plt.title(f'Density vs Distance (Timestep {timestep})')
    plt.legend()
    plt.grid(True)
    plt.show()

# Parse the trajectory data
trajectory_data = parse_lammps_trajectory(trajectory_file_path)

# Compute densities and counts for all timesteps
all_density_results = compute_density_all_timesteps(trajectory_data)

# Save results to a file
save_density_results(all_density_results, output_dir)

# Plot for a specific timestep
specific_timestep = list(trajectory_data.keys())[-1]  # Example: last timestep
plot_density_combined(
    all_density_results[specific_timestep]['midpoints'],
    all_density_results[specific_timestep]['densities_water'],
    all_density_results[specific_timestep]['densities_chitosan'],
    specific_timestep
)

