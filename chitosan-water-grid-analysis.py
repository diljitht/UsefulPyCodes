import os
import sys

# Validate input arguments
if len(sys.argv) != 3:
    print("Usage: python " + sys.argv[0] + " <input_file.lmptrj> <NUM_WATER_MOLECULES>")
    exit()

# Parse input arguments
# Load the trajectory file
trajectory_file_path = sys.argv[1]
try:
    NUM_WATER_MOLECULES = int(sys.argv[2])
    if NUM_WATER_MOLECULES <= 0:
        raise ValueError("NUM_WATER_MOLECULES must be positive.")
except ValueError:
    print("Error: NUM_WATER_MOLECULES input is invalid.")
    exit()

if not os.path.exists(trajectory_file_path):
    print(f"Error: The file '{trajectory_file_path}' was not found.")
    exit()

# Directory for saving analysis files
output_dir = 'grid-analysis/'
os.makedirs(output_dir, exist_ok=True)

# Constants and settings
NUM_POLYMERS = 25
CHITOSAN_MOL_MAX = NUM_POLYMERS
WATER_MOL_MIN = CHITOSAN_MOL_MAX + 1
WATER_MOL_MAX = WATER_MOL_MIN + NUM_WATER_MOLECULES - 1

with open(trajectory_file_path, 'r') as file:
    lines = file.readlines()

# Function to parse box bounds from a timestep block
def parse_box_bounds(lines, start_idx):
    for i in range(start_idx, start_idx + 9):
        if "ITEM: BOX BOUNDS" in lines[i]:
            x_bounds = list(map(float, lines[i + 1].split()))
            y_bounds = list(map(float, lines[i + 2].split()))
            z_bounds = list(map(float, lines[i + 3].split()))
            return (x_bounds, y_bounds, z_bounds)
    return None

def calculate_grid_sizes(bounds, total_grids=1000):
    x_bounds, y_bounds, z_bounds = bounds
    x_range = x_bounds[1] - x_bounds[0]
    y_range = y_bounds[1] - y_bounds[0]
    z_range = z_bounds[1] - z_bounds[0]

    # Divide grids evenly across dimensions (assuming a cube root distribution)
    grids_per_dim = int(round(total_grids ** (1 / 3)))

    # Calculate grid sizes for each dimension
    grid_size_x = x_range / grids_per_dim
    grid_size_y = y_range / grids_per_dim
    grid_size_z = z_range / grids_per_dim

    return (grid_size_x, grid_size_y, grid_size_z)

# Define the grid_index function to convert coordinates to grid indices
def grid_index(coord, bounds, grid_size, grid_count):
    idx = int((coord - bounds[0]) / grid_size)
    return min(max(idx, 0), grid_count - 1)

# Function to parse atom data and update grid
def parse_and_count_atoms_sparse(lines, start_idx, grid, x_bounds, y_bounds, z_bounds, grid_size_x, grid_size_y, grid_size_z):
    atom_data_start = start_idx + 9  # Start after 9 header lines
    for line in lines[atom_data_start:]:
        if line.startswith("ITEM: TIMESTEP"):
            break
        parts = line.split()
        mol = int(parts[1])
        atom_type = int(parts[2])
        x, y, z = map(float, parts[4:7])
        gx = grid_index(x, x_bounds, grid_size_x, int((x_bounds[1] - x_bounds[0]) / grid_size_x))
        gy = grid_index(y, y_bounds, grid_size_y, int((y_bounds[1] - y_bounds[0]) / grid_size_y))
        gz = grid_index(z, z_bounds, grid_size_z, int((z_bounds[1] - z_bounds[0]) / grid_size_z))

        cell_key = (gx, gy, gz)
        if cell_key not in grid:
            grid[cell_key] = {"chitosan": 0, "water": 0}

        if 1 <= mol <= CHITOSAN_MOL_MAX:
            grid[cell_key]["chitosan"] += 1
        elif WATER_MOL_MIN <= mol <= WATER_MOL_MAX and atom_type == 12:
            grid[cell_key]["water"] += 1
    return grid

# Function to check if a cell is an outer cell by examining its neighbors
def is_outer_cell(x, y, z, grid):
    for dx, dy, dz in [(-1, 0, 0), (1, 0, 0), (0, -1, 0), (0, 1, 0), (0, 0, -1), (0, 0, 1)]:
        neighbor_key = (x + dx, y + dy, z + dz)
        if neighbor_key not in grid or grid[neighbor_key]["chitosan"] == 0:
            return True
    return False

# Main loop to process each timestep
start_idx = 0
while start_idx < len(lines):
    if lines[start_idx].startswith("ITEM: TIMESTEP"):
        # Parse the timestep
        timestep = int(lines[start_idx + 1].strip())
        print(f"Processing timestep {timestep}")

        # Parse the box bounds
        box_bounds = parse_box_bounds(lines, start_idx)
        if not box_bounds:
            print(f"Error: Could not parse box bounds at timestep {timestep}.")
            exit()
        x_bounds, y_bounds, z_bounds = box_bounds
        grid_size_x, grid_size_y, grid_size_z = calculate_grid_sizes(box_bounds, total_grids=1000)

        # Create a new grid for this timestep
        grid = {}

        # Parse and count atoms in the grid for the current timestep
        grid_sparse = parse_and_count_atoms_sparse(lines, start_idx, grid, x_bounds, y_bounds, z_bounds, grid_size_x,  grid_size_y, grid_size_z)

        # Identify inner and outer cells
        outer_cells = set()
        inner_cells = set()
        for (gx, gy, gz), cell in grid_sparse.items():
            if cell["chitosan"] > 0:
                if is_outer_cell(gx, gy, gz, grid_sparse):
                    outer_cells.add((gx, gy, gz))
                else:
                    inner_cells.add((gx, gy, gz))

        # Calculate water molecule counts and volumes for current timestep
        outer_water_count = sum(grid_sparse[cell]["water"] for cell in outer_cells)
        inner_water_count = sum(grid_sparse[cell]["water"] for cell in inner_cells)
        outer_volume = len(outer_cells) * (grid_size_x * grid_size_y * grid_size_z)
        inner_volume = len(inner_cells) * (grid_size_x * grid_size_y * grid_size_z)
        total_volume = outer_volume + inner_volume
        total_water_count = outer_water_count + inner_water_count

        # Save results to a text file for the current timestep
        output_file_path = os.path.join(output_dir, trajectory_file_path.rsplit('.', 1)[0] + f'-grid-analysis-timestep-{timestep}.txt')
        with open(output_file_path, 'w') as f:
            f.write(f"Timestep: {timestep}\n")
            f.write(f"Grid Size: Gx = {grid_size_x:0.2f} Å, Gy = {grid_size_y:0.2f} Å, and Gz = {grid_size_z:0.2f} Å,\n")
            f.write(f"Outer Grid Volume: {outer_volume:0.2f} Å³\n")
            f.write(f"Inner Grid Volume: {inner_volume:0.2f} Å³\n")
            f.write(f"Total Grid Volume: {total_volume:0.2f} Å³\n")
            f.write(f"Outer Water Molecules: {outer_water_count}\n")
            f.write(f"Inner Water Molecules: {inner_water_count}\n")
            f.write(f"Total Water Molecules: {total_water_count}\n")

        print(f"Saved analysis for timestep {timestep} to {output_file_path}")

        # Move to the next timestep block
        start_idx += 9 + len(grid_sparse)  # Skip past header and atom data
    else:
        start_idx += 1
