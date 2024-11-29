import os
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

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

# Parse the timestep
timestep = int(next(lines[i + 1].strip() for i, line in enumerate(lines) if line.startswith("ITEM: TIMESTEP")))

# Parse box bounds
def parse_box_bounds(lines):
    for i, line in enumerate(lines):
        if "ITEM: BOX BOUNDS" in line:
            x_bounds = list(map(float, lines[i + 1].split()))
            y_bounds = list(map(float, lines[i + 2].split()))
            z_bounds = list(map(float, lines[i + 3].split()))
            return (x_bounds, y_bounds, z_bounds)
    return None

# Extract box bounds
box_bounds = parse_box_bounds(lines)
x_bounds, y_bounds, z_bounds = box_bounds

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

grid_size_x, grid_size_y, grid_size_z = calculate_grid_sizes(box_bounds, total_grids=1000)

# Helper functions for indexing and parsing
def grid_index(coord, bounds, grid_size, grid_count):
    idx = int((coord - bounds[0]) / grid_size)
    return min(max(idx, 0), grid_count - 1)

# Use a dictionary for sparse grid storage
grid = {}

def parse_and_count_atoms_sparse(lines, grid, x_bounds, y_bounds, z_bounds, grid_size_x, grid_size_y, grid_size_z):
    for line in lines:
        if line.startswith("ITEM: ATOMS"):
            atom_data_start = lines.index(line) + 1
            break
    for line in lines[atom_data_start:]:
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

grid_sparse = parse_and_count_atoms_sparse(lines, grid, x_bounds, y_bounds, z_bounds, grid_size_x,  grid_size_y, grid_size_z)

# Function to determine if a cell is outer by checking neighbors
def is_outer_cell(x, y, z, grid):
    for dx, dy, dz in [(-1, 0, 0), (1, 0, 0), (0, -1, 0), (0, 1, 0), (0, 0, -1), (0, 0, 1)]:
        neighbor_key = (x + dx, y + dy, z + dz)
        if neighbor_key not in grid or grid[neighbor_key]["chitosan"] == 0:
            return True
    return False

# Identify inner and outer cells
outer_cells = set()
inner_cells = set()
for (gx, gy, gz), cell in grid_sparse.items():
    if cell["chitosan"] > 0:
        if is_outer_cell(gx, gy, gz, grid_sparse):
            outer_cells.add((gx, gy, gz))
        else:
            inner_cells.add((gx, gy, gz))

# Calculate water molecule counts and volumes
outer_water_count = sum(grid_sparse[cell]["water"] for cell in outer_cells)
inner_water_count = sum(grid_sparse[cell]["water"] for cell in inner_cells)
outer_volume = len(outer_cells) * (grid_size_x * grid_size_y * grid_size_z)
inner_volume = len(inner_cells) * (grid_size_x * grid_size_y * grid_size_z)
total_volume = outer_volume + inner_volume
total_water_count = outer_water_count + inner_water_count

# Save results to text file
output_file_path = os.path.join(output_dir, trajectory_file_path.rsplit('.', 1)[0] + '-chitosan-water-grid-analysis.txt')
with open(output_file_path, 'w') as f:
    f.write(f"Timestep: {timestep}\n")
    f.write(f"Grid Size: Gx = {grid_size_x:0.2f} Å, Gy = {grid_size_y:0.2f} Å, and Gz = {grid_size_z:0.2f} Å,\n")
    f.write(f"Outer Grid Volume: {outer_volume:0.2f} Å³\n")
    f.write(f"Inner Grid Volume: {inner_volume:0.2f} Å³\n")
    f.write(f"Total Grid Volume: {total_volume:0.2f} Å³\n")
    f.write(f"Outer Water Molecules: {outer_water_count}\n")
    f.write(f"Inner Water Molecules: {inner_water_count}\n")
    f.write(f"Total Water Molecules: {total_water_count}\n")

# Visualization functions
def draw_cube(ax, position, color, face_alpha=0.1, edge_alpha=1.0, edge_width=0.3):
    x, y, z = position
    vertices = [
        [(x, y, z), (x + grid_size_x, y, z), (x + grid_size_x, y + grid_size_y, z), (x, y + grid_size_y, z)],
        [(x, y, z + grid_size_z), (x + grid_size_x, y, z + grid_size_z), (x + grid_size_x, y + grid_size_y, z + grid_size_z), (x, y + grid_size_y, z + grid_size_z)],
        [(x, y, z), (x, y + grid_size_y, z), (x, y + grid_size_y, z + grid_size_z), (x, y, z + grid_size_z)],
        [(x + grid_size_x, y, z), (x + grid_size_x, y + grid_size_y, z), (x + grid_size_x, y + grid_size_y, z + grid_size_z), (x + grid_size_x, y, z + grid_size_z)],
        [(x, y, z), (x + grid_size_x, y, z), (x + grid_size_x, y, z + grid_size_z), (x, y, z + grid_size_z)],
        [(x, y + grid_size_y, z), (x + grid_size_x, y + grid_size_y, z), (x + grid_size_x, y + grid_size_y, z + grid_size_z), (x, y + grid_size_y, z + grid_size_z)]
    ]
    # Set facecolors with low alpha for transparent faces
    cube = Poly3DCollection(vertices, facecolors=color, edgecolor=color, alpha=face_alpha, linewidths=edge_width)
    cube.set_edgecolor((0, 0, 0, edge_alpha))  # Solid edges with specified edge alpha
    ax.add_collection3d(cube)

def plot_water_molecules(ax, grid_cells, grid, color, label, complement=False):
    water_points = [
        (x_bounds[0] + gx * grid_size_x + grid_size_x / 2,
         y_bounds[0] + gy * grid_size_y + grid_size_y / 2,
         z_bounds[0] + gz * grid_size_z + grid_size_z / 2)
        for gx, gy, gz in grid_cells if grid.get((gx, gy, gz), {}).get("water", 0) > 0
    ]

    # Choose complementary color if specified
    water_color = color
    if complement:
        if color == 'blue':  # Outer grid, use orange for water
            water_color = 'orange'
        elif color == 'green':  # Inner grid, use red for water
            water_color = 'red'

    if water_points:  # Only plot if there are water molecules in these cells
        wx, wy, wz = zip(*water_points)
        ax.scatter(wx, wy, wz, color=water_color, marker='o', s=10, label=label)

def plot_grids(inner_cells, outer_cells, grid, title, show_outer, show_inner):
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_title(title)

    # Remove axis spines
    ax.w_xaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
    ax.w_yaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
    ax.w_zaxis.line.set_color((1.0, 1.0, 1.0, 0.0))

    # Remove the axis planes
    ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))

    # Hide axes and labels
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_zlabel("")

    # Get simulation box dimensions
    Lx = x_bounds[1] - x_bounds[0]
    Ly = y_bounds[1] - y_bounds[0]
    Lz = z_bounds[1] - z_bounds[0]

    # Draw the simulation box with transparent faces
    box_vertices = [
        [x_bounds[0], y_bounds[0], z_bounds[0]],
        [x_bounds[1], y_bounds[0], z_bounds[0]],
        [x_bounds[1], y_bounds[1], z_bounds[0]],
        [x_bounds[0], y_bounds[1], z_bounds[0]],
        [x_bounds[0], y_bounds[0], z_bounds[1]],
        [x_bounds[1], y_bounds[0], z_bounds[1]],
        [x_bounds[1], y_bounds[1], z_bounds[1]],
        [x_bounds[0], y_bounds[1], z_bounds[1]]
    ]
    edges = [
        [box_vertices[i] for i in [0, 1, 2, 3]],  # Bottom face
        [box_vertices[i] for i in [4, 5, 6, 7]],  # Top face
        [box_vertices[i] for i in [0, 1, 5, 4]],  # Side edges
        [box_vertices[i] for i in [2, 3, 7, 6]],  # Side edges
        [box_vertices[i] for i in [1, 2, 6, 5]],  # Side edges
        [box_vertices[i] for i in [4, 7, 3, 0]]   # Side edges
    ]
    box = Poly3DCollection(edges, facecolors='white', edgecolor='black', alpha=0.0, linewidths=1)
    ax.add_collection3d(box)

    # Plot outer grids and their water molecules
    if show_outer:
        for gx, gy, gz in outer_cells:
            draw_cube(ax, (x_bounds[0] + gx * grid_size_x, y_bounds[0] + gy * grid_size_y, z_bounds[0] + gz * grid_size_z), color='blue')
        plot_water_molecules(ax, outer_cells, grid, color='blue', label=f'Outer Water Count: {outer_water_count}', complement=True)

    # Plot inner grids and their water molecules
    if show_inner:
        for gx, gy, gz in inner_cells:
            draw_cube(ax, (x_bounds[0] + gx * grid_size_x, y_bounds[0] + gy * grid_size_y, z_bounds[0] + gz * grid_size_z), color='green')
        plot_water_molecules(ax, inner_cells, grid, color='green', label=f'Inner Water Count: {inner_water_count}', complement=True)

    # Annotate the box with dimensions beside each corresponding side in Å
    ax.text((x_bounds[0] + x_bounds[1]) / 2, y_bounds[1], z_bounds[0], f"$L_x$: {Lx:.2f} Å", color='black')
    ax.text(x_bounds[1], (y_bounds[0] + y_bounds[1]) / 2, z_bounds[0], f"$L_y$: {Ly:.2f} Å", color='black')
    ax.text(x_bounds[0], y_bounds[0], (z_bounds[0] + z_bounds[1]) / 2, f"$L_z$: {Lz:.2f} Å", color='black')

    # Set the projection to Ortho
    ax.set_proj_type('ortho')

    # Only add legend if there are labels
    handles, labels = ax.get_legend_handles_labels()
    if labels:
        ax.legend()
    plt.show()

# Display each plot with appropriate labels
plot_grids(inner_cells, outer_cells, grid_sparse, "Outer Grids with Water Molecules", show_outer=True, show_inner=False)
plot_grids(inner_cells, outer_cells, grid_sparse, "Inner Grids with Water Molecules", show_outer=False, show_inner=True)
plot_grids(inner_cells, outer_cells, grid_sparse, "Combined Inner and Outer Grids with Water Molecules", show_outer=True, show_inner=True)
