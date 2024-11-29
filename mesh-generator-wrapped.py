import os
import sys
import numpy as np
import pyvista as pv
from scipy.spatial import KDTree
import alphashape

def load_lammps_dump(file_path):
    """
    Load particle positions (x, y, z), types, and IDs from a LAMMPS dump file across multiple timesteps.
    """
    timesteps_data = []
    current_timestep = None
    positions = []
    types = []
    ids = []
    reading_atoms = False

    with open(file_path, "r") as file:
        for line in file:
            # Detect the start of a new timestep
            if line.startswith("ITEM: TIMESTEP"):
                # Save the data of the previous timestep, if any
                if current_timestep is not None and positions:
                    timesteps_data.append((current_timestep, np.array(positions), np.array(types), np.array(ids)))

                # Read the next line to get the timestep
                current_timestep = int(next(file).strip())
                positions = []
                types = []
                ids = []
                reading_atoms = False

            elif line.startswith("ITEM: ATOMS"):
                reading_atoms = True  # Start reading atom data

            elif reading_atoms and line.strip() and line[0].isdigit():
                data = line.split()

                # Ensure the line has enough elements before accessing data indices
                if len(data) >= 7:
                    atom_id = int(data[0])
                    particle_type = int(data[2])
                    x, y, z = float(data[4]), float(data[5]), float(data[6])
                    positions.append([x, y, z])
                    types.append(particle_type)
                    ids.append(atom_id)
                else:
                    print(f"Skipping line due to insufficient columns: {line.strip()}")

    # Handle the final timestep if the file doesn’t end with a header
    if current_timestep is not None and positions:
        timesteps_data.append((current_timestep, np.array(positions), np.array(types), np.array(ids)))

    return timesteps_data

def create_alpha_shape(positions, alpha=1.0):
    shape = alphashape.alphashape(positions, alpha)
    mesh = pv.wrap(shape).extract_surface()
    mesh.compute_normals(inplace=True)
    return mesh

def offset_mesh_inward(mesh, offset_distance=10.0):
    points = mesh.points
    center = points.mean(axis=0)
    directions = (points - center) / np.linalg.norm(points - center, axis=1)[:, None]
    inner_points = points - directions * offset_distance
    inner_mesh = pv.PolyData(inner_points).delaunay_3d().extract_surface()
    inner_mesh.compute_normals(inplace=True)
    return inner_mesh

def count_oxygen_in_volumes(inner_mesh, outer_mesh, positions, types, inner_threshold=1.0, shell_threshold=10.0, outer_limit_threshold=15.0):
    chitosan_volume = outer_mesh.volume
    inner_volume = inner_mesh.volume
    shell_volume = chitosan_volume - inner_volume
    oxygen_positions = positions[types == 12]

    inner_tree = KDTree(inner_mesh.points)
    outer_tree = KDTree(outer_mesh.points)

    inner_distances, _ = inner_tree.query(oxygen_positions)
    outer_distances, _ = outer_tree.query(oxygen_positions)

    inner_oxygen_count = np.sum(inner_distances < inner_threshold)
    shell_oxygen_count = np.sum((outer_distances < shell_threshold) & (inner_distances >= inner_threshold) & (outer_distances < outer_limit_threshold))

    return chitosan_volume, shell_volume, inner_volume, inner_oxygen_count, shell_oxygen_count

def save_snapshot(mesh, points, color, point_color, title, file_name):
    plotter = pv.Plotter(off_screen=True)
    plotter.add_mesh(mesh, color=color, opacity=0.5, label=title)
    plotter.add_points(points, color=point_color, point_size=5, render_points_as_spheres=True, label="Water Molecules")
    plotter.add_legend()
    plotter.screenshot(file_name)
    plotter.close()

def visualize_and_save_snapshots(inner_mesh, outer_mesh, positions, types, timestep, inner_threshold=1.0, shell_threshold=10.0, outer_limit_threshold=15.0):
    os.makedirs("snapshots", exist_ok=True)

    oxygen_positions = positions[types == 12]
    inner_tree = KDTree(inner_mesh.points)
    outer_tree = KDTree(outer_mesh.points)

    inner_distances, _ = inner_tree.query(oxygen_positions)
    outer_distances, _ = outer_tree.query(oxygen_positions)

    inner_oxygen_positions = oxygen_positions[inner_distances < inner_threshold]
    shell_oxygen_positions = oxygen_positions[(outer_distances < shell_threshold) &
                                              (inner_distances >= inner_threshold) &
                                              (outer_distances < outer_limit_threshold)]

    save_snapshot(inner_mesh, inner_oxygen_positions, color="yellow", point_color="red",
                  title="Inner Volume", file_name=f"snapshots/{file_base}-inner-volume-timestep-{timestep}.png")

    save_snapshot(outer_mesh, shell_oxygen_positions, color="lightblue", point_color="green",
                  title="Outer Shell Volume", file_name=f"snapshots/{file_base}-shell-volume-timestep-{timestep}.png")

    plotter = pv.Plotter(off_screen=True)
    plotter.add_mesh(outer_mesh, color="lightblue", opacity=0.5, label="Outer Shell Volume")
    plotter.add_mesh(inner_mesh, color="yellow", opacity=0.3, label="Inner Volume")
    plotter.add_points(inner_oxygen_positions, color="red", point_size=5, render_points_as_spheres=True, label="Inner Water Molecules")
    plotter.add_points(shell_oxygen_positions, color="green", point_size=5, render_points_as_spheres=True, label="Shell Water Molecules")
    plotter.add_legend()
    plotter.screenshot(f"snapshots/{file_base}-combined-volume-timestep-{timestep}.png")
    plotter.close()

if __name__ == "__main__":
    file_path = sys.argv[1]
    file_base = file_path.rsplit('.', 1)[0]
    alpha_value = 0.1
    inner_threshold = 1.0
    shell_threshold = 10.0
    outer_limit_threshold = 15.0

    timesteps_data = load_lammps_dump(file_path)

    with open(f"{file_base}-volume-counts.dat", "w") as f:
        for timestep, positions, types, ids in timesteps_data:
            chitosan_positions = positions[ids <= 11075]

            if chitosan_positions.size > 0:
                try:
                    outer_mesh = create_alpha_shape(chitosan_positions, alpha=alpha_value)
                    inner_mesh = offset_mesh_inward(outer_mesh, offset_distance=10.0)

                    chitosan_volume, shell_volume, inner_volume, inner_oxygen_count, shell_oxygen_count = count_oxygen_in_volumes(
                        inner_mesh, outer_mesh, positions, types, inner_threshold, shell_threshold, outer_limit_threshold
                    )

                    f.write(f"Timestep {timestep}: Chitosan Volume = {chitosan_volume:.2f}, "
                            f"Shell Volume = {shell_volume:.2f}, "
                            f"Inner Volume = {inner_volume:.2f}, "
                            f"Inner Water Count = {inner_oxygen_count}, "
                            f"Shell Water Count = {shell_oxygen_count}\n")

                    print(f"Timestep {timestep}:")
                    print(f"  Chitosan Volume = {chitosan_volume:.2f}")
                    print(f"  Shell Volume = {shell_volume:.2f}")
                    print(f"  Inner Volume = {inner_volume:.2f}")
                    print(f"  Inner Water Count = {inner_oxygen_count}")
                    print(f"  Shell Water Count = {shell_oxygen_count}")

                    visualize_and_save_snapshots(inner_mesh, outer_mesh, positions, types, timestep,
                                                 inner_threshold, shell_threshold, outer_limit_threshold)

                except ValueError as e:
                    print(f"Error processing timestep {timestep}: {e}")
            else:
                print(f"No valid chitosan particle positions found for timestep {timestep}.")
