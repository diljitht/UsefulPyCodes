import os
import sys
import numpy as np
import pyvista as pv
import alphashape


def load_lammps_dump(file_path, timestep=0):
    """
    Load particle positions (x, y, z) and types from a LAMMPS dump file for a specified timestep.
    """
    positions = []
    types = []
    reading_atoms = False
    current_timestep = None

    with open(file_path, "r") as file:
        for line in file:
            if "ITEM: TIMESTEP" in line:
                current_timestep = int(next(file).strip())
                reading_atoms = current_timestep == timestep
                continue

            if reading_atoms:
                if "ITEM: TIMESTEP" in line:
                    break
                if line.strip() and line[0].isdigit() and len(line.split()) >= 7:
                    data = line.split()
                    try:
                        # Assuming particle type is in the third column
                        particle_type = int(data[2])
                        x, y, z = float(data[4]), float(data[5]), float(data[6])
                        positions.append([x, y, z])
                        types.append(particle_type)
                    except ValueError:
                        continue

    return np.array(positions), np.array(types)


def create_alpha_shape(positions, alpha=1.0):
    """
    Create an alpha shape from the given positions.
    """
    if len(positions.shape) != 2 or positions.shape[1] != 3:
        raise ValueError("Positions must be a 2D array with shape (N, 3).")

    shape = alphashape.alphashape(positions, alpha)
    if shape.is_empty:
        raise ValueError("Alpha shape is empty. Try adjusting the alpha value.")

    mesh = pv.wrap(shape)
    return mesh


def get_color_for_type(particle_type, color_map):
    """
    Get a consistent color for each particle type. Generate a new color if not already in the map.
    """
    if particle_type not in color_map:
        color_map[particle_type] = pv.Color(np.random.rand(3))  # Assign a random color
    return color_map[particle_type]


def visualize_mesh_and_save_snapshot(mesh, positions, types, color_map, timestep):
    """
    Visualize the mesh and save the snapshot with data points colored by type for each timestep.
    """
    # Create the "snapshots" directory if it doesn't exist
    if not os.path.exists("snapshots"):
        os.makedirs("snapshots")

    plotter = pv.Plotter(off_screen=True)
    plotter.add_mesh(mesh, color="lightblue", opacity=0.5)

    # Plot points, colored by type
    for particle_type in np.unique(types):
        points_of_type = positions[types == particle_type]
        color = get_color_for_type(particle_type, color_map)
        plotter.add_points(points_of_type, color=color, point_size=5, render_points_as_spheres=True)

    # Save the screenshot with the timestep in the file name inside "snapshots" folder
    plotter.screenshot(f"snapshots/surface-mesh-timestep-{timestep}.png")
    plotter.close()


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python " + sys.argv[0] + " <input_file.lmptrj>")
    else:
        file_path = sys.argv[1]
        output_file = file_path.rsplit('.', 1)[0] + "-mesh-volume.dat"  # File to store the volume output

        # Clear the output file at the start
        open(output_file, "w").close()

        # Global color map to maintain consistent colors across timesteps
        color_map = {}

		for ts in range(0, 1000000, 10000):
		    desired_timestep = ts
		    alpha_value = 0.1
		
		    try:
		        positions, types = load_lammps_dump(file_path, timestep=desired_timestep)
		
		        if positions.size > 0:
		            try:
		                surface_mesh = create_alpha_shape(positions, alpha=alpha_value)
		
		                # Save volume to output file
		                volume = surface_mesh.volume if surface_mesh.n_points > 0 else 0
		                with open(output_file, "a") as f:
		                    f.write(f"{desired_timestep:>6}: {volume:.2f}\n")
		
		                # Visualize and save the mesh and points as a PNG
		                visualize_mesh_and_save_snapshot(surface_mesh, positions, types, color_map, desired_timestep)
		
		            except ValueError as e:
		                print(e)
		        else:
		            print(f"No valid particle positions found for timestep {desired_timestep}.")
		    except FileNotFoundError:
		        print(f"Error: The file '{file_path}' was not found.")
		    except Exception as e:
		        print(f"An unexpected error occurred: {e}")
