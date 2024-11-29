import os
import sys

def split_timesteps(filename):
    # Create the 'splitted' directory if it doesn't exist
    output_dir = "splitted-trjs"
    os.makedirs(output_dir, exist_ok=True)
    
    with open(filename, 'r') as file:
        lines = file.readlines()
    
    timestep_data = []
    num_atoms = 0
    timestep = None
    filename_base = filename.rsplit('.', 1)[0]  # Base filename without extension
    
    line_idx = 0
    while line_idx < len(lines):
        # Check for new timestep
        if lines[line_idx].startswith("ITEM: TIMESTEP"):
            # Write the previous timestep data to a new file if there is any
            if timestep_data:
                output_filename = os.path.join(output_dir, f"{filename_base}-timestep-{timestep}.lmptrj")
                with open(output_filename, 'w') as out_file:
                    out_file.writelines(timestep_data)
                timestep_data = []
            
            # Begin collecting new timestep data
            timestep_data.append(lines[line_idx])    # ITEM: TIMESTEP
            timestep = lines[line_idx + 1].strip()    # timestep value
            timestep_data.append(lines[line_idx + 1]) # actual timestep number
            
            # Collect header and atom count
            timestep_data.extend(lines[line_idx + 2:line_idx + 9])  # headers (lines 2-9)
            num_atoms = int(lines[line_idx + 3].strip())  # number of atoms

            # Move to the start of atoms data
            line_idx += 9
            # Append the atom data
            timestep_data.extend(lines[line_idx:line_idx + num_atoms])
            line_idx += num_atoms  # Move to the next timestep
        else:
            line_idx += 1

    # Write last timestep data if present
    if timestep_data:
        output_filename = os.path.join(output_dir, f"{filename_base}-timestep-{timestep}.lmptrj")
        with open(output_filename, 'w') as out_file:
            out_file.writelines(timestep_data)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python " + sys.argv[0] + " <input_file.lmptrj>")
    else:
        try:
            filename = sys.argv[1]
        except FileNotFoundError:
            print(f"Error: The file '{filename}' was not found.")
            exit()

        split_timesteps(filename)

