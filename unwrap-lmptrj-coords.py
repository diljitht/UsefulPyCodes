import sys


def convert_wrapped_to_unwrapped(input_filename, output_filename):
    with open(input_filename, "r") as infile, open(output_filename, "w") as outfile:
        box_bounds = []

        while True:
            line = infile.readline()
            if not line:
                break  # End of file

            if line.startswith("ITEM: TIMESTEP"):
                # Write timestep header and value
                outfile.write(line)
                outfile.write(infile.readline())  # Write timestep value

            elif line.startswith("ITEM: NUMBER OF ATOMS"):
                outfile.write(line)
                outfile.write(infile.readline())  # Write atom count

            elif line.startswith("ITEM: BOX BOUNDS"):
                # Capture the box bounds and write to file
                outfile.write(line)
                box_bounds = []
                for _ in range(3):
                    bounds = list(
                        map(float, infile.readline().strip().split()))
                    box_bounds.append(bounds)
                    outfile.write(f"{bounds[0]} {bounds[1]}\n")

            elif line.startswith("ITEM: ATOMS"):
                # Modify the header for atoms with unwrapped coordinates
                outfile.write("ITEM: ATOMS id mol type q xu yu zu\n")
                box_x = box_bounds[0][1] - box_bounds[0][0]
                box_y = box_bounds[1][1] - box_bounds[1][0]
                box_z = box_bounds[2][1] - box_bounds[2][0]

                # Read and process each atom line until we reach a new "ITEM:" header
                while True:
                    line = infile.readline()
                    if not line or line.startswith("ITEM:"):
                        break
                    parts = line.strip().split()
                    atom_id = int(parts[0])
                    mol = int(parts[1])
                    atom_type = int(parts[2])
                    q = float(parts[3])
                    x, y, z = map(float, parts[4:7])
                    ix, iy, iz = map(int, parts[7:10])

                    # Calculate unwrapped coordinates
                    xu = x + ix * box_x
                    yu = y + iy * box_y
                    zu = z + iz * box_z

                    # Write unwrapped data
                    outfile.write(f"{atom_id} {mol} {atom_type} {q} {xu:.6f} {yu:.6f} {zu:.6f}\n")

                # If we broke because of a new "ITEM:" header, put it back for the next loop iteration
                if line.startswith("ITEM:"):
                    infile.seek(infile.tell() - len(line))  # Move back to the start of this "ITEM:" line for reprocessing


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python " + sys.argv[0] + " <input_file.lmptrj>")
    else:
        try:
            input_file = sys.argv[1]
        except FileNotFoundError:
            print(f"Error: The file '{input_file}' was not found.")
            exit()

        output_file = input_file.rsplit('.', 1)[0] + "-unwrapped.lmptrj"
        convert_wrapped_to_unwrapped(input_file, output_file)
