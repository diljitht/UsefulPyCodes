import sys


def sort_and_skip_end_lammps_trajectory(input_file, output_file, atoms_per_frame=11075, header_lines=9, atoms_to_skip=0):
    frame_number = 0  # Track frame number

    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        while True:
            # Read the header lines for the current frame
            header = [infile.readline() for _ in range(header_lines)]
            if not header[0]:  # End of file
                break

            # Modify the 4th line (atoms count) with the adjusted atoms count
            adjusted_atoms_count = atoms_per_frame - atoms_to_skip
            header[3] = f"{adjusted_atoms_count}\n"

            # Write the header lines for the current frame to the output file
            for line in header:
                outfile.write(line)

            # Read atom data for the current frame
            data = []
            for _ in range(atoms_per_frame):
                line = infile.readline()
                if not line:  # End of file
                    break
                data.append(line.strip().split())

            # Sort the atom data by the atom ID (first column)
            data.sort(key=lambda x: int(x[0]))

            # Skip the last `atoms_to_skip` atoms from the sorted data
            data = data[:-atoms_to_skip] if atoms_to_skip > 0 else data

            # Write the sorted and filtered data for the current frame to the output file
            for line in data:
                outfile.write(" ".join(line) + "\n")

            # Increment the frame number after processing the frame
            frame_number += 1
            print(f"Processed frame {frame_number}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python " + sys.argv[0] + " <input_file.lmptrj>")
    else:
        try:
            input_file = sys.argv[1]
        except FileNotFoundError:
            print(f"Error: The file '{input_file}' was not found.")
            exit()

        output_file = input_file.rsplit('.', 1)[0] + "-sorted-skipped-end.lmptrj"

        # Read the atom count from the header of the first frame to determine atoms_per_frame
        with open(input_file, "r") as infile:
            header_lines = [infile.readline() for _ in range(9)]
            atoms_per_frame = int(header_lines[3].strip())  # Assuming atoms_per_frame is on line 4

        atoms_to_skip = atoms_per_frame - 11075  # Number of atoms to skip from the end in each frame

        sort_and_skip_end_lammps_trajectory(input_file, output_file, atoms_per_frame=atoms_per_frame, header_lines=9, atoms_to_skip=atoms_to_skip)
