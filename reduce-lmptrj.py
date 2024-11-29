import sys


def reduce_lammps_trajectory(input_file, output_file, skip_steps):
    """
    Reduces a LAMMPS trajectory file by only including every `skip_steps`-th timestep.

    Args:
        input_file (str): Path to the input LAMMPS trajectory file.
        output_file (str): Path to the output reduced trajectory file.
        skip_steps (int): Number of timesteps to skip between saves.
    """
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        write_block = False
        timestep_count = 0

        for line in infile:
            if "ITEM: TIMESTEP" in line:
                # Check if we should write this block of data
                write_block = timestep_count % skip_steps == 0
                timestep_count += 1

            if write_block:
                outfile.write(line)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python reduce_trajectory.py <input_file> <skip_steps>")
        sys.exit(1)
    else:
        try:
            input_file = sys.argv[1]
        except FileNotFoundError:
            print(f"Error: The file '{input_file}' was not found.")
            exit()

        output_file = input_file.rsplit('.', 1)[0] + "-reduced.lmptrj"

        try:
            skip_steps = int(sys.argv[2])
        except ValueError as e:
            print(e)

        reduce_lammps_trajectory(input_file, output_file, skip_steps)
