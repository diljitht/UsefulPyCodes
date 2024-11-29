import sys
import numpy as np


class LammpsParser:
    def __init__(self, filename):
        self.filename = filename
        self.timesteps = []
        self.box_bounds = None

    def parse(self):
        with open(self.filename, "r") as file:
            timestep_data = None
            atom_section = False
            for line in file:
                if line.startswith("ITEM: TIMESTEP"):
                    timestep_data = {"atoms": []}
                    self.timesteps.append(timestep_data)
                elif line.startswith("ITEM: NUMBER OF ATOMS"):
                    timestep_data["number_of_atoms"] = int(next(file).strip())
                elif line.startswith("ITEM: BOX BOUNDS"):
                    bounds = [
                        list(map(float, next(file).strip().split())) for _ in range(3)
                    ]
                    timestep_data["box_bounds"] = bounds
                    self.box_bounds = bounds
                elif line.startswith("ITEM: ATOMS"):
                    atom_section = True
                    columns = line.strip().split()[
                        2:
                    ]  # Columns like id, mol, type, etc.
                elif atom_section:
                    values = line.strip().split()
                    atom = {
                        col: float(val) if "." in val else int(val)
                        for col, val in zip(columns, values)
                    }
                    timestep_data["atoms"].append(atom)
                else:
                    atom_section = False
        return self.timesteps


def unwrap_atoms(atoms, box_bounds):
    lx, ly, lz = [bound[1] - bound[0] for bound in box_bounds]
    unwrapped_atoms = []
    for atom in atoms:
        # Unwrap coordinates based on the atom's own image flags
        atom["xu"] = atom["x"] + atom["ix"] * lx
        atom["yu"] = atom["y"] + atom["iy"] * ly
        atom["zu"] = atom["z"] + atom["iz"] * lz
        unwrapped_atoms.append(atom)
    return unwrapped_atoms


def replicate_atoms(timesteps, box_bounds, num_waters=1500, num_chitosan_chains=25):
    lx, ly, lz = [bound[1] - bound[0] for bound in box_bounds]
    replicated_data = []
    atom_id_counter = 1  # Start atom ID counter
    chitosan_base_mol_id = (
        num_waters + num_chitosan_chains + 1
    )  # Start unique mol ID for replicated chitosan chains

    for timestep in timesteps:
        original_atoms = timestep["atoms"]

        # Separate chitosan and water atoms
        chitosan_atoms = [
            atom for atom in original_atoms if atom["mol"] <= num_chitosan_chains
        ]
        water_atoms = [
            atom for atom in original_atoms if atom["mol"] > num_chitosan_chains
        ]

        # Unwrap chitosan and water molecules before replication
        chitosan_atoms = unwrap_atoms(chitosan_atoms, box_bounds)
        water_atoms = unwrap_atoms(water_atoms, box_bounds)

        # Track which boxes already contain replicated chitosan
        chitosan_boxes = set()

        # Replicate water molecules and chitosan polymers based on image flags
        for atom in water_atoms:
            ix, iy, iz = atom["ix"], atom["iy"], atom["iz"]
            displacement = np.array([ix * lx, iy * ly, iz * lz])

            # Assign a new unique ID to the water atom while keeping its original mol ID
            atom["id"] = atom_id_counter
            atom_id_counter += 1
            replicated_data.append(atom)

            # Check if the box already has chitosan atoms replicated
            box_key = (ix, iy, iz)
            if box_key not in chitosan_boxes:
                chitosan_boxes.add(box_key)

                # Assign new unique molecular IDs to each chitosan chain in this replicated box
                for chain_index in range(num_chitosan_chains):
                    new_mol_id = (
                        chitosan_base_mol_id + chain_index
                    )  # Unique mol ID for each chain in the replicated box
                    chain_atoms = [
                        atom
                        for atom in chitosan_atoms
                        if atom["mol"] == chain_index + 1
                    ]

                    for chitosan_atom in chain_atoms:
                        displaced_atom = chitosan_atom.copy()
                        displaced_atom["xu"] += displacement[0]
                        displaced_atom["yu"] += displacement[1]
                        displaced_atom["zu"] += displacement[2]
                        displaced_atom['ix'] = ix
                        displaced_atom['iy'] = iy
                        displaced_atom['iz'] = iz
                        # Assign a new unique ID
                        displaced_atom["id"] = atom_id_counter
                        displaced_atom["mol"] = (
                            new_mol_id  # Assign new unique molecular ID for this chain
                        )
                        atom_id_counter += 1
                        replicated_data.append(displaced_atom)

                # Increment base mol ID for the next set of 25 chains in the next replicated box
                chitosan_base_mol_id += num_chitosan_chains

    return replicated_data


def write_replicated_data(filename, replicated_data, box_bounds):
    with open(filename, "w") as file:
        file.write("ITEM: TIMESTEP\n0\n")
        file.write(f"ITEM: NUMBER OF ATOMS\n{len(replicated_data)}\n")
        file.write("ITEM: BOX BOUNDS pp pp pp\n")
        for bound in box_bounds:
            file.write(f"{bound[0]} {bound[1]}\n")
        file.write("ITEM: ATOMS id mol type q xu yu zu ix iy iz\n")  # Include ix, iy, iz here
        for atom in replicated_data:
            file.write(
                f"{atom['id']} {atom['mol']} {atom['type']} {atom['q']} "
                f"{atom['xu']} {atom['yu']} {atom['zu']} "
                f"{atom['ix']} {atom['iy']} {atom['iz']}\n"  # Write ix, iy, iz
            )


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python " + sys.argv[0] + " <input_file.lmptrj>")
    else:
        try:
            filename = sys.argv[1]
        except FileNotFoundError:
            print(f"Error: The file '{filename}' was not found.")
            exit()

        output_filename = filename.rsplit('.', 1)[0] + "-replicated.lmptrj"
        
        parser = LammpsParser(filename)
        timesteps = parser.parse()
        box_bounds = parser.box_bounds
        replicated_data = replicate_atoms(timesteps, box_bounds, num_waters=1500, num_chitosan_chains=25)
        write_replicated_data(output_filename, replicated_data, box_bounds)

        print(f"Replicated data written to {output_filename}")
