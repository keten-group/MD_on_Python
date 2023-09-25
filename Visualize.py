#%%
def read_cg_output(filename):
    with open(filename, "r") as file:
        lines = file.readlines()

    beads = []
    for line in lines:
        data = line.strip().split()
        bead_id = int(data[0])
        atom_type = int(data[1])
        coords = [float(x) for x in data[2:5]]
        beads.append((bead_id, atom_type, coords))
    return beads

def write_xyz(beads, output_filename):
    atom_type_map = {1: "C", 2: "C"}  # Modify this based on your atom types

    with open(output_filename, "w") as out_file:
        out_file.write(f"{len(beads)}\n")
        out_file.write("CG_output.xyz\n")
        for bead in beads:
            atom_symbol = atom_type_map[bead[1]]
            out_file.write(f"{atom_symbol} {bead[2][0]} {bead[2][1]} {bead[2][2]}\n")

if __name__ == "__main__":
    cg_data_file = "CG_output.data"
    beads = read_cg_output(cg_data_file)
    write_xyz(beads, "CG_output.xyz")
# %%
