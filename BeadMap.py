#%%
from pathlib import Path
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis.base import AnalysisFromFunction
from MDAnalysis.analysis import polymer
#%%
# Get current working directory
datapath = Path(Path(__file__).parent.absolute(), "data")

# Set the file paths
dumpfile = Path(datapath, "eq3_last.dump")
datafile = Path(datapath, "eq3_last.data")


# Create MDAnalysis Universe from the equilibrated data file and trajectory file
u = mda.Universe(datafile, dumpfile, format='LAMMPSDUMP')

# Identify alkyl and phenyl groups
alkyl_indices = u.select_atoms("type 1 2 3")  # Update atom types based on your simulation
phenyl_indices = u.select_atoms("type 4 5 6 7")  # Update atom types based on your simulation

# Identify polymer chains by using residue numbers
residue_ids = np.unique(u.residues.resindices)
polymer_chains = [u.select_atoms(f"resid {res_id}") for res_id in residue_ids]

# Check if the chains are selected properly
print(f"Number of chains: {len(polymer_chains)}")
for i, chain in enumerate(polymer_chains):
    print(f"Chain {i+1} has {chain.n_atoms} atoms and {chain.n_residues} residues.")

# Create CG beads for alkyl and phenyl groups
cg_beads = []
for chain in polymer_chains:
    for residue in chain.residues:
        alkyl_bead = residue.atoms.intersection(alkyl_indices)
        phenyl_bead = residue.atoms.intersection(phenyl_indices)
        if alkyl_bead:
            cg_beads.append(alkyl_bead.center_of_mass())
        if phenyl_bead:
            cg_beads.append(phenyl_bead.center_of_mass())

cg_beads = np.array(cg_beads)

# Save the CG bead positions as a new data file (cg_data.data)
with open("cg_data.data", "w") as outfile:
    outfile.write("LAMMPS Data File for Coarse-Grained Polystyrene\n\n")
    outfile.write(f"{len(cg_beads)} atoms\n")
    outfile.write("2 atom types\n")  # 1: alkyl, 2: phenyl
    outfile.write("\n")

    # Include box dimensions from the original data
    box = u.dimensions
    outfile.write(f"0.0 {box[0]} xlo xhi\n")
    outfile.write(f"0.0 {box[1]} ylo yhi\n")
    outfile.write(f"0.0 {box[2]} zlo zhi\n")
    outfile.write("\n")

    outfile.write("Atoms\n\n")
    for i, bead in enumerate(cg_beads):
        atom_type = 1 if i % 2 == 0 else 2  # 1: alkyl, 2: phenyl
        outfile.write(f"{i + 1} {atom_type} {bead[0]} {bead[1]} {bead[2]}\n")

print("Coarse-grained data file 'cg_data.data' has been created.")
# %%
import MDAnalysis as mda
import numpy as np


alkyl_indices = u.select_atoms("type 1 2 3")  # Update atom types based on your simulation
phenyl_indices = u.select_atoms("type 4 5 6 7")  # Update atom types based on your simulation

# Identify polymer chains by using residue numbers
residue_ids = np.unique(u.residues.resindices)
polymer_chains = [u.select_atoms(f"resid {res_id}") for res_id in residue_ids]

# Check if the chains are selected properly
print(f"Number of chains: {len(polymer_chains)}")
for i, chain in enumerate(polymer_chains):
    print(f"Chain {i+1} has {chain.n_atoms} atoms and {chain.n_residues} residues.")

# Create CG beads for alkyl and phenyl groups
cg_beads = []
for chain in polymer_chains:
    alkyl_beads = chain.select_atoms("type 1 2 3").positions.mean(axis=0)  # Update atom types based on your simulation
    phenyl_beads = chain.select_atoms("type 4 5 6 7").positions.mean(axis=0)  # Update atom types based on your simulation
    cg_beads.append(np.vstack((alkyl_beads, phenyl_beads)))

# Write LAMMPS data file with CG beads
with open("cg_data.data", "w") as f:
    f.write("LAMMPS data file for CG polystyrene\n")
    f.write("\n")
    f.write(f"{2*len(polymer_chains)} atoms\n")
    f.write("2 atom types\n")
    f.write("\n")
    f.write("0 100 xlo xhi\n")
    f.write("0 100 ylo yhi\n")
    f.write("0 100 zlo zhi\n")
    f.write("\n")
    f.write("Atoms\n")
    f.write("\n")

    atom_id = 1
    for chain_id, chain_beads in enumerate(cg_beads, start=1):
        for bead_id, bead_pos in enumerate(chain_beads, start=1):
            f.write(f"{atom_id} {chain_id} {bead_id} {bead_pos[0]} {bead_pos[1]} {bead_pos[2]}\n")
            atom_id += 1

print("Coarse-grained data file 'cg_data.data' has been created.")