#%%
import rdkit
import MDAnalysis as mda
import numpy as np
from pathlib import Path
#%%
# Import the strucuture of the repeating unit
# smiles of styrene monomer
monomer_smi = "*CC(*)c1ccccc1"
ter_smi = 'H' # terminal group of the polymers
num_poly = 10
# polymerization on *



monomer_mol = rdkit.Chem.MolFromSmiles(monomer_smi)
# beads definition
# alkyl group
alkyl_group = rdkit.Chem.MolFromSmiles("*CC(*)")
# separate polymer backbone and side chain


# Set the file paths
datapath = Path(Path(__file__).parent.absolute(), "data")
dumpfile = Path(datapath, "eq3_last.dump") # dump file
datafile = Path(datapath, "eq3_last.data") # data file
ter_smi = 'H' # terminal group of the polymers
# Load AA-MD data
u = mda.Universe(datafile, dumpfile, format="LAMMPSDUMP")

# Define atom types
phenyl_hydrogen_type = 1
phenyl_carbon_type = 2
alkyl_carbon_type = 3
alkyl_hydrogen_type_1 = 4
alkyl_hydrogen_type_2 = 5

# Select atoms
phenyl_hydrogens = u.select_atoms(f"type {phenyl_hydrogen_type}")
phenyl_carbons = u.select_atoms(f"type {phenyl_carbon_type}")
alkyl_carbons = u.select_atoms(f"type {alkyl_carbon_type}")
alkyl_hydrogens_1 = u.select_atoms(f"type {alkyl_hydrogen_type_1}")
alkyl_hydrogens_2 = u.select_atoms(f"type {alkyl_hydrogen_type_2}")
alkyl_hydrogens = alkyl_hydrogens_1 + alkyl_hydrogens_2

# Identify polymer chains by using residue numbers (starting from 1)
residue_ids = np.unique(u.residues.resindices) + 1
polymer_chains = [u.select_atoms(f"resid {res_id}") for res_id in residue_ids]

# Check if the chains are selected properly
print(f"Number of chains: {len(polymer_chains)}")
for i, chain in enumerate(polymer_chains):
    print(f"Chain {i+1} has {chain.n_atoms} atoms and {chain.n_residues} residues.")

# Define the bead mapping on each chain
for chain in polymer_chains:
    # Select the alkyl and phenyl groups
    alkyl_groups = chain.select_atoms("type 3 4 5")
    phenyl_groups = chain.select_atoms("type 1 2")

# Initialize lists for CG bead coordinates and types
cg_coords = []
cg_types = []

# Calculate centroids for alkyl beads
for i in range(0, len(alkyl_carbons), 2):
    alkyl_group_carbons = alkyl_carbons[i:i+2]
    alkyl_indices_str = ','.join(map(str, alkyl_group_carbons.atoms.indices))
    alkyl_group_hydrogens = alkyl_hydrogens.select_atoms(f"around 2 (bynum {' '.join(map(str, alkyl_group_carbons.atoms.indices))})")
    bead_coords = np.mean(np.vstack([*[c.position for c in alkyl_group_carbons], *[h.position for h in alkyl_group_hydrogens]]), axis=0)
    cg_coords.append(bead_coords)
    cg_types.append(1)

# Calculate centroids for phenyl beads
for pc in phenyl_carbons:
    bonded_hydrogens = pc.bonded_atoms.intersection(phenyl_hydrogens)
    bead_coords = np.mean(np.vstack([pc.position, *[h.position for h in bonded_hydrogens]]), axis=0)
    cg_coords.append(bead_coords)
    cg_types.append(2)

cg_coords = np.array(cg_coords)

# load size of cell from data file
with open(datafile, "r") as f:
    lines = f.readlines()
    for line in lines:
        if "xlo" in line:
            xlo, xhi = line.split()[0:2]
            xlo, xhi = float(xlo), float(xhi)
        if "ylo" in line:
            ylo, yhi = line.split()[0:2]
            ylo, yhi = float(ylo), float(yhi)
        if "zlo" in line:
            zlo, zhi = line.split()[0:2]
            zlo, zhi = float(zlo), float(zhi)

# Save CG-MD data to LAMMPS data file
with open("cg_data.lammps", "w") as f:
    f.write("LAMMPS CG data file\n\n")
    f.write(f"{len(cg_coords)} atoms\n\n")
    f.write("2 atom types\n\n") # number of bead types
    # box size (xlo, xhi, ylo, yhi, zlo, zhi)
    f.write(f"{xlo} {xhi} xlo xhi\n")
    f.write(f"{ylo} {yhi} ylo yhi\n")
    f.write(f"{zlo} {zhi} zlo zhi\n\n")
    f.write("Atoms\n\n")

    for i, (atype, coord) in enumerate(zip(cg_types, cg_coords), start=1):
        f.write(f"{i} {atype} {coord[0]} {coord[1]} {coord[2]}\n")

print("CG-MD data saved to cg_data.lammps")
# %%
'''
Please write the python code for bead mapping of polystyrene.
The step is as follows:
1. Select the atoms of the polymer chains
2. Identify the polymer chains by using residue numbers (starting from 1)
3. Define the bead mapping on each chain (Note that the beads are mapped on each repeating unit of the polymer)
3-1. Select the atoms of alkyl and phenyl groups.
3-2. 
'''