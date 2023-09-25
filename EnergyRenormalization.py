#%%
import MDAnalysis as mda

from ase.io import read, write
import numpy as np

#%%
from pathlib import Path
import MDAnalysis as mda
import numpy as np
# Get current working directory
datapath = Path(Path(__file__).parent.absolute(), "data")

# Set the file paths
trjfile = Path(datapath, "eq3.lammpstrj")
datafile = Path(datapath, "eq3_last.data")


# Load the LAMMPS trajectory file
u = mda.Universe(datafile, trjfile, format='LAMMPSDUMP')

# Create a new trajectory file for the coarse-grained (CG) simulation
with mda.Writer("cg_trajectory.xtc", n_atoms=u.atoms.n_atoms) as cg_traj:
    for ts in u.trajectory:
        # Assign beads to alkyl and phenyl groups in the polystyrene molecule
        alkyl_indices = u.select_atoms("type 1 2 3")  # Update the atom types based on your simulation
        phenyl_indices = u.select_atoms("type 4 5 6")  # Update the atom types based on your simulation

        # Calculate the center of mass for each group
        alkyl_com = alkyl_indices.center_of_mass()
        phenyl_com = phenyl_indices.center_of_mass()

        # Update the positions of the CG beads
        u.atoms.positions[alkyl_indices.indices] = alkyl_com
        u.atoms.positions[phenyl_indices.indices] = phenyl_com

        # Write the CG frame to the new trajectory
        cg_traj.write(u.atoms)
# %%
