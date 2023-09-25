import os
import sys
import yaml
import subprocess
import os
import platform

# config_path = sys.stdin.readline().strip()
# config = yaml.load(open(config_path, 'r'), Loader=yaml.SafeLoader)
config = yaml.load(open('config.yml', 'r'), Loader=yaml.SafeLoader)

# Input molecule
import re
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
# Conformation Search
from radonpy.core import poly, utils
from radonpy.sim import md, qm
from radonpy.sim.preset import eq 

import time
import pickle
import subprocess
import pathlib
import glob

# forcefields
from radonpy.ff.gaff2_mod import GAFF2_mod
from radonpy.ff.gaff2 import GAFF2

def time_record(command, name_of_process:str):
    time_sta = time.perf_counter()
    exec(command)
    time_end = time.perf_counter()
    tim = time_end- time_sta
    return {name_of_process:tim}


# Meta Data
log_name = config['logname']
work_dir = pathlib.Path(config['workdir'])
try:
    work_dir.mkdir()
except FileExistsError:
    exists = glob.glob(f'{work_dir}*')
    n = len(exists)
    work_dir = pathlib.Path(f'{work_dir}_{n+1}')
    work_dir.mkdir()

# Computational Setting
os.environ['LAMMPS_EXEC'] = config['lammps_path']
omp = config['omp']
mpi = config['mpi']
gpu_count = config['gpu_count']
mem = config['memory']

# Experimental Setting
ff = config['force_field'] # force field
if ff == 'gaff2_mod':
    ff = GAFF2_mod() # force field
elif ff == 'gaff2':
    ff = GAFF2() # force field
else:
    raise Exception('Error: Please select force field from "gaff2_mod" or "gaff2".')

charge = config['charge'] # charge assignment

# Polymer Definition
smiles = config['smiles']
mol = utils.mol_from_smiles(smiles)
ter = None
ter_head = config['head']
ter_tail = ter_head if config['tail'] == None else config['tail']
if ter_head == ter_tail:
    ter = utils.mol_from_smiles(ter_head)
else:
    ter_head = utils.mol_from_smiles(ter_head)
    ter_tail = utils.mol_from_smiles(ter_tail)
n_poly = config['n_poly']

# Conformation Search
time_log = {}
time_sta = time.perf_counter()
# gpu, ompは無しでやる

mol, energy = qm.conformation_search(mol, ff = ff, log_name=log_name,work_dir=work_dir,memory = mem, mpi = mpi)
time_end = time.perf_counter()
tim = time_end- time_sta
print('Conformation Search:',tim)
time_log['Conformation Search'] = tim
#%%
# Charge Assignment for the flagments
if charge is not None:
    time_sta = time.perf_counter()
    qm.assign_charges(mol, charge=charge, opt=False, work_dir=work_dir, omp=omp,  log_name=log_name)
    # terminals
    if ter:
        qm.assign_charges(ter, charge=charge, opt=False, work_dir=work_dir, omp=omp,  log_name=log_name+'_ter')
    else:
        qm.assign_charges(ter_head, charge=charge, opt=False, work_dir=work_dir, omp=omp,  log_name=log_name+'_ter_head')
        qm.assign_charges(ter_tail, charge=charge, opt=False, work_dir=work_dir, omp=omp,  log_name=log_name+'_ter_tail')
    time_end = time.perf_counter()
    tim = time_end- time_sta
    print('Charge Assignment:',tim)
    time_log['Charge Assignment'] = tim
else:
    print('Charge Assignment: Skipped')
    time_log['Charge Assignment'] = 0

# Generate polymer
time_sta = time.perf_counter()
poly_mol = poly.polymerize_rw(mol,n_poly)
if ter:
    poly_mol = poly.terminate_rw(poly_mol,ter)
else:
    poly_mol = poly.terminate_rw(poly_mol,ter_head,ter_tail)
time_end = time.perf_counter()
tim = time_end- time_sta
print('Polymerization:',tim)
time_log['Polymerization'] = tim

# FF assignment
time_sta = time.perf_counter()
result = ff.ff_assign(poly_mol) 
if not result: 
    print('[ERROR: Can not assign force field parameters.]') 
time_end = time.perf_counter()
tim = time_end- time_sta
print('FF Assignment:',tim)
time_log['FF Assignment'] = tim

# Parameters for AA cell
num_of_AA_atoms = config['n_atoms']
num_of_chains = config['n_chains']
density_of_cell = config['cell_density']
if num_of_chains == None and num_of_AA_atoms != None:
    num_of_chains = int(num_of_AA_atoms/Chem.AddHs(poly_mol).GetNumAtoms())
elif num_of_chains != None and num_of_AA_atoms == None:
    pass
else:
    raise Exception('Error: Please pass only one either "num_of_chains" or "num_of_AA_atom".')
# Calculation
# Construct unit cell (Amorphous)
time_sta = time.perf_counter()
ac = poly.amorphous_cell(poly_mol,num_of_chains,density = density_of_cell)
time_end = time.perf_counter()
tim = time_end- time_sta
print(f'Calculation Status: AA cell Generation\n\
Polymerization Degree:{n_poly}\n\
Num of Chains: {num_of_chains}\n\
Num of Atoms: {num_of_AA_atoms}\n\
Density: {density_of_cell}\n\
Calculation Time: {tim}')
time_log['AA Cell Generation'] = tim

time_sta = time.perf_counter()
eqmd = eq.EQ21step(ac, work_dir=work_dir) 

with open(f'{work_dir}/eqmd_in.pkl', 'wb') as f:
    pickle.dump(eqmd, f)

with open(f'{work_dir}/eqmd_in.pkl', 'rb') as f:
    eqmd = pickle.load(f)

ac_eq = eqmd.exec(temp=300, press=1.0, mpi=mpi)
time_end = time.perf_counter()
tim = time_end- time_sta
print(f'Calculation Status: Equilibration\n\
Polymerization Degree:{n_poly}\n\
Num of Chains: {num_of_chains}\n\
Num of Atoms: {num_of_AA_atoms}\n\
Density: {density_of_cell}\n\
Calculation Time: {tim}')
time_log['AA Equilibration'] = tim


with open('Time_log.pkl', 'wb') as f:
    pickle.dump(time_log, f)

# save log as mermaid gannt-chart
textlist = ['gannt',
            f'\ttitle Calculation Workflow for {num_of_chains} PS (Dp = {n_poly}, {num_of_AA_atoms} atoms)',
            '\tdataFormat ss',
            '\taxisFormat %H:%M']

for i, desc, log in enumerate(time_log.items()):
    if i == 0:
        text = f'\t{desc} : 00, {int(log)}s'
    else:
        text = f'\t{desc} : 00, {int(log)}s'
    textlist.append(text)
text = '\n'.join(textlist)
with open('LogGanttChart.mmd', 'w') as f:
    f.write(text)

results = {'BeforeEq':ac,'AfterEq': ac_eq}
with open('result.pkl', 'wb') as f:
    pickle.dump(results,f)

