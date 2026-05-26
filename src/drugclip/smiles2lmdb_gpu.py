import os
from random import seed
import numpy as np
import lmdb
from rdkit import Chem
from rdkit.Chem import AllChem
from tqdm import tqdm
import pickle
from multiprocessing import Pool
import duckdb
from functools import partial
from pathlib import Path
import argparse
import time
from collections import deque
from rdkit.Chem.rdDistGeom import ETKDGv3
from nvmolkit.embedMolecules import EmbedMolecules as nvMolKitEmbed
from nvmolkit.types import HardwareOptions
from nvmolkit.mmffOptimization import MMFFOptimizeMoleculesConfs as nvMolKitMMFFOptimize
from rdkit.Chem import rdmolops

METALS = {
    # Alkali metals
    'Li', 'Na', 'K', 'Rb', 'Cs', 'Fr',
    # Alkaline earth metals
    'Be', 'Mg', 'Ca', 'Sr', 'Ba', 'Ra',
    # Transition metals
    'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
    'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
    'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
    # Post-transition metals
    'Al', 'Ga', 'In', 'Sn', 'Tl', 'Pb', 'Bi', 'Po', 'Fl', 'Lv',
    # Lanthanides
    'La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu',
    # Actinides
    'Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr'
}

def smiles_processing(smi, metals=METALS):
    mol = Chem.MolFromSmiles(smi)

    if mol is None:
        return None
    
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        return None
    
    atoms = mol.GetAtoms()

    for atom in atoms:
        sym = atom.GetSymbol()
        if sym in metals:
            return None
        if abs(atom.GetFormalCharge()) > 3:
            return None
        
    if '.' in smi:
        frags = rdmolops.GetMolFrags(mol, asMols=True)
        if not frags:
            return None
        mol = max(frags, key=lambda m: m.GetNumAtoms())

    try:
        mol = AllChem.AddHs(mol)
    except Exception:
        return None
    
    mol.SetProp("smiles", smi)

    return mol
        
def mol_conformation(mols_list, seed = 42, preprocessingThreads = 2, 
                     batchSize = 25, batchesPerGpu = 2, confsPerMolecule = 1):
    params = ETKDGv3()
    params.randomSeed = seed
    params.useRandomCoords = True

    hardware_opts = HardwareOptions(
        preprocessingThreads=preprocessingThreads,
        batchSize=batchSize,
        batchesPerGpu=batchesPerGpu,
    )

    nvMolKitEmbed(
        molecules=mols_list,
        params=params,
        confsPerMolecule=confsPerMolecule,
        maxIterations=-1,  # Automatic iteration calculation
        hardwareOptions=hardware_opts
    )

    mmff_hardware_opts = HardwareOptions(
        preprocessingThreads=preprocessingThreads,
        batchSize=0)

    #energies = nvMolKitMMFFOptimize(
        #molecules=mols_list,
        #maxIters=200,
        #nonBondedThreshold=100.0,
        #hardwareOptions=mmff_hardware_opts
    #)
    
    return mols_list

def mol2_3Dcoords(mol):
    if mol is not None:
        if mol.GetNumConformers()>0:
            coords = mol.GetConformer().GetPositions()
            smi = mol.GetProp("smiles") if mol.HasProp("smiles") else "UNKNOWN"
            atoms = [atom.GetSymbol() for atom in mol.GetAtoms()]
            data = pickle.dumps({'atoms': atoms, 'coordinates': [coords.astype(np.float32)], 'smi': smi, 'label': 0}, protocol=-1)
            return data
        else:
            smi = mol.GetProp("smiles") if mol.HasProp("smiles") else "UNKNOWN"
            print("FAILED:", smi)
            return None

def process_duckdb_to_lmdb(db_path, output_lmdb, batch_size=5000, 
                           n_molecules=None, seed=42, nthreads=8, gpu_batchsize=25,
                           batchesPerGpu=2, confsPerMolecule=1, commit_interval = 5000):
    """Extract SMILES from our DuckDB chemical DB, generate molecular conformations,
    write 3D coordinates to LMDB in batches"""

    output_path = Path(output_lmdb)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    con = duckdb.connect()

    # Open LMDB once
    env = lmdb.open(
        str(output_path),
        subdir=False,
        readonly=False,
        lock=True,
        readahead=False,
        meminit=False,
        map_size=int(5e11),
        writemap=True,
        map_async=True  
    )

    global_idx = 0
    filtered_count = 0
    times_per_smiles = {}

    if n_molecules is not None:
        result = con.execute(f"""
                             SELECT SMILES
                             FROM read_parquet('{db_path}')
                             LIMIT {n_molecules}""")
    else:
        result = con.execute(f"""
                             SELECT SMILES
                             FROM read_parquet('{db_path}')
                             """)
    try:
        txn = env.begin(write=True)
        pool = Pool(nthreads)
        for batch in result.fetch_record_batch(batch_size):
            batch_mols_list = []
            batch_smiles_list = batch.column(0).to_pylist()
            for mol in pool.imap_unordered(smiles_processing, batch_smiles_list, chunksize=200):
                if mol is not None:
                    batch_mols_list.append(mol)
            if len(batch_mols_list) == 0:
                continue
            batch_mols_list = mol_conformation(mols_list=batch_mols_list, seed=seed, preprocessingThreads=nthreads,
                                               batchSize=gpu_batchsize, batchesPerGpu=batchesPerGpu, confsPerMolecule=confsPerMolecule)
            
            for mol in batch_mols_list:
                data = mol2_3Dcoords(mol=mol)
                if data is not None:
                    key = global_idx.to_bytes(8, "big")
                    txn.put(key, data)
                    global_idx += 1
                else:
                    filtered_count += 1

                if global_idx > 0 and global_idx % commit_interval == 0:
                    txn.commit()    
                    txn = env.begin(write=True)

        pool.close()
        pool.join()

        txn.commit()
        print("Total processed:", global_idx)
        print("Total filtered:", filtered_count)
    finally:
        env.close()
        con.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert SMILES to LMDB"
    )

    parser.add_argument(
        "--output_lmdb",
        "-o",
        required=True,
        help="Output LMDB file path"
    )

    parser.add_argument(
        "--db_path",
        required=True,
        help="Path to parquet database"
    )

    parser.add_argument(
        "--batch_size",
        "-b",
        required=False,
        help="Number of records to process in each batch",
        type=int,
        default=5000
    )

    parser.add_argument(
        "--number_molecules",
        "-n",
        required=False,
        help="Number of molecules to process from the database",
        type=int,
        default=None
    )
    parser.add_argument(
        "--seed",
        required=False,
        help="Random seed for 3D coordinate generation",
        type=int,
        default=42
    )
    parser.add_argument(
        "--nthreads",
        required=False,
        help="Number of threads for parallel processing",
        type=int,
        default=8
    )
    parser.add_argument(
        "--gpu_batchsize",
        required=False,
        help="Batch size for Embed molecule in nvmolkit",
        type=int,
        default=25
    )
    parser.add_argument(
        "--batchesPerGpu",
        required=False,
        help="Number of batches per GPU for Embed molecule in nvmolkit",
        type=int,
        default=2
    )
    parser.add_argument(
        "--confsPerMolecule",
        required=False,
        help="Number of conformations per molecule",
        type=int,
        default=1
    )
    parser.add_argument(
        "--commit_interval",
        required=False,
        help="LMDB commit interval",
        type=int,
        default=5000
    )
    
    args = parser.parse_args()

    process_duckdb_to_lmdb(db_path = args.db_path, output_lmdb = args.output_lmdb, 
                           batch_size=args.batch_size, n_molecules=args.number_molecules, 
                           seed=args.seed, nthreads=args.nthreads, gpu_batchsize=args.gpu_batchsize,
                           batchesPerGpu=args.batchesPerGpu, confsPerMolecule=args.confsPerMolecule,
                           commit_interval=args.commit_interval)