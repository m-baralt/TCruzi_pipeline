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
from collections import defaultdict



def mol2_3Dcoords(smi, mol,seed=42):
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

    if mol is None:
        return None
    
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        return None
    
    atoms = mol.GetAtoms()

    for atom in atoms:
        sym = atom.GetSymbol()
        if sym in METALS:
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
    
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    #params.maxAttempts = 1000
    params.timeout = 1
    #params.pruneRmsThresh = 0.1
    #params.maxIterations = 200
    res = AllChem.EmbedMolecule(mol, params)
    if res != 0:
        return None
    else:
        if AllChem.MMFFHasAllMoleculeParams(mol):
            AllChem.MMFFOptimizeMolecule(mol,maxIters=25)
            coords = mol.GetConformer().GetPositions()
            return [coords.astype(np.float32)]
        else:
            return None
        
def smi2coords(smi, seed=42):
    try:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return None
        else:
            coordinate_list = mol2_3Dcoords(smi, mol, seed=seed)
            if coordinate_list is None:
                return None
            else:
                atoms = [atom.GetSymbol() for atom in mol.GetAtoms()]
                data = pickle.dumps({'atoms': atoms, 'coordinates': coordinate_list, 'smi': smi, 'label': 0}, protocol=-1)
                return data
    except:
        return None

    
def smi2coords_timed(smi):
    t0 = time.time()
    output = smi2coords(smi, seed=42)
    t1 = time.time()
    return output, (t1 - t0), smi

def process_duckdb_to_lmdb(db_path, output_lmdb, batch_size=5000, 
                           commit_interval=5000, n_molecules=None, 
                           seed=42, nthreads=8):
    """Stream from DuckDB and write to LMDB in batches"""
    
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

    pool = Pool(nthreads)
    txn = env.begin(write=True)
    for batch in result.fetch_record_batch(batch_size):
        smiles_list = batch.column(0).to_pylist()
        for result in pool.imap_unordered(smi2coords, smiles_list, chunksize=50):
            if result is not None:
                key = global_idx.to_bytes(8, "big")
                txn.put(key, result)
                global_idx += 1
            else:
                filtered_count += 1
            if global_idx > 0 and global_idx % commit_interval == 0:
                txn.commit()    
                txn = env.begin(write=True)
    txn.commit()
    pool.close()
    pool.join() 

    

    print("Total processed:", global_idx)
    print("Total filtered:", filtered_count)

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
        required=False,
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
    
    args = parser.parse_args()

    process_duckdb_to_lmdb(db_path = args.db_path, output_lmdb = args.output_lmdb, 
                           batch_size=args.batch_size, n_molecules=args.number_molecules, 
                           seed=args.seed, nthreads=args.nthreads)