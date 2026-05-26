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


def smi2_3Dcoords(mol,seed=42):
    mol = AllChem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    #params.maxAttempts = 1000
    #params.timeout = 60
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
    if len(smi) < 300:
        try:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                return None
            else:
                coordinate_list = smi2_3Dcoords(mol, seed=seed)
                if coordinate_list is None:
                    return None
                else:
                    mol = AllChem.AddHs(mol)
                    atoms = [atom.GetSymbol() for atom in mol.GetAtoms()]
                    data = pickle.dumps({'atoms': atoms, 'coordinates': coordinate_list, 'smi': smi, 'label': 0}, protocol=-1)
                    return data
        except:
            return None
    else:
        return None
    
def process_duckdb_to_lmdb(db_path, output_lmdb, smiles_col='SMILES', batch_size=5000, 
                           n_molecules=None, seed=42, nthreads=8):
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
    times_per_smiles = {}

    if n_molecules is not None:
        result = con.execute(f"""
                             SELECT {smiles_col}
                             FROM read_parquet('{db_path}')
                             LIMIT {n_molecules}""")
    else:
        result = con.execute(f"""
                             SELECT {smiles_col}
                             FROM read_parquet('{db_path}')
                             """)

    commit_interval = 50000
    try:
        txn = env.begin(write=True)
        for batch in result.fetch_record_batch(batch_size):
            smiles_list = batch.column(0).to_pylist()
            for smi in smiles_list:
                t_rdkit0 = time.time()
                output = smi2coords(smi, seed=seed)
                t_rdkit1 = time.time()
                rdkit_time = t_rdkit1-t_rdkit0
                t_lmdb = time.time()
                if output is not None:
                    key = global_idx.to_bytes(8, "big")
                    txn.put(key, output)
                    global_idx += 1
                else:
                    filtered_count += 1
                    rdkit_time = 100

                if global_idx > 0 and global_idx % commit_interval == 0:
                    txn.commit()    
                    txn = env.begin(write=True)
                t_lmdb1 = time.time()
                times_per_smiles[smi] = [rdkit_time, t_lmdb1-t_lmdb, len(smi)]


        txn.commit()
        print("Total processed:", global_idx)
        print("Total filtered:", filtered_count)
        with open('/home/mabarr/TCruzi_pipeline/test/smiles_times.pkl', 'wb') as f:
            pickle.dump(times_per_smiles, f)
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