import os
from random import seed

import numpy as np
import lmdb
from rdkit import Chem
from rdkit.Chem import AllChem
from tqdm import tqdm
import pickle
from multiprocessing import Pool, Process, Queue
import duckdb
from functools import partial
from pathlib import Path
import argparse

ETKDG_PARAMS = AllChem.ETKDGv3()
ETKDG_PARAMS.randomSeed = 42
ETKDG_PARAMS.timeout = 1
ETKDG_PARAMS.maxIterations = 200
ETKDG_PARAMS.pruneRmsThresh = 0.1

def smi2_2Dcoords(smi):
    mol = Chem.MolFromSmiles(smi)
    mol = AllChem.AddHs(mol)
    AllChem.Compute2DCoords(mol)
    coordinates = mol.GetConformer().GetPositions().astype(np.float32)
    assert len(mol.GetAtoms()) == len(coordinates), "2D coordinates shape is not align with {}".format(smi)
    return coordinates


def smi2_3Dcoords(smi,seed=42):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return [smi2_2Dcoords(smi).astype(np.float32)]
    mol = AllChem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    res = AllChem.EmbedMolecule(mol, params)
    if res != 0:
        coords = smi2_2Dcoords(smi)
    else:
        if AllChem.MMFFHasAllMoleculeParams(mol):
            AllChem.MMFFOptimizeMolecule(mol)
        coords = mol.GetConformer().GetPositions()
    return [coords.astype(np.float32)]


def smi2coords(smi, seed=42):
    try:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return None
        if len(mol.GetAtoms()) > 400:
            coordinate_list = [smi2_2Dcoords(smi)]
        else:
            coordinate_list = smi2_3Dcoords(smi, seed=seed)
        mol = AllChem.AddHs(mol)
        atoms = [atom.GetSymbol() for atom in mol.GetAtoms()]
        data = pickle.dumps({'atoms': atoms, 'coordinates': coordinate_list, 'smi': smi, 'label': 0}, protocol=-1)
        return data
    except:
        return None

def process_duckdb_to_lmdb(db_path, output_lmdb, smiles_col='SMILES', batch_size=5000, 
                           n_molecules=1000000, seed=42, nthreads=8):
    """Stream from DuckDB and write to LMDB in batches"""
    
    output_path = Path(output_lmdb)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    con = duckdb.connect()
    
    # Open LMDB once
    env = lmdb.open(
        str(output_path),
        subdir=False,
        readonly=False,
        lock=False,
        readahead=False,
        meminit=False,
        max_readers=1,
        map_size=int(5e11),  
    )
    
    global_idx = 0

    #pool_func = partial(smi2coords_with_timeout, seed=seed, timeout=timeout)
    pool_func = partial(smi2coords, seed=seed)

    result = con.execute(f"""
    SELECT {smiles_col}
    FROM read_parquet('{db_path}')
    LIMIT {n_molecules}
""")

    with Pool(nthreads) as pool:
        for batch in result.fetch_record_batch(batch_size):
            smiles_list = batch.column(0).to_pylist()
            txn = env.begin(write=True)
            iterator = pool.imap_unordered(pool_func, smiles_list, chunksize=32)
            for i,output in enumerate(tqdm(iterator, total=len(smiles_list))):
                if output is not None:
                    txn.put(str(global_idx).encode(), output)
                    global_idx += 1
                if i % 1000 == 0:
                    txn.commit()
                    txn = env.begin(write=True)

            txn.commit()
        print("Total processed:", global_idx)

    env.close()
    con.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert SMILES to LMDB"
    )

    parser.add_argument(
        "--db_path",
        "-i",
        required=True,
        help="Path to parquet database"
    )

    parser.add_argument(
        "--output_lmdb",
        "-o",
        required=True,
        help="Output LMDB file path"
    )

    parser.add_argument(
        "--smiles_col",
        "-smi",
        required=False,
        help="Name of the SMILES column in the database",
        default="SMILES"
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
        default=1000000
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

    if not os.path.exists(args.db_path):
        raise ValueError("Input database does not exist")

    process_duckdb_to_lmdb(db_path = args.db_path, output_lmdb=args.output_lmdb, smiles_col=args.smiles_col, 
                           n_molecules=args.number_molecules, batch_size=args.batch_size, seed=args.seed, 
                           nthreads=args.nthreads)
    