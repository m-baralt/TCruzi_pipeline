import os
import argparse
from Bio.PDB import PDBParser
import lmdb
import pickle

def write_lmdb(data, lmdb_path, start_index=0):
    env = lmdb.open(
        lmdb_path,
        subdir=False,
        readonly=False,
        lock=False,
        readahead=False,
        meminit=False,
        map_size=1099511627776
    )
    with env.begin(write=True) as txn:
        for i, d in enumerate(data):
            txn.put(str(start_index + i).encode('ascii'), pickle.dumps(d))
    return start_index + len(data)

def pocket_pdb_to_lmdb(pdb_file, pocket_name=None):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('pocket', pdb_file)
    atoms = [a for a in structure.get_atoms() if a.element != 'H']
    
    pocket_atoms = [a.element for a in atoms]
    pocket_coordinates = [a.coord for a in atoms]

    if pocket_name is None:
        pocket_name = os.path.basename(pdb_file).split('.')[0]

    return {
        'pocket': pocket_name,
        'pocket_atoms': pocket_atoms,
        'pocket_coordinates': pocket_coordinates
    }

def process_pocket_file(pdb_path, output_lmdb, pocket_name=None):
    pocket_dict = pocket_pdb_to_lmdb(pdb_path, pocket_name=pocket_name)
    write_lmdb([pocket_dict], output_lmdb)
    print(f"Wrote pdb file {pdb_path} to {output_lmdb}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert pocket PDB file into LMDB"
    )

    parser.add_argument(
        "--input",
        "-i",
        required=True,
        help="Path to PDB file"
    )

    parser.add_argument(
        "--output",
        "-o",
        required=True,
        help="Output LMDB file path"
    )

    parser.add_argument(
        "--pocket_name",
        "-p",
        required=False,
        help="Name of the pocket (default: derived from PDB filename)"
    )

    args = parser.parse_args()

    if not os.path.exists(args.input):
        raise ValueError("Input file does not exist")

    process_pocket_file(args.input, args.output)