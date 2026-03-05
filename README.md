# Tripanosoma Cruzi AKT-like kinase ligands

In this repository, a workflow to determine potential ligands for the Pleckstrin Homology domain in AKT-like kinase protein in T. Cruzi is made available. Appropriate ligands will fullfill a set of conditions:
- High affinity with the AKT-like kinase protein in T. Cruzi
- Low toxicity in humans
- High solubility

To achieve this goal, several steps will be performed. 

## Setup Instructions

```bash
git clone --recurse-submodules git@github.com:m-baralt/TCruzi_pipeline.git
cd TCruzi_pipeline
```

If already cloned without submodules:

```bash
git submodule update --init --recursive
```

**Instruction for me:** I need to create env files tc_pipeline.yml and drugclip.yml. Instructions for users:

```bash
conda env create -f envs/tc_pipeline.yml
conda env create -f envs/drugclip.yml
```

To download model weights in the data directory run this in your repo directory:

```bash
hf download THU-ATOM/DrugCLIP_data model_weights.zip --repo-type dataset --local-dir external/Drug-The-Whole-Genome/data/
cd external/Drug-The-Whole-Genome/data/
unzip model_weights.zip
```

## Database construction

1. We need to add releases with the files required to built the database, create a bash file for downloading it and change paths of db_creation.py accordingly. 

2. Once the data has been downloaded and placed in the appropriate directory, run:

`python src/db_creation.py` 

This creates a DuckDB database file.

**Note**: The resulting file requires approximately **263** GB of disk space.

3. Convert the database to Parquet:

`python src/duckdb_to_parquet.py` 

This produces a more storage-efficient Parquet version (**~61 GB**).

## Virtual Screening with DrugCLIP

### 1. Download the protein structure

Extract the structure of the *T. cruzi* AKT-like kinase (or another protein/domain of interest) from the Protein Data Bank using the identifier:

**[8OZZ](https://www.rcsb.org/structure/8OZZ)**

Download the structure in **PDB format** and place it in the `data/` directory.

---

### 2. Convert the protein pocket to LMDB format

DrugCLIP requires the protein pocket structure to be stored in **LMDB format**.

Run the following command:

```bash
python src/pocket2lmdb.py -i data/8OZZ.pdb -o data/test_pocket.lmdb
```

To see all available arguments:

```bash
python src/pocket2lmdb.py --help
```

This script converts the raw PDB structure into a cleaned LMDB database compatible with DrugCLIP.

The implementation was modified from the DrugCLIP screening pipeline repository: https://github.com/THU-ATOM/DrugCLIP_screen_pipeline

---

### 3. Convert molecules to LMDB format

DrugCLIP also requires molecules to be stored in LMDB format.

To view all available arguments:

```bash
python src/smiles2lmdb_opt.py --help
```

**Convert molecules from a parquet database**

Example: convert 100 molecules from a parquet database:

```bash
python src/smiles2lmdb_opt.py \
    -o test/test_mols.lmdb \
    --db_path /home/mabarr/db/parquet_db/DB_Source\=ChEMBL/ \
    -n 100
```

**Convert molecules from a SMILES text file**

If you have a `.txt` file containing one SMILES string per line, run:

```bash
python src/smiles2lmdb_opt.py \
    -o test/test_mols.lmdb \
    --txt_path test/smiles_test.txt
```

### 4. Run DrugCLIP virtual screening

Before running the screening, modify the configuration file if needed:

```bash
external/Drug-The-Whole-Genome/retrieval.sh
```

Then launch the DrugCLIP virtual screening pipeline:

```bash
bash scripts/run_drugclip.sh
```

### 5. Output

The pipeline generates a `.txt` file (defined by save_path) containing:

- **SMILES** — molecule representation

- **Cosine similarity** — similarity between protein and molecule embeddings

- **Adjusted robust z-score** — normalized ranking score

- **Label** — extracted from the molecules LMDB database

If the `smiles2lmdb_opt.py` script is used, the label is **always** set to 0.

## Solubility prediction using GROVER and KA-GNN

## Toxicity prediction using GROVER and KA-GNN

## Integrated score for candidates ranking
