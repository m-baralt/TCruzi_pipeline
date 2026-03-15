# Trypanosoma cruzi AKT-like kinase ligand discovery pipeline

This repository provides a workflow to identify potential ligands for the **Pleckstrin Homology (PH) domain** of the AKT-like kinase protein in *Trypanosoma cruzi*.

Candidate ligands are evaluated based on several criteria:

- **High affinity** for the AKT-like kinase PH domain
- **Low predicted toxicity** in humans
- **High predicted solubility**

The workflow integrates structure-based screening and downstream filtering to identify promising compounds.

---

## Setup Instructions

Clone the repository together with its submodules:

```bash
git clone --recurse-submodules git@github.com:m-baralt/TCruzi_pipeline.git
cd TCruzi_pipeline
```

If the repository was cloned **without submodules**, initialize them with:

```bash
git submodule update --init --recursive
```

---

**Environment setup**

Create the required conda environments:

```bash
conda env create -f envs/tc_pipeline.yml
conda env create -f envs/drugclip.yml
```

Activate the main pipeline environment:

```bash
conda activate tc_pipeline
```

---

**Download DrugCLIP model weights**

The DrugCLIP screening pipeline requires pretrained model weights.

Download them using the Hugging Face CLI:

```bash
hf download THU-ATOM/DrugCLIP_data model_weights.zip --repo-type dataset --local-dir external/Drug-The-Whole-Genome/data/
```

Extract the archive:

```bash
cd external/Drug-The-Whole-Genome/data/
unzip model_weights.zip
```

---

## Pipeline overview

1. Create molecular database from publicly available databases. 
2. Extract the protein pocket structure.
3. Convert the pocket structure to LMDB format.
4. Convert molecule databases to LMDB format.
5. Run DrugCLIP virtual screening.
6. Predict toxicity and solubility for each molecule. 
7. Prioritize candidates. 
   
---

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

#### 1. Download the protein structure

Extract the structure of the *T. cruzi* AKT-like kinase (or another protein/domain of interest) from the Protein Data Bank using the identifier:

**[8OZZ](https://www.rcsb.org/structure/8OZZ)**

Download the structure in **PDB format** and place it in the `data/` directory.

---

#### 2. Convert the protein pocket to LMDB format

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

#### 3. Convert molecules to LMDB format

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

#### 4. Run DrugCLIP virtual screening

Before running the screening, modify the configuration file if needed:

```bash
external/Drug-The-Whole-Genome/retrieval.sh
```

Launch the DrugCLIP virtual screening pipeline. The script is designed to run correctly without changing your current environment, so you can stay in the `tc_pipeline` environment:

```bash
bash scripts/run_drugclip.sh
```

#### 5. Output

The pipeline generates a `.txt` file (defined by save_path) containing:

- **SMILES** — molecule representation

- **Cosine similarity** — similarity between protein and molecule embeddings

- **Adjusted robust z-score** — normalized ranking score

- **Label** — extracted from the molecules LMDB database

If the `smiles2lmdb_opt.py` script is used, the label is **always** set to 0.

## Solubility prediction using GROVER and KA-GNN

## Toxicity prediction using GROVER and KA-GNN

## Integrated score for candidates ranking
