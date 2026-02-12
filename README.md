# Tripanosoma Cruzi AKT-like kinase ligands

In this repository, a workflow to determine potential ligands for the Pleckstrin Homology domain in AKT-like kinase protein in T. Cruzi is made available. Appropriate ligands will fullfill a set of conditions:
- High affinity with the AKT-like kinase protein in T. Cruzi
- Low affinity with the AKT-like kinase protein in humans, thus low toxicity for humans
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

1. We need to add releases with the files required to built the database, create an sh for downloading it and change paths of db_creation.py accordingly. 

2. Once the data has been downloaded and placed in the appropriate directory, run:

`python src/db_creation.py` 

This creates a DuckDB database file.

**Note**: The resulting file requires approximately **263** GB of disk space.

3. Convert the database to Parquet:

`python src/duckdb_to_parquet.py` 

This produces a more storage-efficient Parquet version (**~61 GB**).

## T. Cruzi AKT-like Kinase virtual screening with DrugCLIP

1. Extract the structure of the T. Cruzi AKT-like Kinase from Protein Data Bank, accessible through the [8OZZ](https://www.rcsb.org/structure/8OZZ) identifier.
2. Run:

`python pocket2lmdb.py`

This converts the PDB file containing the raw structure to a cleaned lmdb file (required by DrugCLIP). 
Modified from [DrugCLIP screen pipeline repository](https://github.com/THU-ATOM/DrugCLIP_screen_pipeline/tree/main)).

3. Run DrugCLIP virtual screening using a modified file to avoid model weights incompatibilities.

Explanation: model weights were saved with Pytorch<=2.5, but in Pytorch>=2.6, the default call is `weights_only=True`, which only allows tensor objects, so non-tensor objects like `argparse.Namespace` return errors when loaded. `run_retrieval.py` avoids this behaviour. 

## Toxicity and solubility prediction using --


## Integrated score for candidates ranking
