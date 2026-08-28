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
conda env create -f envs/grover.yml
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

**Download GROVER pretrained model weights**

The pretrained weights for the GROVER large model can be downloaded by running:

```bash
bash scripts/download_pretrained_grover.sh
```

This script will automatically download the weights and place them in the appropriate directory (`checkpoints/`).

**Note**: Downloading the pretrained weights is only necessary if you intend to re-run or resume fine-tuning.

---

**Download fine-tuned GROVER weights**

The solubility-finetuned GROVER weights can be downloaded by running:

```bash
bash scripts/download_finetuned_solubility_grover.sh
```

This script automatically downloads and extracts the model weights into the appropriate `model/` directory.

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

## Solubility prediction

If you want to reproduce the GROVER fine-tuning procedure for solubility prediction, follow **Steps 1 and 2**. If you only want to perform solubility inference using a fine-tuned model, skip to **Step 3**.

### 1. Download solubility data for fine-tuning

You can download the solubility datasets using the `download_solubility_data.sh` script. The datasets will be downloaded to `external/grover/solubility_data`.

**Processed data (required for fine-tuning)**

This includes the cleaned and preprocessed dataset ready for direct use in the pipeline.

```bash
bash scripts/download_solubility_data.sh --processed
```

**Raw datasets (optional)**

This contains the original, unprocessed datasets. It is mainly intended for users who want to:

- Inspect the original data
- Reproduce preprocessing steps
- Modify the preprocessing pipeline

```bash
bash scripts/download_solubility_data.sh --raw
```

The jupyter notebook `src/solubility/datasets_preprocessing.ipynb` documents how the original datasets were processed and merged into the final preprocessed dataset.

**All data**

Downloads both raw and processed datasets.

```bash
bash scripts/download_solubility_data.sh --all
```
---

### 2. GROVER fine-tuning for solubility prediction

Run the following script to fine-tune GROVER. There is no need to modify the file paths if the data has not been moved to another location.

```bash
bash scripts/finetune_solubility_grover.sh
```

### 3. Solubility inference using the fine-tuned model

To predict solubility for a new dataset using the fine-tuned GROVER model, run:

```bash
bash scripts/inference_solubility.sh <data_path> <output_path>
```

Where:
- `<data_path>` is the path to the input CSV file containing the molecules SMILES for which solubility will be predicted.
- `<output_path>` is the path to the output CSV file where the solubility predictions will be saved.

The script first generates the RDKit 2D features and then uses the fine-tuned model in `model` directory to predict solubility. The predictions are saved to the specified output path.

### 4. Computation of Solubility-finetuned GROVER embeddings

Compute atom-, bond-, and graph-level embeddings using the solubility-finetuned GROVER model:

```bash
bash scripts/compute_solubility_embeddings.sh <data_path> <features_path>
```

For example:

```bash
bash scripts/compute_solubility_embeddings.sh \
    external/grover/solubility_data/solubility_data.csv \
    external/grover/exampledata/finetune/solubility.npz
```

The resulting embeddings are saved to `grover_embeddings_solubility.pt`.


## Toxicity prediction

## Integrated score for candidates ranking
