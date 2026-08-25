#!/bin/bash
mkdir -p model 

python external/grover/scripts/save_features.py --data_path external/grover/solubility_data/solubility_data.csv \
                                                --save_path external/grover/exampledata/finetune/solubility.npz \
                                                --features_generator rdkit_2d_normalized \
                                                --restart 

accelerate launch external/grover/main.py finetune \
    --data_path external/grover/solubility_data/solubility_data.csv \
    --save_dir model \
    --checkpoint_path external/grover/checkpoints/grover_large.pt \
    --features_path external/grover/exampledata/finetune/solubility.npz \
    --dataset_type regression \
    --split_type scaffold_balanced \
    --split_sizes 0.8 0.1 0.1 \
    --metric rmse \
    --ensemble_size 1 \
    --num_folds 1 \
    --ffn_hidden_size 200 \
    --batch_size 128 \
    --epochs 100 \
    --init_lr 0.00015 \
    --self_attention \
    --save_smiles_splits \
    --use_wandb \
    --wandb_entity mariabar \
    --num_workers 50 \
    2>&1 | tee logs_features.log