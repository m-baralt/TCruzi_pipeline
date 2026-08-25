#!/bin/bash

DATA_PATH="$1"
OUTPUT_PATH="$2"

if [ -z "$DATA_PATH" ] || [ -z "$OUTPUT_PATH" ]; then
    echo "Usage: bash scripts/inference_solubility.sh <data_path> <output_path>"
    exit 1
fi

python external/grover/scripts/save_features.py \
    --data_path "$DATA_PATH" \
    --save_path external/grover/exampledata/finetune/solubility_inference.npz \
    --features_generator rdkit_2d_normalized \
    --restart

python external/grover/main.py predict \
    --data_path "$DATA_PATH" \
    --features_path external/grover/exampledata/finetune/solubility_inference.npz \
    --checkpoint_dir model_features \
    --no_features_scaling \
    --output "$OUTPUT_PATH"