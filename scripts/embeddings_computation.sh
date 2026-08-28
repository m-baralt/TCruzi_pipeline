#!/bin/bash

# Usage:
#   bash scripts/compute_solubility_embeddings.sh <data_path> <features_path>

DATA_PATH="$1"
FEATURES_PATH="$2"

if [ -z "$DATA_PATH" ] || [ -z "$FEATURES_PATH" ]; then
    echo "Usage: bash scripts/compute_solubility_embeddings.sh <data_path> <features_path>"
    exit 1
fi

python external/grover/main.py embeddings \
    --data_path "$DATA_PATH" \
    --checkpoint_path model/fold_0/model_0/model.pt \
    --features_path "$FEATURES_PATH" \
    --no_features_scaling \
    --output_path grover_embeddings_solubility.pt \
    --dataset_type regression \
    --metric rmse

