#!/usr/bin/env bash
set -euo pipefail

mkdir -p checkpoints

FILE="checkpoints/grover_large.pt"
URL="https://github.com/m-baralt/TCruzi_pipeline/releases/download/v1.0.0/grover_large.pt"

if [ ! -f "$FILE" ]; then
  echo "Downloading grover_large.pt..."
  wget -q --show-progress "$URL" -O "$FILE"
  echo "Download complete: $FILE"
else
  echo "File already exists: $FILE (skipping download)"
fi
