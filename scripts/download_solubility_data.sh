#!/usr/bin/env bash
set -euo pipefail

mkdir -p external/grover/solubility_data

TAG="v2.0.0"
BASE_URL="https://github.com/m-baralt/TCruzi_pipeline/releases/download/$TAG"

download() {
  local file=$1

  if [ ! -f "external/grover/solubility_data/$file" ]; then
    echo "Downloading $file..."
    wget -q --show-progress "$BASE_URL/$file" -O "external/grover/solubility_data/$file"
  else
    echo "$file already exists, skipping."
  fi
}

case "${1:-}" in
  --raw)
    download "cui_et_al.csv"
    download "somas_gao_et_al.csv"
    ;;

  --processed)
    download "solubility_data.csv"
    ;;

  --all)
    download "cui_et_al.csv"
    download "somas_gao_et_al.csv"
    download "solubility_data.csv"
    ;;

  *)
    echo "Usage:"
    echo "  bash download_data.sh --raw"
    echo "  bash download_data.sh --processed"
    echo "  bash download_data.sh --all"
    exit 1
    ;;
esac