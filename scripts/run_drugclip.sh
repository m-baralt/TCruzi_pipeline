#!/bin/bash

set -e

conda run -n drugclip bash external/Drug-The-Whole-Genome/retrieval.sh "$@"
