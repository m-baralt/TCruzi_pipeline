#!/bin/bash

# Download the model ZIP from the latest GitHub release,
# extract it into the current directory, and remove the ZIP file.

REPO="m-baralt/TCruzi_pipeline"
ASSET="model.zip"

echo "Downloading $ASSET from the latest release..."

gh release download \
    --repo "$REPO" \
    --pattern "$ASSET" \
    --dir .

if [ $? -ne 0 ]; then
    echo "Error: failed to download $ASSET"
    exit 1
fi

echo "Extracting $ASSET..."

unzip -q "$ASSET"

if [ $? -ne 0 ]; then
    echo "Error: failed to extract $ASSET"
    exit 1
fi

rm "$ASSET"

echo "Model downloaded and extracted successfully."
