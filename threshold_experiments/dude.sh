SCRIPT_DIR="/home/mabarr/TCruzi_pipeline/external/Drug-The-Whole-Genome/"
cd "$SCRIPT_DIR"


echo "First argument: $1"

FOLD_VERSION=6_folds
use_cache=False

for dir in data/DUD-E/*/; do
    target=$(basename "$dir")
    MOL_PATH="${dir}/mols.lmdb" 
    POCKET_PATH="${dir}/pocket.lmdb" 
    save_path="/home/mabarr/TCruzi_pipeline/threshold_experiments/normalizing/${target}.txt"

    CUDA_VISIBLE_DEVICES="1" python ./unimol/run_retrieval.py --user-dir ./unimol $data_path "./dict" --valid-subset test \
        --num-workers 60 --ddp-backend=c10d --batch-size 4 \
        --task drugclip-new --loss in_batch_softmax --arch drugclip  \
        --max-pocket-atoms 511 \
        --fp16 --fp16-init-scale 4 --fp16-scale-window 256  --seed 1 \
        --log-interval 100 --log-format simple \
        --mol-path $MOL_PATH \
        --pocket-path $POCKET_PATH \
        --fold-version $FOLD_VERSION \
        --use-cache $use_cache \
        --save-path $save_path
done