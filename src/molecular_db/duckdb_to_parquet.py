import duckdb
import pyarrow as pa
import pyarrow.parquet as pq
import shutil

def free_gb(path="."):
    usage = shutil.disk_usage(path)
    return usage.free / (1024 ** 3)

con = duckdb.connect("chem_unified_db.db")

chunk_size = 5_000_000
start = 0
min_free_gb = 50
parquet_dir = "parquet_db"

while True:
    free = free_gb(parquet_dir)
    print(f"Free disk space: {free:.1f} GB")

    if free < min_free_gb:
        print("Low disk space, stopping safely.")
        break

    df = con.execute(f"""
        SELECT Standard_InChIKey, Source_Record_ID, SMILES, DB_Source
        FROM compounds
        WHERE rowid BETWEEN {start} AND {start + chunk_size - 1}
    """).fetchdf()

    if len(df) == 0:
        print("All rows processed.")
        break

    table = pa.Table.from_pandas(df)

    pq.write_to_dataset(
        table,
        root_path=parquet_dir,
        partition_cols=["DB_Source"],
        compression="zstd",
        compression_level=9
    )

    start += chunk_size
    print(f"Processed up to rowid {start}")

con.close()