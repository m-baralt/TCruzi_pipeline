import duckdb

# BindingDB downloaded 29/01/2026
# Coconut downloaded 29/01/2026
# Chembl downloaded 30/01/2026
# PubChem downloaded 30/01/2026

con = duckdb.connect("chem_unified_db.db")

# canonical compound table (1 row per InChIKey)
con.execute("""
CREATE OR REPLACE TABLE compounds (
    Standard_InChIKey TEXT PRIMARY KEY,
    Source_Record_ID TEXT,
    SMILES TEXT,
    DB_Source TEXT
);
""")

# statistics table
con.execute("""
CREATE OR REPLACE TABLE db_stats (
    db_name TEXT,
    total_rows BIGINT,
    unique_molecules BIGINT,
    duplicated_rows BIGINT,
    overlap_with_existing BIGINT,
    timestamp TIMESTAMP DEFAULT now()
);
""")

# Binding DB

print("Building Binding DB ...")

con.execute("""
CREATE TEMP TABLE binding_tmp AS
SELECT
    "Ligand InChI Key" AS Standard_InChIKey,
    "Ligand InChI" AS Standard_Inchi,
    "Ligand SMILES" AS SMILES,
    "BindingDB Reactant_set_id" AS Source_Record_ID,
    'BindingDB' AS DB_Source
FROM read_csv_auto('BindingDB_All.tsv')
WHERE "Ligand InChI Key" IS NOT NULL;
""")

con.execute("""
INSERT INTO db_stats
SELECT
    'BindingDB',
    COUNT(*),
    COUNT(DISTINCT Standard_InChIKey),
    COUNT(*) - COUNT(DISTINCT Standard_InChIKey),
    (
        SELECT COUNT(DISTINCT b.Standard_InChIKey)
        FROM binding_tmp b
        JOIN compounds c
        ON b.Standard_InChIKey = c.Standard_InChIKey
    ),
    now()
FROM binding_tmp;
""")

con.execute("""
INSERT INTO compounds
SELECT
    Standard_InChIKey,
    Source_Record_ID,
    SMILES,
    DB_Source
FROM (
    SELECT
        Standard_InChIKey,
        Source_Record_ID,
        SMILES,
        DB_Source,
        ROW_NUMBER() OVER (
            PARTITION BY Standard_InChIKey
            ORDER BY DB_Source, Source_Record_ID
        ) AS rn
    FROM binding_tmp
    WHERE Standard_InChIKey NOT IN (
        SELECT Standard_InChIKey FROM compounds
    )
)
WHERE rn = 1;
""")

con.execute("DROP TABLE binding_tmp;")

print("DONE!")

print("Building Coconut DB ...")

# Coconut
con.execute("""
CREATE TEMP TABLE coconut_tmp AS
SELECT
    standard_inchi_key AS Standard_InChIKey,
    standard_inchi AS Standard_Inchi,
    canonical_smiles AS SMILES,
    identifier AS Source_Record_ID,
    'Coconut' AS DB_Source
FROM read_csv_auto('coconut_csv_lite-01-2026.csv')
WHERE standard_inchi_key IS NOT NULL;
""")

con.execute("""
INSERT INTO db_stats
SELECT
    'Coconut',
    COUNT(*),
    COUNT(DISTINCT Standard_InChIKey),
    COUNT(*) - COUNT(DISTINCT Standard_InChIKey),
    (
        SELECT COUNT(DISTINCT b.Standard_InChIKey)
        FROM coconut_tmp b
        JOIN compounds c
        ON b.Standard_InChIKey = c.Standard_InChIKey
    ),
    now()
FROM coconut_tmp;
""")

con.execute("""
INSERT INTO compounds
SELECT
    Standard_InChIKey,
    Source_Record_ID,
    SMILES,
    DB_Source
FROM (
    SELECT
        Standard_InChIKey,
        Source_Record_ID,
        SMILES,
        DB_Source,
        ROW_NUMBER() OVER (
            PARTITION BY Standard_InChIKey
            ORDER BY DB_Source, Source_Record_ID
        ) AS rn
    FROM coconut_tmp
    WHERE Standard_InChIKey NOT IN (
        SELECT Standard_InChIKey FROM compounds
    )
)
WHERE rn = 1;
""")

con.execute("DROP TABLE coconut_tmp;")

print("DONE!")

print("Building ChEMBL DB...")

# ChEMBL
con.execute("""
CREATE OR REPLACE TABLE chembl_part1 AS
SELECT
    "Inchi Key"    AS Standard_InChIKey,
    "Inchi"        AS Standard_Inchi,
    "Smiles"       AS SMILES,
    "ChEMBL ID"    AS Source_Record_ID,
    'ChEMBL'       AS DB_Source
FROM read_csv(
    'DOWNLOAD-heRwUfDRJj-nAMSjmA31y8RyPmnAG6YgOpVpkx4anHU_eq_.csv',
    delim=';',
    quote='"',
    escape='"',
    header=true,
    strict_mode=false,
    ignore_errors=true,
    max_line_size=10000000
)
WHERE "Inchi Key" IS NOT NULL
""")

con.execute("""
CREATE OR REPLACE TABLE chembl_part2 AS
SELECT
    "column23" AS Standard_InChIKey,
    "column24" AS Standard_Inchi,
    "column22" AS SMILES,
    "column00"  AS Source_Record_ID,
    'ChEMBL' AS DB_Source
FROM read_csv(
    'DOWNLOAD-heRwUfDRJj-nAMSjmA31y8RyPmnAG6YgOpVpkx4anHU_eq__part2.csv',
    delim=';',
    quote='"',
    escape='"',
    header=false,
    strict_mode=false,
    ignore_errors=true,
    max_line_size=10000000
)
WHERE column23 IS NOT NULL
""")

con.execute("""
CREATE TEMP TABLE chembl_tmp AS
SELECT * FROM chembl_part1
UNION ALL
SELECT * FROM chembl_part2
""")

con.execute("""
INSERT INTO db_stats
SELECT
    'ChEMBL',
    COUNT(*),
    COUNT(DISTINCT Standard_InChIKey),
    COUNT(*) - COUNT(DISTINCT Standard_InChIKey),
    (
        SELECT COUNT(DISTINCT b.Standard_InChIKey)
        FROM chembl_tmp b
        JOIN compounds c
        ON b.Standard_InChIKey = c.Standard_InChIKey
    ),
    now()
FROM chembl_tmp;
""")

con.execute("""
INSERT INTO compounds
SELECT
    Standard_InChIKey,
    Source_Record_ID,
    SMILES,
    DB_Source
FROM (
    SELECT
        Standard_InChIKey,
        Source_Record_ID,
        SMILES,
        DB_Source,
        ROW_NUMBER() OVER (
            PARTITION BY Standard_InChIKey
            ORDER BY DB_Source, Source_Record_ID
        ) AS rn
    FROM chembl_tmp
    WHERE Standard_InChIKey NOT IN (
        SELECT Standard_InChIKey FROM compounds
    )
)
WHERE rn = 1;
""")

con.execute("DROP TABLE chembl_tmp;")
con.execute("DROP TABLE chembl_part1;")
con.execute("DROP TABLE chembl_part2;")

print("DONE!")

print("Building PubChem DB ...")

# PubChem

con.execute("""
CREATE OR REPLACE TABLE pubchem_smiles AS
SELECT
    "column0"         AS Source_Record_ID,
    "column1"         AS SMILES
FROM read_csv_auto('CID-SMILES')
""")


con.execute("""
CREATE OR REPLACE TABLE pubchem_inchikey AS
SELECT
    "column0"       AS Source_Record_ID,
    "column1"       AS Standard_Inchi,
    "column2"       AS Standard_InChIKey
FROM read_csv_auto('CID-InChI-Key')
""")

con.execute("""
CREATE TEMP TABLE pubchem_tmp AS
SELECT
    i.Standard_InChIKey,
    i.Standard_Inchi,
    s.SMILES,
    s.Source_Record_ID,
    'PubChem' AS DB_Source
FROM pubchem_smiles s
JOIN pubchem_inchikey i
USING (Source_Record_ID)
""")

con.execute("""
INSERT INTO db_stats
SELECT
    'PubChem',
    COUNT(*),
    COUNT(DISTINCT Standard_InChIKey),
    COUNT(*) - COUNT(DISTINCT Standard_InChIKey),
    (
        SELECT COUNT(DISTINCT b.Standard_InChIKey)
        FROM pubchem_tmp b
        JOIN compounds c
        ON b.Standard_InChIKey = c.Standard_InChIKey
    ),
    now()
FROM pubchem_tmp;
""")

con.execute("""
INSERT INTO compounds
SELECT
    Standard_InChIKey,
    Source_Record_ID,
    SMILES,
    DB_Source
FROM (
    SELECT
        Standard_InChIKey,
        Source_Record_ID,
        SMILES,
        DB_Source,
        ROW_NUMBER() OVER (
            PARTITION BY Standard_InChIKey
            ORDER BY DB_Source, Source_Record_ID
        ) AS rn
    FROM pubchem_tmp
    WHERE Standard_InChIKey NOT IN (
        SELECT Standard_InChIKey FROM compounds
    )
)
WHERE rn = 1;
""")

con.execute("DROP TABLE pubchem_smiles;")
con.execute("DROP TABLE pubchem_inchikey;")
con.execute("DROP TABLE pubchem_tmp;")

print("DONE!")

print("Building Zinc20 DB ...")

import pandas as pd
from bs4 import BeautifulSoup
import requests

url = 'https://cache.docking.org/2D/'
page = requests.get(url)
soup = BeautifulSoup(page.content, 'html.parser')
table = soup.find('table')
headers = [th.get_text(strip=True) for th in table.find_all('th')]

dir_names = []
for tr in table.find_all('tr'):
    cells = tr.find_all('td')
    if cells:
        row = [cell.get_text(strip=True) for cell in cells]
        base_name = row[1]
        if "/" in base_name:
            dir_names.append(url+base_name)

def text_files_scrap(url):
    page = requests.get(url)
    soup = BeautifulSoup(page.content, 'html.parser')
    table = soup.find('table')
    headers = [th.get_text(strip=True) for th in table.find_all('th')]

    rows = []
    for tr in table.find_all('tr'):
        cells = tr.find_all('td')
        if cells:
            row = [cell.get_text(strip=True) for cell in cells]
            base_name = row[1]
            if "txt" in base_name:
                rows.append(url+base_name)
    return rows

tranches_list = []
for name in dir_names:
    urls = text_files_scrap(name)
    tranches_list.append(urls)

flat_tranches_list = [item for sublist in tranches_list for item in sublist]

con.execute("""
CREATE OR REPLACE TABLE zinc20_tmp (
    Standard_InChIKey TEXT,
    Standard_Inchi TEXT,
    SMILES TEXT,
    Source_Record_ID TEXT,
    DB_Source TEXT
)
""")

import os
from zincdl import download_tranches
for url in flat_tranches_list:
    file = download_tranches([url], out_dir="zinc_downloads")
    filename = file[0][0]
    con.execute(f"""
    INSERT INTO zinc20_tmp
    SELECT
        "inchikey" AS Standard_InChIKey,
        '' AS Standard_Inchi,
        "smiles" AS SMILES,
        "zinc_id" AS Source_Record_ID,
        'Zinc20' AS DB_Source
    FROM read_csv_auto(
        '{filename}',
        delim='\t',
        quote='"',
        escape='"',
        header=true,
        strict_mode=false,
        ignore_errors=true,
        max_line_size=10000000
    )
    WHERE "inchikey" IS NOT NULL
    """)
    os.remove(filename)
    print(f"Removed {filename}")

con.execute("""
INSERT INTO db_stats
SELECT
    'Zinc20',
    COUNT(*),
    COUNT(DISTINCT Standard_InChIKey),
    COUNT(*) - COUNT(DISTINCT Standard_InChIKey),
    (
        SELECT COUNT(DISTINCT b.Standard_InChIKey)
        FROM zinc20_tmp b
        JOIN compounds c
        ON b.Standard_InChIKey = c.Standard_InChIKey
    ),
    now()
FROM zinc20_tmp;
""")

con.execute("""
INSERT INTO compounds
SELECT
    Standard_InChIKey,
    Source_Record_ID,
    SMILES,
    DB_Source
FROM (
    SELECT
        Standard_InChIKey,
        Source_Record_ID,
        SMILES,
        DB_Source,
        ROW_NUMBER() OVER (
            PARTITION BY Standard_InChIKey
            ORDER BY DB_Source, Source_Record_ID
        ) AS rn
    FROM zinc20_tmp
    WHERE Standard_InChIKey NOT IN (
        SELECT Standard_InChIKey FROM compounds
    )
)
WHERE rn = 1;
""")

con.execute("DROP TABLE zinc20_tmp;")

print("DONE!")
