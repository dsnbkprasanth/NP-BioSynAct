#!/usr/bin/env bash
# NP-BioSynAct -- download ChEMBL SQLite + extract human bioactivity (the y labels).
# ~5 GB download, ~25 GB unpacked. Runs on your PC's network.
set -euo pipefail
cd "$(dirname "$0")"

BASE="https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest"

# --- 1. discover the real filename (version changes) ------------------------
FILE=$(curl -s "$BASE/checksums.txt" | grep -o 'chembl_[0-9]*_sqlite.tar.gz' | head -1)
[ -n "$FILE" ] || { echo "could not find sqlite tarball name at $BASE"; exit 1; }
echo "file: $FILE"

# --- 2. download (resumable) + verify it's really gzip ----------------------
curl -C - -O "$BASE/$FILE"
if ! tar tzf "$FILE" >/dev/null 2>&1; then
    echo "!! $FILE is not a valid gzip -- download was truncated or blocked."
    echo "   delete it and re-run:  rm -f $FILE && bash 10_get_chembl.sh"
    exit 1
fi
tar -xzf "$FILE"
DB=$(find . -name 'chembl_*.db' | head -1)
echo "db: $DB"

# --- 3. dump all high-confidence human single-protein activities ------------
sqlite3 -header -csv "$DB" > chembl_activities.csv <<'SQL'
SELECT md.chembl_id            AS compound_id,
       cs.canonical_smiles     AS smiles,
       td.chembl_id            AS target_id,
       td.pref_name            AS target_name,
       act.standard_type       AS act_type,
       act.pchembl_value       AS pchembl
FROM activities act
JOIN assays a               ON act.assay_id = a.assay_id
JOIN target_dictionary td   ON a.tid = td.tid
JOIN molecule_dictionary md ON act.molregno = md.molregno
JOIN compound_structures cs ON act.molregno = cs.molregno
WHERE act.pchembl_value >= 6
  AND act.standard_relation = '='
  AND td.target_type = 'SINGLE PROTEIN'
  AND td.organism = 'Homo sapiens';
SQL

echo "rows: $(($(wc -l < chembl_activities.csv) - 1))  ->  chembl_activities.csv"
