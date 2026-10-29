#!/bin/bash
#$ -N download_dolores_RNA
#$ -q all.q
#$ -cwd
#$ -V
#$ -pe smp 1
#$ -j y
#$ -o logs_drive/drive_dolores.$JOB_ID.log

# Download the FASTQ files Dolores shared in 3 Google Drive folders into RAWDATA/.
# - one file at a time, with bandwidth and API limits
# - never overwrites: a file that already exists in RAWDATA/ is skipped and reported
# - checks each download against the MD5 that Google Drive reports
# - records where each file came from in drive_provenance.tsv
#
# Test first (lists what would happen, downloads nothing):
#   bash download_drive_dolores.sh --dry-run
# Then submit from RAWDATA/:
#   qsub download_drive_dolores.sh
# Safe to resubmit: files already downloaded are skipped.

set -uo pipefail

module load miniconda3
conda activate rclone-1.74.4

DRY_RUN=0
[[ "${1:-}" == "--dry-run" ]] && DRY_RUN=1

DEST=/Storage/data1/jorge.munoz/DOLORES/RAWDATA
LOGDIR=${DEST}/logs_drive
PROVENANCE=${DEST}/drive_provenance.tsv
MD5FILE=${DEST}/md5sum_drive_dolores.txt
BWLIMIT=20M          # max download speed (bytes/s)
mkdir -p "$LOGDIR"

# remote name -> Drive folder ID (folders shared by Dolores, 2026-09)
REMOTES=(transcriptome_1 transcriptome_2 transcriptome_3)
declare -A FOLDER=(
    [transcriptome_1]=1ikTM3-lMgVfOTtKEJpG8W644PI7wdc3H
    [transcriptome_2]=12_1GW-8wniw3lmRYo5QrOOPnt-9EkD2Z
    [transcriptome_3]=1nrCjvFUpU0yc7PvC0E-65qLpG7m0Rlfo
)

# team_drive= clears the Shared Drive name saved in the remote config (these are normal shared folders)
RCLONE_OPTS=(--transfers 1 --checkers 1 --bwlimit "$BWLIMIT" --tpslimit 5
             --retries 5 --low-level-retries 10 --ignore-existing --stats 5m)

(( DRY_RUN )) || [[ -f "$PROVENANCE" ]] || echo -e "file\tremote\tfolder_id\tsize_bytes\tdrive_md5\tstatus\tdate" > "$PROVENANCE"

n_ok=0; n_skip=0; n_fail=0
for remote in "${REMOTES[@]}"; do
    src="${remote},team_drive=,root_folder_id=${FOLDER[$remote]}:"
    echo "=== ${remote} (${FOLDER[$remote]})"

    listing=$(rclone lsf --files-only --include "*.fastq.gz" --format "psh" --hash md5 --separator $'\t' "$src")
    if [[ -z "$listing" ]]; then
        echo "ERROR: could not list ${remote} or it has no .fastq.gz files"; n_fail=$((n_fail+1)); continue
    fi

    while IFS=$'\t' read -r file size md5; do
        # never overwrite: skip if a file with this name exists anywhere under RAWDATA/
        existing=$(find "$DEST" -name "$file" -print -quit)
        if [[ -n "$existing" ]]; then
            echo "SKIP (already exists: ${existing#$DEST/}) ${file}"
            n_skip=$((n_skip+1)); continue
        fi

        if (( DRY_RUN )); then
            echo "WOULD DOWNLOAD ${file} (${size} bytes)"; continue
        fi

        echo "DOWNLOAD ${file} (${size} bytes) $(date '+%F %T')"
        if rclone copyto "${src}${file}" "${DEST}/${file}" "${RCLONE_OPTS[@]}" \
                --log-file "${LOGDIR}/rclone_${remote}.log" --log-level INFO; then
            local_md5=$(md5sum "${DEST}/${file}" | cut -d' ' -f1)
            if [[ "$local_md5" == "$md5" ]]; then
                status=OK; n_ok=$((n_ok+1))
                echo "${local_md5}  ${file}" >> "$MD5FILE"
            else
                status="MD5_MISMATCH(local=${local_md5})"; n_fail=$((n_fail+1))
                mv "${DEST}/${file}" "${DEST}/${file}.md5_mismatch"
            fi
        else
            status=RCLONE_FAILED; n_fail=$((n_fail+1))
        fi
        echo "  -> ${status}"
        echo -e "${file}\t${remote}\t${FOLDER[$remote]}\t${size}\t${md5}\t${status}\t$(date '+%F %T')" >> "$PROVENANCE"
    done <<< "$listing"
done

echo "=== done: ${n_ok} downloaded OK, ${n_skip} skipped (already existed), ${n_fail} failed"
(( n_fail == 0 ))
