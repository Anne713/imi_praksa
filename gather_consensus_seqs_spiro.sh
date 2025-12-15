#!/bin/bash

table="../../filogenija_consensus/best_hits_run03.tsv"
output="../../filogenija_consensus/top_spiro_consensus_seqs_run03.fasta"
> "$output"

tail -n +2 "$table" | while IFS=$'\t' read -r sample accession_full _; do

    fasta_file="${sample}.fasta"
    accession="${accession_full%%_*}"

    if [[ -f "$fasta_file" ]]; then

        acc_escaped=$(printf '%s\n' "$accession" | sed 's/[][\.*^$(){}?+|]/\\&/g')

        awk -v acc="$acc_escaped" -v smp="$sample" '
            /^>/ {
                keep=0
                if($0 ~ acc) {
                    keep=1
                    print ">" smp "_" substr($0,2)
                    next
                }
            }
            {if(keep) print}
        ' "$fasta_file" >> "$output"
    else
        echo "Warning: FASTA file '$fasta_file' not found"
    fi

done
