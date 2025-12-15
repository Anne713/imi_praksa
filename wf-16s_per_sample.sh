#!/bin/bash

echo "--- Running Nextflow wf-16s per sample ..."

results_dir="results"
reads_dir="reads_ixodes"

total=$(find "$reads_dir" -maxdepth 1 -name "*.fastq" | wc -l)
count=0

echo "--- Classifying $total samples ..."

for sample in "$reads_dir"/*.fastq; do
    filename=$(basename "$sample")
    prefix="${filename%.fastq}"
    
    ((count++))
    
    if [ -d "$results_dir/${prefix}_results" ]; then
        echo "--- Skipping $sample ($count/$total) - results already exist!"
        continue
    fi
    
    echo "--- Classifying $sample ($count/$total)..."

    nextflow run epi2me-labs/wf-16s \
        --fastq "$sample" \
        --classifier kraken2 \
        --database_set ncbi_16s_18s \
        --bracken_threshold 5 \
        --kraken2_memory_mapping True \
        --out_dir "$results_dir/${prefix}_results" \
        --include_read_assignments True \
        --output_unclassified True
    
    echo "--- Done with $sample ($count/$total)!"
done

if [ -d "work" ]; then
    echo "--- Cleaning up Nextflow work directory to save space ..."
    rm -rf work
fi

echo "--- Finished Nextflow wf-16s! ($count/$total) Results are in $results_dir"

