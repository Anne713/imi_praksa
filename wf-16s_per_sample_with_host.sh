#!/bin/bash

echo "--- Running Nextflow wf-16s per sample ..."

results_dir="results_host_g"

ref_human="reference/GCF_000001405.40_GRCh38.p14_genomic.fna"
ref_ixodes="reference/GCA_964417485.1_Rsp_Iric_1_genomic.fna"
ref_dermacentor="reference/GCA_051549955.1_ASM5154995v1_genomic.fna"

if [[ -f "$ref_human" && -f "$ref_ixodes" && -f "$ref_dermacentor" ]]; then
    echo "--- References already present"
else
    echo "--- Downloading reference genomes ..."
    mkdir -p reference
    cd reference

    wget -nc https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
    wget -nc https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/964/417/485/GCA_964417485.1_Rsp_Iric_1/GCA_964417485.1_Rsp_Iric_1_genomic.fna.gz
    wget -nc https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/051/549/955/GCA_051549955.1_ASM5154995v1/GCA_051549955.1_ASM5154995v1_genomic.fna.gz

    gunzip -f *.gz
    cd ..
    echo "--- Downloaded reference genomes"
fi

total_dermacentor=$(ls reads_dermacentor/*.fastq | wc -l)
total_ixodes=$(ls reads_ixodes/*.fastq | wc -l)
total_human=$(ls reads_human/*.fastq | wc -l)
total_other=$(ls reads_other/*.fastq | wc -l)
total=$(( ${total_dermacentor:-0} + ${total_ixodes:-0} + ${total_human:-0} + ${total_other:-0} ))
tot_count=0

echo "--- Total: $total samples"

echo "--- Classifying $total_dermacentor Dermacentor reticulatus samples ..."
count=0

for sample in reads_dermacentor/*.fastq; do
    filename=$(basename "$sample")
    prefix="${filename%.fastq}"
    
    ((count++))
    ((tot_count++))
    
    if [ -d "$results_dir/${prefix}_results" ]; then
        echo "--- Skipping $sample ($count/$total_dermacentor, $tot_count/$total) - results already exist!"
        continue
    fi
    
    echo "--- Classifying $sample ($count/$total_dermacentor, $tot_count/$total)..."

    nextflow run epi2me-labs/wf-16s \
        --fastq "$sample" \
        --classifier kraken2 \
        --exclude_host "$ref_dermacentor" \
        --database_set ncbi_16s_18s \
        --bracken_threshold 5 \
        --kraken2_memory_mapping True \
        --out_dir "$results_dir/${prefix}_results" \
        --include_read_assignments True \
        --output_unclassified True
    
    echo "--- Done with $sample ($count/$total_dermacentor, $tot_count/$total)!"
done

if [ -d "work" ]; then
    echo "--- Cleaning up Nextflow work directory to save space ..."
    rm -rf work
fi

echo "--- Classifying $total_ixodes Ixodes ricinus samples ..."
count=0

for sample in reads_ixodes/*.fastq; do
    filename=$(basename "$sample")
    prefix="${filename%.fastq}"
    
    ((count++))
    ((tot_count++))
    
    if [ -d "$results_dir/${prefix}_results" ]; then
        echo "--- Skipping $sample ($count/$total_ixodes, $tot_count/$total) - results already exist!"
        continue
    fi
    
    echo "--- Classifying $sample ($count/$total_ixodes, $tot_count/$total)..."

    nextflow run epi2me-labs/wf-16s \
        --fastq "$sample" \
        --classifier kraken2 \
        --exclude_host "$ref_ixodes" \
        --database_set ncbi_16s_18s \
        --bracken_threshold 5 \
        --kraken2_memory_mapping True \
        --out_dir "$results_dir/${prefix}_results" \
        --include_read_assignments True \
        --output_unclassified True
    
    echo "--- Done with $sample ($count/$total_ixodes, $tot_count/$total)!"
done

if [ -d "work" ]; then
    echo "--- Cleaning up Nextflow work directory to save space ..."
    rm -rf work
fi

echo "--- Classifying $total_other other samples ..."
count=0

for sample in reads_other/*.fastq; do
    filename=$(basename "$sample")
    prefix="${filename%.fastq}"
    
    ((count++))
    ((tot_count++))
    
    if [ -d "$results_dir/${prefix}_results" ]; then
        echo "--- Skipping $sample ($count/$total_other, $tot_count/$total) - results already exist!"
        continue
    fi
    
    echo "--- Classifying $sample ($count/$total_other, $tot_count/$total)..."

    nextflow run epi2me-labs/wf-16s \
        --fastq "$sample" \
        --classifier kraken2 \
        --database_set ncbi_16s_18s \
        --bracken_threshold 5 \
        --kraken2_memory_mapping True \
        --out_dir "$results_dir/${prefix}_results" \
        --include_read_assignments True \
        --output_unclassified True
    
    echo "--- Done with $sample ($count/$total_other, $tot_count/$total)!"
done

if [ -d "work" ]; then
    echo "--- Cleaning up Nextflow work directory to save space ..."
    rm -rf work
fi

echo "--- Classifying $total_human human samples ..."
count=0

for sample in reads_human/*.fastq; do
    filename=$(basename "$sample")
    prefix="${filename%.fastq}"
    
    ((count++))
    ((tot_count++))
    
    if [ -d "$results_dir/${prefix}_results" ]; then
        echo "--- Skipping $sample ($count/$total_human, $tot_count/$total) - results already exist!"
        continue
    fi
    
    echo "--- Classifying $sample ($count/$total_human, $tot_count/$total)..."

    nextflow run epi2me-labs/wf-16s \
        --fastq "$sample" \
        --classifier kraken2 \
        --exclude_host "$ref_human" \
        --database_set ncbi_16s_18s \
        --bracken_threshold 5 \
        --kraken2_memory_mapping True \
        --out_dir "$results_dir/${prefix}_results" \
        --include_read_assignments True \
        --output_unclassified True
    
    echo "--- Done with $sample ($count/$total_human, $tot_count/$total)!"
done

if [ -d "work" ]; then
    echo "--- Cleaning up Nextflow work directory to save space ..."
    rm -rf work
fi

echo "--- Finished Nextflow wf-16s! ($tot_count/$total) Results are in $results_dir"


