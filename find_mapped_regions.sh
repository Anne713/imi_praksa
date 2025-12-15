#!/bin/bash

# get reference and make index
echo "--- Getting reference"

wget -nc https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz

if [[ ! -f "GCF_000001405.40_GRCh38.p14_genomic.fna" ]]; then
    echo "--- Unzipping reference ..."
    gunzip -kf GCF_000001405.40_GRCh38.p14_genomic.fna.gz
fi

if [[ ! -f "GCF_000001405.40_GRCh38.p14_genomic.fna.fai" ]]; then
    echo "--- Creating index ..."
    samtools faidx GCF_000001405.40_GRCh38.p14_genomic.fna
fi

ref_genome_idx="GCF_000001405.40_GRCh38.p14_genomic.fna.fai"

# make windows
echo "--- Making windows"

bedtools makewindows -g ${ref_genome_idx} -w 10000 > genome_10kb.bed

# find bams
results_dir="../wf-16s/results_host_g"
samples="to_include.txt"

host_bams=""
while read -r sample; do

    echo "--- $sample ..."

    bam_dir="$results_dir/$sample/host_bam"
    if [[ -d "$bam_dir" ]]; then
        bams=("$bam_dir"/*.bam)
        if [[ -e "${bams[0]}" ]]; then
            for bam in "${bams[@]}"; do
                host_bams="$host_bams $bam"
            done
        else
            echo "Warning: no BAM files in $bam_dir"
        fi
    else
        echo "Warning: no host_bam directory for $sample"
    fi
done < "$samples"

# count overlapping alignments
echo "--- Counting reads ..."

bedtools multicov -bams ${host_bams} -bed genome_10kb.bed > counts_per_window.bed
# help: https://bedtools.readthedocs.io/en/latest/content/tools/multicov.html?highlight=multicov 

echo "--- Summing up ..."

awk 'BEGIN{print "chrom\tstart\tend\tsum_counts"}{
    sum=0; for(i=4;i<=NF;i++) sum+=$i; print $1, $2, $3, sum
}' counts_per_window.bed > counts_per_window.sum.bed

(head -n 1 counts_per_window.sum.bed && tail -n +2 counts_per_window.sum.bed | sort -k4,4nr | head -n 50) > top50_windows.bed

echo "--- Done! Results:"
echo "      * counts_per_window.bed      (per-sample counts per 10kb window)"
echo "      * counts_per_window.sum.bed  (summed across all samples)"
echo "      * top50_windows.bed          (top 50 windows with highest read counts)"

