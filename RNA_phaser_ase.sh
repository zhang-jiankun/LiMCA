#!/bin/bash

# See https://github.com/secastel/phaser
#

cell=$1
array=$2
results="results"
VCF=$3
sample=$4

# export PATH=/share/home/zhangjk/software/htslib:$PATH
# export PYTHONPATH=$PYTHONPATH:/share/home/zhangjk/.local/lib/python2.7/site-packages

echo "[ $(date) ] Submit $cell: $array"
phaser="/share/home/zhangjk/hic_tools/phaser"
scripts="/share/home/zhangjk/scRNAC/MOE/scripts"
GENES="/share/home/zhangjk/database/mouse/gencode.vM25.genes.bed"

sbatch --array=$array <<-EOF
#!/bin/bash
#SBATCH -J $cell
#SBATCH -c 8
#SBATCH --mem 16G
#SBATCH -o ./logs/15_map_transcripts/${cell}_%a-phaser-%A.out

id=\$(printf %3s \${SLURM_ARRAY_TASK_ID} | tr ' ' '0')
mkdir -p $results/${cell}_\${id}/phaser

# Run phASER 
python ${phaser}/phaser/phaser.py --vcf ${VCF} \
    --bam ./results/${cell}_\${id}/${cell}_\${id}_cDNA_report/${cell}_\${id}.gen.sorted.bam \
    --temp_dir  /dshare/xielab/analysisdata/ThreeDGene/zhangjk/scRNAC/Olfactory/temp \
    --paired_end 1 \
    --mapq 1 \
    --baseq 10 \
    --sample $sample \
    --threads 8 \
    --o $results/${cell}_\${id}/phaser/${cell}_\${id}

# Generate haplotype expression quantifications
python ${phaser}/phaser_gene_ae/phaser_gene_ae.py \
    --haplotypic_counts $results/${cell}_\${id}/phaser/${cell}_\${id}.haplotypic_counts.txt \
    --features ${GENES} \
    --o $results/${cell}_\${id}/phaser/${cell}_\${id}_gene_ae.txt

# get a sense for the data
Rscript ${scripts}/phaser_qc.R $results/${cell}_\${id}/phaser/${cell}_\${id}_gene_ae.txt $results/${cell}_\${id}/phaser/${cell}_\${id}_gene_ae.qc.pdf
EOF
