# LiMCA: **L**inking **m**RNA to **C**hromatin **A**rchitecture

<img src="images/limca.png" width="275"/><br>

LiMCA is a new single-cell sequencing method that can simultaneously profile the 3D genome structure and transcriptome in individual cells. On transcriptome side, LiMCA provides full-length mRNA read coverage for thousands of genes. This capacity is crucial for research across various domains, such as olfactory receptor expression. On 3D genome side, LiMCA can accurately identify distinct cell types and reconstruct high-resolution 3D models.

This repository contains data, scripts, pipelines used for LiMCA analysis.

## Workflow

### scRNA-seq analysis 

#### &#x2714; Basic analysis

We followed the standard Smart-seq2 processing workflow documented in the Human Cell Atlas (HCA) Data Portal (https://data.humancellatlas.org/pipelines/smart-seq2-workflow).

#### &#x2714; Allele-specific gene expression from LiMCA 

We can calculate single-cell haplotypic gene expression using [phASER](https://github.com/secastel/phaser). A simple bash script, **RNA_phaser_ase.sh**, is provided.

### scHi-C analysis

&#x2714; For analysis of scHi-C data, we encapsulated **hickit** and **dip-c** commands to creat the bash pipelines, **LiMCA_for_diploid.sh** and **LiMCA_for_olfactory.sh**. See the documentations of [hickit](https://github.com/lh3/hickit) and [dip-c](https://github.com/tanlongzhi/dip-c) for detailed, step-by-step instructions.

```
LiMCA_for_diploid.sh

Usage: LiMCA_for_diploid.sh

Options:
    -g: genome.fa
    -f: R1.fq
    -r: R2.fq
    -t: N, threads
    -v: SNP FILE
    -s: GENOME ASSEMBLY
    -h: help information

Action:
    -A: Reads mapping
    -P: Contacts parsing, haplotype imputation, and scA/B values
    -S: Reconstruction of 3D genome structures, structure filtering
    -D: Typical Dip-C analysis
    -R: Spatial analysis of actively expressed genes in individual cells
        (visualization of actively expressed genes, calculation of spatial relationship between actively expressed genes)
```


```
LiMCA_for_olfactory.sh

Usage: LiMCA_for_olfactory.sh

Options:
    -g: genome.fa
    -f: R1.fq
    -r: R2.fq
    -t: N, threads
    -s: GENOME ASSEMBLY
    -v: SNP FILE
    -h: help information

Action:
    -A: Reads mapping
    -P: Contacts parsing, haplotype imputation
    -S: Reconstruction of 3D genome structures, structure filtering
    -C: Contact analysis of Greek islands and OR genes for olfactory neurons
    -D: Similar to -C, but 3D structure analysis for olfactory neurons

See also:
    GEO README (ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE121nnn/GSE121791/suppl/GSE121791%5F00README%2Emd%2Etxt)
```

