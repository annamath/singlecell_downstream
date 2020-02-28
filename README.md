# scRNAseq_downstream2020

scRNA-seq and scATAC-seq datasets were produced with 10x Genomics (v3) technology and further preprocessed with Cell Ranger.
Since multiple donors were mixed in the same lane, donor demultiplexing based on genotype required. 

*scRNA-seq*

Note: Input fastq file names in the following format:
sample1_S1_L001_R1_001.fastq.gz, sample1_S1_L001_R2_001.fastq.gz

module load cellranger/3.0.1
cellranger count --id=<laneID> \
  --transcriptome=<GRCh38-3.0.0_reference> \
  --sample=<sampleID> \
  --fastqs=<directory containing R1,R2 fastqs>

Variant Calling for single cells using cellSNP
cellSNP v0.1.7
cellSNP -b <filtered_barcodes> \
  -s <possorted_bam> \
  -O <output_directory> \
  -R <region_vcf_file>\
  -p 22 \
  --minMAF 0.1 \
  --minCOUNT 100

vireo v0.1.8
vireo -c <cell_data> -N <n_donors> -o <out_dir>

*scATAC-seq*

Note: Input fastq file names in the following format:
sample1_S1_L001_I1_001.fastq.gz
sample1_S1_L001_R1_001.fastq.gz
sample1_S1_L001_R2_001.fastq.gz
sample1_S1_L001_R3_001.fastq.gz

module load cellranger-atac/1.2.0

cellranger-atac count --id=<laneID> \
  --reference=<GRCh38_reference> \
  --fastqs=<directory containing fastqs and index file> \
  --sample=<sampleID>

