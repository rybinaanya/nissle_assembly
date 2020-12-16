## Bacterial genome assembly and decontamination

*Note*: This repository describes my project performed at Bioinformatics Institute during the fall term - 2020. **something about files located here**
**Results/** - contains results of FastQC analysis, assembly statistics (QUAST output),  completeness and contamination assessment (CheckM report), 
### Background 
Maintaining balanced gut microbiota was shown to be crucial for human health. Probiotics such as *Escgerichia coli* Nissle 1917 could recover beneficial functions in gut microbial communities and prevent it from being populated with pathogenic bacteria. **LINK!!** 

*Escherichia coli* str. Nissle 1917 is one of common probiotics used for maintaining balanced gut microbiota. **LINK!!** Sample of this probiotic obtained by our laboratory at drugstore occured to be contaminated. Sample inoculation led to producing two morphologically distinct types of colonies: of small and big size. We assumed that small colonies (referred further as "Nissle Small", or simply "NS") were formed only by *E. coli* str. Nissle 1917, while big ones ("Nissle Big", or "NB") consisted of both *E. coli* str. Nissle 1917 and some other bacterial contaminants. Both colonies were subjected to NextSeq sequencing in paired-end mode with the read length of 75 nt, generating two sequencing samples: NB and NS. 

**about NCBI - there are 2 complete assembly but they differ**
In this work, we aimed to perform decontamination and *de novo* assembly of two sequencing samples obtained from the pharmaceutical *E. coli* str. Nissle 1917. 

To ,,, we objectives:
.... 
### Required programs used in the study and - check prokka and others - про комп
In this study, the following programs were used:
* FastQC v0.11.9
* SPAdes v3.13.1 
* QUAST v5.1.0rc1
* CONCOCT v1.1.0
* Bowtie2 v2.2.1

### Workflow
#### 1. Quality assessment of raw sequencing data (NS & NB)

Quality of reads was analyzed using FastQC v0.11.9 with default parameters. Variable ${working_dir} specifies the directory where all the output of data analysis would be saved (e.g. /home/rybina/Nissle_project). Variable ${path_to_reads} is the directory where raw reads in fastq format are located (e.g. /home/rybina/Nissle_project/Reads_data).  
```{bash}
fastqc -o ${working_dir}/fastqc_output ${path_to_reads}/*.fastq.gz
```
FastQC report can be found  here **Results/FastQC/**. There was no need to subject reads to trimming.

#### 2. *De novo* assembly (NS & NB)
Bacterial genome assembly was performed via SPAdes v3.13.1 (with option --careful):

```{bash}
# NB sample
spades.py --careful -o ${working_dir}/NB_spades -1 ${path_to_reads}/ARyb_NB3_S54_R1_001.fastq.gz -2 ${path_to_reads}/ARyb_NB3_S54_R2_001.fastq.gz
```
```{bash}
# NS sample
spades.py --careful -o ${working_dir}/NS_spades -1 ${path_to_reads}/ARyb_NS2_S53_R1_001.fastq.gz -2 ${path_to_reads}/ARyb_NS2_S53_R2_001.fastq.gz 
```

#### 3. Assembly statistics (NS & NB)
Statistics on resulting contigs and scaffolds for both samples was obtained using QUAST v5.1.0rc1. Each assembly was compared to complete genome of the *E. coli* Nissle 1917 available under GenBank assembly accession either GCA_003546975.1 (dated on 2018) or GCA_000714595.1 (dated on 2014).  Both QUAST output are located at **Results/QUAST/**

```{bash}
# Download gff files to the working directory ${working_dir}:
wget -P ${working_dir}/ https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/546/975/GCA_003546975.1_ASM354697v1/GCA_003546975.1_ASM354697v1_genomic.gff.gz 
wget -P ${working_dir}/ https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/714/595/GCA_000714595.1_ASM71459v1/GCA_000714595.1_ASM71459v1_genomic.gff.gz

# Download genome fasta files to the working directory ${working_dir}:
wget -P ${working_dir}/ https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/546/975/GCA_003546975.1_ASM354697v1/GCA_003546975.1_ASM354697v1_genomic.fna.gz 
wget -P ${working_dir}/ https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/714/595/GCA_000714595.1_ASM71459v1/GCA_000714595.1_ASM71459v1_genomic.fna.gz

# Decompressing:
gunzip ${working_dir}/GCA_000714595.1_ASM71459v1_genomic.gff.gz
gunzip ${working_dir}/GCA_003546975.1_ASM354697v1_genomic.gff.gz 
gunzip ${working_dir}/GCA_003546975.1_ASM354697v1_genomic.fna.gz
gunzip ${working_dir}/GCA_000714595.1_ASM71459v1_genomic.fna.gz

# Running QUAST: NS & NB assemblies vs GCA_003546975.1 (dated on 2018):
quast.py ${working_dir}/NS_spades/contigs.fasta ${working_dir}/NS_spades/scaffolds.fasta ${working_dir}/NB_spades/contigs.fasta ${working_dir}/NB_spades/scaffolds.fasta -r ${working_dir}/GCA_003546975.1_ASM354697v1_genomic.fna -g ${working_dir}/GCA_003546975.1_ASM354697v1_genomic.gff -o ${working_dir}/quast_output_NS_NB_Nissle2018

# Running QUAST: NS & NB assemblies vs  GCA_000714595.1 (dated on 2014):
quast.py ${working_dir}/NS_spades/contigs.fasta ${working_dir}/NS_spades/scaffolds.fasta ${working_dir}/NB_spades/contigs.fasta ${working_dir}/NB_spades/scaffolds.fasta -r ${working_dir}/GCA_000714595.1_ASM71459v1_genomic.fna -g ${working_dir}/GCA_000714595.1_ASM71459v1_genomic.gff -o ${working_dir}/quast_output_NS_NB_Nissle2014
```

#### 4. Contamination and completeness assessment (NS & NB)
##### 4.1. Binning
Binning was done using CONCOCT v1.1.0. For running CONCOCT, alignment of reads to the contigs should be provided. First, NB sample reads were mapped to the SPAdes-derived contigs via bowtie2 v2.2.1. Resulting alignment was modified using samtools 1.11. 
```{bash} 
# Indexing contigs
bowtie2-build ${working_dir}/NB_spades/contigs.fasta ${working_dir}/NBspades_contigs

# Mapping
bowtie2 -x ${working_dir}/NBspades_contigs -1 ${path_to_reads}/ARyb_NB3_S54_R1_001.fastq.gz -2 ${path_to_reads}/ARyb_NB3_S54_R2_001.fastq.gz -S ${working_dir}/NB3_reads_contigs.sam

# Sam to bam conversion
samtools view -S ${working_dir}/NB3_reads_contigs.sam -b -o ${working_dir}/NB3_reads_contigs.bam

# Samtools sort
samtools sort ${working_dir}/NB3_reads_contigs.bam -o ${working_dir}/NB3_reads_contigs.sorted.bam

# Samtools index
samtools index ${working_dir}/NB3_reads_contigs.sorted.bam
```

Binning:
```{bash}
# Rename a copy of indexed bam file for bining
cp ${working_dir}/NB3_reads_contigs.sorted.bai ${working_dir}/NB3_reads_contigs.sorted.index.bam

# Create a folder for concoct output
mkdir ${working_dir}/concoct_output

# Cut up contigs fasta file in non-overlapping (-o 0) parts of length 10K (chunk size is 10K: -c 10000), the last contig part would be between 10K and 20K long (--merge_last); contig parts would be specified in the BED file (-b contigs_10K.bed):
cut_up_fasta.py ${working_dir}/NB_spades/contigs.fasta -c 10000 -o 0 --merge_last -b ${working_dir}/concoct_output/contigs_10K.bed > ${working_dir}/concoct_output/contigs_10K.fa

# Create coverage table based on  data on contigs parts in BED format and sorted indexed alignment file of reads to contigs:
concoct_coverage_table.py ${working_dir}/concoct_output/contigs_10K.bed ${working_dir}/NB3_reads_contigs.sorted.index.bam > ${working_dir}/concoct_output/coverage_table.tsv 

# Perform unsupervised binning of contigs
concoct --composition_file ${working_dir}/concoct_output/contigs_10K.fa --coverage_file ${working_dir}/concoct_output/coverage_table.tsv -b ${working_dir}/concoct_output/

# Merge subcontig clustering (clustering_gt1000.csv obtained by running concoct command) into original contig clustering (clustering_merged.csv):
merge_cutup_clustering.py ${working_dir}/concoct_output/clustering_gt1000.csv > ${working_dir}/concoct_output/clustering_merged.csv

# Create a folder for bins:
mkdir ${working_dir}/concoct_output/NB_bins

# Extract a fasta file for each cluster specified by concoct
extract_fasta_bins.py ${working_dir}/NB_spades/contigs.fasta ${working_dir}/concoct_output/clustering_merged.csv --output_path ${working_dir}/concoct_output/NB_bins
```
##### 4.2. Bins evaluation using single copy genes
NB sample bins were assessed for completeness and contamination using CheckM v1.1.3. 

```{bash}
# Create a directory for CheckM output:
mkdir ${working_dir}/CheckM_output_NB

# Create a directory for CheckM database:
mkdir ${working_dir}/CheckM_database

# Download CheckM database files
wget -P ${working_dir}/CheckM_database https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz

# Decompressing
tar -xzvf ${working_dir}/CheckM_database/checkm_data_2015_01_16.tar.gz

# Set the directory where CheckM database would be located
checkm data setRoot ${working_dir}/CheckM_database

# Run workflow using lineage-specific marker sets (see https://github.com/Ecogenomics/CheckM/wiki/Workflows#lineage-specific-workflow)
checkm lineage_wf -x fa ${working_dir}/concoct_output/NB_bins ${working_dir}/CheckM_output_NB
```

To evaluate completeness and contamination of NS sample, the following commands were used (memory consuming process):
```{bash}
# Create a directory for CheckM output:
mkdir ${working_dir}/CheckM_output_NS

# Prepare a bin folder:
mkdir ${working_dir}/CheckM_output_NS/NS_bins
cp ${working_dir}/NB_spades/contigs.fasta ${working_dir}/CheckM_output_NS/NS_bins

#  Run workflow using lineage-specific marker sets
checkm lineage_wf ${working_dir}/CheckM_output_NS/NS_bins ${working_dir}/CheckM_output_NS
```
Summary on bins evaluation might be found at **Results/CheckM/** for both samples. 
#### 5. Annotation (NS)
##### 5.1. Prokka
##### 5.2. PGAP
#### 6. Taxonomy identification (NB) 
##### 6.1. 16S rRNA - как назвать ?
##### 6.2. Что-то про то как достали данные
##### 6.3. Kraken preparation
##### 6.4. Kraken running
##### 6.5. Kraken report visualization

