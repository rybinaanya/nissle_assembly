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
### Required programs and - check prokka and others - про комп
FastQC v0.11.9
SPAdes v3.13.1


### Workflow
#### 1. Quality assessment of raw sequencing data (NS & NB)

Quality of reads was analyzed using FastQC v0.11.9 with default parameters. Variable ${working_dir} specifies the directory where all the output of data analysis would be saved (e.g. /home/rybina/Nissle_project). Variable ${path_to_reads} is the directory where raw reads in fastq format are located (e.g. /home/rybina/Nissle_project/Reads_data).  
```{bash}
fastqc -o ${working_dir}/fastqc_output ${path_to_reads}/*.fastq.gz
```
FastQC report can be found  here **Results/FastQC/**. There was no need to subject reads to trimming.

#### 2. *De novo* assembly (NS & NB)
Bacterial genome assembly was performed via SPAdes v3.13.1 

* NB
```
spades.py --careful -o path_to_output -1 path_to_reads/ARyb_NB3_S54_R1_001.fastq.gz -2 path_to_reads/ARyb_NB3_S54_R2_001.fastq.gz
```
* NS
```{bash }spades.py --careful -o path_to_output -1 path_to_reads/ARyb_NS2_S53_R1_001.fastq.gz -2 path_to_reads/ARyb_NS2_S53_R2_001.fastq.gz ```
#### 3. Assembly statistics (NS & NB)
quast

#### 4. Contamination and completeness assessment (NS & NB)
##### 4.1. Binning
##### 4.2. CheckM - как назвать ?
#### 5. Annotation (NS)
##### 5.1. Prokka
##### 5.2. PGAP
#### 6. Taxonomy identification (NB) 
##### 6.1. 16S rRNA - как назвать ?
##### 6.2. Что-то про то как достали данные
##### 6.3. Kraken preparation
##### 6.4. Kraken running
##### 6.5. Kraken report visualization

