## Bacterial genome assembly and decontamination

*Note*: This repository describes my project performed at Bioinformatics Institute during the fall term - 2020. **something about files located here**

### Background 
Maintaining balanced gut microbiota was shown to be crucial for human health. Probiotics such as *Escgerichia coli* Nissle 1917 could recover beneficial functions in gut microbial communities and prevent it from being populated with pathogenic bacteria. **LINK!!** 

*Escherichia coli* str. Nissle 1917 is one of common probiotics used for maintaining balanced gut microbiota. **LINK!!** Sample of this probiotic obtained by our laboratory at drugstore occured to be contaminated. Sample inoculation led to producing two morphologically distinct types of colonies: of small and big size. We assumed that small colonies (referred further as "Nissle Small", or simply "NS") were formed only by *E. coli* str. Nissle 1917, while big ones ("Nissle Big", or "NB") consisted of both *E. coli* str. Nissle 1917 and some other bacterial contaminants.  Two types of colonies were subjected to NextSeq sequencing in paired-end mode with read length of 75 nt, generating two sequencing samples: NB and NS. 

**about NCBI - there are 2 complete assembly but they differ**
In this work, we aimed to perform decontamination and *de novo* assembly of two sequencing samples obtained from the pharmaceutical *E. coli* str. Nissle 1917. 

To ,,, we objectives:
.... 
### Required programs and - check prokka and others - про комп


### Workflow
#### 1. Quality assessment of read ?? (NS & NB)
fastqc
```bash
ls
```
#### 2. *De novo* assembly (NS & NB)
spades
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

