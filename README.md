## Bacterial genome assembly and decontamination

 
### Background 

Our laboratory obtained the sample of *Escherichia coli* str. Nissle 1917 which occurred to be contaminated. Two types of colonies morphologically distinct by size were produced: small and big. We assumed that the first (referred further as "Nissle Small", or simply "NS") were formed only by *E. coli* str. Nissle 1917, while the second ("Nissle Big", or "NB") consisted of *E. coli* str. Nissle 1917 and some other bacterial contaminants. Both colonies were subjected to NextSeq sequencing in paired-end mode with the read length of 75 nt, generating two sequencing samples: NB and NS. 

In our work, we aimed to perform *de novo* assembly and determine taxonomy classification of contaminants for two sequencing samples obtained from the *E. coli* str. Nissle 1917 colonies. 

Therefore, we set the following objectives:
* perform *de novo* assembly of genome; 
* identify taxonomy of bacterial contaminants;
* annotate draft assembly and compare with a reference genome.

### Programs used in the study

In this study, the following programs were used:
* FastQC v0.11.9 (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* SPAdes v3.13.1
* QUAST v5.1.0rc1
* CONCOCT v1.1.0
* CheckM v1.1.3
* Bowtie2 v2.2.1
* Samtools v1.11
* BEDtools v2.27.0
* Prokka v1.12
* PGAP:2020-09-24.build4894
* Barrnap v0.9 (https://github.com/tseemann/barrnap)
* Kraken v1.1.1
* jellyfish  1.1.11 (https://github.com/gmarcais/Jellyfish, for running Kraken)
* Mauve v2.3.1
* Python 3.7
* Biopython

CheckM required:
* HMMER (>=3.1b1)
* prodigal (2.60 or >=2.6.1)
* pplacer (>=1.1)
* python >= 3.6
* numpy >= 1.13.1
* scipy >= 0.19.1
* matplotlib >= 2.1.0
* pysam >= 0.12.0.1
* dendropy >= 4.4.0

For running Prokka, the following dependencies were used:
* BioPerl (mandatory)
* GNU Parallel (mandatory)
* BLAST 2.9.0+ (mandatory)
* Prodigal V2.6.3 (mandatory)
* tbl2asn https://www.ncbi.nlm.nih.gov/genbank/tbl2asn2/ 
* ARAGORN v1.2.36 (recommended)
* barrnap 0.9 (recommended)
* minced 0.4.2 (recommended) https://github.com/ctSkennerton/minced
* signalp-4.1 (recommended)
* infernal-1.1.3 (recommended)

For running PGAP pipeline we also needed:
* ability to run Docker (see https://docs.docker.com/install/) commands without sudo  
* about 100GB of storage for database and working space,
* 2GB-4GB of memory available per CPU used by container.
* The CPU must have SSE 4.2 support (released in 2008).
* cwltool instaled (see https://github.com/common-workflow-language/cwltool)

For running Kraken, HPC cluster should be used. 
All programs were added to the system path. 

### Workflow
#### 1. Quality assessment of raw sequencing data (NS & NB samples)

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

#### 3. Assembly statistics (NS & NB samples)
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

#### 4. Contamination and completeness assessment (NB & NS samples)
##### 4.1. Binning (NB sample)
Binning was done using CONCOCT v1.1.0. For running CONCOCT, alignment of reads to the contigs should be provided. First, NB sample reads were mapped to the SPAdes-derived contigs via bowtie2 v2.2.1. Resulting alignment was modified using samtools 1.11. 
```{bash} 
# Indexing contigs
bowtie2-build ${working_dir}/NB_spades/contigs.fasta ${working_dir}/NBspades_contigs

# Mapping
bowtie2 -x ${working_dir}/NBspades_contigs -1 ${path_to_reads}/ARyb_NB3_S54_R1_001.fastq.gz -2 ${path_to_reads}/ARyb_NB3_S54_R2_001.fastq.gz -S ${working_dir}/NB3_reads_contigs.sam

# Sam to bam conversion
samtools view -S ${working_dir}/NB3_reads_contigs.sam -b -o ${working_dir}/NB3_reads_contigs.bam

# Sorting
samtools sort ${working_dir}/NB3_reads_contigs.bam -o ${working_dir}/NB3_reads_contigs.sorted.bam

# Indexing
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
##### 4.2. Bins evaluation using single copy genes (NB & NS samples)
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
cp ${working_dir}/NS_spades/contigs.fasta ${working_dir}/CheckM_output_NS/NS_bins

# Run workflow using lineage-specific marker sets
checkm lineage_wf -x fasta ${working_dir}/CheckM_output_NS/NS_bins ${working_dir}/CheckM_output_NS
```
Summary on bins evaluation might be found at **Results/CheckM_report/** for both samples. 

#### 5. Annotation (NS sample)
##### 5.1. Prokka

First, prepare output directory and additional input file:
```{bash}
# Create a directory for prokka output:
mkdir ${working_dir}/prokka_output_NS

# Download and decompress last modified curated annotation of E. coli Nissle 1917
wget -P ${working_dir}/ https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/546/975/GCF_003546975.1_ASM354697v1/GCF_003546975.1_ASM354697v1_genomic.gbff.gz
gunzip ${working_dir}/GCF_003546975.1_ASM354697v1_genomic.gbff.gz
```

For running prokka, several options were specified to force overwriting existing output folder (--force), specify genus and species, use filename output prefix "NS" (--prefix NS), add 'gene' features for each 'CDS' feature (--addgenes), use "Nissle" as locus tag prefix (--locustag Nissle), use genus-specific BLAST databases (--usegenus), search for ncRNAs with Infernal+Rfam (--rfam), specify Gram negative (--gram neg),and use RefSeq gbff file to first annotate from (--proteins):
```{bash}
prokka \ 
--outdir ${working_dir}/prokka_output_NS \ 
--force --genus Escherichia \ 
--species coli \ 
--strain Nissle \
--prefix NS \ 
--addgenes \ 
--locustag Nissle \ 
--kingdom Bacteria \ 
--gcode 11 \ 
--usegenus \ 
--rfam \ 
--gram neg \ 
--proteins ${working_dir}/GCF_003546975.1_ASM354697v1_genomic.gbff ${working_dir}/NS_spades/scaffolds.fasta
```

##### 5.2. PGAP
Draft assembly was also annotated using PGAP. First, we need to prepare three input files: 
* input.yaml
```
fasta: 
    class: File
    location: scaffolds.fasta
submol:
  class: File
  location: submol.yaml
```

* submol.yaml
```
organism:
    genus_species: 'Escherichia coli'
    strain: 'Nissle 1917'
```

* scaffolds.fasta

All three files should be put into the folder located in pgap directory:
```{bash}
# Enter the pgap folder git-cloned from the official repository (https://github.com/ncbi/pgap):
cd ${working_dir}/PGAP/pgap

# Create here a folder where three input files would be placed:
mkdir ${working_dir}/PGAP/pgap/Nissle 

# Create input.yaml (see above):
nano ${working_dir}/PGAP/pgap/Nissle/input.yaml

# Create submol.yaml
nano  ${working_dir}/PGAP/pgap/Nissle/submol.yaml

# Copy draft assmebly to ${working_dir}/PGAP/pgap/Nissle/ directory:
cp ${working_dir}/NS_spades/scaffolds.fasta ${working_dir}/PGAP/pgap/Nissle/
```

After that, if PGAP was installed correctrly, we could run the annotation pipeline. We would save output in ${working_dir}/PGAP/Nissle_output directory, skip errors from control analysis to obtain draft annotation (--ignore-all-errors), without reporting to NCBI (-n) and self-updating (--no-self-update):

```{bash}
${working_dir}/PGAP/pgap/scripts/pgap.py --ignore-all-errors -n --no-self-update -o ${working_dir}/PGAP/Nissle_output Nissle/input.yaml
```

#### 6. 16S rRNA gene homology search (NB sample)

Comparing *E. coli* str. Nissle 1917 (GenBank accession number GCA_003546975.1) with NS sample assembly via QUAST (see above) showed an average number of mismatches per 100 kbp aligned bases in assembly equal to 0.78 while assembly alignment to another *E. coli* str. Nissle 1917 genome (GenBank accession number GCA_000714595.1) yielded value of 2.06. Thus, *E. coli* str. Nissle 1917 genome GCA_003546975.1 was chosen as a reference genome for a subsequent analysis (See **Results/QUAST_report**).

We assumed that NB sample reads which failed to be mapped against the reference genome might belong to bacterial contaminant(s). We obtained draft assembly out of NB reads unmapped to the reference genome and extracted 16S rRNA gene. Its homologs were searched and classified in the SILVA database (https://www.arb-silva.de/aligner/ ). Minimum identity with query sequence was sset to 85, other parameters were default. The best hit with 100 % identity belonged to *Bacillus cereus* group.

```{bash}
# Indexing the reference genome (specifying Nissle2018 as prefix):
bowtie2-build ${working_dir}/GCA_003546975.1_ASM354697v1_genomic.fna ${working_dir}/Nissle2018

# Mapping:
bowtie2 -x ${working_dir}/Nissle2018 -1 ${path_to_reads}/ARyb_NB3_S54_R1_001.fastq.gz -2 ${path_to_reads}/ARyb_NB3_S54_R2_001.fastq.gz -S ${working_dir}/NB_Nissle2018.sam

# Sam to bam conversion:
samtools view -b ${working_dir}/NB_Nissle2018.sam  > ${working_dir}/NB_Nissle2018.bam

# Sorting:
samtools sort ${working_dir}/NB_Nissle2018.bam -o ${working_dir}/NB_Nissle2018.sorted.bam

# Indexing:
samtools index ${working_dir}/NB_Nissle2018.sorted.bam 

# Extracting unmapped reads:
samtools view -b -F 2 ${working_dir}/NB_Nissle2018.sorted.bam > ${working_dir}/NB_Nissle2018.unmapped.bam

# Sorting unmapped reads by name (required for bamToFastq v2.27.0):
samtools sort -n ${working_dir}/NB_Nissle2018.unmapped.bam -o ${working_dir}/NB_Nissle2018.unmapped.sortedName.bam

# Bam to fastq conversion:
bamToFastq -i ${working_dir}/NB_Nissle2018.unmapped.sortedName.bam -fq ${working_dir}/NB_Nissle2018.unmappedR1.fastq -fq2 ${working_dir}/NB_Nissle2018.unmappedR2.fastq

# Assembling unmapped reads:
spades.py --careful -o ${working_dir}/NBunmapped_spades -1 ${working_dir}/NB_Nissle2018.unmappedR1.fastq -2 ${working_dir}/NB_Nissle2018.unmappedR2.fastq

# Predicting rRNA genes:
barrnap -o ${working_dir}/rrna_NBunmappedNissle2018.fa ${working_dir}/NBunmapped_spades/scaffolds.fasta

# Extracting sequence(s) of 16S rRNA gene (print line matching the pattern "16S" and the next line):
awk '/16S/{print;getline;print;}' ${working_dir}/rrna_NBunmappedNissle2018.fa > ${working_dir}/16Srrna_NBunmappedNissle2018.fa
```


#### 7. Taxonomy classification using Kraken (NB sample)
Kraken v1.1.1 was used to identify taxonomy of the NB sample sequences (either initial raw reads or contigs obtained in previous step **6**  (```${working_dir}/NBunmapped_spades```)). HPC cluster should be probably used, as a kraken database preparation (```kraken-build --build ``` part) is a computationally intensive process.

##### 7.1. Build a Kraken customized database 
Before runing kraken, we should create custum database which would include only RefSeq complete bacterial genomes

Create a folder with the name of a custom database that would be built (e.g., BacDB):
```{bash}
mkdir ${working_dir}/BacDB
```
Download NCBI taxonomy files (the sequence ID to taxon map, the taxonomic names and tree information). The command would create a folder ```taxonomy/``` within a custom database directory (e.g., ```${working_dir}/BacDB```): 
```{bash}
kraken-build --download-taxonomy --db ${working_dir}/BacDB
```
*Note:  If the error "rsync_from_ncbi.pl: unexpected FTP path (new server?) for na"  was thrown, the commands below might solve the issue (from https://github.com/DerrickWood/kraken2/issues/226):*
```{bash}
awk -v FS='\t' '$20 != "na" {print $0}' ${working_dir}/BacDB/library/bacteria/assembly_summary.txt > ${working_dir}/BacDB/library/bacteria/new_assembly_summary.txt
cp ${working_dir}/BacDB/library/bacteria/new_assembly_summary.txt ${working_dir}/BacDB/library/bacteria/assembly_summary.txt
```
*After that, comment the line ```rm -f assembly_summary.txt```  within the scipt  ```download_genomic_library.sh```,  located at ```$(which download_genomic_library.sh)``` directory.*

Finally, install the library: download all the RefSeq bacterial genomes (about 86Gb size) to a folder ```library/bacteria/``` within your custom database directory. This task is the most computationally expensive. For instance, the following command was run on the computational node with 1Tb RAM during about 2 days:
```{bash}
kraken-build --build  --threads 24 --jellyfish-hash-size 4000000 --db ${working_dir}/BacDB
```

##### 7.2. Taxonomy assignation using customized database
After constructing a custom database, the following commands for taxonomic assignation of DNA reads against the custome database could be run:
```{bash}
# If as an input we would provide scaffolds of NB reads that could not map against the reference genome (see step 6 above):
kraken --db ${working_dir}/BacDB --threads 10 --unclassified-out ${working_dir}/unclassified_NBunmappedScaffolds  --classified-out ${path_out}/classified_NBunmappedScaffolds --output ${working_dir}/kraken_NBunmappedScaffolds.output --fasta-input ${working_dir}/NBunmapped_spades/scaffolds.fasta

# If as input files we would specify raw reads of the NB sample:
kraken --db ${path_db} --threads 24 --unclassified-out ${working_dir}/unclassified_NBreads --classified-out ${working_dir}/classified_NBreads --output ${working_dir}/kraken_NBreads.output  --paired  --fastq-input --gzip-compressed ${path_to_reads}/ARyb_NB3_S54_R1_001.fastq.gz ${path_to_reads}/ARyb_NB3_S54_R2_001.fastq.gz
```

Out of 4632 scaffolds,  4625 sequences were classified (99.85%). Among 6791092 read sequences, 6743684 sequences were assigned taxonomy to (99.30%) while other sequences failed to be classified. The rest sequuences in both cases remained unclassified. 

For displaying kraken results, kraken report could be used. To obtain report file,  run the command:
```{bash}
# in case of scaffolds 
kraken-report --db ${working_dir}/BacDB ${working_dir}/kraken_NBunmappedScaffolds.output > ${working_dir}/kraken_NBunmappedScaffolds.report

# in case of reads
kraken-report --db ${working_dir}/BacDB ${working_dir}/kraken_NBreads.output > ${working_dir}/kraken_NBreads.report
```
Kraken report was visualized at the web server Pavian (https://fbreitwieser.shinyapps.io/pavian/). According to kraken reports, NB sample was metagenome that consisted of bacteria from the species *E. coli* (including *E. coli* str. Nissle 1917) as well as the species *Bacillus cereus*. This outcome corresponded to the result obtained via 16S rRNA gene homology search: predicted 16S rRNA genes of possible contaminants were assigned to *Bacillus cereus* group and *E. coli* with identity varying from 99.88 to 100%. 


#### 8. Genome-wide comparison (NS sample)

We compared draft assembly (annotated via pgap) with *E. coli*  str. K-12 substr. MG1655 (GenBank assembly accesion GCA_000005845) and two complete genomes of *E. coli* str. Nissle 1917 (GenBank assembly accesion GCA_000714595 and GCA_003546975). Files in GBK format were manually retrieved from NCBI and renamed specifying the assembly accession. To download GBK file, go to nucleotide database where complete genome is placed (e.g. https://www.ncbi.nlm.nih.gov/nuccore/U00096.3), then click the button "Send to", and choose "Complete Record", "File" and format "GenBank(full)". For aligning genomes, Mauve (```progressiveMauve```) was used setting a minimum scaled breakpoint penalty to 5000:  
```{bash}
# copy annotation to the working direcotry:
cp ${working_dir}/PGAP/Nissle_output/annot.gbk  ${working_dir}/NS_pgap.gbk
progressiveMauve --min-scaled-penalty=5000 --output= ${working_dir}/Mauve_minscaledpen5000 ${working_dir}/NS_pgap.gbk  ${working_dir}/GCA_000714595_full.gbk ${working_dir}/GCA_003546975_full.gbk ${working_dir}/GCA_000005845_full.gbk > ./Mauve_minscaledpen5000.log
```
Manual inspecting the alignment, we determined location of 28 regions in the draft assembly that were not included into Locally Collinear Blocks (LCB). Respective coordinates were written down into the table and used for extracting sequences from draft assembly (see a piece of biopython code below). Obtained sequences were subjected to online BLASTN search with default parameters.
```{python}
# python3.7

import pandas as pd
import numpy as np
from Bio import SeqIO

# specify working direcotry
working_dir = /home/rybina/Nissle_project/

# import table where coordinates of regions between LCBs are specified. Table contains 2 columns named 'start' and 'end'
df = pd.read_csv(working_dir+'NS_regions_between_LCBs.csv', sep = ',')

# join scaffolds into a single sequence 
with open(working_dir + '/PGAP/Nissle_output/annot_merged.fna','w') as f:
    seq=''
    for record in SeqIO.parse(working_dir + '/PGAP/Nissle_output/annot.fna', 'fasta'):
        seq+=str(record.seq)
    f.write(f">NS_draft_assembly\n{seq}")
    
# Extract  respective sequences from the draft assembly
with open(working_dir+'/NS_regions_between_LCBs.fasta','w') as f_out:
    for i in range(len(df)):
        start=df['start'][i]
        end=df['end'][i]
        for record in SeqIO.parse(working_dir + '/PGAP/Nissle_output/annot_merged.fna'):
            f_out.write(f">start_{start}_end_{end}\n{record.seq[start:end]}\n")
```

BLASTN results demonstrated that all those regions varying by length from 1000 to 10000 bp shared 100% similarity with *E. coli* str. Nissle 1917 (GenBank assembly accession GCA_003546975).

### Conclusion

Sample *Nissle Big* is a metagenome that consists of bacteria from the species *E. coli* (including *E. coli* str. Nissle 1917) and the species *Bacillus cereus*.

From the sample *Nissle Small* we obtained  4,977,707 bp draft assembly which is 4,977,707 bp long and consists of 66 contigs with N50 value of 205,763. At 50.53 %, the GC content is similar to that of the *E. coli* Nissle 1917 genomes available at NCBI. PGAP annotation pipeline revealed 4,773 genes, among them 4,559 are protein-coding, three complete rRNAs (2 5S and 1 23S), 4 partial rRNAs (all 16S) and 58 tRNAs. Comparing NS sample and its assembly with the *E. coli* genomes demonstrated that the draft assembly covered the complete genome of *E. coli* str. Nissle 1917 (GenBank assembly accession: GCA_003546975.1) by 98% and more resembled *E. coli* str. Nissle 1917 genomes than *E. coli* str. K-12 substr. MG1655. These results suggest that the studied small colonies are indeed belong to *E. coli*  str. Nissle 1917.     

### References

_1. Nurk S, Bankevich A, Antipov D, Gurevich AA, Korobeynikov A, Lapidus A, Prjibelski AD, Pyshkin A, Sirotkin A, Sirotkin Y, Stepanauskas R, Clingenpeel SR, Woyke T, McLean JS, Lasken R, Tesler G, Alekseyev MA, Pevzner PA. Assembling single-cell genomes and mini-metagenomes from chimeric MDA products. J Comput Biol. 2013 Oct;20(10):714-37. doi: 10.1089/cmb.2013.0084._

_2. Alexey Gurevich, Vladislav Saveliev, Nikolay Vyahhi and Glenn Tesler. QUAST: quality assessment tool for genome assemblies, Bioinformatics (2013) 29 (8): 1072-1075. doi: 10.1093/bioinformatics/btt086_

_3. Johannes Alneberg, Brynjar Smári Bjarnason, Ino de Bruijn, Melanie Schirmer, Joshua Quick, Umer Z Ijaz, Leo Lahti, Nicholas J Loman, Anders F Andersson & Christopher Quince. 2014. Binning metagenomic contigs by coverage and composition. Nature Methods, doi: 10.1038/nmeth.3103_

_4. Parks DH, Imelfort M, Skennerton CT, Hugenholtz P, Tyson GW. 2015. CheckM: assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes. Genome Research, 25: 1043–1055._

_5. Langmead, B., Salzberg, S. Fast gapped-read alignment with Bowtie 2. Nat Methods 9, 357–359 (2012). doi: 10.1038/nmeth.1923_

_6. Li H. A statistical framework for SNP calling, mutation discovery, association mapping and population genetical parameter estimation from sequencing data. Bioinformatics. 2011 Nov 1;27(21):2987-93. doi: 10.1093/bioinformatics/btr509._

_7. Quinlan AR, Hall IM. BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics. 2010 Mar 15;26(6):841-2. doi: 10.1093/bioinformatics/btq033._

_8. Seemann T. Prokka: rapid prokaryotic genome annotation. Bioinformatics. 2014 Jul 15;30(14):2068-9. PMID:24642063_

_9. Tatusova T, DiCuccio M, Badretdin A, Chetvernin V, Nawrocki EP, Zaslavsky L, Lomsadze A, Pruitt KD, Borodovsky M, Ostell J. NCBI prokaryotic genome annotation pipeline. Nucleic Acids Res. 2016 Aug 19;44(14):6614-24. doi: 10.1093/nar/gkw569._

_10. Seemann T. Barrnap 0.9 : rapid ribosomal RNA prediction._

_11. Wood DE, Salzberg SL: Kraken: ultrafast metagenomic sequence classification using exact alignments. Genome Biology 2014, 15:R46._

_12. Darling AC, Mau B, Blattner FR, Perna NT. Mauve: multiple alignment of conserved genomic sequence with rearrangements. Genome Res. 2004;14(7):1394-1403. doi: 10.1101/gr.2289704_

_13. Cock PA, Antao T, Chang JT, Chapman BA, Cox CJ, Dalke A, Friedberg I, Hamelryck T, Kauff F, Wilczynski B and de Hoon MJL (2009) Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics, 25, 1422-1423_

_14. Finn RD et al. HMMER web server: interactive sequence similarity searching. Nucleic Acids Res. 2011 Jul;39(Web Server issue):W29-37._

_15. Hyatt D et al. Prodigal: prokaryotic gene recognition and translation initiation site identification. BMC Bioinformatics. 2010 Mar 8;11:119._

_16. Stajich et al, The Bioperl toolkit: Perl modules for the life sciences. Genome Res. 2002 Oct;12(10):1611-8._

_17. O. Tange, GNU Parallel - The Command-Line Power Tool, ;login: The USENIX Magazine, Feb 2011:42-47._

_18. Camacho C et al. BLAST+: architecture and applications. BMC Bioinformatics. 2009 Dec 15;10:421._

_19. Hyatt D et al. Prodigal: prokaryotic gene recognition and translation initiation site identification. BMC Bioinformatics. 2010 Mar 8;11:119._

_20 Laslett D, Canback B. ARAGORN, a program to detect tRNA genes and tmRNA genes in nucleotide sequences. Nucleic Acids Res. 2004 Jan 2;32(1):11-6._

_21. Petersen TN et al. SignalP 4.0: discriminating signal peptides from transmembrane regions. Nat Methods. 2011 Sep 29;8(10):785-6._

_22. D. L. Kolbe, S. R. Eddy. Fast Filtering for RNA Homology Search. Bioinformatics, 27:3102-3109, 2011._
