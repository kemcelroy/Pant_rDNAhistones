# Methods for the analyzing rDNA and histone copy number from *Potamopyrgus* Illumina paired-end threads

## Quality control of Illumina reads

make file with list of fastq basenames
```sh
ls *_R1.fastq.gz | sed 's/_R1.fastq.gz//g' > jobs
```
Trim reads with Sickle using array job submission:
```sh
#!/bin/sh
#$ -m ea # send an email when jobs end or abort
#$ -M kyle-mcelroy@uiowa.edu # email address
#$ -cwd # set to current working directory
#$ -N sickle # job name
#$ -t 1-28 # run 28 jobs

# sickle.sh, to run: sickle.sh
# using sickle version 1.33

# pwd=/Users/kemcelroy/reads

# collect the input names from the text file "jobs"
names=($(cat jobs))
echo ${names[${SGE_TASK_ID}]}

/Users/kemcelroy/downloads/sickle/sickle pe \
	-f "${names[${SGE_TASK_ID}]}_R1.fastq.gz" \
	-r "${names[${SGE_TASK_ID}]}_R2.fastq.gz" \
	-t sanger \
	-o "${names[${SGE_TASK_ID}]}_R1.trim.fastq.gz" \
	-p "${names[${SGE_TASK_ID}]}_R2.trim.fastq.gz" \
	-s "${names[${SGE_TASK_ID}]}_S.trim.fastq.gz"
```
Check the read quality with FastQC using array job submission:
```sh
#!/bin/sh
#$ -m ea # send an email when jobs end or abort
#$ -M kyle-mcelroy@uiowa.edu # email address
#$ -cwd # set to current working directory
#$ -N fastqc # job name
#$ -t 1-28 # run 28 jobs

# fastqc.sh, to run: fastqc.sh
# using FastQC v0.11.7

# pwd=/Users/kemcelroy/reads

# collect the input names from the text file "jobs"
names=($(cat jobs))
echo ${names[${SGE_TASK_ID}]}

# make a directory for each set of output
FastQC/"FQC_${names[${SGE_TASK_ID}]}"

/Users/kemcelroy/downloads/FastQC/fastqc "${names[${SGE_TASK_ID}]}_R1.trim.fastq.gz" "${names[${SGE_TASK_ID}]}_R2.trim.fastq.gz" \
	-o FastQC/"FQC_${names[${SGE_TASK_ID}]}"
```

## Generate single copy BUSCO exon sequences
These sequences are used to estimate per haploid copy number of rDNA and histones

Run BUSCO job script
```sh
#!/bin/sh
#$ -m ea
#$ -M kyle-mcelroy@uiowa.edu
#$ -N busco_Pant
#$ -cwd
# BUSCO 3.0.2
export PATH=/Users/kemcelroy/downloads/augustus.2.5.5/bin:/Users/kemcelroy/downloads/augustus.2.5.5/scripts:/Users/kemcelroy/downloads/augustus.2.5.5/config:$PATH
export AUGUSTUS_CONFIG_PATH="/Users/kemcelroy/downloads/augustus.2.5.5/config"
/Users/kemcelroy/downloads/busco/scripts/run_BUSCO.py -i Pant_hapt.fasta -o BUSCO_Pant_hap -l /Users/kemcelroy/downloads/busco/metazoa_odb9/ -m genome
```

First, make a gff with all of the AUGUSTUS predicted exons, numbered for each gene
```bash
for f in *.gff; do awk -v OFS='\t' '{if ($3=="exon") {print $1,$4,$5,var"_exon-",0,$7}}' var="${f%.gff}" $f >> busco_single_copy_exons.gff; done
```
add numbers to each exon
```sh
awk '{ print $0 "\t" ++count[$4] }' busco_single_copy_exons.gff | awk -v OFS='\t' '{print $1,$2,$3,$4$7,$5,$6}' > busco_single_copy_exons-num.gff
```
make a fasta of the exons
```sh
bedtools getfasta -name -s -fi Pant_hap.fasta -bed busco_single_copy_exons-num.gff > busco_single_copy_exons-num.fa
sed -i 's/:.*//' busco_single_copy_exons-num.fa
```
select the longest exon for each BUSCO gene and filter based on length (300 bp)
```bash
samtools faidx busco_single_copy_exons-num.fa
awk -F '[_\t]' '{print $1,$2,$3}' busco_single_copy_exons-num.fa.fai | awk '$3>max[$1]{max[$1]=$3; row[$1]=$0} END{for (i in row) print row[i]}' | awk '{print $1"_"$2,$3}' OFS='\t' > longest_exons_per_busco
awk '{if ($2>300) {print $1}}' longest_exons_per_busco > longest_exons_per_busco_gt300
```
use Joel's script to make a fasta of the selected > 300 bp long exons ... tried a grep loop but it kept not working
```sh
grabContigs_KMedits2.py busco_single_copy_exons-num.fa longest_exons_per_busco_gt300 longest_exons_per_busco_gt300.fa
bwa index longest_exons_per_busco_gt300.fa
```
