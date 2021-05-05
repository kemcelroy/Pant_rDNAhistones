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

Map reads to the single copy exons and measure coverage using array script
```sh
#!/bin/sh
#$ -m ea # send an email when jobs end or abort
#$ -M kyle-mcelroy@uiowa.edu # email address
#$ -cwd # set to current working directory
#$ -pe smp 12 # use 12 cpus
#$ -N map_cov_exon # job name
#$ -t 1-29 # run 29 jobs

# map_cov_exon.sh, to run: map_cov_exon.sh

# pwd=/Users/kemcelroy/rDNA_histones

# load samtools Version: 1.3.1 (using htslib 1.3.2)
module load samtools

# collect the input names from the text file "jobs"
names=($(cat jobs))
echo ${names[${SGE_TASK_ID}]}

# map reads to the set of exons
bwa mem -t 12 0_index/longest_exons_per_busco_gt300.fa "/nfsscratch/Users/kemcelroy/reads/${names[${SGE_TASK_ID}]}_R1.trim.fastq.gz" "/nfsscratch/Users/kemcelroy/reads/${names[${SGE_TASK_ID}]}_R2.trim.fastq.gz" > "/nfsscratch/Users/kemcelroy/rDNA_histones/1_sc-exons/${names[${SGE_TASK_ID}]}_longest_exons_per_busco_gt300.sam"
# convert to sorted bam file
samtools view -bS -q 20 -@ 12 "/nfsscratch/Users/kemcelroy/rDNA_histones/1_sc-exons/${names[${SGE_TASK_ID}]}_longest_exons_per_busco_gt300.sam" | samtools sort -@ 12 > "1_sc-exons/${names[${SGE_TASK_ID}]}_longest_exons_per_busco_gt300.q20.bam"
# index the sorted bam file
samtools index "1_sc-exon/${names[${SGE_TASK_ID}]}_longest_exons_per_busco_gt300.q20.bam"
# calculate per-base coverage
bedtools genomecov -d -ibam "1_sc-exons/${names[${SGE_TASK_ID}]}_longest_exons_per_busco_gt300.q20.bam" -g 0_index/longest_exons_per_busco_gt300.fa > "1_sc-exon/${names[${SGE_TASK_ID}]}_longest_exons_per_busco_gt300.coverage.txt"
# compress the sam file
gzip "/nfsscratch/Users/kemcelroy/rDNA_histones/1_sc-exons/${names[${SGE_TASK_ID}]}_longest_exons_per_busco_gt300.sam"
```

Identify exons that are single copy and consistent coverage for all samples
```sh
cd 1_sc-exons
mkdir split
cp ../jobs ./
# Make a directory for each lineag, then copy the exon coverage file and split it by exon
while read line; do mkdir $line; done < jobs
while read line; do mv $line split; done < jobs
for f in *.txt; do cp $f "split/${f%_longest_exons_per_busco_gt300.coverage.txt}"; done
for d in *; do cd $d && awk -F '\t' '{print>$1}' *.coverage.txt; cd /Users/kemcelroy/rDNA_histones/1_sc-exons/split/; done
# Collect the central tendencies for each exon for each lineag
for d in *; do cd $d && (for f in E*; do sh ../../../stats2.sh $f >> summary.txt; done); cd /Users/kemcelroy/rDNA_histones/1_sc-exons/split; done
for d in *; do cd $d && sed '/,/d' summary.txt | awk '{if ($3 > 0 && $4 > 0 && $5 > 0 && ($3/$4 > 0.9 && $3/$4 < 1.1) && ($3/$5 > 0.9 && $3/$5 < 1.1) && ($4/$5 > 0.9 && $4/$5 < 1.1)) {print $0}}' > central_measures.txt; cd /Users/kemcelroy/rDNA_histones/1_sc-exons/split/; done
cat */central_measures.txt | awk '{print $2}' | awk '{count[$1]++} END {for (word in count) print word, count[word]}' | awk '{if ($2>19) {print $1}}' > consistent
wc -l consistent
12
# The top 12 most evenly covered exons across the 28 lineages
for d in $(cat dirs); do cd $d && (while read line; do grep $line central_measures.txt  >> ../consistent_exon; done < ../consistent); cd /Users/kemcelroy/rDNA_histones/1_sc-exons/split; done
# Collect the per-base coverage for all of the final exons
# 2 exons were clearly multicopy for all lineages, removed those for a final 10 exons
echo Lineage$'\t'Sequence$'\t'Position$'\t'Depth > consistent_coverages
for d in $(cat dirs); do cd $d && (while read line; do grep $line *coverage.txt | awk -v OFS='\t' '{print var,$0}' var="${d}" >> ../consistent_coverages; done < ../consistent_10); cd /Users/kemcelroy/rDNA_histones/1_sc-exons/split; done

# Make a fasta file of the top 10 most consistent evenly covered single copy exons
grabContigs_KMedits2.py longest_exons_per_busco_gt300.fa top10consistent top10consistent.fa
```
