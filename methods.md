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

Identify candidate sequences for single copy use
```sh
# First, make a gff with all of the AUGUSTUS predicted exons, numbered for each gene
for f in *.gff; do awk -v OFS='\t' '{if ($3=="exon") {print $1,$4,$5,var"_exon-",0,$7}}' var="${f%.gff}" $f >> busco_single_copy_exons.gff; done
# add numbers to each exon
awk '{ print $0 "\t" ++count[$4] }' busco_single_copy_exons.gff | awk -v OFS='\t' '{print $1,$2,$3,$4$7,$5,$6}' > busco_single_copy_exons-num.gff
# make a fasta of the exons
bedtools getfasta -name -s -fi /Users/kemcelroy/genome/hap.Dovetail_02012019.racon6.arcs-pilon2.3.gapfiller2.sspace_CCS2_corrected4.racon13.arcs-pilon4.gapfiller2.sspace_CCS1_corrected.6.racon8.pilon7.fasta -bed busco_single_copy_exons-num.gff > busco_single_copy_exons-num.fa
sed -i 's/:.*//' busco_single_copy_exons-num.fa
# select the longest exon for each BUSCO gene and filter based on length (300 bp)
samtools faidx busco_single_copy_exons-num.fa
awk -F '[_\t]' '{print $1,$2,$3}' busco_single_copy_exons-num.fa.fai | awk '$3>max[$1]{max[$1]=$3; row[$1]=$0} END{for (i in row) print row[i]}' | awk '{print $1"_"$2,$3}' OFS='\t' > longest_exons_per_busco
awk '{if ($2>300) {print $1}}' longest_exons_per_busco > longest_exons_per_busco_gt300
# use Joel's script to make a fasta of the selected > 300 bp long exons ... tried a grep loop but it kept not working
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

## Copy number estimation and analysis

Job array script for mapping reads to the single copy exon and rDNA & histone sequences
```sh
#!/bin/sh
#$ -m ea # send an email when jobs end or abort
#$ -M kyle-mcelroy@uiowa.edu # email address
#$ -cwd # set to current working directory
#$ -pe 80cpn 80 # use 80 cpus
#$ -N map_cov_exon # job name
#$ -t 1-29 # run 29 jobs

# map_cov_rDNA_histone_sc-exons.sh, to run: map_cov_rDNA_histone_sc-exons.sh

# pwd=/Users/kemcelroy/rDNA_histones

# load samtools Version: 1.3.1 (using htslib 1.3.2)
module load samtools

# collect the input names from the text file "jobs2"
names=($(cat jobs2))
echo ${names[${SGE_TASK_ID}]}

# map reads to the set of single copy exons along withe genes of interest: 18S, 5.8S, 28S, H2A, H2B, H3, H4, 5S, H1
bwa mem -t 80 0_index/rDNA_histone_sc-exons.fa "/nfsscratch/Users/kemcelroy/reads/${names[${SGE_TASK_ID}]}_R1.trim.fastq.gz" "/nfsscratch/Users/kemcelroy/reads/${names[${SGE_TASK_ID}]}_R2.trim.fastq.gz" > "/nfsscratch/Users/kemcelroy/rDNA_histones/2_cn/${names[${SGE_TASK_ID}]}_rDNA_histone_sc-exons.sam"
# convert to sorted bam file
samtools view -bS -q 20 -@ 80 "/nfsscratch/Users/kemcelroy/rDNA_histones/2_cn/${names[${SGE_TASK_ID}]}_rDNA_histone_sc-exons.sam" | samtools sort -@ 80 > "2_cn/${names[${SGE_TASK_ID}]}_rDNA_histone_sc-exons.q20.bam"
# index the sorted bam file
samtools index "2_cn/${names[${SGE_TASK_ID}]}_rDNA_histone_sc-exons.q20.bam"
# calculate per-base coverage
bedtools genomecov -d -ibam "2_cn/${names[${SGE_TASK_ID}]}_rDNA_histone_sc-exons.q20.bam" -g 0_index/rDNA_histone_sc-exons.fa > "2_cn/${names[${SGE_TASK_ID}]}_rDNA_histone_sc-exons.coverage"
# compress the sam file
gzip "/nfsscratch/Users/kemcelroy/rDNA_histones/2_cn/${names[${SGE_TASK_ID}]}_rDNA_histone_sc-exons.sam"
```

Analyze the coverage results to measure copy number
```sh
# Calculate copy number from the median coverage of rDNA and histones divided by the median coverage of single copy exons
while read line; do sh median_calc.sh $line; done < jobs
rm skip_median_cov.txt
awk '{print $2}' Yellow2_median_cov.txt | sh transpose.awk | awk 'FNR<2 {print "lineage",$0}' OFS='\t' > median_vals.txt
awk -v OFS='\t' '{print $1,$3,$4,$5,$6,$7,$8,$9,$10,$11}' median_vals.txt > copy_numbers.txt
for f in *_cov.txt; do awk '{print $3}' $f | sh transpose.awk | awk -v OFS='\t' '{print var,$0}' var="${f%_median_cov.txt}" >> median_vals.txt; done
awk -v OFS='\t' 'FNR>1 {print $1, $3/$2, $4/$2, $5/$2, $6/$2, $7/$2, $8/$2, $9/$2, $10/$2, $11/$2}' median_vals.txt >> copy_numbers.txt

# generate the file to use in R to make the boxplots
# status is a tab separated file, the first column is the lineage name and the second is "Sex" or "Asex"
# I made a copy number file without Pest
sed '/Pest/d' copy_numbers.txt > copy_numbers_noPest.txt
# merge copy number info with reproductive mode
awk 'NR==FNR{a[NR]=$0;next}{print a[FNR],$0}' status copy_numbers_noPest.txt | awk -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' > copy_number_mode.txt
# format the copy number information into the format I want to make figures in R
# the "tail -n +10" is to just to get rid of the repeated header information, a result of me not quite being able to skip the first row in the format.sh script
echo Mode$'\t'Lineage$'\t'Seq$'\t'CN > rDNA_histone_cpn_mode_h1
sh format.sh copy_number_mode.txt | tail -n +10 >> rDNA_histone_cpn_mode_h1
# make another without H1
sed '/H1/d' rDNA_histone_cpn_mode_h1 > rDNA_histone_cpn_mode
```


## Some scripts I used in calculating median coverage and other stats
```sh
#!/bin/sh
# median_calc.sh
sample=$1
grep exon "${1}_rDNA_histone_sc-exons.coverage" | sort -n -r -k3 | awk -f median.awk | awk -v OFS='\t' '{print var,"exon",$1}' var="${sample}" >> "${1}_median_cov.txt"
grep 5.8S "${1}_rDNA_histone_sc-exons.coverage" | sort -n -r -k3 | awk -f median.awk | awk -v OFS='\t' '{print var,"5.8S",$1}' var="${sample}" >> "${1}_median_cov.txt"
grep 28S "${1}_rDNA_histone_sc-exons.coverage" | sort -n -r -k3 | awk -f median.awk | awk -v OFS='\t' '{print var,"28S",$1}' var="${sample}" >> "${1}_median_cov.txt"
grep 18S "${1}_rDNA_histone_sc-exons.coverage" | sort -n -r -k3 | awk -f median.awk | awk -v OFS='\t' '{print var,"18S",$1}' var="${sample}" >> "${1}_median_cov.txt"
grep 5S "${1}_rDNA_histone_sc-exons.coverage" | sort -n -r -k3 | awk -f median.awk | awk -v OFS='\t' '{print var,"5S",$1}' var="${sample}" >> "${1}_median_cov.txt"
grep H2A "${1}_rDNA_histone_sc-exons.coverage" | sort -n -r -k3 | awk -f median.awk | awk -v OFS='\t' '{print var,"H2A",$1}' var="${sample}" >> "${1}_median_cov.txt"
grep H2B "${1}_rDNA_histone_sc-exons.coverage" | sort -n -r -k3 | awk -f median.awk | awk -v OFS='\t' '{print var,"H2B",$1}' var="${sample}" >> "${1}_median_cov.txt"
grep H3 "${1}_rDNA_histone_sc-exons.coverage" | sort -n -r -k3 | awk -f median.awk | awk -v OFS='\t' '{print var,"H3",$1}' var="${sample}" >> "${1}_median_cov.txt"
grep H4 "${1}_rDNA_histone_sc-exons.coverage" | sort -n -r -k3 | awk -f median.awk | awk -v OFS='\t' '{print var,"H4",$1}' var="${sample}" >> "${1}_median_cov.txt"
grep H1 "${1}_rDNA_histone_sc-exons.coverage" | sort -n -r -k3 | awk -f median.awk | awk -v OFS='\t' '{print var,"H1",$1}' var="${sample}" >> "${1}_median_cov.txt"
```

```sh
#!/bin/sh
# stats2.sh
SEQ=$1
cwd=$(pwd)
LINEAGE="${cwd##*/}"
awk -F'\t' '
{col=$3}{if((col ~  /^-?[0-9]*([.][0-9]+)?$/) && ($0!=""))
{
     sum+=col;
     a[x++]=col;
     b[col]++
     if(b[col]>hf){hf=b[col]}
}
}
END{n = asort(a);idx=int((x+1)/2)
     for (i in b){if(b[i]==hf){(k=="") ? (k=i):(k=k FS i)}{FS=","}}
     print ln,sq,sum/x,((idx==(x+1)/2) ? a[idx] : (a[idx]+a[idx+1])/2),k
}' sq=$SEQ ln=$LINEAGE  OFS='\t' $1
```

```sh
#!/bin/sh
# format.sh
awk -v OFS='\t' '
	{
		print $1,$2,"18S",$5
		print $1,$2,"5.8S",$3
		print $1,$2,"28S",$4
		print $1,$2,"H2A",$7
		print $1,$2,"H2B",$8
		print $1,$2,"H3",$9
		print $1,$2,"H4",$10
		print $1,$2,"5S",$6
		print $1,$2,"H1",$11
	}' $1
  ```
```sh
#transpose.awk, run as: sh transpose.awk FILE
awk '
{
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str" "a[i,j];
        }
        print str
    }
}' $1
```

```sh
#/usr/bin/env awk
#median.awk script, run as: sh median.awk FILE
{
    count[NR] = $3;
}
END {
    if (NR % 2) {
        print count[(NR + 1) / 2];
    } else {
        print (count[(NR / 2)] + count[(NR / 2) + 1]) / 2.0;
    }
}
}
```
