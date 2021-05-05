#Methods for the analyzing rDNA and histone copy number from *Potamopyrgus* Illumina paired-end threads

##Quality control of Illumina reads

make file with list of fastq basenames
```sh
ls *_R1.fastq.gz | sed 's/_R1.fastq.gz//g' > jobs
```
Trim reads with Sickle using array job submission:
```bash
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
```bash
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
