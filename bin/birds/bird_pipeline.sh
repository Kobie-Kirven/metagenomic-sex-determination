################################################################################
# This pipeline will dowonload sequences from "Metagenomic analysis reveals 
# the microbiome and resistome in migratory birds" (PRJNA556790), specifically
# those birds from Tibet.
#
# It then preforms quality trimming on the downloaded sequences and 
# aligns reads to the chicken genome.
################################################################################

#Author - Kobie Kirven
#Date - 10-21-2021

#Tools Needed:
# - esearch
# - efetch
# - trimmomatic
# - Bowtie2
# - GNU parallel

conda activate metagenomicPipeline
module load gcc/8.3.1
module load bowtie2

#Data Directory
DATA=data/10-20-2021
INDEX=bowtie2_databases/GCA_000002315.5_GRCg6a_genomic.fna.gz
mkdir -p reports

seq_depth=10000000

#Run the pipeline
cat ../../${DATA}/SraRunTable.txt | grep "SRR" | grep "Tibet"  \
| cut -d "," -f 1 > tibet_runs.txt

while read run; do
    #Download the FASTQ files at the specified read depth
    fastq-dump -X ${seq_depth} ${run}

    scims -i ../../data/10-20-2021/bowtie2_databases/GCA_000002315.5_GRCg6a_genomic.fna.gz \
    -sca chicken_scaf_ids.txt \
    -r ../../data/10-20-2021/reference_genomes/GCA_000002315.5_GRCg6a_genomic.fna.gz \
    -s ${run}.fastq -t 3 -het CM000121.5 -hom CM000122.5 \
    -o ../../results/tibet_${run}
    
    rm -rf ${run}.fastq

done <tibet_runs.txt

# trimmomatic SE SRR*.fastq ILLUMINACLIP:adapter.fa:2:30:5 SLIDINGWINDOW:4:20

# fastqc SRR*_trimmed.fastq -o reports

#trimmomatic {}.fastq {}_trimmed.fastq 

