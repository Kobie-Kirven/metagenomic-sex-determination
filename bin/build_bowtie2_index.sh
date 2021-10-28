################################################
# Build the index for Bowtie2
###############################################

#Author - Kobie Kirven
#Date - 10-06-2021

#Path to reference genome
REF_PATH=../data/10-06-2021/reference_genome/GRCh38_latest_genomic.fna.gz

#Build the bowtie2 index
bowtie2-build ${REF_PATH} ../data/10-06-2021/bowtie2_index/GRCh38_latest_genomic