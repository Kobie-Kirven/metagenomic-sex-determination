# Download the chicken reference genome and build the bowtie2 reference database

# Author - Kobie Kirven 

#Create the directories 
mkdir -p ../../data/10-20-2021/reference_genomes
mkdir -p ../../data/10-20-2021/bowtie2_databases

#Download the chicken reference genome
cd ../../data/10-20-2021/reference_genomes
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/002/315/GCA_000002315.5_GRCg6a/GCA_000002315.5_GRCg6a_genomic.fna.gz


#Build the bowtie2 database
cd ../../data/10-20-2021/bowtie2_databases
bowtie2-build ../reference_genomes/GCA_000002315.5_GRCg6a_genomic.fna.gz \
GCA_000002315.5_GRCg6a_genomic.fna.gz

