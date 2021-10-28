###################################################################
# Simulate shotgun reads from mock genomes using wgsim
###################################################################

# Author - Kobie Kirven
# Data - 10-07-2021

#Paths and file names
MOCK_PATH=../data/10-06-2021/reference_genome/
mkdir ${MOCK_PATH}reads


# Simulate 10,000 paired-end reads from male genome
wgsim -N 10000 ${MOCK_PATH}male_genome.fa.gz ${MOCK_PATH}reads/male_10000_1.fa \
${MOCK_PATH}reads/male_10000_2.fa

# Simulate 10,000 paired-end reads from female genome
wgsim -N 10000 ${MOCK_PATH}female_genome.fa.gz ${MOCK_PATH}reads/female_10000_1.fa \
${MOCK_PATH}reads/female_10000_2.fa