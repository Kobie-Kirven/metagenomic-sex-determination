###################################################################
# Vary Sequencing depth
###################################################################

# Author - Kobie Kirven
# Data - 10-07-2021


#Paths and file names
MOCK_PATH=../data/10-06-2021/reference_genome/
DATA_PATH=../data/10-06-2021/
mkdir ${MOCK_PATH}reads

for depth in 100 250 500 1000 2500 5000 7500
do
    echo $depth
    # Simulate paired-end reads from male genome
    wgsim -N $depth ${MOCK_PATH}male_genome.fa.gz ${MOCK_PATH}reads/male_${depth}_1.fa \
    ${MOCK_PATH}reads/male_${depth}_2.fa > /dev/null 2>&1
    
    # Simulate paired-end reads from female genome
    wgsim -N $depth ${MOCK_PATH}female_genome.fa.gz ${MOCK_PATH}reads/female_${depth}_1.fa \
    ${MOCK_PATH}reads/female_${depth}_2.fa > /dev/null 2>&1
    
    #Align male reads to reference genome
    bowtie2 --local -x ${DATA_PATH}bowtie2_index/GRCh38_latest_genomic \
    -1 ${DATA_PATH}reference_genome/reads/male_${depth}_1.fa \
    -2 ${DATA_PATH}reference_genome/reads/male_${depth}_2.fa \
    > ../data/10-06-2021/alignments/male_${depth}_align.sam
    
    #Align female reads to reference genome
    bowtie2 --local -x ${DATA_PATH}bowtie2_index/GRCh38_latest_genomic \
    -1 ${DATA_PATH}reference_genome/reads/female_${depth}_1.fa \
    -2 ${DATA_PATH}reference_genome/reads/female_${depth}_2.fa \
    > ../data/10-06-2021/alignments/female_${depth}_align.sam
    
    python3 read_dp.py --depth $depth

done