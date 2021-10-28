##############################################################
# This script downloads the human reference genome from NCBI
##############################################################

# Author - Kobie Kirven
# Date - 10-06-2021

#Link to file
LINK=https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers
FILE_NAME=GRCh38_latest_genomic.fna.gz

#Download the reference genome from NCBI
wget -P ../data/10-06-2021/reference_genome ${LINK}/${FILE_NAME}