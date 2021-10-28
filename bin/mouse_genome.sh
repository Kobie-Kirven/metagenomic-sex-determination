##################################################
# Download mouse genome and build bowtie2 index
##################################################

cd /gpfs/group/exd44/default/kjk6173/metagenomic-sex-pipeline/data/10-06-2021/reference_genome
#Download mouse genome
wget http://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.masked.gz

cd /gpfs/group/exd44/default/kjk6173/metagenomic-sex-pipeline/data/10-06-2021/bowtie2_index
#Build bowtie2 index
bowtie2-build  /gpfs/group/exd44/default/kjk6173/metagenomic-sex-pipeline/data/10-06-2021/reference_genome/mm39.fa.masked.gz \
mm39.fa.masked.gz