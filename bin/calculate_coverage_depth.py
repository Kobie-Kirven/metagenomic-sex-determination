################################################################
# Calculate the coverage for each chromosome
################################################################

#Imports
from Bio import SeqIO
import gzip
import matplotlib
import matplotlib.pyplot as plt

# File paths
small_path = "../data/10-06-2021/"

# Get the lengths of all the chromosomes
chrom_lengths = {}
with gzip.open(small_path + "reference_genome/GRCh38_latest_genomic.fna.gz", "rt") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        if "NW" not in record.id and "NT" not in record.id:
            chrom_lengths[record.id] = len(str(record.seq))
        
chrom_cov = {}
# Calculate total coverage for each chromosome - female
with open(small_path + "alignments/female_10000_align_coverage.txt") as fn:
    lines = fn.readlines()
    for line in lines:
        line = line.strip("\n").split("\t")
        if "NW" not in line[0] and "NT" not in line[0]:
            if line[0][1:] in chrom_cov:
                chrom_cov[line[0]] += int(line[2])
            else:
                chrom_cov[line[0]] = int(line[2])
     
avg_cov = {}
#Compute average coverage
for ele in chrom_cov:
    average = chrom_cov[ele] / chrom_lengths[ele]
    avg_cov[ele] = average
    


id_list = []
avg_list = []
for avg in avg_cov:
    id_list.append(avg)
    avg_list.append(avg_cov["NC_000023.11"]/avg_cov[avg])
    
plt.scatter(id_list, avg_list)
plt.title("X:Autosomal Coverage Ratio - Female")
plt.ylabel('Average Coverage')
plt.xlabel('Chromosome ID')
plt.xticks(rotation = -90)
plt.savefig('../results/figures/x_to_autosome_female_10000.png', bbox_inches='tight')

####################

chrom_cov = {}
# Calculate total coverage for each chromosome - female
with open(small_path + "alignments/male_10000_align_coverage.txt") as fn:
    lines = fn.readlines()
    for line in lines:
        line = line.strip("\n").split("\t")
        if "NW" not in line[0] and "NT" not in line[0]:
            if line[0][1:] in chrom_cov:
                chrom_cov[line[0]] += int(line[2])
            else:
                chrom_cov[line[0]] = int(line[2])
     
avg_cov = {}
#Compute average coverage
for ele in chrom_cov:
    average = chrom_cov[ele] / chrom_lengths[ele]
    avg_cov[ele] = average
    


id_list = []
avg_list = []
for avg in avg_cov:
    id_list.append(avg)
    avg_list.append(avg_cov["NC_000023.11"]/avg_cov[avg])
    
plt.scatter(id_list, avg_list)
plt.title("X:Autosomal Coverage Ratio - male")
plt.ylabel('Average Coverage')
plt.xlabel('Chromosome ID')
plt.xticks(rotation = -90)
plt.savefig('../results/figures/x_to_autosome_male_10000.png', bbox_inches='tight')
