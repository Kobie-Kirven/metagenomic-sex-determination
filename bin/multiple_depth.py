#########################################
# Collect data from multiple runs
#########################################

#Imports
from Bio import SeqIO
import gzip
import matplotlib
import matplotlib.pyplot as plt
import argparse
import os
import statistics as st

parser = argparse.ArgumentParser(description='Zonkey method differnt depths')

parser.add_argument('--depth', dest='depth',help='sequencing depth')

args = parser.parse_args()

class ParseSam:
    """Parse a SAM file and return a SAM record
    for each sequence in the SAM file

    usage example:
    for sam_record in ParseSam("sam_file.sam"):
            print(rec.query)

    """

    def __init__(self, samFile):
        self.samFile = samFile

    def __iter__(self):
        self.fn = open(self.samFile)
        return self

    def __next__(self):
        line = self.fn.__next__()
        while line.startswith("@"):
            line = self.fn.__next__()
        samData = ReadSamLine(line)
        samData.getFeatures()
        return samData


class ReadSamLine:
    """
    Read a line from a SAM file and get each of the fields.
    Naming conventions are based on the SAM format

    """

    def __init__(self, samLine):
        self.samLine = samLine

    def getFeatures(self):
        """
        Takes a line from a SAM file as
        a string collects information
        from the fields
        """

        if not self.samLine.startswith("@"):
            samFields = self.samLine.strip("\n").split("\t")
            self.query = samFields[0]
            self.flag = int(samFields[1])
            self.rnam = samFields[2]
            self.pos = int(samFields[3])
            self.mapq = int(samFields[4])
            self.cigar = samFields[5]
            self.rnext = samFields[6]
            self.pnext = samFields[7]
            self.tlen = int(samFields[8])
            self.seq = samFields[9]
            self.qual = samFields[10]
            self.align_score = None
            self.align_score_best = None
            self.align_score_mate = None
            self.num_ambig = None
            self.mismatches = None
            self.gap_opens = None
            self.gap_ext = None
            self.edit_distance = None
            self.why_filtered = None
            self.yt = None
            self.mismatch_ref_bases = None

            if len(samFields) >= 10:
                # If the file is a bowtie file, get the bowtie-specific
                # fields in the SAM file
                for field in samFields[10:]:
                    if "AS:i" in field:
                        self.align_score = int(field.strip("AS:i"))
                    elif "XS:i" in field:
                        self.align_score_best = int(field.strip("XS:i"))
                    elif "YS:i" in field:
                        self.align_score_mate = int(field.strip("YS:i"))
                    elif "XN:i" in field:
                        self.num_ambig = int(field.strip("XN:i"))
                    elif "XM:i" in field:
                        self.mismatches = int(field.strip("XM:i"))
                    elif "XO:i" in field:
                        self.gap_opens = int(field.strip("XO:i"))
                    elif "XG:i" in field:
                        self.gap_ext = int(field.strip("XG:i"))
                    elif "NM:i" in field:
                        self.edit_distance = int(field.strip("NM:i"))
                    elif "YF:Z" in field:
                        self.why_filtered = int(field.strip("YF:Z"))
                    elif "YT:Z" in field:
                        self.yt = field.strip("YT:Z")
                    elif "MD:Z" in field:
                        self.mismatch_ref_bases = field.strip("MD:Z")

def decompose_sam_flag(flag):
    """
    Decompose SAM flag into its component parts

    Parameters:
        flag(int): Sam flag

    Returns:
        (list): Elements of sam flag
    """
    out_list = []
    flag_list = ["PAIRED", "PROPER_PAIR", "UNMAP", "MUNMAP",
                 "REVERSE", "MREVERSE", "READ1", "READ2",
                 "SECONDARY", "QCFAIL", "DUP", "SUPPLEMENTARY"]

    binary = str(f"{flag:b}"[::-1])
    for i in range(len(binary)):
        if binary[i] == "1":
            out_list.append(flag_list[i])
    return out_list
    
#File paths
SAM_PATH="../data/10-06-2021/"

# Get the lengths of all the chromosomes
chrom_lengths = {}
with open(SAM_PATH + "reference_genome/chrom_lengths.txt") as fn:
    lines = fn.readlines()
    for line in lines:
        line = line.split('\t')
        chrom_lengths[line[0]] = line[1]
            
###################################################################################

def not_in_list(not_list, input_list):
    flag = True
    for ele in not_list:
        if ele in input_list:
            return False
    return True

fn = open("multiple_depth_rep_results.txt","w")
depths = [100, 250, 1000, 2500, 3500, 5000, 7500, 10000]
for depth in depths:
    fn.write(">" + str(depth) + "\n")
    runs = 100
    MOCK_PATH = "../data/10-06-2021/reference_genome/"
    DATA_PATH = "../data/10-06-2021/"

    female_hit_counts = {"NC_000001.11":[], "NC_000002.12":[], "NC_000003.12":[], "NC_000004.12":[],
    "NC_000005.10":[],"NC_000006.12":[],"NC_000007.14":[],"NC_000008.11":[], "NC_000009.12":[],
    "NC_000010.11":[], "NC_000011.10":[], "NC_000012.12":[], "NC_000013.11":[], "NC_000014.9":[],
    "NC_000015.10":[], "NC_000016.10":[], "NC_000017.11":[], "NC_000018.10":[], "NC_000019.10":[],
    "NC_000020.11":[], "NC_000021.9":[], "NC_000022.11":[], "NC_000023.11":[]}

    # male_hit_counts = {"NC_000001.11":[], "NC_000002.12":[], "NC_000003.12":[], "NC_000004.12":[],
    # "NC_000005.10":[],"NC_000006.12":[],"NC_000007.14":[],"NC_000008.11":[], "NC_000009.12":[],
    # "NC_000010.11":[], "NC_000011.10":[], "NC_000012.12":[], "NC_000013.11":[], "NC_000014.9":[],
    # "NC_000015.10":[], "NC_000016.10":[], "NC_000017.11":[], "NC_000018.10":[], "NC_000019.10":[],
    # "NC_000020.11":[], "NC_000021.9":[], "NC_000022.11":[], "NC_000023.11":[]}

    for i in range(runs):
        # # Simulate paired-end reads from male genome
        # os.system("wgsim -N " + str(depth) + " " +  MOCK_PATH + 
        # "male_genome.fa.gz "  + MOCK_PATH + "reads/male_" + str(depth) + 
        # "_1.fa " + MOCK_PATH + "reads/male_"  + str(depth)  + "_2.fa > /dev/null 2>&1")
    
        # Simulate paired-end reads from female genome
        os.system("wgsim -N " + str(depth) + " " +  MOCK_PATH + 
        "female_genome.fa.gz "  + MOCK_PATH + "reads/female_" + str(depth) + 
        "_1.fa " + MOCK_PATH + "reads/female_" + str(depth)  + "_2.fa > /dev/null 2>&1")
    
        # #Align male reads to reference genome
        # os.system("bowtie2 --local -x " + DATA_PATH + "bowtie2_index/GRCh38_latest_genomic \
        # -1 " + DATA_PATH + "reference_genome/reads/male_" + str(depth) + "_1.fa \
        # -2 " + DATA_PATH + "reference_genome/reads/male_" + str(depth) + "_2.fa \
        # > ../data/10-06-2021/alignments/male_" + str(depth) + "_align.sam")
    
        #Align female reads to reference genome
        os.system("bowtie2 --local -x " + DATA_PATH + "bowtie2_index/GRCh38_latest_genomic \
        -1 " + DATA_PATH + "reference_genome/reads/female_" + str(depth) + "_1.fa \
        -2 " + DATA_PATH + "reference_genome/reads/female_" + str(depth) + "_2.fa \
        > ../data/10-06-2021/alignments/female_" + str(depth) + "_align.sam")
    
        # ################################################################################
    
    
    
        for rec in ParseSam(SAM_PATH + "alignments/female_" + str(depth) + "_align.sam"):
            if not_in_list(["SECONDARY", "UNMAP", "MUNMAP", "SUPPLEMENTARY"], decompose_sam_flag(rec.flag)) == True:
                if "NW" not in rec.rnam and "NT" not in rec.rnam and rec.rnam in female_hit_counts:
                    if len(female_hit_counts[rec.rnam]) == i:
                        female_hit_counts[rec.rnam].append(1)
                    else:
                        female_hit_counts[rec.rnam][-1] += 1
    
    
    
        # for rec in ParseSam(SAM_PATH + "alignments/male_" + str(depth) + "_align.sam"):
        #     if not_in_list(["SECONDARY", "UNMAP", "MUNMAP", "SUPPLEMENTARY"], decompose_sam_flag(rec.flag)) == True:
        #         if "NW" not in rec.rnam and "NT" not in rec.rnam and rec.rnam in male_hit_counts:
        #             if len(male_hit_counts[rec.rnam]) == i:
        #                 male_hit_counts[rec.rnam].append(1)
        #             else:
        #                 male_hit_counts[rec.rnam][-1] += 1
    

    for ele in female_hit_counts:
        for g in range(len(female_hit_counts[ele])):
            female_hit_counts[ele][g] = female_hit_counts[ele][g] / float(chrom_lengths[ele])
    
    
    # # for ele in male_hit_counts:
    #     for g in range(len(male_hit_counts[ele])):
    #         male_hit_counts[ele][g] = male_hit_counts[ele][g] / float(chrom_lengths[ele])
    

    these_ids = ["NC_000001.11", "NC_000002.12", "NC_000003.12", "NC_000004.12",
    "NC_000005.10","NC_000006.12","NC_000007.14","NC_000008.11", "NC_000009.12",
    "NC_000010.11", "NC_000011.10", "NC_000012.12", "NC_000013.11", "NC_000014.9",
    "NC_000015.10", "NC_000016.10", "NC_000017.11", "NC_000018.10", "NC_000019.10",
    "NC_000020.11", "NC_000021.9", "NC_000022.11", "NC_000023.11"]


    female_avg_list = []
    for ids in these_ids:
        for n in range(len(female_hit_counts[ids])):
            female_hit_counts[ids][n] = female_hit_counts["NC_000023.11"][n] / female_hit_counts[ids][n]
        female_avg_list.append(female_hit_counts[ids])
    fn.write(str(female_avg_list) + "\n")
    
    # male_avg_list = []
    # for ids in these_ids:
    #     for g in range(len(male_hit_counts[ids])):
    #         male_hit_counts[ids][g] = male_hit_counts["NC_000023.11"][g] / male_hit_counts[ids][g]
    #     male_avg_list.append(male_hit_counts[ids])
    # fn.write(str(male_avg_list) + "\n")

fn.close()
# plt.figure(0)   
# plt.boxplot(female_avg_list)
# plt.xticks([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
# 20, 21, 22, 23], these_ids)
# plt.ylim([0,2])
# plt.title("X:Autosomal Coverage Ratio - Female")
# plt.ylabel('Normalized Coverage')
# plt.xlabel('Chromosome ID')
# plt.xticks(rotation = -90)
# plt.savefig('../results/figures/x_to_autosome_female_zonkey_5000_100_reps.png', bbox_inches='tight')
###

# male_avg_list = []

# for ids in these_ids:
    # for g in range(len(male_hit_counts[ids])):
        # male_hit_counts[ids][g] = male_hit_counts["NC_000023.11"][g] / male_hit_counts[ids][g]
    # male_avg_list.append(male_hit_counts[ids])

# plt.figure(1)   
# plt.boxplot(male_avg_list)
# plt.xticks([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
# 20, 21, 22, 23], these_ids)
# plt.ylim([0,2])
# plt.title("X:Autosomal Coverage Ratio - Male")
# plt.ylabel('Normalized Coverage')
# plt.xlabel('Chromosome ID')
# plt.xticks(rotation = -90)
# plt.savefig('../results/figures/x_to_autosome_male_zonkey_5000_100_reps.png', bbox_inches='tight')