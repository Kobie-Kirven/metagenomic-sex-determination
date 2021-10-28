###############################################################################
# This scrip will generate mock diploid genomes for both the male and 
# female sex.
# Male genome - creates a duplication of each of the chromosomes
# in the reference genome except for the X and Y chromosome
# Female genome - removes the Y chromosome then duplicates all of the
# chromosomes in the reference genome.
###############################################################################

#Author - Kobie Kirven
#Date - 10-06-2021

#imports
from Bio import SeqIO
import os
import gzip

#File paths
ref_path = "../data/10-06-2021/reference_genome/"
ref_name = "GRCh38_latest_genomic.fna.gz"
male_output = "male_genome.fa"
female_output = "female_genome.fa"


def duplicate_chroms_female(reference, exclude_list, output_file):
    """
    Duplicate the chromosomes that are not in exclude_list
    
    Parameters:
        reference(str): Path to reference genome in FASTA format
        exclude_list(list): Chromosomes or keywords in FASTA IDs to exclude
        output_file(str): Path to output_file
    """
    fn = open(output_file, "w")
    with gzip.open(reference, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            for ele in exclude_list:
                if not ele in str(record.description):
                    #Duplicate each chromosome
                    for i in range(2):
                        fn.write(">" + str(record.description) + "\n")
                        fn.write(str(record.seq) + "\n")
    fn.close()

def ele_not_in_str(ele_list, description):
    flag = False
    for ele in ele_list:
        if ele in description:
            flag = True
    return flag

def duplicate_chroms_male(reference, exclude_list, output_file):
    """
    Duplicate the chromosomes that are not in exclude_list
    
    Parameters:
        reference(str): Path to reference genome in FASTA format
        exclude_list(list): Chromosomes or keywords in FASTA IDs to exclude
        output_file(str): Path to output_file
    """
    fn = open(output_file, "w")
    with gzip.open(reference, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if ele_not_in_str(exclude_list, str(record.description)) == True:
                fn.write(">" + str(record.description) + "\n")
                fn.write(str(record.seq) + "\n")
            else:
                for i in range(2):
                    #Duplicate the autosomes
                    fn.write(">" + str(record.description) + "\n")
                    fn.write(str(record.seq) + "\n")
    fn.close()


####################################
# Generate female genome
####################################

# Generate the FASTA file with female genome
duplicate_chroms_female((ref_path + ref_name), ["chromosome Y"], (ref_path + female_output))

# Gzip the file to save space
os.system("gzip {}{}".format(ref_path, female_output))


####################################
# Generate male genome
####################################
    
# Generate the FASTA file with the male genome    
duplicate_chroms_male((ref_path + ref_name), ["chromosome Y","chromosome X"], (ref_path + male_output))

# Gzip the file to save space
os.system("gzip {}{}".format(ref_path, male_output))


