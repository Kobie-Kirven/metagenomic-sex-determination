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
                        self.why_filtered = field.strip("YF:Z")
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
    
    
import matplotlib
import matplotlib.pyplot as plt
from Bio import SeqIO
import gzip    
    
#File paths
SAM_PATH="../data/10-06-2021/"

# Get the lengths of all the chromosomes
chrom_lengths = {}
with gzip.open(SAM_PATH + "reference_genome/mm39.fa.masked.gz", "rt") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        if "_" not in record.id and "GL" not in record.id:
            chrom_lengths[record.id] = len(str(record.seq))
            
#####
hit_counts = {}

def not_in_list(not_list, input_list):
    flag = True
    for ele in not_list:
        if ele in input_list:
            return False
    return True

for rec in ParseSam(SAM_PATH + "reads/mouse_filtered.sam"):
    if not_in_list(["SECONDARY", "UNMAP", "MUNMAP", "SUPPLEMENTARY"], decompose_sam_flag(rec.flag)) == True:
        if "_" not in rec.rnam and "NT" not in rec.rnam:
            if rec.rnam in hit_counts:
                hit_counts[rec.rnam] += 1
            else:
                hit_counts[rec.rnam] = 1
print(hit_counts)           
            
avg_cov = {}
#Compute average coverage
for ele in hit_counts:
    average = hit_counts[ele] / chrom_lengths[ele]
    avg_cov[ele] = average
    


id_list = []
avg_list = []
for avg in avg_cov:
    if avg != "chrY":
        id_list.append(avg)
        avg_list.append(avg_cov["chrX"]/avg_cov[avg])
plt.figure(0)    
plt.scatter(id_list, avg_list)
plt.ylim([0,2])
plt.title("X:Autosomal Coverage Ratio Mouse")
plt.ylabel('Average Coverage')
plt.xlabel('Chromosome ID')
plt.xticks(rotation = -90)
plt.savefig('../results/figures/x_to_autosome_mouse.png', bbox_inches='tight')
