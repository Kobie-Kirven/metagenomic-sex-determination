import gzip
from Bio import SeqIO

SAM_PATH="../data/10-06-2021/reference_genome/"

chrom_lengths = {}
with gzip.open(SAM_PATH + "GRCh38_latest_genomic.fna.gz", "rt") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        if "NW" not in record.id and "NT" not in record.id:
            chrom_lengths[record.id] = len(str(record.seq))
            
with open(SAM_PATH + "chrom_lengths.txt", "w") as fn:
    for chrom in chrom_lengths:
        fn.write("{}\t{}\n".format(chrom, chrom_lengths[chrom]))