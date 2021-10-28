## bin/ - Scripts for data analysis:

contents:
- download_reference_genome.sh - Download the human reference genome from NCBI (10-06-2021)
- build_bowtie2_index.sh - Build bowtie2 index for human reference genome (10-06-2021)
- make_simulated_chromosomes.py - Create mock male and female genomes (10-07-2021)
- simulate_shotgun_reads.sh - Simulate 10,000 paired-end reads from male and female mock genomes respectively (10-07-2021)
- calculate_coverage_depth.py - Calculate coverage by summing coverage for each 
base from the "samtools depth" output (10-07-2021)
- wgsim_validation.sh - Generate 10 reads from male and female mock genomes respectively using wgsim (10-07-2021)
- try_zonkey.py - Test the method for sex determination from the Zonkey pipeline (10-07-2021)
- align_with_bowtie2.sh - Align reads and generate coverage files using samtools depth (10-08-2021)
- read_depth_pipeline.sh - Use the Zonkey method for simulations at various sequencing depths (10-08-2021)
- read_dp.py - Calculates X:Autosomal ratios from a bam file (used in read_depth_pipeline.sh) (10-08-2021)
- get_chrom_lengths.py - Calculate scaffold lengths from a reference genome
- mouse_genome.sh -  Download mouse genome and build bowtie2 index (10-14-2021)
- scims_cluster.pbs - Submit a scims run from one sample to run on Roar (10-14-2021)
- mouse_test.py - Test the Zonkey method on mouse cecal metagenomic data
- pbs_script.pbs - Submission script for "multiple_depth.py"
- multiple_depth.py - Test the Zonkey method at multiple sequencing depths (10-14-2021)
- birds/ - Scripts for analyzing metagenomes from birds 
- twinsuk/ - Scripts for analyzing metagenomes from TwinsUK study