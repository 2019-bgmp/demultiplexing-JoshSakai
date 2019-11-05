#!/bin/bash

#SBATCH --account=bgmp          ### SLURM account which will be charged for the job
#SBATCH --partition=bgmp        ### Partition (like a queue in PBS)
#SBATCH --job-name=demultiplexing     ### Job Name
#SBATCH --output=demul.out         ### File in which to store job output
#SBATCH --error=demul.err          ### File in which to store job error messages
#SBATCH --time=24:00:00       ### Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=1               ### Node count required for the job
#SBATCH --ntasks-per-node=1     ### Nuber of tasks to be launched per Node
#SBATCH --cpus-per-task=1       ### Number of cpus (cores) per task

conda activate bgmp_py3

#./demultiplex.py -r1 test_R1.fastq.gz -r2 test_R4.fastq.gz -i1 test_R2.fastq.gz -i2 test_R2.fastq.gz -q 25 -b barcodes.tsv

/usr/bin/time -v ./demultiplex.py -r1 1294_S1_L008_R1_001.fastq.gz -r2 1294_S1_L008_R4_001.fastq.gz -i1 1294_S1_L008_R2_001.fastq.gz -i2 1294_S1_L008_R3_001.fastq.gz -q 25 -b barcodes.tsv
