#!/usr/bin/env bash
#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --job-name=deduper_pairedend
#SBATCH --output=deduper_pairedend.out
#SBATCH --error=deduper_pairedend.err
#SBATCH --time=0-12:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=pberry@uoregon.edu
#SBATCH --mail-type=ALL

conda deactivate
conda deactivate
conda deactivate
conda activate bgmp_py3
module load samtools/1.5

/usr/bin/time -v samtools sort /projects/bgmp/shared/deduper/paired_end.sam -o paired_end.sorted.sam
/usr/bin/time -v ./berry_deduper.py -f paired_end.sorted.sam \
-u STL96.txt \
-w \
-q \
-p
