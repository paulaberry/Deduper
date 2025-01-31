#!/usr/bin/env bash
#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --job-name=deduper_dataset2
#SBATCH --output=deduper_dataset2.out
#SBATCH --error=deduper_dataset2.err
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

/usr/bin/time -v samtools sort /projects/bgmp/shared/deduper/Dataset2.sam -o Dataset2.sorted.sam
/usr/bin/time -v ./berry_deduper.py -f Dataset2.sorted.sam \
-u STL96.txt \
-w \
