#!/usr/bin/env bash
#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --job-name=deduper_dataset1
#SBATCH --output=deduper_dataset1.out
#SBATCH --error=deduper_dataset1.err
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

/usr/bin/time -v samtools sort /projects/bgmp/shared/deduper/Dataset1.sam -o Dataset1.sorted.sam
/usr/bin/time -v ./paulaberry_deduper.py -f Dataset1.sorted.sam \
-u STL96.txt \
-w \
-q
