#!/bin/sh
#SBATCH --array=0-42
#SBATCH --verbose
#SBATCH --job-name=BEDgraph2CSV
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err
#SBATCH --time=20:00:00
#SBATCH --nodes=1
#SBATCH --mem=16GB


module load r/intel/3.3.2


# capture the output of the command line and store in an environmental variable
# OUTSIDE () STORES IT AS SEPARATE ELEMENTS
FILES=($(ls *.bedgraph))

# array mode
INPUT=${FILES[$SLURM_ARRAY_TASK_ID]}

OUTPUT=${INPUT}.MYC.bedgraph

# run script with 1 input
Rscript createInputMatricesChipSeq_argsVersion.R $INPUT