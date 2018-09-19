#!/bin/sh
#SBATCH --array=0-
#SBATCH --verbose
#SBATCH --job-name=BAM2Bedgraph
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err
#SBATCH --time=20:00:00
#SBATCH --nodes=1
#SBATCH --mem=16GB

# ARRAY IS 0 TO 3 BECAUSE THERE ARE 4 files



module load deeptools/intel/2.4.2
module load pysam/intel/0.10.0
module load pybigwig/intel/0.3.3
module load samtools/intel/1.3.1
module load bedtools/intel/2.26.0


# capture the output of the command line and store in an environmental variable
# OUTSIDE () STORES IT AS SEPARATE ELEMENTS
FILES=($(ls *.bam))

# array mode
INPUT=${FILES[$SLURM_ARRAY_TASK_ID]}

BAMSUBSET=${FILES}.subset.bam

OUTPUT=${FILES}.MYC.bedgraph



# need to index the bam file first!
samtools index $INPUT

# take subset to genomic position, then take per bp coverage through the subsetted bam file -- much faster
# subset full bam to only desired coordinates
samtools view -b CBX2-human.A549.ENCFF573RUL.bam "chr8:125000000-135000000" | bedtools genomecov -ibam stdin -g /scratch/cgsb/sanjana/commons/hg19_chrom_sizes.genome -d > test.d.bedgraph 


# optionally divide the raw read counts in this output file by the 
awk '{printf "%s %.2f\n",$1,$2/1024/1024}' $OUTPUT
# change the .2f to gain precision in decimal points



# take the per bp coverage through the full bam file, then subset to genomic position
# find coverage per bp
# bedtools genomecov -ibam $INPUT -g /scratch/cgsb/sanjana/commons/hg19_chrom_sizes.genome -d | awk '$1 == "chr8" && $2 > 126000000 && $2 < 133000000 {print $1,$2,$3}' > $OUTPUT

# bedtools genomecov -ibam $FILES -g /scratch/cgsb/sanjana/commons/hg19_chrom_sizes.genome -bga | awk '$1 == "chr8" && $2 > 126000000 && $2 < 133000000 {print $1,$2,$3,$4}' > $OUTPUT