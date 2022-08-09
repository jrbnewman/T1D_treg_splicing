#!/bin/sh
#SBATCH --mail-user=jrbnewman@ufl.edu
#SBATCH --job-name=run_miso_star
#SBATCH --account=concannon
#SBATCH --qos=concannon-b
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /blue/concannon/share/jnewman/bri_case_control/scripts/SLURM_LOGS/out.aln_bwa.%j.%A.%a.out
#SBATCH -t 6:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=30gb
#SBATCH --array=1-104

## 104
###1-247%30
#

## Modules

module load miso/0.5.3
module load samtools/1.3.1
module load bedtools
#module load bwa/0.7.15
module load python/2.7.6
module load star

## Set directories
    PROJ=/blue/concannon/share/jnewman
    ORIG=$PROJ/bri_case_control/fastq_by_sample

    OUTPUT=$PROJ/bri_case_control/star_aln_se_cc
    if [ ! -e $OUTPUT ]; then mkdir -p $OUTPUT; fi

 DESIGN_FILE=$PROJ/bri_case_control/design_files/case_control_sample_list_425_426.txt
 DESIGN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DESIGN_FILE)
 IFS=',' read -ra ARRAY <<< "$DESIGN"

 SAMPLE=${ARRAY[0]}
 NAME=${SAMPLE}


## set references
REFDIR=/blue/concannon/share/jnewman/references/hg19_aceview/STAR_index_aceview
#
NUMPROC=2

ROZ=$(pwd)/star_roz/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
#mkdir -p $ROZ

#Align to genome with STAR

STAR --runThreadN $NUMPROC \
     --genomeDir $REFDIR \
     --readFilesIn $PROJ/bri_case_control/fastq_unzip/${NAME}.fastq \
     --outSAMtype BAM SortedByCoordinate \
     --outFileNamePrefix ${OUTPUT}/${NAME} \
     --outTmpDir $ROZ \
     --alignSJDBoverhangMin 3 \
     --alignSJoverhangMin 8
