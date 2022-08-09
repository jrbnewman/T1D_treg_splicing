#!/bin/sh
#SBATCH --mail-user=jrbnewman@ufl.edu
#SBATCH --account=concannon
#SBATCH --qos=concannon-b
#SBATCH --job-name=rsem_calc
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /ufrc/concannon/share/jnewman/scripts/SLURM_LOGS/out.%j.%A.%a.out
#SBATCH -t 24:00:00
#SBATCH --ntasks=1
#SBATCH --mem=16gb
#SBATCH --cpus-per-task=4
#SBATCH --array=182,329
#

#2-413

mkdir -p tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
export ROZ=$(pwd)/tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

# Load modules
module load R bowtie/0.12.9 rsem/1.2.28

## Set directories

# Project folder
PROJ=/ufrc/concannon/share/jnewman

# FASTQ folder
FASTQ=$PROJ/fastq_by_sample

# References folder
REFDIR=$PROJ/references/event_analysis

# Output
OUTDIR=$PROJ/rsem_output_bri_cc
  if [ ! -e $OUTDIR ]; then mkdir -p $OUTDIR; fi

## Design file

DESIGN_FILE=$PROJ/design_files/case_control_design_file_v2.csv
DESIGN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DESIGN_FILE)
IFS=',' read -ra ARRAY <<< "$DESIGN"

SAMPLE=${ARRAY[0]}
CELL=${ARRAY[1]}
STATUS=${ARRAY[2]}
AGE=${ARRAY[3]}
SEX=${ARRAY[4]}
PLATE=${ARRAY[5]}
WELL=${ARRAY[6]}
NAME=${ARRAY[7]}

## Set sample FQ
READS=$FASTQ/${SAMPLE}.fastq

## Set cell type-specific reference
REF=${REFDIR}/hg19_aceview_${CELL}_RSEM_v3

## Set sample-specific output name
OUTPUT=${OUTDIR}/hg19_aceview_apn2_100prc_cell_${CELL}_${SAMPLE}_

## Run RSEM

echo "Scripted started at:" $(date +"%m/%d/%Y %H:%M:%S $HOSTNAME")

rsem-calculate-expression -p 4 \
                          --temporary-folder ${ROZ} \
                          --time \
                          --calc-ci \
                          --estimate-rspd \
                          --fragment-length-mean 174 \
                          --fragment-length-sd 49 \
                          ${READS} ${REF} ${OUTPUT}

echo "Scripted ended at:" $(date +"%m/%d/%Y %H:%M:%S $HOSTNAME")

