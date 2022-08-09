#!/bin/sh
#SBATCH --mail-user=jrbnewman@ufl.edu
#SBATCH --job-name=aln_genome_bwa
#SBATCH --account=concannon
#SBATCH --qos=concannon-b
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /ufrc/concannon/share/jnewman/scripts/SLURM_LOGS/out.aln_bwa_genome.%j.%A.%a.out
#SBATCH -t 12:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=48gb
#SBATCH --array=1-247%20
#

##442

mkdir -p tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
export TMPDIR=$(pwd)/tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

module load bwa/0.7.15
module load python/2.7.6
module load samtools/1.3.1

## Set directories
    PROJ=/ufrc/concannon/share/jnewman
    

    ORIG=$PROJ/fastq_by_sample
#    ORIG=$PROJ/aln_junctions_cc/unaln_reads
    
    OUTPUT=$PROJ/bwa_mem_aln_se_cc
    if [ ! -e $OUTPUT ]; then mkdir -p $OUTPUT; fi

    ## Output folder for split SAM files
    SPLITSAM=$OUTPUT/bwa_split
    if [ ! -e $SPLITSAM ]; then mkdir -p $SPLITSAM; fi

    ## Design file
    
    
 DESIGN_FILE=$PROJ/design_files/case_control_sample_list_425_426.txt
 DESIGN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DESIGN_FILE)
 IFS=',' read -ra ARRAY <<< "$DESIGN"

 SAMPLE=${ARRAY[0]}
 NAME=${SAMPLE}

## set references
    REF=$PROJ/references/event_analysis/hg19_all_chr_for_bedtools_plus_mito.fa

#### create bwa-mem reference

#    bwa index -p $REF -a bwtsw $FASTA

#### run bwa-mem

    bwa mem -t 2 $REF $ORIG/${NAME}.fastq > $OUTPUT/${NAME}.sam

#### parse bwa alignment file
    python $PROJ/scripts/BWASplitSAM_07mai.py \
           -fq1 $ORIG/${NAME}.fastq \
            -s $OUTPUT/${NAME}.sam \
            --outdir $SPLITSAM

### Clean up

rm $OUTPUT/${NAME}.sam
rm $SPLITSAM/${NAME}_*.fq
rm $SPLITSAM/${NAME}_ambiguous.sam
rm $SPLITSAM/${NAME}_chimeric.sam

