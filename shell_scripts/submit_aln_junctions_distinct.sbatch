#!/bin/sh
#SBATCH --mail-user=jrbnewman@ufl.edu
#SBATCH --job-name=aln_splicing
#SBATCH --account=concannon
#SBATCH --qos=concannon
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /ufrc/concannon/share/jnewman/scripts/SLURM_LOGS/out.aln_splicing.%j.%A.%a.out
#SBATCH -t 8:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=12gb
#SBATCH --array=243
#

##442

mkdir -p tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
export TMPDIR=$(pwd)/tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}


module load bowtie/1.2.2
module load python/2.7.6

## Store number of processors that will be used, i.e. the number of files we will split into to run LAST on
    NUMPROCS=4

## Set Directories
    PROJ=/ufrc/concannon/share/jnewman

## Data
    ORIG=$PROJ/fastq_by_sample

## Get info from design file. Need a design file with unaligned reads from step 1 (aligning to junctions)

 DESIGN_FILE=$PROJ/design_files/case_control_sample_list.txt
 DESIGN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DESIGN_FILE)
 IFS=',' read -ra ARRAY <<< "$DESIGN"

 SAMPLE=${ARRAY[0]}
 NAME=${SAMPLE}

    
## Set Different References
    REF=$PROJ/references/event_analysis/hg19_aceview_junctions_BT1

## Create Necessary Folders
    OUTPUT=$PROJ/aln_junctions_cc
    if [ ! -e $OUTPUT ] ; then mkdir -p $OUTPUT; fi

    # Create JOB LOG directory and start log
        LOGS=$OUTPUT/job_logs 
        ALN_ERROR_LOG=$LOGS/size_errors.txt
        if [ ! -d $LOGS ]; then mkdir -p $LOGS; fi

        MYLOG=$LOGS/${NAME}.log
        printf "`date` $NAME PBS_ARRAYID:$PBS_ARRAYID HOSTNAME:$HOSTNAME\n" > $MYLOG
    # Create ALN LOG directory
        ALNLOGS=$OUTPUT/aln_logs 
        if [ ! -d $ALNLOGS ]; then mkdir -p $ALNLOGS; fi

    # Create UNALN READ directory
        UNALNDIR=$OUTPUT/unaln_reads
        if [ ! -d $UNALNDIR ]; then mkdir -p $UNALNDIR; fi

### FUNCTIONS FOR ALIGNMENT PIPELINE ###
# I have created a separate file that hold the alignment functions
source $PROJ/scripts/alignment_functions.sh

## Start Alignment Pipeline
    printf "<------------------- STARTING SE alignment process for $NAME [`date`] ------------------->\n" >> $MYLOG
    READS=$ORIG/${NAME}.fastq

    qual=`python /ufrc/mcintyre/share/python.git/identify_quality.py -i $READS`
    if [ $qual == 'phred64' ]; 
    then 
        # set to old illumina quality scores phred64/solexa 1.3
        btqual='--phred64-quals'
        lastqual='3'
    else
        # change to sanger format which is what all new illumina data is
        btqual='--phred33-quals'
        lastqual='1'
    fi

    bowtie_se_uniq

## Combine all Sam files 
    echo "START Combine SAM files">>$MYLOG
    cat *.sam >$OUTPUT/${NAME}.sam 2>>$MYLOG
    echo "FINISH Combining SAM files [`date`]" >>$MYLOG

## Combine all Unaln FASTQ 
   echo "START Combine Unaln FQ">>$MYLOG
   cat *_unaln_bt.fq >$UNALNDIR/${NAME}_unaln.fq 2>>$MYLOG
   echo "FINISH Combining SAM files [`date`]" >>$MYLOG

echo "Script complete, [`date`]" >>$MYLOG

