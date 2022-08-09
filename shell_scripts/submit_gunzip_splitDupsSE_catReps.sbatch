#!/bin/sh
#SBATCH --mail-user=jrbnewman@ufl.edu
#SBATCH --job-name=cc_dups
#SBATCH --mail-type=FAIL
#SBATCH --account=concannon
#SBATCH --qos=concannon
#SBATCH --no-requeue
#SBATCH -o /ufrc/concannon/share/jnewman/scripts/SLURM_LOGS/out.cc_dups_splitDupsSE_catReps.%j.%A.%a.out
#SBATCH -t 24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16gb
#SBATCH --array=25,36,49,53,54,59,71,100,141,160,264,353,358


### MAX: 442

mkdir -p roz/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
export TMPDIR=$(pwd)/roz/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

module load python/2.7.6

# Set directories
    PROJ=/ufrc/concannon/share/jnewman
    ORIG=$PROJ/original_data/PC2014
    FASTQ=$PROJ/fastq_by_sample
    if [ ! -e $FASTQ ]; then mkdir -p $FASTQ; fi


MCPYTHON=/ufrc/mcintyre/share/python.git

OUTPUT=$PROJ/fastq_split_dups
if [ ! -e $OUTPUT ]; then mkdir -p $OUTPUT; fi

OUTCOUNTS=$OUTPUT/duplicate_counts
if [ ! -e $OUTCOUNTS ]; then mkdir -p $OUTCOUNTS; fi

OUTLOGS=$OUTPUT/duplicate_counts/logs
if [ ! -e $OUTLOGS ]; then mkdir -p $OUTLOGS; fi

### Set FASTQ file

     DESIGN_FILE=$PROJ/design_files/case_control_sample_list.txt

     DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
     IFS=',' read -ra ARRAY <<< "$DESIGN"

     SAMPLE=${ARRAY[0]}
     NAME=${SAMPLE}

     ## If output FASTQ file already exists, then delete and we will remake
     if [ -f $FASTQ/${SAMPLE}.fastq ]; then
        rm $FASTQ/${SAMPLE}.fastq
     fi

     DESIGN_FILE_SUB=$PROJ/design_files/case_control_design_file_02jrbn.csv
     
     ## Iterate over all reps for sample and cat into one FASTQ file for alignments
     for ENTRY in {2..3526}
     do
          SUBDESIGN=$(cat ${DESIGN_FILE_SUB} | head -n ${ENTRY} | tail -n1 )
          IFS=',' read -ra SUBARRAY <<< "$SUBDESIGN"

          FLOWCELL=${SUBARRAY[0]}
          PROJECT=${SUBARRAY[1]}
          PLATE=${SUBARRAY[2]}
          WELL=${SUBARRAY[3]}
          INDEX=${SUBARRAY[4]}
          BARCODE=${SUBARRAY[5]}
          LANE=${SUBARRAY[6]}
          REP=${SUBARRAY[7]}
          NUM=${SUBARRAY[8]}

          SUBNAME=${PLATE}_${WELL}_${INDEX}_${BARCODE}_${LANE}_${REP}_${NUM}
          SUBSAMPLE=Sample_${PLATE}_${WELL}_${INDEX}
          FILEPATH=${ORIG}/${FLOWCELL}/${PROJECT}/Sample_${PLATE}_${WELL}_${INDEX}
      #    echo $SUBSAMPLE

          if [ ${SUBSAMPLE} == ${SAMPLE} ];
            then echo ${SUBNAME}
            echo $(less ${FILEPATH}/${SUBNAME}.fastq.gz | wc -l) 
            gunzip -c ${FILEPATH}/${SUBNAME}.fastq.gz > $FASTQ/${SAMPLE}.${SUBNAME}.tmp.fastq
            $MCPYTHON/fastqSplitDups.py -r1 $FASTQ/${SAMPLE}.${SUBNAME}.tmp.fastq --outdir $OUTPUT -o $OUTCOUNTS/${NAME}_duplicate_counts_summary.csv -t $OUTCOUNTS/${NAME}_001_duplicate_sequences.tsv -g $OUTLOGS/${NAME}_001_fastqSplitDups.log

            cat $OUTPUT/${SAMPLE}.${SUBNAME}.tmp_distinct.fq >> $FASTQ/${SAMPLE}.fastq
            rm $FASTQ/${SAMPLE}.${SUBNAME}.tmp.fastq
            rm $OUTPUT/${SAMPLE}.${SUBNAME}.tmp_*.fq
          fi
       done

