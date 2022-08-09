#!/bin/sh
#SBATCH --mail-user=jrbnewman@ufl.edu
#SBATCH --job-name=fusion_counts
#SBATCH --account=concannon
#SBATCH --qos=concannon
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /ufrc/concannon/share/jnewman/scripts/SLURM_LOGS/out.cc_fusions.%j.%A.%a.out
#SBATCH -t 8:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=12gb
#SBATCH --array=243
#

### 442

mkdir -p tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
export TMPDIR=$(pwd)/tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

module load python/2.7.6

    PROJ=/ufrc/concannon/share/jnewman

#### Make Output Directory
        #OUTPUT=$PROJ/coverage_counts_fusions_cc
        #OUTPUT=$PROJ/coverage_counts_fragments_cc
        OUTPUT=$PROJ/coverage_counts_introns_cc

        if [ ! -e $OUTPUT ]; then mkdir -p $OUTPUT; fi

## Design file
     DESIGN_FILE=$PROJ/design_files/case_control_sample_list.txt

     DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
     IFS=',' read -ra ARRAY <<< "$DESIGN"

     SAMPLE=${ARRAY[0]}

     NAME=${SAMPLE}

# Create LOG directory and start log
        LOGS=$OUTPUT/logs
        if [ ! -e $LOGS ]; then mkdir -p $LOGS; fi
        MYLOG=$LOGS/${FILE}.log
        printf "`date` $FILE SLURM_ID:$SLURM_ARRAY_JOB_ID HOSTNAME:$HOSTNAME \n\n" > "${MYLOG}"

#### COVERAGE COUNTS
    #BED=$PROJ/references/event_analysis/hg19_aceview_fusions_coverage
    #BED=$PROJ/references/event_analysis/hg19_aceview_exon_fragments_coverage
    BED=$PROJ/references/event_analysis/hg19_aceview_introns_from_fusions

    SAM=$PROJ/bwa_mem_aln_se_cc/bwa_split/${NAME}_unique
    MPILEUP=$PROJ/genome_mpileups_cc/${NAME}.mpileup

## List of chromosomes including X and Y.
CHROM=/ufrc/concannon/share/jnewman/design_files/chrom_list.txt

COUNT=0
for chrom in $(cat $CHROM)
do
    awk -v chr="$chrom" '{if ($3 == chr) print $0}' $SAM.sam > $TMPDIR/${NAME}_${chrom}.sam

    awk -v chr="$chrom" '{if ($1 == chr) print $0}' $BED.bed > $TMPDIR/hg19_fragments_${chrom}.bed

    echo "Starting Coverage Counts for $NAME `date`" >> "${MYLOG}"
    python /ufrc/concannon/share/jnewman/scripts/rpkm_calculate.py \
        -b $TMPDIR/hg19_fragments_${chrom}.bed \
        -m $MPILEUP \
        -s $TMPDIR/${NAME}_${chrom}.sam \
        -n ${NAME} \
        --cv \
        -g "${MYLOG}" \
        -o $TMPDIR/${NAME}_${chrom}.csv
    echo "Finished Coverage Counts for $NAME `date`" >> "${MYLOG}"

    if [ $COUNT == 0 ]
    then
        cat $TMPDIR/${NAME}_${chrom}.csv > $OUTPUT/${NAME}.csv
        COUNT=1
    else
        tail -n+2 $TMPDIR/${NAME}_${chrom}.csv >> $OUTPUT/${NAME}.csv
    fi
done


echo "Script complete [`date`]" >> "${MYLOG}"

rm $TMPDIR/*


