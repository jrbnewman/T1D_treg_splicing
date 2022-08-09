#!/bin/sh
#SBATCH --mail-user=jrbnewman@ufl.edu
#SBATCH --job-name=mpileup_genome
#SBATCH --account=concannon
#SBATCH --qos=concannon
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /ufrc/concannon/share/jnewman/scripts/SLURM_LOGS/out.genome_mpileups.%j.%A.%a.out
#SBATCH -t 4:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8gb
#SBATCH --array=243


#51-100
#

##442


mkdir -p tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
export TMPDIR=$(pwd)/tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

module load samtools/1.3.1

#### Set Directories
    PROJ=/ufrc/concannon/share/jnewman
    INPUT=$PROJ/bwa_mem_aln_se_cc/bwa_split


## Design file


     DESIGN_FILE=$PROJ/design_files/case_control_sample_list.txt

     DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
     IFS=',' read -ra ARRAY <<< "$DESIGN"

     SAMPLE=${ARRAY[0]}

     NAME=${SAMPLE}

     FILE=${SAMPLE}



    # Create OUTPUT directory if needed.
        OUTPUT=$PROJ/genome_mpileups_cc
        if [ ! -e $OUTPUT ]; then mkdir -p $OUTPUT; fi

    # Create LOG directory and start log
        LOGS=$OUTPUT/logs
        if [ ! -e $LOGS ]; then mkdir -p $LOGS; fi
        MYLOG=$LOGS/${FILE}_mpileip.log
        printf "`date` $LINE $MATE PBS_ID:$PBS_ARRAYID HOSTNAME:$HOSTNAME \n\n" > $MYLOG

#### References
    REF=$PROJ/references/event_analysis/hg19_all_chr_for_bedtools_plus_mito.fa

#### Convert SAM to BAM -- Read 1
    SAM=$INPUT/${NAME}_unique.sam
    SAMO=$INPUT/${NAME}_oposite.sam
    SAMU=$INPUT/${NAME}_mapped.sam
    
    cat $SAMO $SAMU > $SAM
    
    BAM=$INPUT/temp2/${NAME}

if [ ! -e $INPUT/temp2 ]; then mkdir $INPUT/temp2; fi


    printf "<-------------------- Convert SAM to BAM -------------------->" >> $MYLOG
    echo `date`": Starting SAM to BAM conversion" >> $MYLOG
    samtools view -ut $REF.fai -o $BAM.bam $SAM 2>> $MYLOG
    samtools sort -T ${BAM}.tmp.sorted -o $BAM.sorted.bam $BAM.bam 2>> "${MYLOG}"
    samtools index $BAM.sorted.bam >> $MYLOG
    echo `date`": Finished SAM to BAM conversion" >> $MYLOG
#### Make mpielup
    PILEUP=$OUTPUT/${NAME}.mpileup

    printf "<-------------------- Convert BAM to MPILEUP -------------------->" >> $MYLOG
    echo `date`": Generating pileup" >> $MYLOG
    samtools mpileup -d 1000000000 -f $REF $BAM.sorted.bam  > $PILEUP 2>> $MYLOG

echo `date`": Finished Script complete" >> $MYLOG

rm $SAMU
rm $SAMO
rm $BAM.*

