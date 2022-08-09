#!/bin/sh
#SBATCH --mail-user=jrbnewman@ufl.edu
#SBATCH --account=concannon
#SBATCH --qos=concannon
#SBATCH --job-name=star_refs
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /blue/concannon/share/jnewman/out.%j.%A.%a.out
#SBATCH -t 8:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=20gb

module load star

PROJ=/blue/concannon/share/jnewman/references/hg19_aceview
BRI=/blue/concannon/share/jnewman/bri_case_control
NUMPROC=4

REFDIR=$PROJ/STAR_index_aceview_425_all
if [ ! -e $REFDIR ]; then mkdir -p $REFDIR; fi

FASTA=$PROJ/hg19_all_chr_for_bedtools_plus_mito.fa
GTF=$BRI/AceView.subset.GRCh37.hg19.gtf
JUNCLENGTH=101


STAR --runThreadN ${NUMPROC} \
     --runMode genomeGenerate \
     --genomeDir ${REFDIR} \
     --genomeFastaFiles ${FASTA} \
     --sjdbGTFfile ${GTF} \
     --sjdbGTFtagExonParentTranscript transcript_id \
     --sjdbOverhang ${JUNCLENGTH}