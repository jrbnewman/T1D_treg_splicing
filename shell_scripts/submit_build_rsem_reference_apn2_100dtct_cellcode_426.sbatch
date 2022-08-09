#!/bin/sh
#SBATCH --mail-user=jrbnewman@ufl.edu
#SBATCH --account=concannon
#SBATCH --qos=concannon-b
#SBATCH --job-name=rsem_refs
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /ufrc/concannon/share/jnewman/scripts/SLURM_LOGS/out.%j.%A.%a.out
#SBATCH -t 6:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=12gb
#

module load R
module load bowtie/0.12.9
module load rsem/1.2.28

PROJ=/ufrc/concannon/share/jnewman
REF=$PROJ/references/event_analysis

    GENE2XS=hg19_aceview_transcripts_cellcode_426_100perc_apn2_gene2xs_v3.txt
    FASTA=hg19_aceview_transcripts_cellcode_426_100perc_apn2_v3.fasta
    RSEMREF=hg19_aceview_426_RSEM_v3


rsem-prepare-reference --bowtie -p 2 \
                       --transcript-to-gene-map $REF/${GENE2XS} \
                       $REF/${FASTA} \
                       $REF/${RSEMREF}

