#!/bin/sh
#SBATCH --mail-user=jrbnewman@ufl.edu
#SBATCH --job-name=ea_refs
#SBATCH --account=concannon
#SBATCH --qos=concannon
#SBATCH --partition=hpg1-compute
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /ufrc/concannon/share/jnewman/scripts/SLURM_LOGS/ase_aln.%j.%A.%a.out
#SBATCH -t 6:00:00
#SBATCH --mem=12gb

module load bowtie/1.2.2
module load bwa/0.7.17
module load samtools/1.3.1

PROJ=/ufrc/concannon/share/jnewman
REFDIR=$PROJ/references/event_analysis

### Bowtie database for junctions

bowtie-build ${REFDIR}/hg19_aceview_junctions.fa  ${REFDIR}/hg19_aceview_junctions_BT1

samtools faidx ${REFDIR}/hg19_aceview_junctions.fa

### BWA-MEM database for genome

#bwa index -a bwtsw ${REFDIR}/hg19_all_chr_for_bedtools_plus_mito.fa

