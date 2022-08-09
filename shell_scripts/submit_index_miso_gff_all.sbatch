#!/bin/sh
#SBATCH --mail-user=jrbnewman@ufl.edu
#SBATCH --job-name=miso_star
#SBATCH --account=concannon
#SBATCH --qos=concannon-b
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /blue/concannon/share/jnewman/bri_case_control/scripts/SLURM_LOGS/out.miso.%j.%A.%a.out
#SBATCH -t 4:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=12gb


module load miso/0.5.3
module load samtools/1.3.1
module load bedtools
#module load bwa/0.7.15
module load python/2.7.6
module load star


## Set directories
    PROJ=/blue/concannon/share/jnewman

     echo "Indexing GFF file"

     GFF=$PROJ/bri_case_control/miso_treg_genes/commonshortest/A3SS.treg.gff3
     GFFDIR=$PROJ/bri_case_control/miso_treg_genes/A3SS_gff_indexed
          if [ ! -e $GFFDIR ]; then mkdir -p $GFFDIR; fi

     index_gff --index ${GFF} ${GFFDIR}

     GFF=$PROJ/bri_case_control/miso_treg_genes/commonshortest/A5SS.treg.gff3
     GFFDIR=$PROJ/bri_case_control/miso_treg_genes/A5SS_gff_indexed
          if [ ! -e $GFFDIR ]; then mkdir -p $GFFDIR; fi

     index_gff --index ${GFF} ${GFFDIR}

     GFF=$PROJ/bri_case_control/miso_treg_genes/commonshortest/MXE.treg.gff3
     GFFDIR=$PROJ/bri_case_control/miso_treg_genes/MXE_gff_indexed
          if [ ! -e $GFFDIR ]; then mkdir -p $GFFDIR; fi
     index_gff --index ${GFF} ${GFFDIR}


     GFF=$PROJ/bri_case_control/miso_treg_genes/commonshortest/SE.treg.gff3
     GFFDIR=$PROJ/bri_case_control/miso_treg_genes/SE_gff_indexed
          if [ ! -e $GFFDIR ]; then mkdir -p $GFFDIR; fi

     index_gff --index ${GFF} ${GFFDIR}


     GFF=$PROJ/bri_case_control/miso_treg_genes/commonshortest/RI.treg.gff3
     GFFDIR=$PROJ/bri_case_control/miso_treg_genes/RI_gff_indexed
          if [ ! -e $GFFDIR ]; then mkdir -p $GFFDIR; fi

     index_gff --index ${GFF} ${GFFDIR}
