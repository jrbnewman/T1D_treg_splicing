#!/bin/sh
#SBATCH --mail-user=jrbnewman@ufl.edu
#SBATCH --job-name=run_miso_star
#SBATCH --account=concannon
#SBATCH --qos=concannon-b
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /blue/concannon/share/jnewman/bri_case_control/scripts/SLURM_LOGS/out.aln_bwa.%j.%A.%a.out
#SBATCH -t 6:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=30gb
#SBATCH --array=11-104

## 104
###1-247%30
#

## Modules

module load miso/0.5.3
module load samtools/1.3.1
module load bedtools
#module load bwa/0.7.15
module load python/2.7.6
module load star

## Set directories
    PROJ=/blue/concannon/share/jnewman
    ORIG=$PROJ/bri_case_control/fastq_by_sample

    OUTPUT=$PROJ/bri_case_control/star_aln_se_cc
    if [ ! -e $OUTPUT ]; then mkdir -p $OUTPUT; fi


 DESIGN_FILE=$PROJ/bri_case_control/design_files/case_control_sample_list_425_426.txt
 DESIGN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DESIGN_FILE)
 IFS=',' read -ra ARRAY <<< "$DESIGN"

 SAMPLE=${ARRAY[0]}
 NAME=${SAMPLE}

REFDIR=/blue/concannon/share/jnewman/references/hg19_aceview/STAR_index_aceview



NUMPROC=4

ROZ=$(pwd)/star_roz/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
#mkdir -p $ROZ

#Align to genome with STAR

STAR --runThreadN $NUMPROC \
     --genomeDir $REFDIR \
     --readFilesIn $PROJ/bri_case_control/fastq_unzip/${NAME}.fastq \
     --outSAMtype BAM SortedByCoordinate \
     --outFileNamePrefix ${OUTPUT}/${NAME} \
     --outTmpDir $ROZ \
     --alignSJDBoverhangMin 3 \
     --alignSJoverhangMin 8



### (3) Convert SAM to BAM, sort and index

    BAMDIR=$PROJ/bri_case_control/miso_bam_files/${NAME}
    if [ ! -e $BAMDIR ]; then mkdir -p $BAMDIR; fi

    ## Fasta reference
    FASTA=$PROJ/references/hg19_aceview/hg19_all_chr_for_bedtools_plus_mito.fa

    echo "Converting SAM to BAM"

    sam_to_bam --convert=$OUTPUT/${NAME}.sam ${BAMDIR} --ref=${FASTA}.fai

### RI

     echo "Indexing GFF file"

     GFF=$PROJ/bri_case_control/miso_DS_genes/commonshortest/RI.DS425.gff3
     GFFDIR=$PROJ/bri_case_control/miso_DS_genes/RI_gff_indexed
          if [ ! -e $GFFDIR ]; then mkdir -p $GFFDIR; fi

     index_gff --index ${GFF} ${GFFDIR}
     echo "Running MISO"

     MISO_OUTDIR=$PROJ/bri_case_control/miso_output_RI_STAR/${NAME}
     if [ ! -e ${MISO_OUTDIR} ]; then mkdir -p ${MISO_OUTDIR}; fi

samtools index ${OUTPUT}/${NAME}Aligned.sortedByCoord.out.bam
miso --run ${GFFDIR} ${OUTPUT}/${NAME}Aligned.sortedByCoord.out.bam --output-dir ${MISO_OUTDIR} --read-len 101

     echo "Summarizing MISO output"

     SUMOUT=$PROJ/bri_case_control/miso_summary_RI_STAR/${NAME}
     if [ ! -e ${SUMOUT} ]; then mkdir -p ${SUMOUT}; fi

    summarize_miso --summarize-samples ${MISO_OUTDIR} ${SUMOUT}
     echo "all done!"



### A3SS

     echo "Indexing GFF file"

     GFF=$PROJ/bri_case_control/miso_DS_genes/commonshortest/A3SS.DS425.gff3
     GFFDIR=$PROJ/bri_case_control/miso_DS_genes/A3SS_gff_indexed
          if [ ! -e $GFFDIR ]; then mkdir -p $GFFDIR; fi

     index_gff --index ${GFF} ${GFFDIR}

     echo "Running MISO"

     MISO_OUTDIR=$PROJ/bri_case_control/miso_output_A3SS_STAR/${NAME}
     if [ ! -e ${MISO_OUTDIR} ]; then mkdir -p ${MISO_OUTDIR}; fi

     #samtools index ${OUTPUT}/${NAME}Aligned.sortedByCoord.out.bam
     miso --run ${GFFDIR} ${OUTPUT}/${NAME}Aligned.sortedByCoord.out.bam --output-dir ${MISO_OUTDIR} --read-len 101


     echo "Summarizing MISO output"

     SUMOUT=$PROJ/bri_case_control/miso_summary_A3SS_STAR/${NAME}
     if [ ! -e ${SUMOUT} ]; then mkdir -p ${SUMOUT}; fi

    summarize_miso --summarize-samples ${MISO_OUTDIR} ${SUMOUT}


### A5SS

     echo "Indexing GFF file"

     GFF=$PROJ/bri_case_control/miso_DS_genes/commonshortest/A5SS.DS425.gff3
     GFFDIR=$PROJ/bri_case_control/miso_DS_genes/A5SS_gff_indexed
          if [ ! -e $GFFDIR ]; then mkdir -p $GFFDIR; fi

     index_gff --index ${GFF} ${GFFDIR}

     echo "Running MISO"

     MISO_OUTDIR=$PROJ/bri_case_control/miso_output_A5SS_STAR/${NAME}
     if [ ! -e ${MISO_OUTDIR} ]; then mkdir -p ${MISO_OUTDIR}; fi

     #samtools index ${OUTPUT}/${NAME}Aligned.sortedByCoord.out.bam
     miso --run ${GFFDIR} ${OUTPUT}/${NAME}Aligned.sortedByCoord.out.bam --output-dir ${MISO_OUTDIR} --read-len 101


     echo "Summarizing MISO output"

     SUMOUT=$PROJ/bri_case_control/miso_summary_A5SS_STAR/${NAME}
     if [ ! -e ${SUMOUT} ]; then mkdir -p ${SUMOUT}; fi

    summarize_miso --summarize-samples ${MISO_OUTDIR} ${SUMOUT}




### MXE

     echo "Indexing GFF file"

     GFF=$PROJ/bri_case_control/miso_DS_genes/commonshortest/MXE.DS425.gff3
     GFFDIR=$PROJ/bri_case_control/miso_DS_genes/MXE_gff_indexed
          if [ ! -e $GFFDIR ]; then mkdir -p $GFFDIR; fi

     index_gff --index ${GFF} ${GFFDIR}

     echo "Running MISO"

     MISO_OUTDIR=$PROJ/bri_case_control/miso_output_MXE_STAR/${NAME}
     if [ ! -e ${MISO_OUTDIR} ]; then mkdir -p ${MISO_OUTDIR}; fi

     #samtools index ${OUTPUT}/${NAME}Aligned.sortedByCoord.out.bam
     miso --run ${GFFDIR} ${OUTPUT}/${NAME}Aligned.sortedByCoord.out.bam --output-dir ${MISO_OUTDIR} --read-len 101


     echo "Summarizing MISO output"

     SUMOUT=$PROJ/bri_case_control/miso_summary_MXE_STAR/${NAME}
     if [ ! -e ${SUMOUT} ]; then mkdir -p ${SUMOUT}; fi

    summarize_miso --summarize-samples ${MISO_OUTDIR} ${SUMOUT}




### SE

     echo "Indexing GFF file"

     GFF=$PROJ/bri_case_control/miso_DS_genes/commonshortest/SE.DS425.gff3
     GFFDIR=$PROJ/bri_case_control/miso_DS_genes/SE_gff_indexed
          if [ ! -e $GFFDIR ]; then mkdir -p $GFFDIR; fi

     index_gff --index ${GFF} ${GFFDIR}

     echo "Running MISO"

     MISO_OUTDIR=$PROJ/bri_case_control/miso_output_SE_STAR/${NAME}
     if [ ! -e ${MISO_OUTDIR} ]; then mkdir -p ${MISO_OUTDIR}; fi

     #samtools index ${OUTPUT}/${NAME}Aligned.sortedByCoord.out.bam
     miso --run ${GFFDIR} ${OUTPUT}/${NAME}Aligned.sortedByCoord.out.bam --output-dir ${MISO_OUTDIR} --read-len 101


     echo "Summarizing MISO output"

     SUMOUT=$PROJ/bri_case_control/miso_summary_SE_STAR/${NAME}
     if [ ! -e ${SUMOUT} ]; then mkdir -p ${SUMOUT}; fi

    summarize_miso --summarize-samples ${MISO_OUTDIR} ${SUMOUT}
