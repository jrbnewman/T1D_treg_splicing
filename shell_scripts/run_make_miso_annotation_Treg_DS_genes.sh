#!/bin/bash


## subset GTF file

PROJ=/TB14/TB14/sandbox/miso_sandbox


GTFIN=$PROJ/gtf/AceView.sorted.GRCh37.hg19.gtf
GTFOUT=$PROJ/gtf/AceView.subset.GRCh37.hg19.gtf

XSLIST=$PROJ/xscript_list.txt

for XS in $(cat $XSLIST); do
           grep $XS $GTFIN >> $GTFOUT
   done



## Create a MISO-like annotation for the PacBio transcripts

PROJ=/TB14/TB14/sandbox/miso_sandbox

# Output directory
ANNOTOUT=$PROJ/miso_treg_genes
if [ ! -e $ANNOTOUT ]; then mkdir -p $ANNOTOUT ; fi

# Create a directory for creating UCSC-style table
# Note that UCSC tables are essentially genePred files, so we generate one

TABLEDIR=$PROJ/genepred_treg_genes
if [ ! -e $TABLEDIR ]; then mkdir -p $TABLEDIR ; fi

# Location of GTF
GTF=$PROJ/gtf/AceView.subset.GRCh37.hg19.gtf


## Create genePred/UCSC-tyle table
## Note this MUST be called knownGene.txt, otherwise the annotation script doesn't work
$PROJ/scripts/gtfToGenePred -genePredExt ${GTF} ${TABLEDIR}/knownGene.txt

# Create MISO-like annotation
python2.7 $PROJ/scripts/rnaseqlib-clip/rnaseqlib/gff/gff_make_annotation.py \
                  $TABLEDIR $ANNOTOUT --genome-label treg


