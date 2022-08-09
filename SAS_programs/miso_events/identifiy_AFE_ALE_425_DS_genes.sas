/* Import MISO annotations and assign to events */

ods listing; ods html close;
libname cc '!PATCON/case_control/sas_data/fixed';
libname hg19 '!PATCON/useful_human_data/aceview_hg19/sas_data';
libname cclocal '/TB14/TB14/immigrans/store/cc_sandbox/sas_data';


proc datasets lib=work kill noprint;
run;
quit;



proc import datafile="!PATCON/case_control/references/hg19_aceview/annotations/hg19_aceview_exon_annotations.csv"
    out=exon_annot dbms=csv replace;
    guessingrows=max;
run;

proc import datafile="!PATCON/case_control/references/hg19_aceview/annotations/hg19_aceview_fusion_annotations.csv"
    out=fusion_annot dbms=csv replace;
    guessingrows=max;
run;




/* Compare first and last fusions between transcripts of the same gene (DS genes) -- AFE/ALE */

data ds_genes;
   set cc.m37_anova_fdr_ds_case4;
   where cell_Type=425 and flag_fdr05=1;
   keep gene_id;
run;


data xs2keep;
   set cc.flag_xscript_tpm_gt0_ge05_425;
   where flag_425_CTL_tpm_gt0=1 or flag_425_T1D_tpm_gt0=1;
   keep transcript_id;
run;

data xs2gene;
  set hg19.hg19_xscript2gene;
run;

proc sort data=xs2keep;
  by transcript_id;
proc sort data=xs2gene;
  by transcript_id;
run;

data xs2keep2;
  merge xs2gene (in=in1) xs2keep (in=in2);
  by transcript_id;
  if in1 and in2;
run;

proc sort data=xs2keep2;
   by gene_id;
proc sort data=ds_genes;
    by gene_id;
run;

data ds_genes2xs;
  merge ds_genes (in=in1) xs2keep2 (in=in2);
  by gene_id;
  if in1 and in2;
run;

/* Get first and last exons per transcript */

data exon_annot2;
  set exon_annot;
  length transcript_id2 $100.;
  do i = 1 by 1 while(scan(transcript_id,i,"|") ^="");
     transcript_id2=compress(scan(transcript_id,i,"|"));
     output;
     end;
   keep gene_id exon_id start stop strand transcript_id2;
   rename transcript_id2=transcript_id;
run;

proc sort data=exon_annot2;
  by gene_id transcript_id;
proc sort data=ds_genes2xs;
  by gene_id transcript_id;
run;

data exon_annot3;
    merge exon_annot2 (in=in1) ds_genes2xs (in=in2);
    by gene_id transcript_id;
    if in1 and in2;
run;



proc sort data=exon_annot3;
  by gene_id transcript_id start stop;
run;

data first_exon last_exon non_first_exon non_last_exon;
   set exon_annot3;
   by gene_id transcript_id;
   if first.transcript_id then output first_exon;
   else output non_first_exon;
   if last.transcript_id then output last_exon;
   else output non_last_exon;
run;

data fus2exon;
   set fusion_annot;
   length exon_id2 $100.;
  do i = 1 by 1 while(scan(exon_id,i,"|") ^="");
     exon_id2=compress(scan(exon_id,i,"|"));
     output;
     end;
   keep fusion_id exon_id2;
   rename exon_id2=exon_id;
run;

proc sort data=fus2exon;
  by exon_id;
proc sort data=first_exon;
   by exon_id;
proc sort data=last_exon;
   by exon_id;
proc sort data=non_first_exon;
   by exon_id;
proc sort data=non_last_exon;
   by exon_id;
run;


data first_exon2fus;
   merge first_exon (in=in1) fus2exon (in=in2);
   by exon_id;
   if in1 and in2;
run;

data last_exon2fus;
   merge last_exon (in=in1) fus2exon (in=in2);
   by exon_id;
   if in1 and in2;
run;


data non_first_exon2fus;
   merge non_first_exon (in=in1) fus2exon (in=in2);
   by exon_id;
   if in1 and in2;
run;

data non_last_exon2fus;
   merge non_last_exon (in=in1) fus2exon (in=in2);
   by exon_id;
   if in1 and in2;
run;

proc freq data=first_exon2fus noprint;
  tables gene_id*fusion_id / out=first_fus_count;
proc freq data=last_exon2fus noprint;
  tables gene_id*fusion_id / out=last_fus_count;
run;


/*  remove internal exonic regions that are also AFE/ALEs: these are too hard to resolve as the inclusion counts are
    dependent on exclusion counts 


   get list of first and last fusions for each gene
   eliminate first/last fusion if it's not first/last of another transcript (i.e. it's internal)
*/


data alt_first_fus;
   set first_exon2fus;
   keep gene_id fusion_id ;
run;

data non_first_fus;
   set non_first_exon2fus;
   keep gene_id fusion_id ;
run;

proc sort data=alt_first_fus nodup;
  by gene_id fusion_id;
proc sort data=non_first_fus nodup;
  by gene_id fusion_id;
run;

data alt_first_fus_set;
  merge alt_first_fus (in=in1) non_First_fus (in=in2);
  by gene_id fusion_id;
  if in2 then delete;
run;

data alt_last_fus;
   set last_exon2fus;
   keep gene_id fusion_id ;
run;

data non_last_fus;
   set non_last_exon2fus;
   keep gene_id fusion_id ;
run;

proc sort data=alt_last_fus nodup;
  by gene_id fusion_id;
proc sort data=non_last_fus nodup;
  by gene_id fusion_id;
run;

data alt_last_fus_set;
  merge alt_last_fus (in=in1) non_last_fus (in=in2);
  by gene_id fusion_id;
  if in2 then delete;
run;



proc freq data=alt_first_fus_set noprint;
   tables gene_id / out=num_first_fus;
run;

proc freq data=alt_last_fus_set noprint;
   tables gene_id / out=num_last_fus;
run;


data alt_first_fus_set2keep;
    set num_first_fus;
    where count > 1;
    keep gene_id;
run;

data alt_last_fus_set2keep;
    set num_last_fus;
    where count > 1;
    keep gene_id;
run;


proc sort data=alt_first_fus_set2keep;
   by gene_id;
proc sort data=alt_last_fus_set2keep;
   by gene_id;
proc sort data=alt_first_fus_set;
   by gene_id;
proc sort data=alt_last_fus_set;
   by gene_id;
run;

data alt_first_fus_set2;
  merge alt_first_fus_set (in=in1) alt_first_fus_set2keep (in=in2);
  by gene_id;
  if in1 and in2;
run;

data alt_last_fus_set2;
  merge alt_last_fus_set (in=in1) alt_last_fus_set2keep (in=in2);
  by gene_id;
  if in1 and in2;
run;

data gene2strand;
    set exon_annot;
    keep gene_id strand;
run;

data fus_start;
   set fusion_annot;
   keep fusion_id fusion_start;
run;

proc sort data=gene2strand nodup;
   by gene_id;
proc sort data=alt_first_fus_set2;
  by gene_id;
proc sort data=alt_last_fus_set2;
  by gene_id;
run;

data alt_first_fus_set3;
   merge alt_first_fus_set2 (in=in1) gene2strand (in=in2);
   by gene_id;
   if in1 and in2;
run;

data alt_last_fus_set3;
   merge alt_last_fus_set2 (in=in1) gene2strand (in=in2);
   by gene_id;
   if in1 and in2;
run;

proc sort data=alt_first_fus_set3;
   by fusion_id;
proc sort data=alt_last_fus_set3;
   by fusion_id;
proc sort data=fus_start;
   by fusion_id;
run;



data alt_first_fus_set4;
   merge alt_first_fus_set3 (in=in1) fus_start (in=in2);
   by fusion_id;
   if in1 and in2;
run;

data alt_last_fus_set4;
   merge alt_last_fus_set3 (in=in1) fus_start (in=in2);
   by fusion_id;
   if in1 and in2;
run;

data alt_fusion_set_all;
    set alt_first_fus_set4 (in=in1) alt_last_fus_set4 (in=in2) ;
    if in1 and strand="+" then flag_first_last="F";
    else if in1 and strand="-" then flag_first_last="L";
    else if in2 and strand="+" then flag_first_last="L";
    else if in2 and strand="-" then flag_first_last="F";
run;

proc sort data=alt_fusion_set_all;
   by gene_id strand flag_first_last fusion_start;
run;

data ref_fusion alt_fusion;
  set alt_fusion_set_all;
  by gene_id strand flag_first_last;
  if first.flag_first_last then output ref_fusion;
  else output alt_fusion;
run;



proc import datafile="/TB14/TB14/immigrans/store/cc_sandbox/event_analysis_output/hpc/bri_cc_fusion_counts_wide.tsv"
     out=fus_counts_sbys dbms=tab replace;
     guessingrows=max;
run;

proc transpose data=fus_counts_sbys out=fus_counts;
   by event_id;
   var Sample_: ;
run;

data fus_counts2;
  set fus_counts;
   if _NAME_="Sample_P1_C1_97288" then delete;
   if _NAME_="Sample_P2_H2_11333425" then delete;
   if _NAME_="Sample_P3_H8_22378426" then delete;
   if _NAME_="Sample_P4_E5_81333426" then delete;
   if _NAME_="Sample_P4_E7_1641425" then delete;
   if _NAME_="Sample_P5_G7_97721127" then delete;
   if _NAME_="Sample_P6_A12_49925400" then delete;
   if _NAME_="Sample_P6_A9_28992400" then delete;
   if _NAME_="Sample_P6_B11_65680400" then delete;
   if _NAME_="Sample_P6_B12_29014400" then delete;
   if _NAME_="Sample_P6_B7_29941400" then delete;
   if _NAME_="Sample_P6_C11_47601400" then delete;
   if _NAME_="Sample_P6_C12_97721400" then delete;
   if _NAME_="Sample_P6_D3_24509425" then delete;
   if _NAME_="Sample_P6_E2_48902425" then delete;
   if _NAME_="Sample_P6_F11_8913040" then delete;
   if _NAME_="Sample_P6_B7_29941400" then delete;
   if _NAME_="Sample_P3_G6_85413426" then delete;
   if _NAME_="Sample_P3_H4_72730426" then delete;
   if _NAME_="Sample_P6_H2_25560426" then delete;
   if _NAME_="Sample_P1_C1_97288" then delete;
   if _NAME_="Sample_P2_E4_21239425" then delete;

  rename event_id=fusion_id _NAME_=sample_id col1=apn;
run;

data samples;
  set cc.design_file_new;
  keep sample_id ;
run;



proc sort data=fus_counts2;
  by sample_id;
proc sort data=samples;
  by sample_id;
run;

data counts2;
  merge samples (in=in1) fus_counts2 (in=in2);
  by sample_id;
  if in1 and in2;
run;


data counts3;
  set counts2;
  if apn > 0;
run;



/* quartiles */

proc sort data=counts3;
   by sample_id;
proc means data=counts3 noprint;
   by sample_id;
   var apn;
   output out=quartiles_apn0 q3=q3;
   run;


proc means data=quartiles_apn0 noprint;
            var q3;
            output out=q3_q3 q3=q3;
            run;

data _null_;
   set q3_q3;
   call symput('q3q3',q3);
run;

%put &q3q3.;

/*
Analysis Variable : q3

            Upper
         Quartile
     ------------
     36.264232673
     ------------
*/

    proc sort data=counts2;
        by sample_id;        run;

    proc sort data=quartiles_apn0;
        by sample_id;
        run;

/* merge and oops test */


    data counts_w_q3 oops_apn0;
        merge counts2 (in=in1) quartiles_apn0 (in=in2);
        by sample_id;
        if in1 and in2 then output counts_w_q3;
        else output oops_apn0; 
        run;

data fusion_q3_norm_counts;
        set counts_w_q3;
        log_apn=log(apn + 1);
        * Q3 normalization;
        q3_apn=(apn/q3) * &q3q3.;
        q3_ff=&q3q3./q3;
        log_q3_apn=log(q3_apn + 1);
run;


data ref_fusion_first ref_fusion_last;
  set ref_fusion;
  if flag_first_last="F" then output ref_fusion_First;
  else if flag_first_last="L" then output ref_fusion_last;
run;

data alt_fusion_first alt_fusion_last;
  set alt_fusion;
  if flag_first_last="F" then output alt_fusion_First;
  else if flag_first_last="L" then output alt_fusion_last;
run;

 


proc sort data=fusion_q3_norm_counts;
   by fusion_id ;
proc sort data=ref_fusion_first;
   by fusion_id ;
proc sort data=ref_fusion_last;
   by fusion_id ;
proc sort data=alt_fusion_first;
   by fusion_id ;
proc sort data=alt_fusion_last;
   by fusion_id ;
run;

data ref_fusion_first_counts; 
   merge ref_fusion_first (in=in1) fusion_q3_norm_counts (in=in2);
  by fusion_id;
  if in1 and in2;
run;

data alt_fusion_first_counts; 
   merge alt_fusion_first (in=in1) fusion_q3_norm_counts (in=in2);
  by fusion_id;
  if in1 and in2;
run;
data ref_fusion_last_counts; 
   merge ref_fusion_last (in=in1) fusion_q3_norm_counts (in=in2);
  by fusion_id;
  if in1 and in2;
run;

data alt_fusion_last_counts; 
   merge alt_fusion_last (in=in1) fusion_q3_norm_counts (in=in2);
  by fusion_id;
  if in1 and in2;
run;

data ref_fusion_first_counts2;
   set ref_fusion_first_counts;
   keep gene_id fusion_id strand sample_id q3_apn;
   rename fusion_id=ref_fusion_id q3_apn=ref_q3_apn;
run;

data alt_fusion_first_counts2;
   set alt_fusion_first_counts;
   keep gene_id fusion_id strand sample_id q3_apn;
   rename fusion_id=alt_fusion_id q3_apn=alt_q3_apn;
run;

proc sort data=ref_fusion_first_counts2;
   by gene_id strand sample_id;
proc sort data=alt_fusion_first_counts2;
   by gene_id strand sample_id;
run;

data afe_dataset;
  merge ref_fusion_first_counts2 (in=in1)  alt_fusion_first_counts2 (in=in2);
   by gene_id strand sample_id;
  if in1 and in2;
run;




data ref_fusion_last_counts2;
   set ref_fusion_last_counts;
   keep gene_id fusion_id strand sample_id q3_apn;
   rename fusion_id=ref_fusion_id q3_apn=ref_q3_apn;
run;

data alt_fusion_last_counts2;
   set alt_fusion_last_counts;
   keep gene_id fusion_id strand sample_id q3_apn;
   rename fusion_id=alt_fusion_id q3_apn=alt_q3_apn;
run;

proc sort data=ref_fusion_last_counts2;
   by gene_id strand sample_id;
proc sort data=alt_fusion_last_counts2;
   by gene_id strand sample_id;
run;

data ale_dataset;
  merge ref_fusion_last_counts2 (in=in1)  alt_fusion_last_counts2 (in=in2);
   by gene_id strand sample_id;
  if in1 and in2;
run;

data afe_psi;
  set afe_dataset;
  q3_apn_psi = (alt_q3_apn / (alt_q3_apn + ref_q3_apn)) * 100;
run;

data ale_psi;
  set ale_dataset;
  q3_apn_psi = (alt_q3_apn / (alt_q3_apn + ref_q3_apn)) * 100;
run;

data sample_key;
  set cc.design_file_new;
  where cell_type=425;
  keep sample_id status cell_type;
run; 

proc sort data=sample_key;
  by sample_id;
proc sort data=afe_psi;
  by sample_id;
proc sort data=ale_psi;
  by sample_id;
run;

data afe_psi_2;
  merge afe_psi (in=in1) sample_key (in=in2);
  by sample_id;
  if  in1 and in2;
run;


data ale_psi_2;
  merge ale_psi (in=in1) sample_key (in=in2);
  by sample_id;
  if  in1 and in2;
run;


proc sort data=afe_psi_2;
  by gene_id ref_fusion_id alt_fusion_id status;
proc means data=afe_psi_2 noprint;
  by gene_id ref_fusion_id alt_fusion_id status;
  var q3_apn_psi;
  output out=afe_psi_mean mean=;
run;

proc transpose data=afe_psi_mean out=afe_psi_mean_sbys prefix=psi_;
   by gene_id ref_fusion_id alt_fusion_id;
   id status;
   var q3_apn_psi;
run;

data afe_psi_mean_sbys2;
  set afe_psi_mean_sbys;
  psi_diff=psi_T1D-psi_CTL;
run;




proc sort data=ale_psi_2;
  by gene_id ref_fusion_id alt_fusion_id status;
proc means data=ale_psi_2 noprint;
  by gene_id ref_fusion_id alt_fusion_id status;
  var q3_apn_psi;
  output out=ale_psi_mean mean=;
run;

proc means data=afe_psi_2 noprint;
  by gene_id ref_fusion_id alt_fusion_id status;
  var q3_apn_psi;
  output out=afe_psi_sd stddev=;
run;

proc means data=ale_psi_2 noprint;
  by gene_id ref_fusion_id alt_fusion_id status;
  var q3_apn_psi;
  output out=ale_psi_sd stddev=;
run;



proc transpose data=ale_psi_mean out=ale_psi_mean_sbys prefix=psi_;
   by gene_id ref_fusion_id alt_fusion_id;
   id status;
   var q3_apn_psi;
run;

data ale_psi_mean_sbys2;
  set ale_psi_mean_sbys;
  psi_diff=psi_T1D-psi_CTL;
run;

proc transpose data=ale_psi_sd out=ale_psi_sd_sbys prefix=psi_sd_;
   by gene_id ref_fusion_id alt_fusion_id;
   id status;
   var q3_apn_psi;
run;


proc transpose data=afe_psi_sd out=afe_psi_sd_sbys prefix=psi_sd_;
   by gene_id ref_fusion_id alt_fusion_id;
   id status;
   var q3_apn_psi;
run;


data cc.afe_psi_mean_sbys;
   set afe_psi_mean_sbys2;
run;

data cc.ale_psi_mean_sbys;
   set ale_psi_mean_sbys2;
run;



data cc.afe_psi_sd_sbys;
   set afe_psi_sd_sbys;
run;

data cc.ale_psi_sd_sbys;
   set ale_psi_sd_sbys;
run;

/* now I shoale_psi_mean_sbys2uld add this to my splicing table
   then also check fragments in common to determine if mutually exclusive isoforms or complex splicing patterns 


OR leave it alone and finish the manuscript !

*/

/* Compare FRAGMENTS in only DS425 genes!! */


ods html close; ods listing;
libname cc '!PATCON/case_control/sas_data/fixed';
libname hg19 '!PATCON/useful_human_data/aceview_hg19/sas_data';

/* Event-to-transcript-to-gene index */

    data WORK.EVENT2XS    ;
    %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
    infile
    '!PATCON/case_control/references/hg19_aceview/annotations/hg19_aceview_event2transcript2gene_index.csv'
    delimiter = ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
       informat event_id $24. ;
       informat transcript_id $53. ;
       informat gene_id $36. ;
       informat annotation_frequency $12. ;
       informat flag_multigene best32. ;
       format event_id $24. ;
       format transcript_id $53. ;
       format gene_id $36. ;
       format annotation_frequency $12. ;
       format flag_multigene best12. ;
    input
                event_id $
                transcript_id $
                gene_id $
                annotation_frequency $
                flag_multigene
    ;
    if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
    run;

/* For each cell type, I first need to get a list of genes
   For each gene, I then need to get the list of detected transcripts
		-> If gene has 1 transcript, then event_concordance=1. Move onto the next
   For each transcript, count:
			the number events per transcript
			the number events in total between the two transcripts
			the number of events common
			the proportion of total events in common (n_XS1,2 / U_XS1,2)
   No need to do ALL pairwise permutations, just the unique combinations (these can stack by flipping XS1 and XS2)
*/

data xs2gene;
  set event2xs;
  keep transcript_id gene_id;
run;

proc sort data=xs2gene nodup;
  by transcript_id gene_id;
run;

data event2xs2;
   set event2xs;
   where event_id ^? "junction_";
run;


proc sort data=event2xs2 ;
  by transcript_id gene_id ;
run;

%macro loopGenes(cell);

 /* First determine what genes and transcripts I need to look at */

  data xs_dtct;
     set cc.flag_xscript_tpm_gt0_ge05_&cell.;
     if flag_&cell._CTL_tpm_gt0 = 1 or flag_&cell._T1D_tpm_gt0 = 1;
     keep transcript_id;
  run;
  
  proc sort data=xs_dtct;
    by transcript_id;
  run;
  
  data xs2gene_exp;
    merge xs_dtct (in=in1) xs2gene (in=in2);
    by transcript_id;
    if in1 and in2;
  run;
  
  data event2xs_exp;
    merge xs_dtct (in=in1) event2xs2 (in=in2);
    by transcript_id;
    if in1 and in2;
  run;
  
  data gene_list;
    set xs2gene_exp;
    keep gene_id;
  run;
  
  proc sort data=gene_list nodup;
    by gene_id;
  run;

  /* Set up dataset to append to */
  
  data xs_compare_cell_&cell.;
    format cell_type best12.;
    length gene_id $36.;
    length transcript_id_A $53.;
    length transcript_id_B $53.;
    format in_both best12.;
    format in_A_only best12.;
    format in_B_only best12.;
    format in_none best12.;
    format num_frag_xs_A best12.;
    format num_frag_xs_B best12.;
    format num_frag_total best12.;
    format num_frag_common best12.;
    format perc_frag_common best12.;
    format perc_frag_common_refA best12.;
    if _n_ > 0 then delete;
  run;


  /* Pull out transcripts from each gene to iterate over */
  %macro compareXSpairs();

    /* declare macro variables to use */
    %local iterGene iterXS queryXS refXS iterXSref countGene countXS inGene xsPerGene;
  
    /* Get the number of genes to iterate over */

    data _null_;
      set gene_list nobs=n;
      call symputx('countGene', n);
      stop;
    run;

    /* initiate loop for genes */
    %let iterGene=1;
    %do %while (&iterGene. <= &countGene. );

        /* Extract gene info */
        data _null_ ;
            set gene_list (firstobs=&iterGene. obs=&iterGene.); * read 1 record;
            call symput("inGene", strip(gene_id)); *get gene ID;
        run;
    
        /* Subset transcripts within selected gene */
        data xscripts;
            set xs2gene_exp;
            where gene_id="&inGene.";
        run;

        /* Get number of transcripts to iterate over */
        data _null_;
            set xscripts nobs=n;
            call symputx('countXS', n);
            stop;
        run;

        /* Initiate loop for transcripts */
        %let iterXS=1;
        %do %while (&iterXS. <= &countXS.);

            /* Initiate loop for transcript A */
            %let iterXSref=1;
             %do %while (&iterXSref. <= &countXS.);

               /* Get transcript to reference (transcript A) */
               data _null_;
                  set xscripts (firstobs=&iterXSref. obs=&iterXSref.); * read 1 record;
                  call symput("refXS", strip(transcript_id));
               run;
              
               /* Get transcript to query (transcript B) */
               data _null_;
                  set xscripts (firstobs=&iterXS. obs=&iterXS.); * read 1 record;
                  call symput("queryXS", strip(transcript_id));
               run;

               data event_xsA;
                   set event2xs_exp;
                   where transcript_id="&refXS.";
                   keep event_id;
               run;

               data event_xsB;
                   set event2xs_exp;
                   where transcript_id="&queryXS.";
                   keep event_id;
               run;

               proc sort data=event_xsA;
                  by event_id;
               proc sort data=event_xsB;
                  by event_id;
               run;

               data event_xsAvB;
                  merge event_xsA (in=in1) event_xsB (in=in2);
                  by event_id;
                  if in1 then flag_in_A=1; else flag_in_A=0;
                  if in2 then flag_in_B=1; else flag_in_B=0;
               run;

               proc freq data=event_xsAvB noprint;
                    tables flag_in_A*flag_in_B / out=event_ctab;
               run;

               data event_ctab2;
                   set event_ctab;
                   length varname $32.;
                   if flag_in_A=1 and flag_in_B=1 then varname="in_both";
                   else if flag_in_A=1 and flag_in_B=0 then varname="in_A_only";
                   else if flag_in_A=0 and flag_in_B=1 then varname="in_B_only";
                   else varname="in_none";
                   keep varname count;
               run;
             
               proc transpose data=event_ctab2 out=event_ctab_sbys(drop=_NAME_ _LABEL_);
                    var count;
                    id varname;
               run;

               /* Assemble comparison data */


               data xsPair_roz;
                 format cell_type best12.;
                 length gene_id $36.;
                 length transcript_id_A $53.;
                 length transcript_id_B $53.;
                 format in_both best12.;
                 format in_A_only best12.;
                 format in_B_only best12.;
                 format in_none best12.;
                 format num_frag_xs_A best12.;
                 format num_frag_xs_B best12.;
                 format num_frag_total best12.;
                 format num_frag_common best12.;
                 format perc_events_common best12.;
                 format perc_events_common_refA best12.; 
                 set event_ctab_sbys;
                 cell_type=&cell.;
                 gene_id="&inGene.";
                 transcript_id_A="&refXS.";
                 transcript_id_B="&queryXS.";
                 if in_both=. then in_both=0;
                 if in_A_only=. then in_A_only=0;
                 if in_B_only=. then in_B_only=0;
                 if in_none=. then in_none=0;
                 num_events_xs_A=in_both + in_A_only;
                 num_events_xs_B=in_both + in_B_only;
                 num_events_total=in_both + in_A_only + in_B_only;
                 num_events_common=in_both;
                 perc_events_common=num_events_common/num_events_total;
                 perc_events_common_refA=num_events_common/num_events_xs_A;
               run;
  
               /* Append to dataset to save */

               proc append base=xs_compare_cell_&cell. data=xsPair_roz force;
               run;

               /* increment ref transcript loop */
               %let iterXSref=%eval(&iterXSref. + 1);
             %end;

             /* increment query transcript loop */
             %let iterXS=%eval(&iterXS. + 1);

           %end;

           /* increment query transcript loop */
           %let iterGene=%eval(&iterGene. + 1);

        %end;
  %mend;

  %compareXSpairs();

%mend;


proc printto log="/home/jrbnewman/saslog.txt";
run;

%loopGenes(127);
%loopGenes(400);
%loopGenes(425);
%loopGenes(426);

proc printto ;
run;


/* Stack and make permenant */

data cc.xscript_pairwise_frag_compare;
   set xs_compare_cell_127 xs_compare_cell_400
       xs_compare_cell_425 xs_compare_cell_426;
run;


