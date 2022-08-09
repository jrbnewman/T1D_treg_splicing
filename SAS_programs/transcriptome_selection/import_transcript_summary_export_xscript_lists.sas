ods listing; ods html close;

libname cc '!PATCON/case_control/sas_data/fixed';
libname hg19 '!PATCON/useful_human_data/aceview_hg19/sas_data';

/* Import transcript summary and export transcripts with 100% of features detected for each cell type */

    data WORK.XSCRIPT_SUMMARY    ;
    %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
    infile
'!PATCON/case_control/event_analysis_output/bri_cc_summary_of_transcripts_exp_genes_apn2_min_10bp.tsv'
   delimiter='09'x MISSOVER DSD lrecl=32767 firstobs=2 ;
       informat gene_id $32. ;
       informat transcript_id $49. ;
       informat treatment_group $7. ;
       informat flag_gene_has_multigene_exon best32. ;
       informat num_fragments best32. ;
       informat num_detected_fragments best32. ;
       informat perc_fragments_detected best32. ;
       informat num_unique_fragments best32. ;
       informat num_detected_unique_fragments best32. ;
       informat perc_unique_fragments_detected best32. ;
       informat num_junctions best32. ;
       informat num_detected_junctions best32. ;
       informat perc_junctions_detected best32. ;
       informat num_unique_junctions best32. ;
       informat num_detected_unique_junctions best32. ;
       informat perc_unique_junctions_detected best32. ;
       informat num_events best32. ;
       informat num_detected_events best32. ;
       informat perc_events_detected best32. ;
       informat num_unique_events best32. ;
       informat num_detected_unique_events best32. ;
       informat perc_unique_events_detected best32. ;
          format gene_id $12. ;
          format transcript_id $49. ;
          format treatment_group $7. ;
          format flag_gene_has_multigene_exon best12. ;
          format num_fragments best12. ;
          format num_detected_fragments best12. ;
          format perc_fragments_detected best12. ;
          format num_unique_fragments best12. ;
          format num_detected_unique_fragments best12. ;
          format perc_unique_fragments_detected best12. ;
          format num_junctions best12. ;
          format num_detected_junctions best12. ;
          format perc_junctions_detected best12. ;
          format num_unique_junctions best12. ;
          format num_detected_unique_junctions best12. ;
          format perc_unique_junctions_detected best12. ;
          format num_events best12. ;
          format num_detected_events best12. ;
          format perc_events_detected best12. ;
          format num_unique_events best12. ;
          format num_detected_unique_events best12. ;
          format perc_unique_events_detected best12. ;
      input
                  gene_id $
                  transcript_id $
                  treatment_group $
                  flag_gene_has_multigene_exon
                  num_fragments
                  num_detected_fragments
                  perc_fragments_detected
                  num_unique_fragments
                  num_detected_unique_fragments
                  perc_unique_fragments_detected
                  num_junctions
                  num_detected_junctions
                  perc_junctions_detected
                  num_unique_junctions
                  num_detected_unique_junctions
                  perc_unique_junctions_detected
                  num_events
                  num_detected_events
                  perc_events_detected
                  num_unique_events
                  num_detected_unique_events
                  perc_unique_events_detected
      ;
      if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
      run;

data av_xscripts;
  set hg19.hg19_aceview_xs2gene_fasta_index;
run;

/* make permenant */

data cc.events_xscript_summary_apn2;
   set xscript_summary;
run;


/* We only want to look at transcripts with 100% of events detected per cell type,
   so let's subset these by condition */

data con_127 case_127 con_400 case_400 con_425 case_425 con_426 case_426 oops;
   set xscript_summary;
   where perc_events_detected = 100;
        if treatment_group="127_CTL" then output con_127;
   else if treatment_group="127_T1D" then output case_127;
   else if treatment_group="400_CTL" then output con_400;
   else if treatment_group="400_T1D" then output case_400;
   else if treatment_group="425_CTL" then output con_425;
   else if treatment_group="425_T1D" then output case_425;
   else if treatment_group="426_CTL" then output con_426;
   else if treatment_group="426_T1D" then output case_426;
   else output oops;
   keep transcript_id gene_id;
run;

/* Summary:
127_CTL:	13632 transcripts
127_T1D:	15271 transcripts
400_CTL:	25002 transcripts
400_T1D:	25628 transcripts
425_CTL:	11779 transcripts
425_T1D:	14927 transcripts
426_CTL:	13048 transcripts
426_T1D:	12349 transcripts

Check crosstabs within cell type
*/

proc sort data=con_127 nodup;
   by transcript_id;
proc sort data=case_127 nodup;
   by transcript_id;
proc sort data=con_400 nodup;
   by transcript_id;
proc sort data=case_400 nodup;
   by transcript_id;
proc sort data=con_425 nodup;
   by transcript_id;
proc sort data=case_425 nodup;
   by transcript_id;
proc sort data=con_426 nodup;
   by transcript_id;
proc sort data=case_426 nodup;
   by transcript_id;
run;

data xscript_list;
  merge con_127 (in=in1) case_127 (in=in2) 
        con_400 (in=in3) case_400 (in=in4) 
        con_425 (in=in5) case_425 (in=in6) 
        con_426 (in=in7) case_426 (in=in8) ;
   by transcript_id;
   if in1 then flag_in_con127=1; else flag_in_con127=0;
   if in2 then flag_in_case127=1; else flag_in_case127=0;

   if in3 then flag_in_con400=1; else flag_in_con400=0;
   if in4 then flag_in_case400=1; else flag_in_case400=0;

   if in5 then flag_in_con425=1; else flag_in_con425=0;
   if in6 then flag_in_case425=1; else flag_in_case425=0;

   if in7 then flag_in_con426=1; else flag_in_con426=0;
   if in8 then flag_in_case426=1; else flag_in_case426=0;
run;

/* 30844 distinct transcripts: so there is probably a core set of genes expressed,
    plus some cell-type specific ones */

proc freq data=xscript_list;
  tables flag_in_con127*flag_in_case127
         flag_in_con400*flag_in_case400
         flag_in_con425*flag_in_case425
         flag_in_con426*flag_in_case426;
run;


/* OVERLAP:

CELL	CONTROL	BOTH	CASE	TOTAL
127		321		13311	1960	15592
400		1817	23185	2443	27445
425		184		11595	3332	15111
426		1052	11996	353		13401


Export these lists and gene2xs indices for RSEM
*/

data av_xscripts;
  set hg19.hg19_aceview_xs2gene_fasta_index;
run;


%macro extractXS(cell);

data xs_&cell.;
   set xscript_list;
   where flag_in_con&cell.=1 or flag_in_case&cell.=1;
   keep transcript_id;
run;

proc sort data=xs_&cell.;
   by transcript_id;
proc sort data=av_xscripts;
   by transcript_id;
run;

data xs2gene_&cell.;
  merge av_xscripts (in=in1) xs_&cell. (in=in2);
  by transcript_id;
  if in1 and in2;
  keep gene_id transcript_fasta_id;
run;

data xs_&cell._2;
  set xs2gene_&cell.;
  drop gene_id;
run;

proc export data=xs2gene_&cell.
     outfile="!PATCON/case_control/references/hg19_aceview_transcripts_cellcode_&cell._100perc_apn2_gene2xs_v3.txt"
     dbms=tab replace; putnames=no;
run;

proc export data=xs_&cell._2
     outfile="!PATCON/case_control/references/hg19_aceview_transcripts_cellcode_&cell._100perc_apn2_list_v3.txt"
     dbms=tab replace; putnames=no;
run;

%mend;

%extractXS(127);
%extractXS(400);
%extractXS(425);
%extractXS(426);



