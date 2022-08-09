/* Count genes that have no "complete" transcripts by condition */

ods listing; ods html close;
libname cc '!PATCON/case_control/sas_data/fixed';
libname con '!PATCON/sas_data';

data flag_100;
  set cc.events_xscript_summary_apn2;
  if perc_events_detected=100 then flag_dtct_100=1;
  else flag_dtct_100=0;
  keep gene_id transcript_id treatment_group flag_dtct_100;
run;

proc sort data=flag_100 nodup;
  by gene_id transcript_id treatment_group;
run;

proc freq data=flag_100 noprint;
   tables transcript_id*treatment_group /out=check;
run;

data check2;
  set check;
  where count > 1;
run;

/* Okay no dups, good */

proc sort data=flag_100;
  by gene_id treatment_group;
proc means data=flag_100 noprint;
  by gene_id treatment_group;
  var flag_dtct_100;
  output out=xs_per_gene_by_trt(rename=(_FREQ_=num_transcripts)) sum=num_transcripts_100dtct;
run;

data flag_gene;
  set xs_per_gene_by_trt;
   if num_transcripts_100dtct = num_transcripts then flag_all_xs=1;
   else flag_all_xs=0;
   if num_transcripts_100dtct = 0 then do;
        flag_no_xs=1;
        flag_some_xs=0;
        end;
   else do;
        flag_some_xs=1; 
        flag_no_xs=0;
        end;
run;

proc sort data=flag_gene;
   by treatment_group;
proc freq data=flag_gene noprint;
   tables treatment_group*flag_all_xs*flag_some_xs*flag_no_xs / out=xs_count;
run;

proc print data=xs_count;
run;

/*
Group	N_all	N_some	N_none
127_CTL	2905	5067	11394
127_T1D	3166	5423	11168
400_CTL	4834	6970	10692
400_T1D	5271	7033	10663
425_CTL	2624	4732	11842
425_T1D	3121	5447	11543
426_CTL	2802	4997	11597
426_T1D	2711	4750	11562


About 50-60% of all expressed genes do so without complete transcripts
according to the event data (at APN>2)

Make density plots and decide if I should also look at APN>0 and APN>5 */

data genes2keep;
  set con.immunogene_flags;
  where flag_pseudogene=0;
  keep gene_id;
run;

data xs2gene;
   set flag_100;
   keep gene_id transcript_id;
run;


ods graphics;


%macro plotDensity(cell);

data con_xs;
   set flag_100;
   where flag_dtct_100=1 and treatment_group="&cell._CTL";
   keep transcript_id;
run;

data t1d_xs;
   set flag_100;
   where flag_dtct_100=1 and treatment_group="&cell._T1D";
   keep transcript_id;
run;

proc sort data=con_xs  nodup;
   by transcript_id;
proc sort data=t1d_xs  nodup;
   by transcript_id;
run;

data xs_counts;
  merge con_xs (in=in1) t1d_xs (in=in2);
  by transcript_id;
run;

proc sort data=xs_counts nodup;
  by transcript_id;
proc sort data=xs2gene nodup;
  by transcript_id gene_id;
run;

data xs_counts2;
  merge xs_counts (in=in1) xs2gene (in=in2) ;
  by transcript_id;
  if in1 and in2;
run;

proc sort data=xs_counts2;
   by gene_id;

proc sort data=genes2keep nodup;
   by gene_id;
run;

data xs_counts3;
  merge xs_counts2 (in=in1) genes2keep (in=in2);
  by gene_id;
  if in1 and in2;
run;


proc freq data=xs_counts3 noprint;
  tables gene_id / out=xs_per_gene;
run;

proc sgplot data=xs_per_gene;
   where count < 13;
   histogram count;
   xaxis label = "Transcripts per gene";
   title "Distribution of transcripts per gene for cell_type &cell., APN >= 2.";
run;

%mend;

%plotDensity(127);
%plotDensity(400);
%plotDensity(425);
%plotDensity(426);

/*

mv ~/SGPlot.png $PATCON/case_control/data_visualization/retained_xscript_per_gene_cell127_apn2_100dtct_fixed.png
mv ~/SGPlot1.png $PATCON/case_control/data_visualization/retained_xscript_per_gene_cell400_apn2_100dtct_fixed.png
mv ~/SGPlot2.png $PATCON/case_control/data_visualization/retained_xscript_per_gene_cell425_apn2_100dtct_fixed.png
mv ~/SGPlot3.png $PATCON/case_control/data_visualization/retained_xscript_per_gene_cell426_apn2_100dtct_fixed.png

Distributions aren't changing too much, so my transcript selection is probably okay
I will run the APN0 and APN5 in the background for completeness, but I suspect APN2 all the way

*/




