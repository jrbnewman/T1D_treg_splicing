ods listing; ods html close;

libname cc '!PATCON/case_control/sas_data/fixed';
libname hg19 '!PATCON/useful_human_data/aceview_hg19/sas_data';

/* Import fragment counts */

proc import datafile="/mnt/store/cc_sandbox/event_analysis_output/hpc/bri_cc_fragment_counts_wide.tsv"
   out=frag_counts dbms=tab replace;
   guessingrows=5000;
run;

proc sort data=frag_counts;
   by event_id;
proc transpose data=frag_counts out=frag_counts_tall;
  by event_id;
run;

/* Make permenant */

data cc.counts_by_fragment;
  set frag_counts_tall;
  rename event_id=fragment_id
         _NAME_=sample_id
         COL1=apn;
run;

