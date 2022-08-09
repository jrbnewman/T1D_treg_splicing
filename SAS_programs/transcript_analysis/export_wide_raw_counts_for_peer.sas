/* Libraries */

ods listing; ods html close;
libname cc '!HOME/concannon/case_control/sas_data/fixed';
libname hg19 '!HOME/concannon/useful_human_data/aceview_hg19/sas_data';

/* Export normalized and centered counts for running PEER
   These need to be in a wide format: columns=samples, rows=transcripts
   No headers, no row names.

   I want to make PEER sets for each the test set and the "common transcripts" set
   All samples together, rather than by cell type */

/* "complete common" set */

%macro xsCell(cell);
data xs&cell.;
   set cc.flag_xscript_tpm_gt0_ge05_&cell.;
   where flag_&cell._CTL_tpm_gt0 = 1 and flag_&cell._T1D_tpm_gt0 = 1;
   keep transcript_id;
run;

proc sort data=xs&cell.;
  by transcript_id;
run;

%mend;

%xsCell(127);
%xsCell(400);
%xsCell(425);
%xsCell(426);

data xs2keep;
  merge xs127 (in=in1) xs400 (in=in2) xs425 (in=in3) xs426 (in=in4);
  by transcript_id;
  if in1 and in2 and in3 and in4;
run;

data counts;
   set cc.tpm_counts_q3_norm_cent_all2;
   keep sample_id transcript_id log_tpm; 
run;

proc sort data=xs2keep;
   by transcript_id;
proc sort data=counts;
  by transcript_id;
run;

data counts2;
  merge xs2keep (in=in1) counts (in=in2); 
  by transcript_id;
  if in1 and in2;
run;

proc sort data=counts2;
   by transcript_id sample_id;
proc transpose data=counts2 out=counts_sbys_norm;
   by transcript_id;
   id sample_id;
   var log_tpm;
run;


* Export a "with headers" and "without headers" version. The "without headers" version is used by PEER, the "with headers"
  version is for tracking the order of samples and transcripts;

data counts_sbys_norm2;
  set counts_sbys_norm;
  drop _NAME_;
run;


proc export data=counts_sbys_norm2
     outfile="!PATCON/case_control/text_data/counts_sbys_raw_complete_common_counts_fixed.csv"
     dbms=csv replace;
run;
data counts_sbys_norm2;
  set counts_sbys_norm;
  drop transcript_id _NAME_;
run;

proc export data=counts_sbys_norm2
     outfile="!PATCON/case_control/text_data/counts_sbys_raw_complete_common_counts_noheaders_fixed.csv"
     dbms=csv replace; putnames=no;
run;

