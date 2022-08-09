/* FDR on differential splicing model and count significance at:

   P<0.05
   P<0.01
   P<0.001
   FDR < 0.05
   FDR < 0.10
   FDR < 0.20

   per cell type, for transcript*status and for status

CASES:

Case 0: Gene has no detected transcripts -> gene is off, no model, no need to worry about this

Case 1: Gene has detected transcripts in only T1D or CTL -> gene is DD, no model
Case 2: Gene has 1 transcript detected in both T1D and CTL -> use DE model
		-> check "status" effect in DS model ONLY to be sure
Case 3: Gene has 1 or more transcripts but these are mutually exclusive between T1D and CTL
		-> gene is DS by default, no model needed
		-> check DS model to be sure
Case 4: Gene has multiple transcripts, at least one expressed in both T1D and CTL -> DS model
*/

ods html close; ods listing;
libname cc '!PATCON/case_control/sas_data';
libname hg19 '!PATCON/useful_human_data/aceview_hg19/sas_data';
libname cc 'D:\concannon\case_control\sas_data\fixed';
libname hg19 'D:\concannon\aceview_hg19\sas_data';
/* First I need to bin genes into cases: ONLY doing FDR on things I can test, based on above */

data xs2gene;
   set hg19.hg19_xscript2gene;
run;

%macro assignCase(cell);

/* Get number of detected EA transcripts per gene (overall), and for cases, controls */

data xs_&cell._t1d xs_&cell._ctl xs_&cell._both ;
   set cc.flag_xscript_tpm_gt0_ge05_&cell.;
   if flag_&cell._CTL_tpm_gt0=1 then output xs_&cell._ctl;
   if flag_&cell._T1D_tpm_gt0=1 then output xs_&cell._t1d;
   if flag_&cell._T1D_tpm_gt0=1 and flag_&cell._CTL_tpm_gt0=1 then output xs_&cell._both;
   keep transcript_id;
run;

data xs_&cell._all;
  set xs_&cell._t1d xs_&cell._ctl xs_&cell._both;
run;

proc sort data=xs_&cell._all nodup;
   by transcript_id;
proc sort data=xs2gene nodup;
  by transcript_id gene_id;
run;

data xs_&cell._all2;
  merge xs2gene (in=in1) xs_&cell._all (in=in2);
  by transcript_id;
  if in1 and in2;
run;

proc sort data=xs_&cell._all2;
   by gene_id;
proc freq data=xs_&cell._all2 noprint;
   tables gene_id / out=xs_per_gene_all_&cell.;
run;

data xs_per_gene_all_&cell._2;
  set xs_per_gene_all_&cell.;
  keep gene_id count;
  rename count=xs_per_gene_total;
run;

/* Get number of detected case transcripts */
proc sort data=xs_&cell._t1d nodup;
   by transcript_id;
proc sort data=xs2gene nodup;
  by transcript_id gene_id;
run;

data xs_&cell._t1d2;
  merge xs2gene (in=in1) xs_&cell._t1d (in=in2);
  by transcript_id;
  if in1 and in2;
run;

proc sort data=xs_&cell._t1d2;
   by gene_id;
proc freq data=xs_&cell._t1d2 noprint;
   tables gene_id / out=xs_per_gene_t1d_&cell.;
run;

data xs_per_gene_t1d_&cell._2;
  set xs_per_gene_t1d_&cell.;
  keep gene_id count;
  rename count=xs_per_gene_t1d;
run;

/* Get number of detected control transcripts */
proc sort data=xs_&cell._ctl nodup;
   by transcript_id;
proc sort data=xs2gene nodup;
  by transcript_id gene_id;
run;

data xs_&cell._ctl2;
  merge xs2gene (in=in1) xs_&cell._ctl (in=in2);
  by transcript_id;
  if in1 and in2;
run;

proc sort data=xs_&cell._ctl2;
   by gene_id;
proc freq data=xs_&cell._ctl2 noprint;
   tables gene_id / out=xs_per_gene_ctl_&cell.;
run;

data xs_per_gene_ctl_&cell._2;
  set xs_per_gene_ctl_&cell.;
  keep gene_id count;
  rename count=xs_per_gene_ctl;
run;

/* Get number of transcripts detected in both cases and controls */
proc sort data=xs_&cell._both nodup;
   by transcript_id;
proc sort data=xs2gene nodup;
  by transcript_id gene_id;
run;

data xs_&cell._both2;
  merge xs2gene (in=in1) xs_&cell._both (in=in2);
  by transcript_id;
  if in1 and in2;
run;

proc sort data=xs_&cell._both2;
   by gene_id;
proc freq data=xs_&cell._both2 noprint;
   tables gene_id / out=xs_per_gene_both_&cell.;
run;

data xs_per_gene_both_&cell._2;
  set xs_per_gene_both_&cell.;
  keep gene_id count;
  rename count=xs_per_gene_both;
run;

proc sort data=xs_per_gene_all_&cell._2;
   by gene_id;
proc sort data=xs_per_gene_both_&cell._2;
   by gene_id;
proc sort data=xs_per_gene_ctl_&cell._2;
   by gene_id;
proc sort data=xs_per_gene_t1d_&cell._2;
   by gene_id;
run;

data xs_per_gene_&cell.;
   merge xs_per_gene_all_&cell._2 (in=in1) xs_per_gene_both_&cell._2 (in=in2)
         xs_per_gene_ctl_&cell._2 (in=in3) xs_per_gene_t1d_&cell._2 (in=in4); 
   by gene_id;
   if not in1 then xs_per_gene_total=0;
   if not in2 then xs_per_gene_both=0;
   if not in3 then xs_per_gene_ctl=0;
   if not in4 then xs_per_gene_t1d=0;
   cell_type=&cell.;
run;

%mend;

%assignCase(127);
%assignCase(400);
%assignCase(425);
%assignCase(426);

data xs_per_gene_all;
   set xs_per_gene_127 xs_per_gene_400 xs_per_gene_425 xs_per_gene_426;
run;

/* Set cases */

data xs_per_gene_cases;
   set xs_per_gene_all;
   /* Case 0: Gene is off */
   if xs_per_gene_ctl=0 and xs_per_gene_t1d = 0 then case=0;
  
   /* Case 1: Gene is differentially detected */
   else if xs_per_gene_ctl=0 and xs_per_gene_t1d > 0 then case=1;
   else if xs_per_gene_ctl > 0 and xs_per_gene_t1d = 0 then case=1;

   /* Case 2: DE model */
   else if xs_per_gene_total=1 and xs_per_gene_ctl=1 and xs_per_gene_t1d = 1 then case=2;

   /* Case 3: All transcripts mutually exclusive */ 
   else if xs_per_gene_total > 1 and xs_per_gene_both = 0 then case=3;

   /* Case 4: DS model */
   else if xs_per_gene_total > 1 and xs_per_gene_both > 0 then case=4;

   /* Case -9: Catch all for logic check */
   else case=-9;
run;

/* PROC FREQ to check cell_type*case */

proc sort data=xs_per_gene_cases;
  by cell_type case;
proc freq data=xs_per_gene_cases noprint;
  tables cell_type*case / out=cases_per_cell;
run;

proc print data=cases_per_cell;
run;

/*

         cell_
  Obs     type    case    COUNT    PERCENT

    1     127       1        17     0.0457
    2     127       2      5883    15.7997
    3     127       4      2612     7.0149
    4     400       1        46     0.1235
    5     400       2      7913    21.2515
    6     400       4      4570    12.2734
    7     425       1        18     0.0483
    8     425       2      5909    15.8695
    9     425       4      2534     6.8054
   10     426       1        14     0.0376
   11     426       2      5535    14.8650
   12     426       4      2184     5.8654



Case 1 (Gene is DD): 14-46 genes
Case 2 (DE model): 5535-7913 genes
Case 3 (Mutually exclusive): 0 genes
Case 4 (DS model): 2184-4570 genes

*/

/* Make permenant */

data cc.xscripts_per_gene_set_model_case;
   set xs_per_gene_cases;
run;



/* How do I want to summarize these data?

(1) Transcript DE/DD summary -- easy
		- up/down regulation? estimates

(2) Gene DS/DE summary:
    Case 1: Flag gene DD
    Case 2: DE model (status), DS model (status) --> check that isoform is DE
			-> monotranscript genes, 2-way Venn with Case 1 (all Case 1 have a single transcript)
	Case 4: DS model (status) for DE, DS model (transcript*status) for splicing
	Need to look at up/down reg of gene?

Compare: Gene DS vs Gene DE vs Gene has DE/DD isoform (3-way Venn diagram)
	-> case 4 only

After, then look at event-level data (for Monday, pull out IKZFs and any other sig T1D gene)
*/

