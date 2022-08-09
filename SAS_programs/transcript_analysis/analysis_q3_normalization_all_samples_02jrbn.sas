/* Import transcript counts using iterdataset */

ods listing; ods html close;
libname cc '!PATCON/case_control/sas_data/fixed';
libname hg19 '!PATCON/useful_human_data/aceview_hg19/sas_data';

/* Normalizations : should this be done across the entire experiment or per cell type (given different transcriptomes?

Will use Q3/Q3 normalizations, as this has typically served me well in the past. Include line centering as well (try cell_type, cell_type*treatment, sample)

Going to try both ways and look to see if one way is better than the other (want to plot everything together)
First do across whole experiment in this program, and a second program to do by cell type

*/

data counts;
   set cc.hg19_rsem_100perc_apn2_xscripts;
   keep sample_id tpm transcript_id;
run;


data counts2;
   set counts;
   if tpm > 0;
run;

/* quartiles */

proc sort data=counts2;
   by sample_id;
proc means data=counts2 noprint;
   by sample_id;
   var tpm;
   output out=quartiles_tpm0 q3=q3;
   run;

/* calculate the q3 of q3 for normalizations */
proc means data=quartiles_tpm0 q3;
   var q3;
   run;

/*
Analysis Variable : q3

            Upper
         Quartile
     ------------
       50.6275000
     ------------
*/


/* Now we need to merge the stats back into the original dataset */

    proc sort data=counts;
        by sample_id;
        run;

    proc sort data=quartiles_tpm0;
        by sample_id;
        run;

/* merge and oops test */


    data counts_w_q3 oops_tpm0;
        merge counts (in=in1) quartiles_tpm0 (in=in2);
        by sample_id;
        if in1 and in2 then output counts_w_q3;
        else output oops_tpm0; 
        run;

/* Calculate the adjustments per sample */
/* Divide the across-sample UQ by the sample UQ and then output */
/* Then we want to export this as plot the distributions */

data cc.tpm_counts_q3_norm_all;
        set counts_w_q3;
        log_tpm=log(tpm + 1);
        * Q3 normalization;
        q3_tpm=(tpm/q3) * 50.6275000;
        q3_ff=50.6275000/q3;
        log_q3_tpm=log(q3_tpm + 1);
run;

/* Line centering : different options: Mean/median/q3, cell/condition/sample */

*Get design file, as I need treatment, genotype and timepoint;

data sample_key;
  set cc.design_file_new;
  keep sample_id cell_type status;
run;

proc sort data=sample_key nodup;
   by sample_id;
proc sort data=cc.tpm_counts_q3_norm_all;
   by sample_id;
run;


data norm_q3_w_key;
   merge sample_key (in=in1) cc.tpm_counts_q3_norm_all (in=in2);
   by sample_id;
   if in1 and in2;
run;

/* First center by cell type */

proc sort data=norm_q3_w_key;
   by cell_type;
proc means data=norm_q3_w_key noprint;
   where log_q3_tpm > 0;
   by cell_type;
   var log_q3_tpm;
   output out=means mean=mean_log_q3_tpm_cell median=median_log_q3_tpm_cell q3=q3_log_q3_tpm_cell;
run;

proc sort data=means;
   by cell_type;
run;

data norm_q3_w_key2;
   merge norm_q3_w_key (in=in1) means (in=in2);
   by cell_type;
   if in1;
   drop _type_ _freq_;
run;

*Calculate centered values;

data tpm_q3_norm_cent_cell;
   set norm_q3_w_key2;
   meancent_log_q3_tpm_cell=log_q3_tpm - mean_log_q3_tpm_cell;
   medcent_log_q3_tpm_cell=log_q3_tpm - median_log_q3_tpm_cell;
   q3cent_log_q3_tpm_cell=log_q3_tpm - q3_log_q3_tpm_cell;
run;

/* Now cell*status */

proc sort data=norm_q3_w_key;
   by cell_type status;
proc means data=norm_q3_w_key noprint;
   where log_q3_tpm > 0;
   by cell_type status;
   var log_q3_tpm;
   output out=means mean=mean_log_q3_tpm_condit median=median_log_q3_tpm_condit q3=q3_log_q3_tpm_condit;
run;

proc sort data=means;
   by cell_type status;
run;

data norm_q3_w_key2;
   merge norm_q3_w_key (in=in1) means (in=in2);
   by cell_type status;
   if in1;
   drop _type_ _freq_;
run;

*Calculate centered values;

data tpm_q3_norm_cent_condit;
   set norm_q3_w_key2;
   meancent_log_q3_tpm_condit=log_q3_tpm - mean_log_q3_tpm_condit;
   medcent_log_q3_tpm_condit=log_q3_tpm - median_log_q3_tpm_condit;
   q3cent_log_q3_tpm_condit=log_q3_tpm - q3_log_q3_tpm_condit;
run;

/* Now by individual */

data sample2subject;
   set cc.sample2subject_key;
   keep sample_id new_subject_id;
run;

proc sort data=sample2subject;
   by sample_id;
proc sort data=norm_q3_w_key;
   by sample_id;
run;

data norm_q3_w_key2;
  merge sample2subject (in=in1) norm_q3_w_key (in=in2);
  by sample_id;
  if in1 and in2;
run;


proc sort data=norm_q3_w_key2;
   by new_subject_id;
proc means data=norm_q3_w_key2 noprint;
   where log_q3_tpm > 0;
   by new_subject_id;
   var log_q3_tpm;
   output out=means mean=mean_log_q3_tpm_subj median=median_log_q3_tpm_subj q3=q3_log_q3_tpm_subj;
run;

proc sort data=means;
   by new_subject_id;
run;

data norm_q3_w_key3;
   merge norm_q3_w_key2 (in=in1) means (in=in2);
   by new_subject_id;
   if in1;
   drop _type_ _freq_;
run;

*Calculate centered values;

data tpm_q3_norm_cent_subj;
   set norm_q3_w_key3;
   meancent_log_q3_tpm_subj=log_q3_tpm - mean_log_q3_tpm_subj;
   medcent_log_q3_tpm_subj=log_q3_tpm - median_log_q3_tpm_subj;
   q3cent_log_q3_tpm_subj=log_q3_tpm - q3_log_q3_tpm_subj;
run;


/* Merge all together, make permenant and export for plots */

data tpm_q3_norm_cent_subj2;
   set tpm_q3_norm_cent_subj;
   keep sample_id transcript_id mean_log_q3_tpm_subj median_log_q3_tpm_subj q3_log_q3_tpm_subj
        meancent_log_q3_tpm_subj medcent_log_q3_tpm_subj q3cent_log_q3_tpm_subj;
run;

data tpm_q3_norm_cent_condit2;
   set tpm_q3_norm_cent_condit;
   keep sample_id transcript_id mean_log_q3_tpm_condit median_log_q3_tpm_condit q3_log_q3_tpm_condit
        meancent_log_q3_tpm_condit medcent_log_q3_tpm_condit q3cent_log_q3_tpm_condit;
run;

proc sort data=tpm_q3_norm_cent_cell;
   by sample_id transcript_id;
proc sort data=tpm_q3_norm_cent_subj2;
   by sample_id transcript_id;
proc sort data=tpm_q3_norm_cent_condit2;
   by sample_id transcript_id;
run;

data tpm_q3_norm_cent_all;
  merge tpm_q3_norm_cent_cell (in=in1) tpm_q3_norm_cent_subj2 (in=in2) tpm_q3_norm_cent_condit2 (in=in3);
  by sample_id transcript_id;
  if in1 and in2 and in3;
run;



data tpm_q3_norm_cent_all_gt0;
   set  tpm_q3_norm_cent_all;
   where tpm>0;
   length condition $8.;
   condition=catx("_",cell_type,status);
run;

proc export data=tpm_q3_norm_cent_all_gt0
     outfile="!PATCON/case_control/text_data/counts_by_transcript_all_tpm_norm_cent_gt0_fixed_02jrbn.csv"
     dbms=csv replace;
run;


data cc.tpm_counts_q3_norm_cent_all2;
   set tpm_q3_norm_cent_all;
run;



