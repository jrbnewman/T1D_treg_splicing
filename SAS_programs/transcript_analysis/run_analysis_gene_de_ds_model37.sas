/* Set libraries */

*ods listing; ods html close;

libname oldcc '!HOME/concannon/case_control/sas_data';
libname cc "!HOME/concannon/BRI/RNAseq/sas_data";


/* get counts for test set */

%macro prepData(cell);

data xs2keep_&cell.;
  set oldcc.flag_xscript_tpm_gt0_ge05_&cell.;
  where flag_&cell._CTL_tpm_gt0=1 or flag_&cell._T1D_tpm_gt0=1;
  keep transcript_id;
run;

data counts;
  set oldcc.tpm_counts_q3_norm_cent_all2;
  where cell_type=&cell.;
if sample_id="Sample_P1_C1_97288" then delete;
if sample_id="Sample_P2_H2_11333425" then delete;
if sample_id="Sample_P3_H8_22378426" then delete;
if sample_id="Sample_P4_E5_81333426" then delete;
if sample_id="Sample_P4_E7_1641425" then delete;
if sample_id="Sample_P5_G7_97721127" then delete;
if sample_id="Sample_P6_A12_49925400" then delete;
if sample_id="Sample_P6_A9_28992400" then delete;
if sample_id="Sample_P6_B11_65680400" then delete;
if sample_id="Sample_P6_B12_29014400" then delete;
if sample_id="Sample_P6_B7_29941400" then delete;
if sample_id="Sample_P6_C11_47601400" then delete;
if sample_id="Sample_P6_C12_97721400" then delete;
if sample_id="Sample_P6_D3_24509425" then delete;
if sample_id="Sample_P6_E2_48902425" then delete;
if sample_id="Sample_P6_F11_8913040" then delete;
if sample_id="Sample_P6_B7_29941400" then delete;
if sample_id="Sample_P3_G6_85413426" then delete;
if sample_id="Sample_P3_H4_72730426" then delete;
if sample_id="Sample_P6_H2_25560426" then delete;
if sample_id="Sample_P1_C1_97288" then delete;
if sample_id="Sample_P2_E4_21239425" then delete;
  keep sample_id transcript_id log_q3_tpm;
run;

proc sort data=counts;
  by transcript_id;
proc sort data=xs2keep_&cell.;
  by transcript_id;
run;

data counts_&cell.;
   merge xs2keep_&cell. (in=in1) counts (in=in2);
   by transcript_id;
   if in1 and in2;
run;

%mend;

%prepData(425); %prepData(426);

data counts2;
  set counts_425 counts_426;
run;

data xs2gene;
   set oldcc.hg19_xscript2gene;
   *set hg19.hg19_xscript2gene;
run;

proc sort data=counts2;
  by transcript_id;
proc sort data=xs2gene;
  by transcript_id gene_id;
run;

data counts3 no_gene no_counts;
  merge xs2gene (in=in1) counts2 (in=in2);
  by transcript_id;
  if in1 and in2 then output counts3;
  else if in1 then output no_counts;
  else output no_gene;
run;


/* Get design file and covariates */

data design;
   set oldcc.design_file_final;
   keep sample_id cell_type status age sex plate;
run;

proc sort data=counts3;
  by sample_id;
proc sort data=design;
  by sample_id;
run;

data counts_w_key;
  merge design (in=in1) counts3 (in=in2);
  by sample_id;
  if in1 and in2;
run;

ods listing close;

data peer;
  set oldcc.peer_factors_q3norm_common; 
run;

proc sort data=counts_w_key;
  by sample_id;
proc sort data=peer;
  by sample_id;
run;

data counts_w_key2;
  merge peer (in=in1) counts_w_key (in=in2);
  by sample_id;
  if in1 and in2;
run;


proc sort data=counts_w_key2;
  by cell_type gene_id transcript_id status sample_id;
run;

*ods graphics;

/* Run Model 37 */

proc glimmix data=counts_w_key2 plots=studentpanel;
   by cell_type gene_id;
   class transcript_id status sex plate;
   model log_q3_tpm = transcript_id status transcript_id*status sex status*sex plate p1|p2 / ddfm=kr;
   random plate;
   lsmeans status transcript_id*status / diff ;
   /* CTL T1D */
   contrast 'CTL-T1D' status -1 1 ;
   estimate 'CTL-T1D' status -1 1 ;
   output out=resid_by_cell resid=resid pred=pred student=stu;
   ods output tests3=anova
              lsmeans=lsmeans_test
              diffs = diffs_test
              contrasts = model_contrasts
              estimates = model_estimates
              FitStatistics=fit_stats;

run;
quit;


* Flag residuals;
proc univariate data = resid_by_cell normal noprint;
  by cell_type gene_id ;
  var Resid;
  output out = normtest probn=pnorm;
  run;

data flag_resids_by_cell;
  set normtest;
  if pnorm = . then flag_fail_norm = .;
        else if pnorm le 0.05 then flag_fail_norm = 1;
        else flag_fail_norm = 0;
  run;

proc freq data = flag_resids_by_cell  noprint ;
  tables cell_type*flag_fail_norm / out=resid_fail_count;
  run;


/* Make permenant */

data cc.m37_anova_results_all_ds2;
  set anova;
run;

data cc.m37_anova_residuals_all_ds2;
  set resid_by_cell;
run;

data cc.m37_anova_resid_failrate_all_ds2;
  set resid_fail_count;
run;

data cc.m37_model_fit_stats_all_ds2;
 set fit_stats;
run;

data cc.m37_model_lsmeans_all_ds2;
 set lsmeans_test;
run;

data cc.m37_model_lsmeans_diff_all_ds2;
 set diffs_test;
run;

data cc.m37_model_constrasts_all_ds2;
 set  model_contrasts;
run;

data cc.m37_model_estimates_all_ds2;
 set model_estimates;
run;



/* FDR and counts -- gene DE */

data genes2keep;
  set oldcc.xscripts_per_gene_set_model_case;
  where case in (2,4);
  keep gene_id cell_type;
run;

data ds_model_gene_de;
  set cc.m37_anova_results_all_ds2;
  where effect in ("status");
  keep gene_id cell_type fvalue probf ;
run;

proc sort data=genes2keep;
  by cell_type gene_id;
proc sort data=ds_model_gene_de;
  by cell_type gene_id;
run;

data ds_model_gene_de2;
   merge genes2keep (in=in1) ds_model_gene_de (in=in2);
   by cell_type gene_id;
   if in1 and in2;
run;

/* FDR by cell type and effect */

proc sort data=ds_model_gene_de2;
  by cell_type ProbF;
run;

proc multtest inpvalues(ProbF)=ds_model_gene_de2 fdr
              out=ds_model_gene_de_fdr noprint;
			  by cell_type;
run; quit;

/* Flag P and FDR at different levels */

data flag_gene_de_fdr;
  set ds_model_gene_de_fdr;
  if ProbF < 0.05 then flag_p05 = 1; else flag_p05=0;
  if ProbF < 0.01 then flag_p01 = 1; else flag_p01=0;
  if ProbF < 0.001 then flag_p001 = 1; else flag_p001=0;
  if fdr_p < 0.05 then flag_fdr05 = 1; else flag_fdr05=0;
  if fdr_p < 0.1 then flag_fdr10 = 1; else flag_fdr10=0;
  if fdr_p < 0.2 then flag_fdr20 = 1; else flag_fdr20=0;
run;

proc freq data=flag_gene_de_fdr;
  by cell_type ;
  tables flag_p05 flag_p01 flag_p001
         flag_fdr05 flag_fdr10 flag_fdr20;
run;

/* Make permenant */

data cc.m37_anova_fdr_gene_de_cases_2_4;
   set flag_gene_de_fdr;
run;



/* FDR and counts -- gene DS */


data genes2keep;
  set oldcc.xscripts_per_gene_set_model_case;
  where case=4;
  keep gene_id cell_type;
run;

data ds_model;
  set cc.m37_anova_results_all_ds2;
  where effect in ("transcript_id*status");
  keep gene_id cell_type fvalue probf ;
run;

proc sort data=genes2keep;
  by cell_type gene_id;
proc sort data=ds_model;
  by cell_type gene_id;
run;

data ds_model2;
   merge genes2keep (in=in1) ds_model (in=in2);
   by cell_type gene_id;
   if in1 and in2;
run;

/* FDR by cell type and effect */

proc sort data=ds_model2;
  by cell_type ProbF;
run;

proc multtest inpvalues(ProbF)=ds_model2 fdr
              out=ds_model_fdr noprint;
	by cell_type;
run; quit;

/* Flag P and FDR at different levels */

data flag_ds_fdr;
  set ds_model_fdr;
  if ProbF < 0.05 then flag_p05 = 1; else flag_p05=0;
  if ProbF < 0.01 then flag_p01 = 1; else flag_p01=0;
  if ProbF < 0.001 then flag_p001 = 1; else flag_p001=0;
  if fdr_p < 0.05 then flag_fdr05 = 1; else flag_fdr05=0;
  if fdr_p < 0.1 then flag_fdr10 = 1; else flag_fdr10=0;
  if fdr_p < 0.2 then flag_fdr20 = 1; else flag_fdr20=0;
run;
ods listing; ods html close;
proc freq data=flag_ds_fdr;
  by cell_type ;
  tables flag_p05 flag_p01 flag_p001
         flag_fdr05 flag_fdr10 flag_fdr20;
run;

/* Make permenant */

data cc.m37_anova_fdr_ds_case4;
   set flag_ds_fdr;
run;



