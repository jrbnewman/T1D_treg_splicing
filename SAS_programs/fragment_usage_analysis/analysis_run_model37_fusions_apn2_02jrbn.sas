*ods listing; ods html close;
libname oldcc '!PATCON/case_control/sas_data/fixed';
libname cc '!PATCON/BRI/RNAseq/sas_data';
libname cclocal '/TB14/TB14/immigrans/store/cc_sandbox/sas_data';

proc datasets lib=work kill noprint;
run;
quit;

%macro prepData(cell);

data xs2keep_&cell.;
  set oldcc.flag_fusion_apn_gt0_ge2_&cell.;
  where flag_&cell._CTL_apn_ge2=1 and flag_&cell._T1D_apn_ge2=1;
  keep event_id;
run;

data counts;
  set cclocal.fusion_q3_norm_counts;
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
  keep sample_id event_id log_q3_apn;
run;

data sample2keep;
  set cc.design_file_final;
  where cell_type=&cell.;
  keep sample_id;
run;

proc sort data=sample2keep nodup;
   by sample_id;
proc sort data=counts;
  by sample_id;
run;

data counts1;
  merge sample2keep (in=in1) counts (in=in2);
  by sample_id;
  if in1 and in2;
run;


proc sort data=counts1;
  by event_id;
proc sort data=xs2keep_&cell.;
  by event_id;
run;

data counts2;
   merge xs2keep_&cell. (in=in1) counts1 (in=in2);
   by event_id;
   if in1 and in2;
run;

/* Transpose and fill in missing data */

proc sort data=counts2;
  by event_id sample_id;
proc transpose data=counts2 out=counts_sbys(drop=_NAME_);
  by event_id;
  id sample_id;
  var log_q3_apn;
run;

data counts_sbys2;
  set counts_sbys;
  array change _numeric_;
  do over change;
     if change=. then change=0;
     end;
run;

proc transpose name=sample_id data=counts_sbys2 out=counts3(rename=(col1=log_q3_apn));
   by event_id;
   var Sample_:;
run;



data counts_&cell.;
   set counts3;
run;

%mend;

%prepData(425); %prepData(426);

data counts2;
  set counts_425 counts_426;
run;


/* Get fragment to gene annotation */

data frag2gene;
  set oldcc.events_summary_by_fusion_apn2;
  where flag_multigene=0;
  keep event_id gene_id;
run;

proc sort data=frag2gene;
  by event_id;
proc sort data=counts2;
  by event_id;
run;

data counts3;
  merge frag2gene (in=in1) counts2 (in=in2);
  by event_id;
  if in1 and in2;
run;




/* Get design file and covariates */

data design;
        set cc.design_file_final;
        keep sample_id cell_type status_king01_d sex_king01_a plate Flag_for_removal_ambig3;
        rename status_king01_d=status sex_king01_a=sex;
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
  by cell_type gene_id status sample_id event_id ;
run;




*ods graphics;

/* Run Model 37 */
ods listing close;
proc glimmix data=counts_w_key2 plots=studentpanel;
   by cell_type gene_id ;
   where Flag_for_removal_ambig3=0;
   class event_id status sex plate;
   model log_q3_apn = status|event_id sex status*sex plate p1|p2 / ddfm=kr;
   random plate;
   lsmeans status / diff ;
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


data samples;
  set counts_w_key;
   where Flag_for_removal_ambig3=0;
  keep sample_id cell_type status;
run;

proc sort data=samples nodup;
  by sample_id;
run;

ods listing;
proc freq data=samples;
  tables cell_type*status;
run;


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

data cc.m37_anova_results_fus_ds;
  set anova;
run;

data cc.m37_anova_residuals_fus_ds;
  set resid_by_cell;
run;

data cc.m37_anova_resid_failrate_fus_ds;
  set resid_fail_count;
run;

data cc.m37_model_fit_stats_fus_ds;
 set fit_stats;
run;

data cc.m37_model_lsmeans_fus_ds;
 set lsmeans_test;
run;

data cc.m37_model_lsmeans_diff_fus_ds;
 set diffs_test;
run;

data cc.m37_model_constrasts_fus_ds;
 set  model_contrasts;
run;

data cc.m37_model_estimates_fus_ds;
 set model_estimates;
run;


/* FDR on splicing test */




/* 425 vs 426 gene DE and DS FDRs and counts */


data ds_model_gene_de;
  set cc.m37_anova_results_fus_ds;
  where effect in ("status") and probf ne .;
  keep gene_id cell_type fvalue probf ;
run;

proc sort data=ds_model_gene_de;
  by cell_type ProbF;
run;

proc multtest inpvalues(ProbF)=ds_model_gene_de fdr
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
  tables  cell_type*flag_p05 cell_type*flag_fdr05;
run;


data cc.m37_anova_fdr_gene_fus_de;
   set flag_gene_de_fdr;
run;


/* DS model */

data ds_model;
  set cc.m37_anova_results_fus_ds;
  where effect in ("event_id*status") and probf ne .;
  keep gene_id cell_type fvalue probf ;
run;
/* FDR by cell type and effect */

proc sort data=ds_model;
  by cell_type ProbF;
run;

proc multtest inpvalues(ProbF)=ds_model fdr
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

data cc.m37_anova_fdr_gene_fus_ds;
   set flag_ds_fdr;
run;


/* To do:
1. compare between cell types
2. compare with weighted rank results
3. RRM1/splicing factor overrepresentation
*/


/* what genes were analyzed in both? */
data ds425 ds426;
  set cc.m37_anova_fdr_gene_fus_ds;
  if cell_type=425 then output ds425;
  if cell_type=426 then output ds426;
  keep gene_id;
run;

proc sort data=ds425;
  by gene_id;
proc sort data=ds426;
  by gene_id;
run;

data ds425_426;
  merge ds425 (in=in1) ds426 (in=in2);
  by gene_id;
  if in1 then flag_425=1; else flag_425=0;
  if in2 then flag_426=1; else flag_426=0;
run;

proc freq data=ds425_426;
 tables flag_425*flag_426;
run;

/* What sig genes are in both? */

data ds425 ds426;
  set cc.m37_anova_fdr_gene_fus_ds;
  where flag_fdr05=1;
  if cell_type=425 then output ds425;
  if cell_type=426 then output ds426;
  keep gene_id;
run;

proc sort data=ds425;
  by gene_id;
proc sort data=ds426;
  by gene_id;
run;

data ds425_426;
  merge ds425 (in=in1) ds426 (in=in2);
  by gene_id;
  if in1 then flag_425=1; else flag_425=0;
  if in2 then flag_426=1; else flag_426=0;
run;

proc freq data=ds425_426;
 tables flag_425*flag_426;
run;



/* How does transcript and fragment DS compare? */



data ds425 ds426;
  set cc.m37_anova_fdr_gene_fus_ds;
  where flag_fdr05=1;
  if cell_type=425 then output ds425;
  if cell_type=426 then output ds426;
  keep gene_id;
run;

data ds425_xs ds426_xs;
  set cc.m37_anova_fdr_ds_case4;
  where flag_fdr05=1;
  if cell_type=425 then output ds425_xs;
  if cell_type=426 then output ds426_xs;
  keep gene_id;
run;

proc sort data=ds425 nodup; by gene_id;
proc sort data=ds426 nodup; by gene_id;
proc sort data=ds425_xs nodup; by gene_id;
proc sort data=ds426_xs nodup; by gene_id;
run;

data ds425_compare;
  merge ds425 (in=in1) ds425_xs (in=in2);
  by gene_id;
  if in1 then flag_425_DEFU=1; else flag_425_DEFU=0;
  if in2 then flag_425_DSXS=1; else flag_425_DSXS=0;
run;

data ds426_compare;
  merge ds426 (in=in1) ds426_xs (in=in2);
  by gene_id;
  if in1 then flag_426_DEFU=1; else flag_426_DEFU=0;
  if in2 then flag_426_DSXS=1; else flag_426_DSXS=0;
run;

proc freq data=ds425_compare;
  tables flag_425_DEFU * flag_425_DSXS ;
run;

proc freq data=ds426_compare;
  tables flag_426_DEFU * flag_426_DSXS ;
run;





proc import datafile="!PATCON/useful_human_data/aceview_hg19/SpliceAid_factors_aceview_id.csv"
   out=splicing_factor_list dbms=csv replace;
   guessingrows=max;
run;

proc import datafile="!PATCON/useful_human_data/aceview_hg19/genes_with_rrm1_domains.txt"
   out=rrm1_genes dbms=csv replace;
   guessingrows=max;
run;

data frag_ds;
   set cc.m37_anova_fdr_gene_fus_ds;
   keep cell_type gene_id flag_fdr05;
run;

proc sort data=frag_ds;
  by gene_id;
proc sort data=splicing_factor_list;
  by gene_id;
proc sort data=rrm1_genes;
  by gene_id;
run;

data frag_ds_sf_rrm1;
  merge frag_ds (in=in1) splicing_factor_list (in=in2) rrm1_genes (in=in3);
  by gene_id;
  if in2 then flag_sf=1; else flag_sf=0;
  if in3 then flag_rrm1=1; else flag_rrm1=0;
  if in1 then output;
run;

proc sort data=frag_ds_sf_rrm1;
  by cell_type;
proc freq data=frag_ds_sf_rrm1;
  by cell_type;
  tables flag_fdr05*flag_sf / chisq ;
  tables flag_fdr05*flag_rrm1 / chisq ;
run;


proc freq data=frag_ds_sf_rrm1;
  by cell_type;
  where flag_rrm1=1;
  tables flag_fdr05*flag_sf / chisq ;
run;
