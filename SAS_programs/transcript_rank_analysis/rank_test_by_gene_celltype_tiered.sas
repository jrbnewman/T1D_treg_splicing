

ods listing; ods HOME/concannon close;
libname oldcc '!PATCON/case_control/sas_data/fixed';
libname cc '!HOME/concannon/BRI/RNAseq/sas_data';
libname hg19 '!HOME/concannon/useful_human_data/aceview_hg19/sas_data';


/* Rank tests on transcripts:

Try base model (no covariates)
And adaptation of Model 37 */

data xs_rank;
set cc.xscript_rank_by_cell_subj_all;
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
run;




data design;
set cc.design_file_final;
keep sample_id cell_type status_king01_d sex_king01_a plate Flag_for_removal_ambig3;
rename status_king01_d=status sex_king01_a=sex;
run;



data peer;
set oldcc.peer_factors_q3norm_common;
keep sample_id p1 p2;
run;

proc sort data=xs_rank;
by sample_id;
proc sort data=design;
by sample_id;
proc sort data=peer;
by sample_id;
run;

data xs_rank_w_key;
merge design (in=in1) peer (in=in2) xs_rank (in=in3);
by sample_id;
if in1 and in2 and in3;
run;

proc sort data=xs_rank_w_key;
by cell_type gene_id;
run;

ods listing close;
proc glm data=xs_rank_w_key ;
by cell_type gene_id;
where Flag_for_removal_ambig3=0;
class transcript_id status sex plate;
model xscript_3lvl = transcript_id status transcript_id*status sex status*sex plate p1|p2 / ss3;
lsmeans transcript_id*status / pdiff cl;
output out=exon_rank_resid_m37 predicted=pred residual=resid student=stu;
ods output modelanova=xs_rank_glm_m37 diff=pdiff_m37 ;
run;

proc glm data=xs_rank_w_key ;
by cell_type gene_id;
where Flag_for_removal_ambig3=0;
class transcript_id status sex plate;
model xscript_3lvl = transcript_id status transcript_id*status sex status*sex plate p1|p2 / ss3;
lsmeans transcript_id*status / slice=transcript_id;
output out=exon_rank_resid_m37_2 predicted=pred residual=resid student=stu;
ods output modelanova=xs_rank_glm_m37_2 slicedanova=slice_m37 ;
run;


/* Make permenant */


data cc.xscript_3rank_test_m37_anova_v3;
set xs_rank_glm_m37;
run;

data cc.xscript_3rank_test_m37_resid_v3;
set exon_rank_resid_m37;
run;

data cc.xscript_3rank_test_m37_pdiffs;
set pdiff_m37;
run;

data cc.xscript_3rank_test_m37_anova_v3b;
set xs_rank_glm_m37;
run;

data cc.xscript_3rank_test_m37_resid_v3b;
set exon_rank_resid_m37;
run;

data cc.xscript_3rank_test_m37_slices;
set slice_m37;
run;


/* Calculate FDR */


data ranktest;
  set cc.xscript_3rank_test_m37_anova_v3;
  where source ="transcript_id*status";
  keep gene_id cell_type fvalue probf ;
run;

/* FDR by cell type and effect */

proc sort data=ranktest;
  by cell_type ProbF;
run;

proc multtest inpvalues(ProbF)=ranktest fdr
              out=ranktest_fdr noprint;
run; quit;

/* Flag P and FDR at different levels */

data flag_ranktest_fdr;
  set ranktest_fdr;
  if probF = . then do;
      flag_p05=.;
      flag_p01=.;
      flag_p001=.;
      flag_fdr05=.;
      flag_fdr10=.;
      flag_fdr20=.;
      end;
  else do;
  if ProbF < 0.05 and ProbF > 0 then flag_p05 = 1; else flag_p05=0;
  if ProbF < 0.01 and ProbF > 0  then flag_p01 = 1; else flag_p01=0;
  if ProbF < 0.001 and ProbF > 0  then flag_p001 = 1; else flag_p001=0;
  if fdr_p < 0.05 and ProbF > 0  then flag_fdr05 = 1; else flag_fdr05=0;
  if fdr_p < 0.1 and ProbF > 0  then flag_fdr10 = 1; else flag_fdr10=0;
  if fdr_p < 0.2 and ProbF > 0  then flag_fdr20 = 1; else flag_fdr20=0;
  end;
run;

proc freq data=flag_ranktest_fdr;
  by cell_type ;
  tables flag_p05 flag_p01 flag_p001
         flag_fdr05 flag_fdr10 flag_fdr20;
run;

/* Make permenant */

data cc.m37_anova_fdr_ranktest_3lvl;
   set flag_ranktest_fdr;
run;

