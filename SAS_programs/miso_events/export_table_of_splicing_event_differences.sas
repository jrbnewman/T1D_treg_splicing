/* Import MISO annotations and assign to events */

ods listing; ods html close;
libname cc '!PATCON/case_control/sas_data/fixed';
libname hg19 '!PATCON/useful_human_data/aceview_hg19/sas_data';
libname cclocal '/TB14/TB14/immigrans/store/cc_sandbox/sas_data';


proc datasets lib=work kill noprint;
run;
quit;


/* Create big splicing table

to miso_annot_w_psi, I need to append AFE and ALE events

need columns:
event type
gene_id
spliced_in_coord
miso_id
percent_spliced_in_CTL
percent_spliced_in_T1D
diff_percent_spliced_in
*/

proc import datafile="!PATCON/case_control/references/hg19_aceview/annotations/hg19_aceview_fusion_annotations.csv"
    out=fusion_annot dbms=csv replace;
    guessingrows=max;
run;


data ref_fus_annot;
   set fusion_annot;
   length ref_coord $100.;
   ref_coord=catx(":",chr,fusion_start,fusion_stop,strand);
   keep fusion_id ref_coord;
   rename fusion_id=ref_fusion_id;
run;

data alt_fus_annot;
   set fusion_annot;
   length alt_coord $100.;
   alt_coord=catx(":",chr,fusion_start,fusion_stop,strand);
   keep fusion_id alt_coord;
   rename fusion_id=alt_fusion_id;
run;



data stack_afe_ale;
     set cc.afe_psi_mean_sbys (in=in1)
         cc.ale_psi_mean_sbys (in=in2)  ;
     length event_type $5.;
     if in1 then event_type="AFE";
     if in2 then event_type="ALE";
     drop _NAME_;
     rename psi_CTL=percent_spliced_in_CTL
 psi_T1D=percent_spliced_in_T1d
 psi_diff=diff_percent_spliced_in;
run;

proc sort data=stack_afe_ale;
   by ref_fusion_id;
proc sort data=ref_fus_annot;
   by ref_fusion_id;
run;

data afe_ale_fus1;
  merge stack_afe_ale (in=in1) ref_fus_annot (in=in2);
  by ref_fusion_id;
  if in1 and in2;
run;

proc sort data=afe_ale_fus1;
   by alt_fusion_id;
proc sort data=alt_fus_annot;
   by alt_fusion_id;
run;

data afe_ale_fus2;
  merge afe_ale_fus1 (in=in1) alt_fus_annot (in=in2);
  by alt_fusion_id;
  if in1 and in2;
run;

data afe_ale_fus3;
   set afe_ale_fus2;
   length spliced_in_coord $100.;
   spliced_in_coord=catx("/",ref_coord, alt_coord);
   drop ref_coord alt_coord alt_fusion_id ref_fusion_id;
   run;

/* Stack with other events */

data miso_events;
   set cc.miso_annot_w_psi;
run;

data all_events;
   set miso_events afe_ale_fus3;
   drop p_psi_q3_apn_avg miso_id;
run;


/* Now merge in DS results -- anova, ordinal rank, tiered rank */

data ds_test;
  set cc.m37_anova_fdr_ds_case4;
  where cell_type=425 and flag_fdr05=1;
  keep gene_id Probf fdr_p;
  rename probf=nominal_P_DS_test fdr_p=FDR_P_DS_test;
run;


data orank_test;
  set cc.m37_anova_fdr_ranktest_1_n;
  where cell_type=425;
  keep gene_id Probf fdr_p;
  rename probf=nominal_P_ordinal_rank fdr_p=FDR_P_ordinal_rank;
run;

data trank_test;
  set cc.m37_anova_fdr_ranktest_wght;
  where cell_type=425;
  keep gene_id Probf fdr_p;
  rename probf=nominal_P_tiered_rank fdr_p=FDR_P_tiered_rank;
run;

proc sort data=all_events;
  by gene_id;
proc sort data=ds_test;
  by gene_id;
proc sort data=orank_test;
  by gene_id;
proc sort data=trank_test;
  by gene_id;
run;


data DS_table;
  merge ds_test (in=in1) orank_test (in=in2) trank_test (in=in3) all_events (in=in4);
  by gene_id;
  if not in4 then event_type="n/a";
  if in1 then output;
run;


proc import datafile="!PATCON/useful_human_data/aceview_hg19/SpliceAid_factors_aceview_id.csv"
   out=splicing_factor_list dbms=csv replace;
   guessingrows=max;
run;

proc sort data=ds_table;
  by gene_id;
proc sort data=splicing_factor_list;
   by gene_id;
run;

data ds_table2;
  merge ds_table (in=in1) splicing_factor_list (in=in2);
   by gene_id;
  if in2 then flag_splicing_Factor=1;
   else flag_splicing_Factor=0;
  if in1;
run;


proc export data=ds_table2
   outfile="!PATCON/case_control/DS_table_425_sig_genes_all_events.csv"
   dbms=csv replace;
run;
