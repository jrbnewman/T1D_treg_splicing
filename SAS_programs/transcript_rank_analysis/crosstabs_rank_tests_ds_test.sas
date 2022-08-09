ods listing; ods html close;
libname cc '!PATCON/case_control/sas_data/fixed';
libname hg19 '!PATCON/useful_human_data/aceview_hg19/sas_data';

proc datasets lib=work kill noprint;
run;
quit;

/* Compare rank differences with sig DS */

%macro compareRank(rankDiff,rankTest,sigTest);

data rank_diff;
   set cc.&rankDiff.;
run;

data rank_test;
   set cc.&rankTest.;
   keep cell_type gene_id flag_fdr05;
   rename flag_fdr05=flag_Rank_sig;
run;


data ds_test;
   set cc.&sigTest.;
   keep cell_type gene_id flag_fdr05;
   rename flag_fdr05=flag_ds_sig;
run;

proc sort data=rank_diff;
  by cell_type gene_id;
proc sort data=rank_test;
  by cell_type gene_id;
proc sort data=ds_test;
  by cell_type gene_id;
run;

data compare;
   merge rank_diff (in=in1) rank_test (in=in2) ds_Test (in=in3);
   by cell_type gene_id;
   if in1 then flag_has_rankdiff=1; else flag_has_rankdiff=0;
   if in2 then flag_has_ranksig=1; else flag_has_ranksig=0;
   if in3 then flag_has_dssig=1; else flag_has_dssig=0;
   if flag_Rank_sig=. then flag_rank_sig=0;
run;


proc freq data=compare noprint;
   where cell_type=425 or cell_type=426;
   tables cell_type*flag_rank_freq_diff_gt_10perc*flag_Rank_sig*flag_ds_sig / out=perc10;
   tables cell_type*flag_rank_freq_diff_gt_20perc*flag_Rank_sig*flag_ds_sig / out=perc20;
   tables cell_type*flag_rank_freq_diff_gt_30perc*flag_Rank_sig*flag_ds_sig / out=perc30;
run;

proc print data=perc10;
proc print data=perc20;
proc print data=perc30;
run;


%mend;

%compareRank(rank_freq_diff_quant_by_gene,m37_anova_fdr_ranktest_1_n,m37_anova_fdr_ds_case4);
%compareRank(rank_freq_diff_quant_by_gn_rnk1,m37_anova_fdr_ranktest_1_n,m37_anova_fdr_ds_case4);
%compareRank(rank_freq_diff_wght_by_gene,m37_anova_fdr_ranktest_wght,m37_anova_fdr_ds_case4);
%compareRank(rank_freq_diff_wght_by_gn_rnk1,m37_anova_fdr_ranktest_wght,m37_anova_fdr_ds_case4);

