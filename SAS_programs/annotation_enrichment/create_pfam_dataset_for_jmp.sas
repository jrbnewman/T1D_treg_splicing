libname cc '!PATCON/case_control/sas_data/fixed';

ods listing;
ods html close;


data ds425 ds426;
    set cc.m37_anova_fdr_ds_case4;
    if cell_type=425 then output ds425;
    if cell_type=426 then output ds426;
    keep gene_id flag_fdr05;
run;

data de425 de426;
    set cc.m37_anova_fdr_gene_de_cases_2_4;
    if cell_type=425 then output de425;
    if cell_type=426 then output de426;
    keep gene_id flag_fdr05;
run;

data ds_frag425 ds_frag426;
    set cc.m37_anova_fdr_gene_frag_ds;
    if cell_type=425 then output ds_frag425;
    if cell_type=426 then output ds_frag426;
    keep gene_id flag_fdr05;
run;

data de_frag425 de_frag426;
    set cc.m37_anova_fdr_gene_frag_de;
    if cell_type=425 then output de_frag425;
    if cell_type=426 then output de_frag426;
    keep gene_id flag_fdr05;
run;

data q425 q426;
  set cc.m37_anova_fdr_ranktest_1_n;
     if cell_type=425 then output q425;
    if cell_type=426 then output q426;
    keep gene_id flag_fdr05;
run;


data w425 w426;
  set cc.m37_anova_fdr_ranktest_wght;
     if cell_type=425 then output w425;
    if cell_type=426 then output w426;
    keep gene_id flag_fdr05;
run;
 

data ds425_2; set ds425; rename flag_fdr05=flag_425_XS_DS; run;
data de425_2; set de425; rename flag_fdr05=flag_425_XS_DE; run;
data ds426_2; set ds426; rename flag_fdr05=flag_426_XS_DS; run;
data de426_2; set de426; rename flag_fdr05=flag_426_XS_DE; run;

data ds_frag425_2; set ds_frag425; rename flag_fdr05=flag_425_frag_DS; run;
data de_frag425_2; set de_frag425; rename flag_fdr05=flag_425_frag_DE; run;
data ds_frag426_2; set ds_frag426; rename flag_fdr05=flag_426_frag_DS; run;
data de_frag426_2; set de_frag426; rename flag_fdr05=flag_426_frag_DE; run;



data q425_2; set q425; rename flag_fdr05=flag_425_qrank; run;
data q426_2; set q426; rename flag_fdr05=flag_426_qrank; run;
data w425_2; set w425; rename flag_fdr05=flag_425_wrank; run;
data w426_2; set w426; rename flag_fdr05=flag_426_wrank; run;


proc sort data=ds425_2;  by gene_id;
proc sort data=de425_2;  by gene_id;
proc sort data=ds426_2;  by gene_id;
proc sort data=de426_2;  by gene_id;

proc sort data=ds_frag425_2;  by gene_id;
proc sort data=de_frag425_2;  by gene_id;
proc sort data=ds_frag426_2;  by gene_id;
proc sort data=de_frag426_2;  by gene_id;

proc sort data=q425_2;  by gene_id;
proc sort data=q426_2;  by gene_id;
proc sort data=w425_2;  by gene_id;
proc sort data=w426_2;  by gene_id;
run;

data all_flags;
  merge ds425_2 de425_2 q425_2 w425_2 ds_frag425_2 de_frag425_2
        ds426_2 de426_2 q426_2 w426_2 ds_frag426_2 de_frag426_2;
   by gene_id;
run;

data pfam;
   set cc.pfam_dataset_for_jmp_v2;
   keep gene_id pfam_id_cat pfam_accession_cat  ;
run;

proc sort data=pfam nodup;
  by gene_id;
run;



proc import datafile="!PATCON/useful_human_data/aceview_hg19/downloaded_files/AceView.ncbi_37.pfamhits.txt"
   out=pfam_annot dbms=tab replace;
   guessingrows=195388;
run;

/* one domain per gene */

data pfam_annot2;
  set pfam_annot;
  keep __Gene pfam accession;
  rename __Gene=gene_id;
run;

proc sort data=pfam_annot2 nodup;
  by gene_id pfam;
run;

proc freq data=pfam_annot2 noprint;
  tables gene_id / out=pfam_per_gene;
proc sort data=pfam_per_gene;
  by descending count;
run;


data pfam_cat;
   array pf[13] $500.;
   array access[13] $500.;
   retain pf1-pf13;
   retain access1-access13;
   set pfam_annot2;
   by gene_id;
   if first.gene_id then do;
       call missing(of pf1-pf13); 
       call missing(of access1-access13); 
       records=0;
       end;
   records + 1;
   pf[records]=pfam;
   access[records]=accession;
   if last.gene_id then output;
run;

data pfam_cat2;
  set pfam_cat;
  length pfam_id_cat $3000.;
  length pfam_accession_cat $3000.;
  pfam_id_cat = catx("|", OF pf1-pf13);
  pfam_accession_cat = catx("|", OF access1-access13);
  keep pfam_id_cat pfam_accession_cat gene_id;
run;

proc sort data=pfam_cat2;
  by gene_id;
proc sort data=all_flags;
  by gene_id;
run;

data all_flags2;
  merge all_flags (in=in1) pfam_cat2 (in=in2);
  by gene_id;
  if in1;
run;

data all_flags3;
   set all_flags2;
   flag_425_any_xs=sum(flag_425_XS_DS,flag_425_XS_DE,flag_425_qrank,flag_425_wrank);
   flag_425_any_xs_ds=sum(flag_425_XS_DS,flag_425_qrank,flag_425_wrank);
   flag_425_any=sum(flag_425_XS_DS,flag_425_XS_DE,flag_425_qrank,flag_425_wrank,flag_425_frag_DS,flag_425_frag_DE);
   flag_425_any_ds=sum(flag_425_XS_DS,flag_425_qrank,flag_425_wrank,flag_425_frag_DS);
   flag_425_xs_wght_ds=sum(flag_425_XS_DS,flag_425_wrank);
   flag_425_xs_quant_ds=sum(flag_425_XS_DS,flag_425_qrank);
   flag_425_xs_rank_ds=sum(flag_425_XS_DS,flag_425_qrank,flag_425_wrank);

   flag_426_any_xs=sum(flag_426_XS_DS,flag_426_XS_DE,flag_426_qrank,flag_426_wrank);
   flag_426_any_xs_ds=sum(flag_426_XS_DS,flag_426_qrank,flag_426_wrank);
   flag_426_any=sum(flag_426_XS_DS,flag_426_XS_DE,flag_426_qrank,flag_426_wrank,flag_426_frag_DS,flag_426_frag_DE);
   flag_426_any_ds=sum(flag_426_XS_DS,flag_426_qrank,flag_426_wrank,flag_426_frag_DS);
   flag_426_xs_wght_ds=sum(flag_426_XS_DS,flag_426_wrank);
   flag_426_xs_quant_ds=sum(flag_426_XS_DS,flag_426_qrank);
   flag_426_xs_rank_ds=sum(flag_426_XS_DS,flag_426_qrank,flag_426_wrank);
run;
   

data cc.pfam_dataset_uniq_any;
  set all_flags3;
run;

