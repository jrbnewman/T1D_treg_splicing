/* Idea: In multi-transcript genes, there is (usually) only one major, functional isoform
   The rest of the transcripts are involved in regulating how much of that major isoform
   constitutes total gene expression. How to test?

   Case/control or control/treatment experiments: does the isoform switch?
   Multi-tissue expression: is the major isoform always the same?

   Look at case-control data first, and then pull out any interesting genes to make plots

 */

ods listing; ods html close;
libname cc '!PATCON/case_control/sas_data/fixed';
libname hg19 '!PATCON/useful_human_data/aceview_hg19/sas_data';

/* Rank transcripts within gene, within subject */

%macro rankXS(cell);

data xs2keep_&cell.;
   set cc.xscripts_multi_xs_genes_&cell.;
run;

data counts;
   set cc.tpm_counts_q3_norm_cent_all2;
   where cell_type=&cell.;
   /* Remove samples not used in DE/DS models */
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

   keep sample_id status transcript_id tpm;
run;

proc sort data=counts;
  by transcript_id;
proc sort data=xs2keep_&cell.;
  by transcript_id;
run;

data counts2;
  merge xs2keep_&cell. (in=in1) counts (in=in2);
  by transcript_id;
  if in1 and in2;
run;

data counts_tpm;
  set counts2;
  keep gene_id sample_id tpm;
run;

proc sort data=counts_tpm nodup;
  by sample_id gene_id descending tpm;
run;

data rank_tpm;
  retain xscript_rank;
  retain xscript_3lvl;
  set counts_tpm;
  by sample_id gene_id;
  if first.gene_id then do;
      xscript_rank=1;
      xscript_3lvl=1;
      flag_met=1;
      end;
  else do;
      if last.gene_id then xscript_3lvl=3;
      else xscript_3lvl=2;
      xscript_rank=xscript_rank+1;
      flag_met=0;
      end;
run;

proc sort data=rank_tpm;
  by sample_id gene_id tpm;
proc sort data=counts2;
  by sample_id gene_id tpm;
run;

data counts_w_rank_&cell.;
  merge counts2 (in=in1) rank_tpm (in=in2);
  by sample_id gene_id tpm;
  if in1 and in2;
run;

%mend;

%rankXS(425);
%rankXS(426);

/* Stack, make permenant */

data xs_counts_w_rank;
  set counts_w_rank_127 counts_w_rank_400
      counts_w_rank_425 counts_w_rank_426;
run;

data cc.xscript_rank_by_cell_subj_all;
   set xs_counts_w_rank;
run;


/* 3-level rank is confounded by lots of low-expressed transcripts "competing" to be lowest ranked transcript
   i.e. noise.


   Here I am implementing a weighed 3-level rank, where:
   Rank 1: Highest TPM per gene/sample, or within 1 SD of highest TPM
   Rank 3: Lowest TPM per gene/sample, or within 1 SD of lowest TPM
   Rank 2: Everything else

 */

ods listing; ods html close;
libname cc '!PATCON/case_control/sas_data/fixed';
libname hg19 '!PATCON/useful_human_data/aceview_hg19/sas_data';

/* Rank transcripts within gene, within subject */

%macro rankXS(cell);

data xs2keep_&cell.;
   set cc.xscripts_multi_xs_genes_&cell.;
run;

data counts;
   set cc.tpm_counts_q3_norm_cent_all2;
   where cell_type=&cell.;
   /* Remove samples not used in DE/DS models */
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

   keep sample_id status transcript_id tpm;
run;

proc sort data=counts;
  by transcript_id;
proc sort data=xs2keep_&cell.;
  by transcript_id;
run;

data counts2;
  merge xs2keep_&cell. (in=in1) counts (in=in2);
  by transcript_id;
  if in1 and in2;
run;

/* Calculate the standard deviation of TPM values within a subject */

proc sort data=counts2;
   by sample_id gene_id;
proc means data=counts2 noprint;
   by sample_id gene_id;
   var tpm;
   output out=tpm_sd_by_gene_subj(drop=_TYPE_ _FREQ_) stddev=tpm_sd max=tpm_max min=tpm_min;
run;


data counts3;
  merge counts2 (in=in1) tpm_sd_by_gene_subj (in=in2);
  by sample_id gene_id;
  if in1 and in2;
run;

/* Assign rank */

data counts_w_rank;
  set counts3;
  if tpm = tpm_max then xscript_rank=1;
  else if tpm > (tpm_max - tpm_sd) then xscript_rank=1;
  else if tpm = tpm_min then xscript_rank=3;
  else if tpm < (tpm_min + tpm_sd) then xscript_rank=3;
  else xscript_rank=2;

  if tpm = tpm_max then xscript_rank2=1;
  else if tpm > (tpm_max - tpm_sd/2) then xscript_rank2=1;
  else if tpm = tpm_min then xscript_rank2=3;
  else if tpm < (tpm_min + tpm_sd/2) then xscript_rank2=3;
  else xscript_rank2=2;

run;

ods listing; 
proc freq data=counts_w_rank;
  tables xscript_rank xscript_rank2 xscript_rank*xscript_rank2;
run;

data counts_w_rank_&cell.;
  set counts_w_rank;
run;

%mend;


%rankXS(425);
%rankXS(426);

/* Stack, make permenant */

data xs_counts_w_rank;
  set counts_w_rank_127 counts_w_rank_400
      counts_w_rank_425 counts_w_rank_426;
run;

data cc.xscript_wghtd_rank_by_cell_subj;
   set xs_counts_w_rank;
run;


