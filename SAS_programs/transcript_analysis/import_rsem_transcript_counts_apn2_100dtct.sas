/* Import transcript counts using iterdataset */

ods listing; ods html close;
libname cc '!PATCON/case_control/sas_data/fixed';
libname hg19 '!PATCON/useful_human_data/aceview_hg19/sas_data';


%include "!PATCON/case_control/sas_programs/event_analysis/fixed/iterdataset.sas";


data libraries;
  set cc.design_file_new;
  keep sample_id cell_type;
run;

proc datasets noprint;
  delete rsem_: ;
run; quit;



%macro importRSEM(sample,cell);

           data WORK.&sample._est    ;
           %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
           infile "/home/jrbnewman/concannon/case_control/alignment_output/rsem_output_bri_cc/hg19_aceview_apn2_100prc_cell_&cell._&sample._.isoforms.results" delimiter='09'x MISSOVER DSD
       lrecl=32767 firstobs=2 ;
              informat transcript_id $105. ;
              informat gene_id $36. ;
              informat length best32. ;
              informat effective_length best32. ;
              informat expected_count best32. ;
              informat TPM best32. ;
              informat FPKM best32. ;
              informat IsoPct best32. ;
              informat posterior_mean_count best32. ;
              informat posterior_SD_count best32. ;
              informat pme_TPM best32. ;
              informat pme_FPKM best32. ;
              informat IsoPct_from_pme_TPM best32. ;
              informat TPM_ci_lower_bound best32. ;
              informat TPM_ci_upper_bound best32. ;
              informat TPM_CVquartile best32. ;
              informat FPKM_ci_lower_bound best32. ;
              informat FPKM_ci_upper_bound best32. ;
              informat FPKM_CVquartile best32. ;
             format transcript_id $105. ;
             format gene_id $36. ;
             format length best12. ;
             format effective_length best12. ;
             format expected_count best12. ;
             format TPM best12. ;
             format FPKM best12. ;
             format IsoPct best12. ;
             format posterior_mean_count best12. ;
             format posterior_SD_count best12. ;
             format pme_TPM best12. ;
             format pme_FPKM best12. ;
             format IsoPct_from_pme_TPM best12. ;
             format TPM_ci_lower_bound best12. ;
             format TPM_ci_upper_bound best12. ;
             format TPM_CVquartile best12. ;
             format FPKM_ci_lower_bound best12. ;
             format FPKM_ci_upper_bound best12. ;
             format FPKM_CVquartile best12. ;
          input
                      transcript_id $
                      gene_id $
                      length
                      effective_length
                      expected_count
                      TPM
                      FPKM
                      IsoPct
                      posterior_mean_count
                      posterior_SD_count
                      pme_TPM
                      pme_FPKM
                      IsoPct_from_pme_TPM
                      TPM_ci_lower_bound
                      TPM_ci_upper_bound
                      TPM_CVquartile
                      FPKM_ci_lower_bound
                      FPKM_ci_upper_bound
                      FPKM_CVquartile
       ;
       if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
       run;




*trim data -- only want TPM estimate and transcript_id;

data rsem_&sample.;
  length sample_id $32.;
  set &sample._est;
  sample_id="&sample.";
  cell_type=&cell.;
  keep transcript_id sample_id TPM cell_type;
run;



%mend;


   %iterdataset(dataset=libraries, function=%nrstr(%importRSEM(&sample_id.,&cell_type.);));

data all_xscript_est;
   set rsem_Sample_: ;
run;

/* Update IDs */

data xs2fasta;
   set hg19.hg19_aceview_xs2gene_fasta_index;
   keep transcript_id transcript_fasta_id;
   rename transcript_id=aceview_xs_id transcript_fasta_id=transcript_id;
run;

proc sort data=all_xscript_est;
   by transcript_id;
proc sort data=xs2fasta;
   by transcript_id;
run;

data all_xscript_est2;
  merge all_xscript_est (in=in1) xs2fasta (in=in2);
   by transcript_id;
  if in1 and in2;
run;

/* Make permenant */

data cc.hg19_rsem_100perc_apn2_xscripts;
   set all_xscript_est2;
   drop transcript_id;
   rename aceview_xs_id=transcript_id;
run;


