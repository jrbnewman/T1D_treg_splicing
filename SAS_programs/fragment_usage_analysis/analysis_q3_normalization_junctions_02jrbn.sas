ods html close; ods listing;
libname cc '!PATCON/case_control/sas_data/fixed';
libname cclocal '/TB14/TB14/immigrans/store/cc_sandbox/sas_data';

data samples;
  set cc.design_file_new;
  keep sample_id ;
run;

data counts;
   set cclocal.junction_counts_by_sample;
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

   keep sample_id apn event_id;
run;

proc sort data=counts;
  by sample_id;
proc sort data=samples;
  by sample_id;
run;

data counts2;
  merge samples (in=in1) counts (in=in2);
  by sample_id;
  if in1 and in2;
run;


data counts3;
  set counts2;
  if apn > 0;
run;

proc datasets noprint;
  delete counts;
run;

/* quartiles */

proc sort data=counts3;
   by sample_id;
proc means data=counts3 noprint;
   by sample_id;
   var apn;
   output out=quartiles_apn0 q3=q3;
   run;

proc means data=quartiles_apn0 noprint;
            var q3;
            output out=q3_q3 q3=q3;
            run;

data _null_;
   set q3_q3;
   call symput('q3q3',q3);
run;

%put %q3q3.;

/*
Analysis Variable : q3

            Upper
         Quartile
     ------------
       33.8900000
     ------------
*/


proc datasets noprint;
  delete counts3;
run;


    proc sort data=counts2;
        by sample_id;
        run;

    proc sort data=quartiles_apn0;
        by sample_id;
        run;

/* merge and oops test */


    data counts_w_q3 oops_apn0;
        merge counts2 (in=in1) quartiles_apn0 (in=in2);
        by sample_id;
        if in1 and in2 then output counts_w_q3;
        else output oops_apn0; 
        run;


data missing;
  set oops_apn0;
  keep sample_id;
run;

proc sort data=missing nodup;
   by sample_id;
run;

data fragment_q3_norm_counts;
        set counts_w_q3;
        log_apn=log(apn + 1);
        * Q3 normalization;
        q3_apn=(apn/q3) * &q3q3.;
        q3_ff=&q3q3./q3;
        log_q3_apn=log(q3_apn + 1);
run;

proc datasets noprint;
   delete counts_w_q3 counts2;
run;


/* Export data for making plots */

*Get design file, as I need add in some covariates;

data sample_key;
  set cc.design_file_new;
  keep sample_id cell_type status plate;
run;

proc sort data= sample_key;
  by sample_id;
proc sort data=fragment_q3_norm_counts;
  by sample_id;
run;

data frag_counts_all2;
  merge sample_key (in=in1) fragment_q3_norm_counts (in=in2);
  by sample_id;
  if in1 and in2;
  if apn > 0 then output;
run;

proc export data=frag_counts_all2
     outfile="!PATCON/case_control/text_data/counts_by_junction_all_q3_norm_gt0_02jrbn.csv"
     dbms=csv replace;
run;

/* Make counts permanent */

data cclocal.junction_q3_norm_counts;
  set fragment_q3_norm_counts;
run;

data cclocal.junction_missing_counts;
  set missing;
run;


