/* Import transcript counts using iterdataset */

ods listing; ods html close;
libname cc '!PATCON/case_control/sas_data/fixed';
libname hg19 '!PATCON/useful_human_data/aceview_hg19/sas_data';

/* For each cell type, flag transcripts as on/off at TPM > 0 at >=0.5 */

/* Get counts */

data counts;
  set cc.hg19_rsem_100perc_apn2_xscripts;
run;

/* Get design file to relate samples to conditions */

data design;
  set cc.design_file_new;
  keep sample_id status cell_type;
run;


%macro flagTPM(cell);
 
   data counts_roz;
      set counts;
      where cell_type = &cell.;
   run;

   data design_roz;
      set design;
      where cell_type = &cell.;
      drop cell_type;
   run;

   proc sort data=counts_roz;
      by sample_id;
   proc sort data=design_roz;
      by sample_id;
   run;
   
   data counts_w_key oops1_&cell. oops2_&cell. ;
     merge design_roz (in=in1) counts_roz (in=in2);
     by sample_id;
     if in1 and in2 then output counts_w_key;
     else if in1 then output oops1_&cell.;
     else output oops2_&cell.;
   run;
   
   %macro flagTPM2(tpm,tpmvar);

      data flag_tpm;
         set counts_w_key;
         %if &tpm. = 0 %then %do;
            if tpm > &tpm. then flag_tpm_gt&tpmvar. = 1; else flag_tpm_gt&tpmvar. = 0;
            %end;
         %else %do;
            if tpm >= &tpm. then flag_tpm_ge&tpmvar. = 1; else flag_tpm_ge&tpmvar. = 0;
            %end;
      run;

      /* Calculate the proportion of samples where TPM is gt/ge specified threshold */

      proc sort data=flag_tpm;
         by transcript_id status;
      proc means data=flag_tpm noprint;
         by transcript_id status;
         %if &tpm. = 0 %then %do;
            var flag_tpm_gt&tpmvar. ;
            %end;
         %else %do;
            var flag_tpm_ge&tpmvar. ;
            %end;
      output out=flag_xs_by_trt mean=trt_perc;

      /* If proportion is 50% or greater then flag transcript as on */

     data flag_trt_perc;
        set flag_xs_by_trt;
        if trt_perc >= 0.5 then flag_xscript_on=1;
        else if trt_perc = 0 then flag_xscript_on=0;
        else flag_xscript_on=.;
        keep transcript_id status flag_xscript_on;
     run;

     /* Put side-by-side */
     proc sort data=flag_trt_perc;
       by transcript_id status;
     proc transpose data=flag_trt_perc out=flag_trt_perc_sbys;
       by transcript_id;
       var flag_xscript_on;
       id status;
     run;

     data flag_xs_by_trt_tpm&tpmvar.;
        set flag_trt_perc_sbys;
        keep transcript_id T1D CTL;
        %if &tpm. = 0 %then %do;
           rename CTL=flag_&cell._CTL_tpm_gt0 T1D=flag_&cell._T1D_tpm_gt0;
           %end;
        %else %do;
           rename CTL=flag_&cell._CTL_tpm_ge&tpmvar. T1D=flag_&cell._T1D_tpm_ge&tpmvar.;
           %end;
     run;

   %mend;
   
   %flagTPM2(0,0);
   %flagTPM2(0.5,05);

   proc sort data=flag_xs_by_trt_tpm0;
      by transcript_id;
   proc sort data=flag_xs_by_trt_tpm05;
      by transcript_id;
   run;

   data flag_xs_by_trt_&cell.;
     merge flag_xs_by_trt_tpm0 (in=in1) flag_xs_by_trt_tpm05 (in=in2);
     by transcript_id;
     if in1 and in2;
   run;

%mend;

%flagTPM(127);
%flagTPM(400);
%flagTPM(425);
%flagTPM(426);

/* Make permenant */

data cc.flag_xscript_tpm_gt0_ge05_127;
   set flag_xs_by_trt_127;
run;


data cc.flag_xscript_tpm_gt0_ge05_400;
   set flag_xs_by_trt_400;
run;


data cc.flag_xscript_tpm_gt0_ge05_425;
   set flag_xs_by_trt_425;
run;


data cc.flag_xscript_tpm_gt0_ge05_426;
   set flag_xs_by_trt_426;
run;

/* Counts */

proc freq data=flag_xs_by_trt_127 noprint;
    tables flag_127_CTL_tpm_gt0 * flag_127_T1D_tpm_gt0 / out=xs_cnt_127;
    tables flag_127_CTL_tpm_ge05 * flag_127_T1D_tpm_ge05 / out=xs_cnt_127_05;
run;
proc freq data=flag_xs_by_trt_400 noprint;
    tables flag_400_CTL_tpm_gt0 * flag_400_T1D_tpm_gt0 / out=xs_cnt_400;
    tables flag_400_CTL_tpm_ge05 * flag_400_T1D_tpm_ge05 / out=xs_cnt_400_05;
run;
proc freq data=flag_xs_by_trt_425 noprint;
    tables flag_425_CTL_tpm_gt0 * flag_425_T1D_tpm_gt0 / out=xs_cnt_425;
    tables flag_425_CTL_tpm_ge05 * flag_425_T1D_tpm_ge05 / out=xs_cnt_425_05;
run;
proc freq data=flag_xs_by_trt_426 noprint;
    tables flag_426_CTL_tpm_gt0 * flag_426_T1D_tpm_gt0 / out=xs_cnt_426;
    tables flag_426_CTL_tpm_ge05 * flag_426_T1D_tpm_ge05 / out=xs_cnt_426_05;
run;

proc print data=xs_cnt_127;
proc print data=xs_cnt_127_05;
proc print data=xs_cnt_400;
proc print data=xs_cnt_400_05;
proc print data=xs_cnt_425;
proc print data=xs_cnt_425_05;
proc print data=xs_cnt_426;
proc print data=xs_cnt_426_05;
run; quit;

/*
127:
flag_127_    flag_127_
 CTL_tpm_     T1D_tpm_
   gt0          gt0       COUNT

    .            .           94
    .            0           28
    .            1          137	(case-specific)
    0            .           29
    0            0          271
    1            .            3	(control-specific)
    1            1        21319

 flag_127_    flag_127_
  CTL_tpm_     T1D_tpm_
    ge05         ge05      COUNT

     .            .          455
     .            0            2
     .            1          213	(case-specific)
     0            .           35
     0            0          731
     1            .           62	(control-specific)
     1            1        20481

400:
flag_400_    flag_400_
 CTL_tpm_     T1D_tpm_
   gt0          gt0       COUNT

    .            .          241
    .            0          106
    .            1            3	(case-specific)
    0            .           48
    0            0          450
    1            .            1	(control-specific)
    1            1        32495

 flag_400_    flag_400_
  CTL_tpm_     T1D_tpm_
    ge05         ge05      COUNT

     .            .         1213
     .            0            3
     .            1          263	(case-specific)
     0            .           12
     0            0          946
     1            .          105	(control-specific)
     1            1        31707


425:
flag_425_    flag_425_
 CTL_tpm_     T1D_tpm_
   gt0          gt0       COUNT

    .            .           68
    .            0           12
    .            1            5	(case-specific)
    0            .           71
    0            0          246
    1            .            6	(control-specific)
    1            1        20848

flag_425_    flag_425_
 CTL_tpm_     T1D_tpm_
   ge05         ge05      COUNT

    .            .         1652
    .            0            4
    .            1          417	(case-specific)
    0            .            5
    0            0          636
    1            .           56	(control-specific)
    1            1        19842



426:
flag_426_    flag_426_
 CTL_tpm_     T1D_tpm_
   gt0          gt0       COUNT

    .            .           90
    .            0           69
    .            1           63	(case-specific)
    0            .           27
    0            0          218
    1            .            9	(control-specific)
    1            1        19027

flag_426_    flag_426_
 CTL_tpm_     T1D_tpm_
   ge05         ge05      COUNT

    .            .          444
    .            0            1
    .            1           89	(case-specific)
    0            .           33
    0            0          593
    1            .           72	(control-specific)
    1            1        18336

*/


