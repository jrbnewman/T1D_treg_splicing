/* Flag junctions on/off */

ods listing; ods html close;

libname cc '!PATCON/case_control/sas_data/fixed';
libname cclocal '/TB14/TB14/immigrans/store/cc_sandbox/sas_data';
libname hg19 '!PATCON/useful_human_data/aceview_hg19/sas_data';

data junc_counts_tall;
  set cclocal.junction_counts_by_sample;
run;


data design;
  set cc.design_file_new;
  keep sample_id status cell_type;
run;

%macro flagapn(cell);
 

   data design_roz;
      set design;
      where cell_type = &cell.;
      drop cell_type;
   run;

   proc sort data=junc_counts_tall;
      by sample_id;
   proc sort data=design_roz;
      by sample_id;
   run;

   data counts_w_key ;
     merge design_roz (in=in1) junc_counts_tall (in=in2);
     by sample_id;
     if in1 and in2 then output counts_w_key;
   run;
   
   %macro flagapn2(apn,apnvar);

      data flag_apn;
         set counts_w_key;
         %if &apn. = 0 %then %do;
            if apn > &apn. then flag_apn_gt&apnvar. = 1; else flag_apn_gt&apnvar. = 0;
            %end;
         %else %do;
            if apn >= &apn. then flag_apn_ge&apnvar. = 1; else flag_apn_ge&apnvar. = 0;
            %end;
      run;

      /* Calculate the proportion of samples where apn is gt/ge specified threshold */

      proc sort data=flag_apn;
         by event_id status;
      proc means data=flag_apn noprint;
         by event_id status;
         %if &apn. = 0 %then %do;
            var flag_apn_gt&apnvar. ;
            %end;
         %else %do;
            var flag_apn_ge&apnvar. ;
            %end;
      output out=flag_junc_by_trt mean=trt_perc;

      /* If proportion is 50% or greater then flag transcript as on */

     data flag_trt_perc;
        set flag_junc_by_trt;
        if trt_perc >= 0.5 then flag_junc_on=1;
        else if trt_perc = 0 then flag_junc_on=0;
        else flag_junc_on=.;
        keep event_id status flag_junc_on;
     run;

     /* Put side-by-side */
     proc sort data=flag_trt_perc;
       by event_id status;
     proc transpose data=flag_trt_perc out=flag_trt_perc_sbys;
       by event_id;
       var flag_junc_on;
       id status;
     run;

     data flag_junc_by_trt_apn&apnvar.;
        set flag_trt_perc_sbys;
        keep event_id T1D CTL;
        %if &apn. = 0 %then %do;
           rename CTL=flag_&cell._CTL_apn_gt0 T1D=flag_&cell._T1D_apn_gt0;
           %end;
        %else %do;
           rename CTL=flag_&cell._CTL_apn_ge&apnvar. T1D=flag_&cell._T1D_apn_ge&apnvar.;
           %end;
     run;

   %mend;
   
   %flagapn2(0,0);
   %flagapn2(2,2);
   %flagapn2(5,5);

   proc sort data=flag_junc_by_trt_apn0;
      by event_id;
   proc sort data=flag_junc_by_trt_apn2;
      by event_id;
   proc sort data=flag_junc_by_trt_apn5;
      by event_id;
   run;

   data flag_junc_by_trt_&cell.;
     merge flag_junc_by_trt_apn0 (in=in1) flag_junc_by_trt_apn2 (in=in2) flag_junc_by_trt_apn5 (in=in3);
     by event_id;
     if in1 and in2;
   run;

%mend;

%flagapn(127);
%flagapn(400);
%flagapn(425);
%flagapn(426);

/* Make permenant */

data cc.flag_junction_apn_gt0_ge2_127;
   set flag_junc_by_trt_127;
run;

data cc.flag_junction_apn_gt0_ge2_400;
   set flag_junc_by_trt_400;
run;

data cc.flag_junction_apn_gt0_ge2_425;
   set flag_junc_by_trt_425;
run;

data cc.flag_junction_apn_gt0_ge2_426;
   set flag_junc_by_trt_426;
run;
