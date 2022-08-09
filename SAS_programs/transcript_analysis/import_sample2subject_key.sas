/* Import fixed subject IDs */

ods listing; ods html close;
libname cc '!PATCON/case_control/sas_data/fixed';

proc import datafile="!PATCON/case_control/design_files/case_control_sample2subject_key_fixed.csv"
     out=sample2subject_key dbms=csv replace;
     guessingrows=max;
run;

/* Make permenant */

data cc.sample2subject_key;
   set sample2subject_key;
run;

