libname cc "!HOME/concannon/BRI/RNAseq/sas_data";


proc import datafile="!HOME/concannon/BRI/RNAseq/design_files/final_design_file_super_duper.csv"
   out=design_file dbms=csv replace;
   guessingrows=max;
run;

data cc.design_file_final;
   length condition $10.;
   set design_file;
   condition=catx("_",cell_type,status);
run;

