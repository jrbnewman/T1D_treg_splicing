/* Import MISO annotations and assign to events */

ods listing; ods html close;
libname cc '!HOME/concannon/case_control/sas_data/fixed';
libname hg19 '!HOME/concannon/useful_human_data/aceview_hg19/sas_data';
libname cclocal '/TB14/TB14/immigrans/store/cc_sandbox/sas_data';


proc datasets lib=work kill noprint;
run;
quit;



proc import datafile="/TB14/TB14/sandbox/miso_sandbox/miso_treg_genes/commonshortest/A3SS.treg.gff3"
     out=miso_a3ss dbms=tab replace;
     guessingrows=max;
     getnames=no;
run;

proc import datafile="/TB14/TB14/sandbox/miso_sandbox/miso_treg_genes/commonshortest/A5SS.treg.gff3"
     out=miso_a5ss dbms=tab replace;
     guessingrows=max;
     getnames=no;
run;

proc import datafile="/TB14/TB14/sandbox/miso_sandbox/miso_treg_genes/commonshortest/MXE.treg.gff3"
     out=miso_mxe dbms=tab replace;
     guessingrows=max;
     getnames=no;
run;

proc import datafile="/TB14/TB14/sandbox/miso_sandbox/miso_treg_genes/commonshortest/RI.treg.gff3"
     out=miso_ri dbms=tab replace;
     guessingrows=max;
     getnames=no;
run;
proc import datafile="/TB14/TB14/sandbox/miso_sandbox/miso_treg_genes/commonshortest/SE.treg.gff3"
     out=miso_se dbms=tab replace;
     guessingrows=max;
     getnames=no;
run;




/* Get junction annotations,
       Exon skipping annotations
       Exon annotations
       exon fragment annotations */


proc import datafile="!HOME/concannon/case_control/references/hg19_aceview/annotations/hg19_aceview_exon_annotations.csv"
    out=exon_annot dbms=csv replace;
    guessingrows=max;
run;

proc import datafile="!HOME/concannon/case_control/references/hg19_aceview/annotations/hg19_aceview_exon_fragment_annotations.csv"
    out=fragment_annot dbms=csv replace;
    guessingrows=max;
run;


proc import datafile="!HOME/concannon/case_control/references/hg19_aceview/annotations/hg19_aceview_junctions_full_annotation.csv"
    out=junc_annot dbms=csv replace;
    guessingrows=max;
run;



/* Reformat A3SS entres */

data index_A3SS;
   set miso_a3ss;
   where var3="gene";
   length annotation $200.;
   length miso_exon_id $200.;
   length exon_assignment $20.;
   annotation=compress(tranwrd(scan(VAR9,1,";"),"ID=",""));
   do i = 1 by 1 while(scan(annotation,i,"@") ^= "");
      miso_exon_id=compress(scan(annotation,i,"@"));
      if i = 1 then exon_assignment="common";
      else exon_assignment="AB";
      output;
      end;
   keep annotation miso_exon_id exon_assignment;
run;

data index_A3SS_2;
  set index_A3SS;
  format chr $4.;
  format exon_start best12.;
  format exon_stop best12.;
  format strand $1.;
  format exon_assignment2 $20.;
  if exon_assignment="common" then do;
      chr=compress(scan(miso_exon_id,1,":"));
      exon_start=scan(miso_exon_id,2,":")  - 1 ; *exon data is stored as 0-based;
      exon_stop=scan(miso_exon_id,3,":") + 0;
      strand=compress(scan(miso_exon_id,4,":"));
      exon_assignment2="common";
          output;
      end;
  else do;
      chr=compress(scan(miso_exon_id,1,":"));
      exon_stop=scan(miso_exon_id,3,":") + 0;
      strand=compress(scan(miso_exon_id,4,":"));
      do i = 1 by 1 while (scan(scan(miso_exon_id,2,":"),i,"|") ^= "" );
          exon_start=scan(scan(miso_exon_id,2,":"),i,"|")  - 1 ; *exon data is stored as 0-based;
          if i = 1 then exon_assignment2="A";
          else exon_assignment2="B";
          output;
          end;
       end;
run;


/* Reformat A5SS entres */

data index_A5SS;
   set miso_a5ss;
   where var3="gene";
   length annotation $200.;
   length miso_exon_id $200.;
   length exon_assignment $20.;
   annotation=compress(tranwrd(scan(VAR9,1,";"),"ID=",""));
   do i = 1 by 1 while(scan(annotation,i,"@") ^= "");
      miso_exon_id=compress(scan(annotation,i,"@"));
      if i = 1 then exon_assignment="AB";
      else exon_assignment="common";
      output;
      end;
   keep annotation miso_exon_id exon_assignment;
run;

data index_A5SS_2;
  set index_A5SS;
  format chr $4.;
  format exon_start best12.;
  format exon_stop best12.;
  format strand $1.;
  format exon_assignment2 $20.;
  if exon_assignment="common" then do;
      chr=compress(scan(miso_exon_id,1,":"));
      exon_start=scan(miso_exon_id,2,":")  - 1 ; *exon data is stored as 0-based;
      exon_stop=scan(miso_exon_id,3,":") + 0;
      strand=compress(scan(miso_exon_id,4,":"));
      exon_assignment2="common";
          output;
      end;
  else do;
      chr=compress(scan(miso_exon_id,1,":"));
      exon_start=scan(miso_exon_id,2,":")  - 1 ; *exon data is stored as 0-based;
      strand=compress(scan(miso_exon_id,4,":"));
      do i = 1 by 1 while (scan(scan(miso_exon_id,3,":"),i,"|") ^= "" );
          exon_stop=scan(scan(miso_exon_id,3,":"),i,"|") + 0;
          if i = 1 then exon_assignment2="A";
          else exon_assignment2="B";
          output;
          end;
       end;
run;


data index_MXE;
   set miso_MXE;
   where var3="gene";
   length annotation $200.;
   length miso_exon_id $200.;
   length exon_assignment $20.;
   annotation=compress(tranwrd(scan(VAR9,1,";"),"ID=",""));
   do i = 1 by 1 while(scan(annotation,i,"@") ^= "");
      miso_exon_id=compress(scan(annotation,i,"@"));
      if i = 1 then exon_assignment="common1";
      else if i = 2 then exon_assignment="A";
      else if i = 3 then exon_assignment="B";
      else exon_assignment="common2";
      output;
      end;
   keep annotation miso_exon_id exon_assignment;
run;

data index_MXE_2;
  set index_MXE;
  format chr $4.;
  format exon_start best12.;
  format exon_stop best12.;
  format strand $1.;
  format exon_assignment2 $20.;
  chr=compress(scan(miso_exon_id,1,":"));
  exon_start=scan(miso_exon_id,2,":") - 1 ; *exon data is stored as 0-based;
  exon_stop=scan(miso_exon_id,3,":") + 0;
  strand=compress(scan(miso_exon_id,4,":"));
  exon_assignment2=exon_assignment;
run;

data index_SE;
   set miso_SE;
   where var3="gene";
   length annotation $200.;
   length miso_exon_id $200.;
   length exon_assignment $20.;
   annotation=compress(tranwrd(scan(VAR9,1,";"),"ID=",""));
   do i = 1 by 1 while(scan(annotation,i,"@") ^= "");
      miso_exon_id=compress(scan(annotation,i,"@"));
      if i = 1 then exon_assignment="common1";
      else if i = 2 then exon_assignment="A";
      else exon_assignment="common2";
      output;
      end;
   keep annotation miso_exon_id exon_assignment;
run;

data index_SE_2;
  set index_SE;
  format chr $4.;
  format exon_start best12.;
  format exon_stop best12.;
  format strand $1.;
  format exon_assignment2 $20.;
  chr=compress(scan(miso_exon_id,1,":"));
  exon_start=scan(miso_exon_id,2,":")  - 1 ; *exon data is stored as 0-based;
  exon_stop=scan(miso_exon_id,3,":") + 0;
  strand=compress(scan(miso_exon_id,4,":"));
  exon_assignment2=exon_assignment;
run;



data index_RI;
   set miso_RI;
   where var3="gene";
   length annotation $200.;
   length miso_exon_id $200.;
   length exon_assignment $20.;
   annotation=compress(tranwrd(scan(VAR9,1,";"),"ID=",""));
   do i = 1 by 1 while(scan(annotation,i,"@") ^= "");
      miso_exon_id=compress(scan(annotation,i,"@"));
      if i = 1 then exon_assignment="common1";
      else exon_assignment="common2";
      output;
      end;
   keep annotation miso_exon_id exon_assignment;
run;

data index_RI_2;
  set index_RI;
  format chr $4.;
  format exon_start best12.;
  format exon_stop best12.;
  format strand $1.;
  format exon_assignment2 $20.;
  chr=compress(scan(miso_exon_id,1,":"));
  exon_start=scan(scan(miso_exon_id,2,":"),1,"-")  - 1 ; *exon data is stored as 0-based;
  exon_stop=scan(scan(miso_exon_id,2,":"),2,"-") + 0;
  strand=compress(scan(miso_exon_id,3,":"));
  exon_assignment2=exon_assignment;
run;

/* Make a master index */

data miso_index;
    set index_a3ss_2 (in=in1) index_a5ss_2 (in=in2)
        index_mxe_2 (in=in3) index_se_2 (in=in4) index_ri_2 (in=in5);
    length event_type $4.;
    if in1 then event_type="A3SS";
    if in2 then event_type="A5SS";
    if in3 then event_type="MXE";
    if in4 then event_type="SE";
    if in5 then event_type="RI";
    if exon_start > exon_stop then do;
       exon_stop2=exon_start+1;
       exon_start2=exon_stop-1;
       end;
   else do;
       exon_stop2=exon_stop;
       exon_start2=exon_start;
       end;
   drop exon_start exon_stop;
   rename exon_start2=exon_start exon_stop2=exon_stop;
run;

/* take the first exon and map to a gene */

proc sort data=miso_index;
   by event_type annotation chr exon_start exon_stop strand;
run;

data exon2gene;
  set exon_annot;
  keep gene_id chrom start stop strand;
  rename chrom=chr start=exon_start stop=exon_stop;
run;

proc sort data=miso_index;
   by chr exon_start exon_stop strand;
proc sort data=exon2gene;
   by chr exon_start exon_stop strand;
run;


data miso_index_w_gene no_gene;
  merge miso_index (in=in1) exon2gene (in=in2);
  by chr exon_start exon_stop strand;
  if in1 and in2 then output miso_index_w_gene;
  else if in1 then output no_gene;
run;

    
data annot2gene;
  set miso_index_w_gene;
   keep gene_id annotation;
run;

proc sort data=annot2gene nodup;
  by annotation gene_id;
proc freq data=annot2gene noprint;
   tables annotation / out=gene_per_annot;
run;
data check;
  set gene_per_annot;
  where count > 1;
run;


proc sort data=annot2gene;
  by annotation;
proc sort data=miso_index;
  by annotation;
run;

data miso_index2;
  merge annot2gene (in=in1) miso_index (in=in2);
  by annotation;
  if in1 and in2;
  drop exon_assignment i;
run;


/*  Inclusion/exclusion events:

A3SS:
Inclusion (A) shorter acceptor exon: annotated junction from exon_common_stop to exon_B_start
Exclusion (B)  longer acceptor exon:  fragment between exon_A_start and exon_B_start and annotated junction from exon_common_stop to exon_A_start

*/


data A3SS_A_donor;
    set miso_index2;
    where ((event_type="A3SS" and strand="+") or (event_type="A5SS" and strand="-")) and exon_assignment2="common";
    drop exon_assignment2 exon_start  miso_exon_id;
    rename exon_stop=donor_stop;
run;

data A3SS_A_acceptor;
    set miso_index2;
    where ((event_type="A3SS" and strand="+") or (event_type="A5SS" and strand="-")) and exon_assignment2="A";
    drop exon_assignment2 exon_stop miso_exon_id;
    rename exon_start=acceptor_start ;
run;

data A3SS_B_donor;
    set miso_index2;
    where ((event_type="A3SS" and strand="+") or (event_type="A5SS" and strand="-")) and exon_assignment2="common";
    drop exon_assignment2 exon_start  miso_exon_id;
    rename exon_stop=donor_stop;
run;

data A3SS_B_acceptor;
    set miso_index2;
    where ((event_type="A3SS" and strand="+") or (event_type="A5SS" and strand="-")) and exon_assignment2="B";
    drop exon_assignment2 exon_stop  miso_exon_id;
    rename exon_start=acceptor_start;
run;

data A3SS_B_frag_start;
    set miso_index2;
    where ((event_type="A3SS" and strand="+") or (event_type="A5SS" and strand="-")) and exon_assignment2="A";
    drop exon_assignment2 exon_stop  miso_exon_id;
    rename exon_start=fragment_range_start;
run;

data A3SS_B_frag_stop;
    set miso_index2;
    where((event_type="A3SS" and strand="+") or (event_type="A5SS" and strand="-")) and exon_assignment2="B";
    drop exon_assignment2 exon_stop  miso_exon_id;
    rename exon_start=fragment_range_stop;
run;

proc sort data=A3SS_A_donor;
   by annotation gene_id chr strand event_type;
proc sort data=A3SS_A_acceptor;
   by annotation gene_id chr strand event_type;

proc sort data=A3SS_B_donor;
   by annotation gene_id chr strand event_type;
proc sort data=A3SS_B_acceptor;
   by annotation gene_id chr strand event_type;

proc sort data=A3SS_B_frag_start;
   by annotation gene_id chr strand event_type;
proc sort data=A3SS_B_frag_stop;
   by annotation gene_id chr strand event_type;
run;


data A3SS_A_junction;
   merge A3SS_A_donor (in=in1) A3SS_A_acceptor (in=in2);
   by annotation gene_id chr strand event_type;
   if in1 and in2;
run;

data A3SS_B_junction;
   merge A3SS_B_donor (in=in1) A3SS_B_acceptor (in=in2);
   by annotation gene_id chr strand event_type;
   if in1 and in2;
run;

data A3SS_B_frag2get;
   merge A3SS_B_frag_start (in=in1) A3SS_B_frag_stop (in=in2);
   by annotation gene_id chr strand event_type;
   if in1 and in2;
run;


/*

A5SS:
Inclusion (A) longer donor exon:  fragment between exon_A_stop and exon_B_stop and annotated junction from exon_B to exon_common
Exclusion (B) shorter donor exon:  annotated junction from exon_A to exon_common

*/

data A5SS_A_donor;
    set miso_index2;
    where ((event_type="A3SS" and strand="-") or (event_type="A5SS" and strand="+")) and exon_assignment2="A";
    drop exon_assignment2 exon_start  miso_exon_id;
    rename exon_stop=donor_stop;
run;

data A5SS_A_acceptor;
    set miso_index2;
    where ((event_type="A3SS" and strand="-") or (event_type="A5SS" and strand="+")) and exon_assignment2="common";
    drop exon_assignment2 exon_stop miso_exon_id;
    rename exon_start=acceptor_start ;
run;

data A5SS_B_donor;
    set miso_index2;
    where ((event_type="A3SS" and strand="-") or (event_type="A5SS" and strand="+")) and exon_assignment2="B";
    drop exon_assignment2 exon_start  miso_exon_id;
    rename exon_stop=donor_stop;
run;

data A5SS_B_acceptor;
    set miso_index2;
    where ((event_type="A3SS" and strand="-") or (event_type="A5SS" and strand="+")) and exon_assignment2="common";
    drop exon_assignment2 exon_stop  miso_exon_id;
    rename exon_start=acceptor_start;
run;

data A5SS_B_frag_start;
    set miso_index2;
    where ((event_type="A3SS" and strand="-") or (event_type="A5SS" and strand="+")) and exon_assignment2="A";
    drop exon_assignment2 exon_start  miso_exon_id;
    rename exon_stop=fragment_range_start;
run;

data A5SS_B_frag_stop;
    set miso_index2;
    where((event_type="A3SS" and strand="-") or (event_type="A5SS" and strand="+")) and exon_assignment2="B";
    drop exon_assignment2 exon_start  miso_exon_id;
    rename exon_stop=fragment_range_stop;
run;

proc sort data=A5SS_A_donor;
   by annotation gene_id chr strand event_type;
proc sort data=A5SS_A_acceptor;
   by annotation gene_id chr strand event_type;

proc sort data=A5SS_B_donor;
   by annotation gene_id chr strand event_type;
proc sort data=A5SS_B_acceptor;
   by annotation gene_id chr strand event_type;

proc sort data=A5SS_B_frag_start;
   by annotation gene_id chr strand event_type;
proc sort data=A5SS_B_frag_stop;
   by annotation gene_id chr strand event_type;
run;


data A5SS_A_junction;
   merge A5SS_A_donor (in=in1) A5SS_A_acceptor (in=in2);
   by annotation gene_id chr strand event_type;
   if in1 and in2;
run;

data A5SS_B_junction;
   merge A5SS_B_donor (in=in1) A5SS_B_acceptor (in=in2);
   by annotation gene_id chr strand event_type;
   if in1 and in2;
run;

data A5SS_B_frag2get;
   merge A5SS_B_frag_start (in=in1) A5SS_B_frag_stop (in=in2);
   by annotation gene_id chr strand event_type;
   if in1 and in2;
run;

/*

MXE:
Inclusion (A) reads: fragments between exon_B_start and exon_B_stop + annotated junctions between common1/2 and exon_B
Exclusion (B) reads: fragments between exon_A_start and exon_A_stop + annotated junctions between common1/2 and exon_A

*/

data MXE_common_donor_plus;
    set miso_index2;
    where event_type="MXE" and strand="+" and exon_assignment2="common1" ;
    drop exon_assignment2 exon_start miso_exon_id;
    rename exon_stop=donor_stop;
run;

data MXE_A1_acceptor_plus;
    set miso_index2;
    where event_type="MXE" and strand="+" and exon_assignment2="A" ;
    drop exon_assignment2 exon_stop miso_exon_id;
    rename exon_start=acceptor_start;
run;

data MXE_A2_donor_plus;
    set miso_index2;
    where event_type="MXE" and strand="+" and exon_assignment2="A" ;
    drop exon_assignment2 exon_start miso_exon_id;
    rename exon_stop=donor_stop;
run;

data MXE_common_acceptor_plus;
    set miso_index2;
    where event_type="MXE" and strand="+" and exon_assignment2="common2" ;
    drop exon_assignment2 exon_stop miso_exon_id;
    rename exon_start=acceptor_start;
run;

data MXE_common_donor_minus;
    set miso_index2;
    where event_type="MXE" and strand="-" and exon_assignment2="common2" ;
    drop exon_assignment2 exon_start miso_exon_id;
    rename exon_stop=donor_stop;
run;

data MXE_A1_acceptor_minus;
    set miso_index2;
    where event_type="MXE" and strand="-" and exon_assignment2="A" ;
    drop exon_assignment2 exon_stop miso_exon_id;
    rename exon_start=acceptor_start;
run;

data MXE_A2_donor_minus;
    set miso_index2;
    where event_type="MXE" and strand="-" and exon_assignment2="A" ;
    drop exon_assignment2 exon_start miso_exon_id;
    rename exon_stop=donor_stop;
run;

data MXE_common_acceptor_minus;
    set miso_index2;
    where event_type="MXE" and strand="-" and exon_assignment2="common1" ;
    drop exon_assignment2 exon_stop miso_exon_id;
    rename exon_start=acceptor_start;
run;



data MXE_B1_acceptor_plus;
    set miso_index2;
    where event_type="MXE" and strand="+" and exon_assignment2="B" ;
    drop exon_assignment2 exon_stop miso_exon_id;
    rename exon_start=acceptor_start;
run;

data MXE_B2_donor_plus;
    set miso_index2;
    where event_type="MXE" and strand="+" and exon_assignment2="B" ;
    drop exon_assignment2 exon_start miso_exon_id;
    rename exon_stop=donor_stop;
run;

data MXE_B1_acceptor_minus;
    set miso_index2;
    where event_type="MXE" and strand="-" and exon_assignment2="B" ;
    drop exon_assignment2 exon_stop miso_exon_id;
    rename exon_start=acceptor_start;
run;

data MXE_B2_donor_minus;
    set miso_index2;
    where event_type="MXE" and strand="-" and exon_assignment2="B" ;
    drop exon_assignment2 exon_start miso_exon_id;
    rename exon_stop=donor_stop;
run;

proc sort data=MXE_common_donor_plus;
   by annotation gene_id chr strand event_type;
proc sort data=MXE_common_acceptor_plus;
   by annotation gene_id chr strand event_type;
proc sort data=MXE_A2_donor_plus;
   by annotation gene_id chr strand event_type;
proc sort data=MXE_A1_acceptor_plus;
   by annotation gene_id chr strand event_type;
proc sort data=MXE_B2_donor_plus;
   by annotation gene_id chr strand event_type;
proc sort data=MXE_B1_acceptor_plus;
   by annotation gene_id chr strand event_type;

proc sort data=MXE_common_donor_minus;
   by annotation gene_id chr strand event_type;
proc sort data=MXE_common_acceptor_minus;
   by annotation gene_id chr strand event_type;
proc sort data=MXE_A2_donor_minus;
   by annotation gene_id chr strand event_type;
proc sort data=MXE_A1_acceptor_minus;
   by annotation gene_id chr strand event_type;
proc sort data=MXE_B2_donor_minus;
   by annotation gene_id chr strand event_type;
proc sort data=MXE_B1_acceptor_minus;
   by annotation gene_id chr strand event_type;
run;

data MXE_A1_junction_plus;
  merge MXE_common_donor_plus (in=in1) MXE_A1_acceptor_plus (in=in2);
  by annotation gene_id chr strand event_type;
  if in1 and in2;
run;

data MXE_B1_junction_plus;
  merge MXE_common_donor_plus (in=in1) MXE_B1_acceptor_plus (in=in2);
  by annotation gene_id chr strand event_type;
  if in1 and in2;
run;

data MXE_A2_junction_plus;
  merge MXE_A2_donor_plus (in=in1) MXE_common_acceptor_plus (in=in2);
  by annotation gene_id chr strand event_type;
  if in1 and in2;
run;

data MXE_B2_junction_plus;
  merge MXE_B2_donor_plus (in=in1) MXE_common_acceptor_plus (in=in2);
  by annotation gene_id chr strand event_type;
  if in1 and in2;
run;

data MXE_A1_junction_minus;
  merge MXE_common_donor_minus (in=in1) MXE_A1_acceptor_minus (in=in2);
  by annotation gene_id chr strand event_type;
  if in1 and in2;
run;

data MXE_B1_junction_minus;
  merge MXE_common_donor_minus (in=in1) MXE_B1_acceptor_minus (in=in2);
  by annotation gene_id chr strand event_type;
  if in1 and in2;
run;

data MXE_A2_junction_minus;
  merge MXE_A2_donor_minus (in=in1) MXE_common_acceptor_minus (in=in2);
  by annotation gene_id chr strand event_type;
  if in1 and in2;
run;

data MXE_B2_junction_minus;
  merge MXE_B2_donor_minus (in=in1) MXE_common_acceptor_minus (in=in2);
  by annotation gene_id chr strand event_type;
  if in1 and in2;
run;


data MXE_A_frag2get;
    set miso_index2;
    where event_type="MXE" and exon_assignment2="A" ;
    drop exon_assignment2  miso_exon_id;
    rename exon_start=fragment_range_start
           exon_stop=fragment_range_stop;
run;


data MXE_B_frag2get;
    set miso_index2;
    where event_type="MXE" and exon_assignment2="B" ;
    drop exon_assignment2  miso_exon_id;
    rename exon_start=fragment_range_start
           exon_stop=fragment_range_stop;
run;

/*

SE:
Inclusion reads: fragments from exon_B_start to exon_A_stop + annotated junctions from common1 to exon A and common2 to exon A
Exclusion reads: annotated junctions excluding exon_A
*/



data SE_common_donor_plus;
    set miso_index2;
    where event_type="SE" and strand="+" and exon_assignment2="common1" ;
    drop exon_assignment2 exon_start miso_exon_id;
    rename exon_stop=donor_stop;
run;

data SE_A1_acceptor_plus;
    set miso_index2;
    where event_type="SE" and strand="+" and exon_assignment2="A" ;
    drop exon_assignment2 exon_stop miso_exon_id;
    rename exon_start=acceptor_start;
run;

data SE_A2_donor_plus;
    set miso_index2;
    where event_type="SE" and strand="+" and exon_assignment2="A" ;
    drop exon_assignment2 exon_start miso_exon_id;
    rename exon_stop=donor_stop;
run;

data SE_common_acceptor_plus;
    set miso_index2;
    where event_type="SE" and strand="+" and exon_assignment2="common2" ;
    drop exon_assignment2 exon_stop miso_exon_id;
    rename exon_start=acceptor_start;
run;


data SE_common_donor_minus;
    set miso_index2;
    where event_type="SE" and strand="-" and exon_assignment2="common2" ;
    drop exon_assignment2 exon_start miso_exon_id;
    rename exon_stop=donor_stop;
run;

data SE_A1_acceptor_minus;
    set miso_index2;
    where event_type="SE" and strand="-" and exon_assignment2="A" ;
    drop exon_assignment2 exon_stop miso_exon_id;
    rename exon_start=acceptor_start;
run;

data SE_A2_donor_minus;
    set miso_index2;
    where event_type="SE" and strand="-" and exon_assignment2="A" ;
    drop exon_assignment2 exon_start miso_exon_id;
    rename exon_stop=donor_stop;
run;

data SE_common_acceptor_minus;
    set miso_index2;
    where event_type="SE" and strand="-" and exon_assignment2="common1" ;
    drop exon_assignment2 exon_stop miso_exon_id;
    rename exon_start=acceptor_start;
run;

proc sort data=SE_common_donor_plus;
   by annotation gene_id chr strand event_type;
proc sort data=SE_common_acceptor_plus;
   by annotation gene_id chr strand event_type;
proc sort data=SE_A2_donor_plus;
   by annotation gene_id chr strand event_type;
proc sort data=SE_A1_acceptor_plus;
   by annotation gene_id chr strand event_type;
proc sort data=SE_common_donor_minus;
   by annotation gene_id chr strand event_type;
proc sort data=SE_common_acceptor_minus;
   by annotation gene_id chr strand event_type;
proc sort data=SE_A2_donor_minus;
   by annotation gene_id chr strand event_type;
proc sort data=SE_A1_acceptor_minus;
   by annotation gene_id chr strand event_type;
run;

data SE_A1_junction_plus;
  merge SE_common_donor_plus (in=in1) SE_A1_acceptor_plus (in=in2);
  by annotation gene_id chr strand event_type;
   if in1 and in2;
run;

data SE_A2_junction_plus;
  merge SE_A2_donor_plus (in=in1) SE_common_acceptor_plus (in=in2);
  by annotation gene_id chr strand event_type;
   if in1 and in2;
run;

data SE_B_junction_plus;
  merge SE_common_donor_plus (in=in1) SE_common_acceptor_plus (in=in2);
  by annotation gene_id chr strand event_type;
   if in1 and in2;
run;

data SE_A1_junction_minus;
  merge SE_common_donor_minus (in=in1) SE_A1_acceptor_minus (in=in2);
  by annotation gene_id chr strand event_type;
   if in1 and in2;
run;

data SE_A2_junction_minus;
  merge SE_A2_donor_minus (in=in1) SE_common_acceptor_minus (in=in2);
  by annotation gene_id chr strand event_type;
   if in1 and in2;
run;

data SE_B_junction_minus;
  merge SE_common_donor_minus (in=in1) SE_common_acceptor_minus (in=in2);
  by annotation gene_id chr strand event_type;
   if in1 and in2;
run;


data SE_A_frag2get;
    set miso_index2;
    where event_type="SE" and exon_assignment2="A" ;
    drop exon_assignment2  miso_exon_id;
    rename exon_start=fragment_range_start
           exon_stop=fragment_range_stop;
run;

/*
RI:
Inclusion events: fragments between exon_A1_stop and exon_A2_start
Exclusion events: junction from exon_A1_stop to exon_A2_start
*/

data RI_B_donor_plus;
   set miso_index2;
   where event_type="RI" and strand="+" and exon_assignment2="common1";
   drop exon_assignment2 miso_exon_id exon_start;
   rename exon_stop=donor_stop;
run;

data RI_B_donor_minus;
   set miso_index2;
   where event_type="RI" and strand="-" and exon_assignment2="common2";
   drop exon_assignment2 miso_exon_id exon_start;
   rename exon_stop=donor_stop;
run;


data RI_B_acceptor_plus;
   set miso_index2;
   where event_type="RI" and strand="+" and exon_assignment2="common2";
   drop exon_assignment2 miso_exon_id exon_stop;
   rename exon_start=acceptor_start;
run;

data RI_B_acceptor_minus;
   set miso_index2;
   where event_type="RI" and strand="-" and exon_assignment2="common1";
   drop exon_assignment2 miso_exon_id exon_stop;
   rename exon_start=acceptor_start;
run;

proc sort data=RI_B_donor_plus;
  by annotation gene_id chr strand;
proc sort data=RI_B_acceptor_plus;
  by annotation gene_id chr strand;
proc sort data=RI_B_donor_minus;
  by annotation gene_id chr strand;
proc sort data=RI_B_acceptor_minus;
  by annotation gene_id chr strand;
run;

data RI_B_junction_plus;
   merge RI_B_donor_plus (in=in1) RI_B_acceptor_plus (in=in2);
   by annotation gene_id chr strand;
   if in1 and in2;
run;

data RI_B_junction_minus;
   merge RI_B_donor_minus (in=in1) RI_B_acceptor_minus (in=in2);
   by annotation gene_id chr strand;
   if in1 and in2;
run;

data RI_A_frag2get;
   set RI_B_junction_plus RI_B_junction_minus;
   rename donor_stop=fragment_range_start
          acceptor_start=fragment_range_stop;
run;

/* Miso to Event indices
If it includes FRAGMENTS/EXONS then it it INCLUSION

if MXE, INCLUSION is the first exon, EXCLUSION is the second exon

*/

data miso2ea_frag;
   set a3ss_b_frag2get (in=in1) a5ss_b_frag2get (in=in1)
       mxe_a_frag2get (in=in1) mxe_b_frag2get (in=in2)
       RI_A_Frag2get (in=in1) se_a_frag2get (in=in1);
       length exc_inc $1.;
       if in1 then exc_inc="I";
       if in2 then exc_inc="E";
run;

data miso2ea_junc;
   set a3ss_a_junction (in=in2) a3ss_b_junction (in=in1)
       a5ss_a_junction (in=in2) a5ss_b_junction (in=in1)
       mxe_a1_junction_plus (in=in1) mxe_a1_junction_minus (in=in1)
       mxe_a2_junction_plus (in=in1) mxe_a2_junction_minus (in=in1)
       mxe_b1_junction_plus (in=in2) mxe_b1_junction_minus (in=in2)
       mxe_b2_junction_plus (in=in2) mxe_b2_junction_minus (in=in2) 
       ri_b_junction_plus (in=in2) ri_b_junction_minus (in=in2)
       se_a1_junction_plus (in=in1) se_a1_junction_minus (in=in1)
       se_a2_junction_plus (in=in1) se_a2_junction_minus (in=in1)
       se_b_junction_plus (in=in2) se_b_junction_minus (in=in2);
       length exc_inc $1.;
       if in1 then exc_inc="I";
       if in2 then exc_inc="E";
run;  

proc sort data=miso2ea_junc;
   by chr donor_stop acceptor_start;
proc sort data=junc_annot;
  by chr donor_stop acceptor_start;
run;

data junc_annot2;
  set junc_annot;
  keep chr donor_stop acceptor_Start;
run;

proc sort data=junc_annot2 nodup;
  by chr donor_stop acceptor_start;
run;


data miso_bad;
  merge miso2ea_junc (in=in1) junc_annot2 (in=in2);
   by chr donor_stop acceptor_start;
  if in1 and not in2 ;
  keep annotation ;
run;

proc sort data=miso2ea_junc;
   by annotation;
proc sort data=miso2ea_frag;
   by annotation;
proc sort data=miso_bad nodup;
   by annotation;
run;

data miso2ea_junc2;
  merge miso2ea_junc (in=in1) miso_bad (in=in2);
  by annotation;
  if in1 and not in2;
run;

data miso2ea_frag2;
  merge miso2ea_frag (in=in1) miso_bad (in=in2);
  by annotation;
  if in1 and not in2;
run;

/* For junctions: merge in junction ID and transcript ID */

data junc_annot3;
  set junc_annot;
   keep junction_id chr donor_stop acceptor_start transcript_id;
run;


proc sort data=miso2ea_junc2;
   by chr donor_stop acceptor_start;
proc sort data=junc_annot3;
   by chr donor_stop acceptor_start;
run;

data miso2ea_junc_w_iso;
  merge miso2ea_junc2 (in=in1) junc_annot3 (in=in2);
  by chr donor_stop acceptor_start;
  if in1 and in2;
run;


/* For fragments: iterate through miso events and take fragments withint fragment range */

data miso2ea_frag_w_iso;
   length annotation $200.;
   length gene_id $36.;
   length chr $4.;
   length strand $1.;
   length event_type $4.;
   format fragment_range_start best12.;
   format fragment_range_stop best12.;
   length exc_inc $1.;
   length transcript_id $1612.;
   length fragment_id $13.;
   if _n_ > 0 then delete;
run;

%macro mergeFrag2MISO();

   
    /* declare macro variables to use */
    %local iterAnnot countAnnot misoAnnot misoGene misoChr misoStrand misoType misoStart misoStop misoEI;
  
    /* Get the number of genes to iterate over */

    data _null_;
      set miso2ea_frag2 nobs=n;
      call symputx('countAnnot', n);
      stop;
    run;

    /* iterAnnot loop for genes */
    %let iterAnnot=1;
    %do %while (&iterAnnot. <= &countAnnot.  );

        /* Extract MISO event info */
        data _null_ ;
            set miso2ea_frag2 (firstobs=&iterAnnot. obs=&iterAnnot.); * read 1 record;
            call symput("misoAnnot", strip(annotation)); 
            call symput("misoGene", strip(gene_id)); 
            call symput("misoChr", strip(chr)); 
            call symput("misoStrand", strip(strand)); 
            call symput("misoType", strip(event_type)); 
            call symput("misoStart", strip(fragment_range_start)); 
            call symput("misoStop", strip(fragment_range_stop)); 
            call symput("misoEI", strip(exc_inc)); 
        run;
    
        /* Extract fragments */

        data frag2keep;
            set fragment_annot;
            if chr="&misoChr." and fragment_start >="&misoStart." and fragment_stop <="&misoStop." ;
            keep fragment_id transcript_id;
        run;

        data frag_roz;
         length annotation $200.;
         length gene_id $36.;
         length chr $4.;
         length strand $1.;
         length event_type $4.;
         format fragment_range_start best12.;
         format fragment_range_stop best12.;
         length exc_inc $1.;
         set frag2keep;
         annotation="&misoAnnot.";
         gene_id="&misoGene.";
         chr="&misoChr.";
         strand="&misoStrand.";
         event_type="&misoType.";
         fragment_range_start="&misoStart.";
         fragment_range_stop="&misoStop.";
         exc_inc="&misoEI.";
      run;
  
               /* Append to dataset to save */

               proc append base=miso2ea_frag_w_iso data=frag_roz force;
               run;

           /* increment query transcript loop */
           %let iterAnnot=%eval(&iterAnnot. + 1);

        %end;

%mend;


%mergeFrag2MISO();


/* Short miso ID */

data miso_events_all;
   set miso_index;
   keep annotation event_Type;
run;

proc sort data=miso_events_all nodup;
   by event_type annotation;
run;

data miso_shortID;
  retain event_num;
  set miso_events_all;
  length miso_id $20.;
  by event_type;
  if first.event_type then event_num=1;
  else event_num+1;
  miso_id = catx("_",event_type,"425",event_num);
  keep annotation event_type miso_id;
run;

proc sort data=miso_shortID;
   by annotation event_type;
proc sort data=miso2ea_frag_w_iso;
   by annotation event_type;
proc sort data=miso2ea_junc_w_iso;
   by annotation event_type;
run;


data miso2ea_frag2short;
  merge miso2ea_frag_w_iso (in=in1) miso_shortID (in=in2);
  by annotation event_type;
  if in1 and in2;
run;

data miso2ea_junc2short;
  merge miso2ea_junc_w_iso (in=in1) miso_shortID (in=in2);
  by annotation event_type;
  if in1 and in2;
run;

data frag2short;
  set miso2ea_frag2short;
  keep fragment_id miso_id;
  rename fragment_id=event_id;
run;



proc sort data=frag2short;
   by event_id miso_id;
proc freq data=frag2short noprint;
   tables event_id / out=miso_count;
proc sort data=miso_count;
  by descending count;
run; * 6 max miso events per fragment;


data frag2short_cat;
   array miso[6] $20.;
   retain miso1-miso6;
   set frag2short ;
   by event_id;
   if first.event_id then do;
       call missing(of miso1-miso6);
       records=0;
       end;
   records + 1;
   miso[records]=miso_id;
   if last.event_id then output;
run;

data frag2short_cat2;
  set frag2short_cat;
  length miso_id_cat $200.;
  miso_id_cat = catx("|", OF miso1-miso6);
  keep event_id miso_id_cat;
run;


proc import datafile="!HOME/concannon/case_control/references/hg19_aceview/annotations/hg19_aceview_junction_to_sequence_index.csv"
    out=junc2seq dbms=csv replace;
    guessingrows=max;
run;

data junc2seq2;
  set junc2seq;
  keep junction_id sequence_id;
  rename sequence_id=event_id;
run;

data junc2short;
  set miso2ea_junc2short;
  keep miso_id junction_id;
run;

proc sort data=junc2seq2;
  by junction_id;
proc sort data=junc2short;
  by junction_id;
run;

data junc2short2;
  merge junc2short (in=in1) junc2seq2 (in=in2);
  by junction_id;
  if in1 and in2;
run;
proc sort data=junc2short2;
   by event_id miso_id;
proc freq data=junc2short2 noprint;
   tables event_id / out=miso_count;
proc sort data=miso_count;
  by descending count;
run; * 6 max miso events per fragment;


data junc2short_cat;
   array miso[6] $20.;
   retain miso1-miso6;
   set junc2short2 ;
   by event_id;
   if first.event_id then do;
       call missing(of miso1-miso6);
       records=0;
       end;
   records + 1;
   miso[records]=miso_id;
   if last.event_id then output;
run;

data junc2short_cat2;
  set junc2short_cat;
  length miso_id_cat $200.;
  miso_id_cat = catx("|", OF miso1-miso6);
  keep event_id miso_id_cat;
run;


data frag_counts;
   set cclocal.fragment_q3_norm_counts;
   keep sample_id event_id apn q3_apn;
run;

data junc_counts;
   set cclocal.junction_q3_norm_counts;
   keep sample_id event_id apn q3_apn;
run;

proc sort data=junc2short_cat2;
   by event_id;
proc sort data=frag2short_cat2;
   by event_id;
proc sort data=frag_counts;
   by event_id;
proc sort data=junc_counts;
   by event_id;
run;

data frag2short_cat_counts;
  merge frag2short_cat2 (in=in1) frag_counts (in=in2);
  by event_id;
  if in1;
run;

data junc2short_cat_counts;
  merge junc2short_cat2 (in=in1) junc_counts (in=in2);
  by event_id;
  if in1;
run;

data samples2keep;
  set cc.design_file_new;
  where cell_type=425;
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

run;

proc sort data=frag2short_cat_counts;
  by sample_id;
proc sort data=junc2short_cat_counts;
  by sample_id;
proc sort data=samples2keep;
  by sample_id;
run;

data frag2short_cat_counts2;
   merge samples2keep (in=in1) frag2short_cat_counts (in=in2);
   by sample_id;
   if in1;
run;

data junc2short_cat_counts2;
   merge samples2keep (in=in1) junc2short_cat_counts (in=in2);
   by sample_id;
   if in1;
run;

/* Back-calculate region_depth */

data frag_length;
  set fragment_annot;
  event_length = fragment_stop - fragment_start + 1;
  if event_length < 10 then flag_short_event=1; else flag_short_event=0;
  keep fragment_id event_length flag_short_event;
  rename fragment_id=event_id;
run;

proc sort data=frag_length;
   by event_id;
proc sort data=frag2short_cat_counts2;
   by event_id;
run;

data frag2short_cat_counts3;
  merge frag2short_cat_counts2 (in=in1) frag_length (in=in2);
   by event_id;
  if in1 and in2 ;
  region_depth = event_length * apn ;
  q3_region_depth = event_length * q3_apn ;
run;


data junc_length;
  set junc2seq;
  event_length = length(sequence);
  if event_length < 100 then flag_short_event=1; else flag_short_event=0;
  keep sequence_id event_length flag_short_event;
  rename sequence_id=event_id;
run;

proc sort data=junc_length nodup;
   by event_id;
proc sort data=junc2short_cat_counts2;
   by event_id;
run;

data junc2short_cat_counts3;
  merge junc2short_cat_counts2 (in=in1) junc_length (in=in2);
   by event_id;
  if in1 and in2 ;
  region_depth = event_length * apn ;
  q3_region_depth = event_length * q3_apn ;
run;

/* expand catted MISO IDs, I need to calculate E and I counts for fragments and junctions */

data frag2short_counts3;
  set frag2short_cat_counts3;
  length miso_id $20.;
  do i = 1 by 1 while (scan(miso_id_cat,i,"|") ^= "");
       miso_id=compress(scan(miso_id_cat,i,"|"));   
      output;
      end;
  run;

data junc2short_counts3;
  set junc2short_cat_counts3;
  length miso_id $20.;
  do i = 1 by 1 while (scan(miso_id_cat,i,"|") ^= "");
       miso_id=compress(scan(miso_id_cat,i,"|"));   
      output;
      end;
  run;


data frag2short_1;
   set miso2ea_frag2short ;
   keep miso_id fragment_id exc_inc;
   rename fragment_id=event_id;
run;

proc sort data=frag2short_counts3;
   by miso_id event_id;
proc sort data=frag2short_1;
   by miso_id event_id;
run;

data frag2short_counts4;
  merge frag2short_counts3 (in=in1) frag2short_1 (in=in2);
  by miso_id event_id;
  if in1 and in2;
run;



data junc2short_1;
   set miso2ea_junc2short ;
   keep miso_id junction_id exc_inc;
run;

proc sort data=junc2short_1;
   by junction_id;
proc sort data=junc2seq2;
   by junction_id;
run;

data junc2short_2;
  merge junc2short_1 (in=in1) junc2seq2 (in=in2);
  by junction_id;
  if in1 and in2;
run;


proc sort data=junc2short_counts3;
   by miso_id event_id;
proc sort data=junc2short_2;
   by miso_id event_id;
run;

data junc2short_counts4;
  merge junc2short_counts3 (in=in1) junc2short_2 (in=in2);
  by miso_id event_id;
  if in1 and in2;
run;


/* ID miso events that won't be calculated due to small event size */

data frag2short_event_size;
   set frag2short_counts4;
   keep exc_inc miso_id flag_short_event event_id;
run;

data junc2short_event_size;
   set junc2short_counts4;
   keep exc_inc miso_id flag_short_event event_id;
run;

proc sort data=frag2short_event_size nodup;
  by miso_id exc_inc event_id flag_short_event;
proc sort data=junc2short_event_size nodup;
  by miso_id exc_inc event_id flag_short_event;
run;

proc means data=frag2short_event_size noprint; 
  by miso_id exc_inc;
  var flag_short_event ;
  output out=cnt_frag_len_ok N=n_frag_total sum=n_frag_short;
run;

proc means data=junc2short_event_size noprint; 
  by miso_id exc_inc;
  var flag_short_event ;
  output out=cnt_junc_len_ok N=n_junc_total sum=n_junc_short;
run;



/* flag ok if at least 1 event usable */

data flag_frag_len_ok;
   set cnt_frag_len_ok;
   if n_frag_short < n_frag_total then flag_frag_ok=1; 
   else flag_frag_ok=0;
run;

data flag_junc_len_ok;
   set cnt_junc_len_ok;
   if n_junc_short < n_junc_total then flag_junc_ok=1; 
   else flag_junc_ok=0;
run;

proc sort data=flag_frag_len_ok;
   by miso_id exc_inc;
proc transpose data=flag_frag_len_ok out=flag_frag_len_ok_sbys(drop=_NAME_) prefix=flag_frag_ok_;
   by miso_id;
   id exc_inc;
   var flag_frag_ok;
run;

proc sort data=flag_junc_len_ok;
   by miso_id exc_inc;
proc transpose data=flag_junc_len_ok out=flag_junc_len_ok_sbys(drop=_NAME_) prefix=flag_junc_ok_;
   by miso_id;
   id exc_inc;
   var flag_junc_ok;
run;

data flag_event_len_ok_sbys;
  merge flag_frag_len_ok_sbys flag_junc_len_ok_sbys ;
  by miso_id;
run;

data flag_event_len_ok_sbys2;
   set flag_event_len_ok_sbys;
   length event_type $4.;
   event_type=compress(scan(miso_id,1,"_"));
   if flag_frag_ok_I = 1 or flag_junc_ok_I = 1 then flag_I_ok=1; else flag_I_ok=0;
   if flag_frag_ok_E = 1 or flag_junc_ok_E = 1 then flag_E_ok=1; else flag_E_ok=0;
   if flag_I_ok=1 and flag_E_ok=1 then flag_psi_estimable=1; else flag_psi_estimable=0;
run;

proc freq data=flag_event_len_ok_sbys2;
  tables flag_I_ok*flag_E_ok flag_psi_estimable;
run;


proc freq data=flag_event_len_ok_sbys2 noprint;
  tables event_type*flag_I_ok*flag_E_ok / out=event_e_i_check;
  tables event_type*flag_psi_estimable / out=event_est_check;
run;

proc print data=event_e_i_check;
proc print data=event_est_check;
run;


/*


       Table of flag_I_ok by flag_E_ok

     flag_I_ok     flag_E_ok

     Frequency|
     Percent  |
     Row Pct  |
     Col Pct  |       0|       1|  Total
     ---------+--------+--------+
            0 |      2 |      0 |      2
              |   0.25 |   0.00 |   0.25
              | 100.00 |   0.00 |
              |   3.33 |   0.00 |
     ---------+--------+--------+
            1 |     58 |    756 |    814
              |   7.11 |  92.65 |  99.75
              |   7.13 |  92.87 |
              |  96.67 | 100.00 |
     ---------+--------+--------+
     Total          60      756      816
                  7.35    92.65   100.00


       flag_psi_                             Cumulative    Cumulative
       estimable    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0          60        7.35            60         7.35
               1         756       92.65           816       100.00


 event_    flag_    flag_
  type      I_ok     E_ok    COUNT

  A3SS       0        0         1
  A3SS       1        0        13
  A3SS       1        1       117

  A5SS       0        0         1
  A5SS       1        0         6
  A5SS       1        1        83

  MXE        1        1        25

  RI         1        0        24
  RI         1        1       262

  SE         1        0        15
  SE         1        1       269


 event_    flag_psi_
  type     estimable    COUNT

  A3SS         0          14
  A3SS         1         117

  A5SS         0           7
  A5SS         1          83

  MXE          1          25

  RI           0          24
  RI           1         262

  SE           0          15
  SE           1         269

*/


proc sort data=frag2short_counts4;
   by sample_id miso_id exc_inc;
proc means data=frag2short_counts4 noprint;
  by sample_id miso_id exc_inc;
  var apn q3_apn event_length region_depth q3_region_depth;
  output out=frag2short_counts_by_miso sum(region_depth)=region_depth sum(q3_region_depth)=q3_region_depth sum(event_length)=event_length 
         mean(apn)=apn mean(q3_apn)=q3_apn;
run;

proc sort data=frag2short_counts_by_miso;
   by miso_id sample_id exc_inc;
proc transpose data=frag2short_counts_by_miso out=frag_miso_counts_depth(drop=_NAME_) prefix=depth_frag_;
   by miso_id sample_id;
   id exc_inc;
   var region_depth;
run;

proc transpose data=frag2short_counts_by_miso out=frag_miso_counts_q3_depth(drop=_NAME_) prefix=q3_depth_frag_;
   by miso_id sample_id;
   id exc_inc;
   var q3_region_depth;
run;

proc transpose data=frag2short_counts_by_miso out=frag_miso_counts_length(drop=_NAME_) prefix=length_frag_;
   by miso_id sample_id;
   id exc_inc;
   var event_length;
run;


proc transpose data=frag2short_counts_by_miso out=frag_miso_counts_apn(drop=_NAME_) prefix=apn_frag_;
   by miso_id sample_id;
   id exc_inc;
   var apn;
run;

proc transpose data=frag2short_counts_by_miso out=frag_miso_counts_q3_apn(drop=_NAME_) prefix=q3_apn_frag_;
   by miso_id sample_id;
   id exc_inc;
   var q3_apn;
run;





proc sort data=junc2short_counts4;
   by sample_id miso_id exc_inc;
proc means data=junc2short_counts4 noprint;
  by sample_id miso_id exc_inc;
  var apn q3_apn event_length region_depth q3_region_depth;
  output out=junc2short_counts_by_miso sum(region_depth)=region_depth sum(q3_region_depth)=q3_region_depth sum(event_length)=event_length 
         mean(apn)=apn mean(q3_apn)=q3_apn;
run;

proc sort data=junc2short_counts_by_miso;
   by miso_id sample_id exc_inc;
proc transpose data=junc2short_counts_by_miso out=junc_miso_counts_depth(drop=_NAME_) prefix=depth_junc_;
   by miso_id sample_id;
   id exc_inc;
   var region_depth;
run;

proc transpose data=junc2short_counts_by_miso out=junc_miso_counts_q3_depth(drop=_NAME_) prefix=q3_depth_junc_;
   by miso_id sample_id;
   id exc_inc;
   var q3_region_depth;
run;

proc transpose data=junc2short_counts_by_miso out=junc_miso_counts_length(drop=_NAME_) prefix=length_junc_;
   by miso_id sample_id;
   id exc_inc;
   var event_length;
run;


proc transpose data=junc2short_counts_by_miso out=junc_miso_counts_apn(drop=_NAME_) prefix=apn_junc_;
   by miso_id sample_id;
   id exc_inc;
   var apn;
run;

proc transpose data=junc2short_counts_by_miso out=junc_miso_counts_q3_apn(drop=_NAME_) prefix=q3_apn_junc_;
   by miso_id sample_id;
   id exc_inc;
   var q3_apn;
run;

data miso_all_counts_for_psi;
     merge frag_miso_counts_depth frag_miso_counts_q3_depth frag_miso_counts_length frag_miso_counts_apn frag_miso_counts_q3_apn
     junc_miso_counts_depth junc_miso_counts_q3_depth junc_miso_counts_length junc_miso_counts_apn junc_miso_counts_q3_apn;
     by miso_id sample_id;
run;


data miso_all_counts_for_psi2;
   set miso_all_counts_for_psi;
   array change _numeric_;
   do over change;
    if change=. then change=0;
    end;
run;

data miso_all_counts_for_psi3;
   set miso_all_counts_for_psi2;
   psi_depth_sum=(depth_frag_I + depth_junc_I) /  ((depth_frag_E + depth_junc_E) + (depth_frag_I + depth_junc_I));
   psi_q3_depth_sum=(q3_depth_frag_I + q3_depth_junc_I) /  ((q3_depth_frag_E + q3_depth_junc_E) + (q3_depth_frag_I + q3_depth_junc_I));

   psi_depth_avg = ((depth_frag_I + depth_junc_I)/2) / (((depth_frag_I + depth_junc_I)/2) + ((depth_frag_E + depth_junc_E)/2));

   psi_q3_depth_avg = ((q3_depth_frag_I + q3_depth_junc_I)/2) / (((q3_depth_frag_E + q3_depth_junc_E)/2) + ((q3_depth_frag_I + q3_depth_junc_I)/2));

   psi_apn_avg = ((apn_frag_I + apn_junc_I)/2) / (((apn_frag_I + apn_junc_I)/2) + ((apn_frag_E + apn_junc_E)/2));

   psi_q3_apn_avg = ((q3_apn_frag_I + q3_apn_junc_I)/2) / (((q3_apn_frag_I + q3_apn_junc_I)/2) + ((q3_apn_frag_E + q3_apn_junc_E)/2));

   psi_recalc_apn_avg = ((apn_frag_I/(length_frag_I+1) + apn_junc_I/(length_junc_I+1))/2) /
                        (((apn_frag_I/(length_frag_I+1) + apn_junc_I/(length_junc_I+1))/2) + 
                         ((apn_frag_E/(length_frag_E+1) + apn_junc_E/(length_junc_E+1))/2));
   psi_recalc_q3_apn_avg = ((q3_apn_frag_I/(length_frag_I+1) + q3_apn_junc_I/(length_junc_I+1))/2) /
                        (((q3_apn_frag_I/(length_frag_I+1) + q3_apn_junc_I/(length_junc_I+1))/2) + 
                         ((q3_apn_frag_E/(length_frag_E+1) + q3_apn_junc_E/(length_junc_E+1))/2));
run;


/* Import MISO PSIs */


%include "!PATCON/case_control/sas_programs/event_analysis/fixed/iterdataset.sas";


%macro impMISOpsi(sampleID);

proc import datafile="/TB14/TB14/sandbox/miso_sandbox/miso_summary_A3SS/&sampleID./summary/&sampleID..miso_summary"
   out=A3SS_&sampleID. dbms=tab replace;
   guessingrows=max;
run;

data A3SS2_&sampleID.;
  set A3SS_&sampleID.;
  length sample_id $22.;
  sample_id="&sampleID.";
  keep sample_id event_name miso_posterior_mean;
  rename event_name=annotation miso_posterior_mean=miso_PSI;
run;


proc import datafile="/TB14/TB14/sandbox/miso_sandbox/miso_summary_A5SS/&sampleID./summary/&sampleID..miso_summary"
   out=A5SS_&sampleID. dbms=tab replace;
   guessingrows=max;
run;

data A5SS2_&sampleID.;
  set A5SS_&sampleID.;
  length sample_id $22.;
  sample_id="&sampleID.";
  keep sample_id event_name miso_posterior_mean;
  rename event_name=annotation miso_posterior_mean=miso_PSI;
run;


proc import datafile="/TB14/TB14/sandbox/miso_sandbox/miso_summary_MXE/&sampleID./summary/&sampleID..miso_summary"
   out=MXE_&sampleID. dbms=tab replace;
   guessingrows=max;
run;

data MXE2_&sampleID.;
  set MXE_&sampleID.;
  length sample_id $22.;
  sample_id="&sampleID.";
  keep sample_id event_name miso_posterior_mean;
  rename event_name=annotation miso_posterior_mean=miso_PSI;
run;


proc import datafile="/TB14/TB14/sandbox/miso_sandbox/miso_summary_RI/&sampleID./summary/&sampleID..miso_summary"
   out=RI_&sampleID. dbms=tab replace;
   guessingrows=max;
run;

data RI2_&sampleID.;
  set RI_&sampleID.;
  length sample_id $22.;
  sample_id="&sampleID.";
  keep sample_id event_name miso_posterior_mean;
  rename event_name=annotation miso_posterior_mean=miso_PSI;
run;


proc import datafile="/TB14/TB14/sandbox/miso_sandbox/miso_summary_SE/&sampleID./summary/&sampleID..miso_summary"
   out=SE_&sampleID. dbms=tab replace;
   guessingrows=max;
run;

data SE2_&sampleID.;
  set SE_&sampleID.;
  length sample_id $22.;
  sample_id="&sampleID.";
  keep sample_id event_name miso_posterior_mean;
  rename event_name=annotation miso_posterior_mean=miso_PSI;
run;

%mend;

data samples;
  set miso_all_counts_for_psi3;
  keep sample_id;
run;

proc sort data=samples nodup;
  by sample_id;
run;

   %iterdataset(dataset=samples, function=%nrstr(%impMISOpsi(&sample_id.);));


data miso_psi_all;
   length annotation $200.;
   set MXE2_: A3SS2_: A5SS2_: RI2_: SE2_:;
run;


proc sort data=miso_psi_all;
   by annotation;
proc sort data=miso_shortid;
   by annotation;
run;

proc freq data=miso_shortid noprint;
   tables annotation / out=check;
run;


data miso_psi_all2;
  merge miso_psi_all (in=in1) miso_shortid (in=in2);
  by annotation;
  if in1 and in2;
run;


proc sort data=miso_psi_all2;
   by miso_id sample_id;
proc sort data=miso_all_counts_for_psi3;
   by miso_id sample_id;
run;

data miso_psi_vs_ea;
  merge miso_all_counts_for_psi3 (in=in1) miso_psi_all2 (in=in2);
  by miso_id sample_id;
  *if in1 and in2;
run;


proc sort data=miso_psi_vs_ea;
   by miso_id;
proc sort data=flag_event_len_ok_sbys2;
   by miso_id;
run;

data miso_psi_vs_ea2;
  merge miso_psi_vs_ea (in=in1) flag_event_len_ok_sbys2 (in=in2);
  by miso_id;
run;

ods graphics / antialiasmax=200000;

proc sgplot data=miso_psi_vs_ea2;
  where flag_psi_estimable=1 and event_type="SE";
  scatter x=psi_depth_sum y=miso_psi;
run;

proc sgplot data=miso_psi_vs_ea2;
  where flag_psi_estimable=1;
  scatter x=psi_q3_depth_sum y=miso_psi;
run;


proc sgplot data=miso_psi_vs_ea2;
  where flag_psi_estimable=1;
  scatter x=psi_apn_avg y=miso_psi;
run;

proc sgplot data=miso_psi_vs_ea2;
  where flag_psi_estimable=1;
  scatter x=psi_recalc_q3_apn_avg y=miso_psi;
run;


proc corr data=miso_psi_vs_ea2 pearson out=miso_corr noprint;
   by miso_id;
   var psi_depth_sum psi_q3_depth_sum psi_apn_avg psi_recalc_q3_apn_avg miso_psi;
run;

data miso_corr2;
  set miso_corr;
  where _NAME_ = "miso_PSI";
run;

proc sort data=miso_corr2;
   by descending psi_depth_sum;
run;



data design;
  set  cc.design_file_new;
run;
proc sort data=miso_psi_vs_ea ;
  by sample_id;
proc sort data=design;
  by sample_id;
run;

data miso_psi_vs_ea2;
  merge miso_psi_vs_ea (in=in1) design (in=in2);
  by sample_id;
  if in1 and in2;
run;

proc sort data=miso_psi_vs_ea2;
  by  miso_id cell_type;
run;

proc glimmix data=miso_psi_vs_ea2 ;
  by  miso_id  cell_type;
  class status ;
  model psi_apn_avg = status / ddfm=kr;
  output out=pred_apn_avg pred=pred resid=resid student=stu;
  ods output tests3=anova_apn_avg;
run;


proc glimmix data=miso_psi_vs_ea2 ;
  by  miso_id  cell_type;
  class status ;
  model psi_q3_apn_avg = status / ddfm=kr;
  output out=pred_q3_apn_avg pred=pred resid=resid student=stu;
  ods output tests3=anova_q3_apn_avg;
run;


proc glimmix data=miso_psi_vs_ea2 ;
  by  miso_id  cell_type;
  class status ;
  model psi_q3_depth_sum = status / ddfm=kr;
  output out=pred_q3_depth_sum pred=pred resid=resid student=stu;
  ods output tests3=anova_q3_depth_sum;
run;


proc glimmix data=miso_psi_vs_ea2 ;
  by  miso_id  cell_type;
  class status ;
  model miso_psi = status / ddfm=kr;
  output out=pred_miso pred=pred resid=resid student=stu;
  ods output tests3=anova_miso;
run;


data anova_apn_avg2;
   set anova_apn_avg;
   keep miso_id cell_type fvalue;
   rename fvalue=F_psi_apn_avg;
run;

data anova_q3_apn_avg2;
   set anova_q3_apn_avg;
   keep miso_id cell_type fvalue;
   rename fvalue=F_psi_q3_apn_avg;
run;
data anova_q3_depth_sum2;
   set anova_q3_depth_sum;
   keep miso_id cell_type fvalue;
   rename fvalue=F_psi_q3_depth_Sum;
run;

data anova_miso2;
   set anova_miso;
   keep miso_id cell_type fvalue;
   rename fvalue=F_psi_miso;
run;

proc sort data=anova_apn_avg2;
   by miso_id cell_type;
proc sort data=anova_q3_apn_avg2;
   by miso_id cell_type;
proc sort data=anova_q3_depth_sum2;
   by miso_id cell_type;
proc sort data=anova_miso2;
   by miso_id cell_type;
run;

data anvoa_compare;
  merge anova_apn_avg2 anova_q3_apn_avg2 anova_q3_depth_sum2 anova_miso2;
   by miso_id cell_type;
run;

proc corr data=anvoa_compare pearson;
   var F_: ;
run;




data anova_apn_avg2;
   set anova_apn_avg;
   keep miso_id cell_type probf;
   rename probf=p_psi_apn_avg;
run;

data anova_q3_apn_avg2;
   set anova_q3_apn_avg;
   keep miso_id cell_type probf;
   rename probf=p_psi_q3_apn_avg;
run;
data anova_q3_depth_sum2;
   set anova_q3_depth_sum;
   keep miso_id cell_type probf;
   rename probf=p_psi_q3_depth_Sum;
run;

data anova_miso2;
   set anova_miso;
   keep miso_id cell_type probf;
   rename probf=p_psi_miso;
run;

proc sort data=anova_apn_avg2;
   by miso_id cell_type;
proc sort data=anova_q3_apn_avg2;
   by miso_id cell_type;
proc sort data=anova_q3_depth_sum2;
   by miso_id cell_type;
proc sort data=anova_miso2;
   by miso_id cell_type;
run;

data anova_compare;
  merge anova_apn_avg2 anova_q3_apn_avg2 anova_q3_depth_sum2 anova_miso2;
   by miso_id cell_type;
run;

data anova_compare2;
   set anova_compare;
   if p_psi_apn_avg = . then flag_apn=.;
   else if p_psi_apn_avg < 0.05 then flag_apn=1; else flag_apn=0;

   if p_psi_q3_apn_avg = . then flag_q3_apn=.;
   else if p_psi_q3_apn_avg < 0.05 then flag_q3_apn=1; else flag_q3_apn=0;

   if p_psi_q3_depth_Sum = . then flag_q3_depth=.;
   else if p_psi_q3_depth_Sum < 0.05 then flag_q3_depth=1; else flag_q3_depth=0;

   if p_psi_miso = . then flag_miso=.;
   else if p_psi_miso < 0.05 then flag_miso=1; else flag_miso=0;

run;

proc freq data=anova_compare2 noprint;
  tables flag_apn*flag_q3_apn*flag_q3_Depth*flag_miso / out=compare;
proc print data=compare;
run;


/*
             flag_q3_    flag_q3_    flag_
 flag_apn       apn        depth      miso    COUNT    PERCENT

     .           .           .         .         8       .
     .           .           .         0        59       .
     .           .           .         1         4       .
     0           0           0         .         4       .
     0           0           0         0       441     59.4340
     0           0           0         1        62      8.3558
     0           0           1         0        14      1.8868
     0           0           1         1        21      2.8302
     0           1           0         0         3      0.4043
     0           1           0         1         3      0.4043
     0           1           1         0         3      0.4043
     0           1           1         1         6      0.8086
     1           0           0         0        19      2.5606
     1           0           0         1        13      1.7520
     1           0           1         0         3      0.4043
     1           0           1         1        14      1.8868
     1           1           0         0         2      0.2695
     1           1           0         1         7      0.9434
     1           1           1         0        26      3.5040
     1           1           1         1       105     14.1509


*/

proc sort data=miso_psi_vs_ea2;
   by miso_id cell_type status;
proc means data=miso_psi_vs_ea2 noprint;
   by miso_id cell_type status;
   var depth_frag_I depth_frag_E depth_junc_I depth_junc_E ;
   output out=sum_inc_exc_depth_by_miso sum=;
run;


proc means data=miso_psi_vs_ea2 noprint;
   by miso_id cell_type status;
   var depth_frag_I depth_frag_E depth_junc_I depth_junc_E miso_psi psi_depth_sum;
   output out=mean_inc_exc_depth_by_miso mean=;
run;


proc sort data=miso_psi_vs_ea2;
   by miso_id cell_type status;
proc means data=miso_psi_vs_ea2 noprint;
  by miso_id cell_type status;
  var   psi_apn_avg   psi_q3_apn_avg   psi_q3_depth_sum psi_depth_sum   miso_psi ;
  output out=mean_psi_by_event mean=;
run;

/* Annotate transcripts 

A3SS, A5SS, MXE : take transcripts from junctions
SE: take E transcript from skipping junction, I from including BOTH junctions (proc freq on transcritps, take count = 2)
RI: take E transcript from junction, I from ALL fragments (proc freq on transcript, take max count) 

*/


data miso_junc2xs;
  set miso2ea_junc_w_iso;
  keep annotation gene_id junction_id transcript_id exc_inc;
run;

data miso_frag2xs;
  set miso2ea_frag_w_iso;
  keep annotation gene_id fragment_id transcript_id exc_inc;
run;

proc sort data=miso_junc2xs;
   by annotation;
proc sort data=miso_frag2xs;
   by annotation;
proc sort data=miso_shortid;
   by annotation;
run;
  

data miso_junc2xs2;
  merge miso_junc2xs (in=in1) miso_shortid (in=in2);
  by annotation;
  if in1 and in2;
run;

data miso_frag2xs2;
  merge miso_frag2xs (in=in1) miso_shortid (in=in2);
  by annotation;
  if in1 and in2;
run;

/* Get E and I transcripts from junctions */

data miso_xs_a5ss_a3ss;
   set miso_junc2xs2;
   where event_type in ("A5SS","A3SS");
   length transcript_id2 $100.;
   do i = 1 by 1 while(scan(transcript_id,i,"|") ^= "");
     transcript_id2=compress(scan(transcript_id,i,"|"));
     output;
     end;
   keep miso_id event_type gene_id exc_inc transcript_id2;
run;

proc sort data=miso_xs_a5ss_a3ss;
   by miso_id event_type gene_id exc_inc transcript_id2;
proc freq data=miso_xs_a5ss_a3ss noprint;
   by miso_id event_type gene_id exc_inc;
   tables transcript_id2 / out=xs_count_per_junc;
run;

proc sort data=xs_count_per_junc;
   by miso_id gene_id exc_inc descending count;
run;

data xs_count_per_junc2;
  set xs_count_per_junc;
   by miso_id gene_id exc_inc;
   if first.exc_inc then output;
   drop PERCENT transcript_id2;
run;

proc sort data=xs_count_per_junc2;
  by miso_id gene_id exc_inc count;
proc sort data=xs_count_per_junc;
  by miso_id gene_id exc_inc count;
run;

data miso_xs_a5ss_a3ss2;
  merge xs_count_per_junc (in=in1) xs_count_per_junc2 (in=in2);
  by miso_id gene_id exc_inc count;
  if in1 and in2;
run;


data miso_xs_mxe;
   set miso_junc2xs2;
   where event_type in ("MXE");
   length transcript_id2 $100.;
   do i = 1 by 1 while(scan(transcript_id,i,"|") ^= "");
     transcript_id2=compress(scan(transcript_id,i,"|"));
     output;
     end;
   keep miso_id event_type gene_id exc_inc transcript_id2;
run;

proc sort data=miso_xs_mxe;
   by miso_id event_type gene_id exc_inc transcript_id2;
proc freq data=miso_xs_mxe noprint;
   by miso_id event_type gene_id exc_inc;
   tables transcript_id2 / out=xs_count_per_junc;
run;

proc sort data=xs_count_per_junc;
   by miso_id gene_id exc_inc descending count;
run;

data miso_xs_mxe2;
  set xs_count_per_junc;
   by miso_id gene_id exc_inc;
   if count = 2;
   drop PERCENT ;
run;



  

data miso_xs_se_exc;
   set miso_junc2xs2;
   where event_type in ("SE") and exc_inc="E";
   length transcript_id2 $100.;
   do i = 1 by 1 while(scan(transcript_id,i,"|") ^= "");
     transcript_id2=compress(scan(transcript_id,i,"|"));
     output;
     end;
   keep miso_id event_type gene_id exc_inc transcript_id2;
run;


data miso_xs_se_inc;
   set miso_junc2xs2;
   where event_type in ("SE") and exc_inc="I";
   length transcript_id2 $100.;
   do i = 1 by 1 while(scan(transcript_id,i,"|") ^= "");
     transcript_id2=compress(scan(transcript_id,i,"|"));
     output;
     end;
   keep miso_id event_type gene_id exc_inc transcript_id2;
run;


proc sort data=miso_xs_se_inc;
   by miso_id event_type gene_id exc_inc transcript_id2;
proc freq data=miso_xs_se_inc noprint;
   by miso_id event_type gene_id exc_inc;
   tables transcript_id2 / out=xs_count_per_junc;
run;

proc sort data=xs_count_per_junc;
   by miso_id gene_id exc_inc descending count;
run;

data miso_xs_se_inc2;
  set xs_count_per_junc;
   by miso_id gene_id exc_inc;
   if count = 2;
   drop PERCENT ;
run;



data miso_xs_se_inc_frag;
   set miso_frag2xs2;
   where event_type in ("SE") and exc_inc="I";
   length transcript_id2 $100.;
   do i = 1 by 1 while(scan(transcript_id,i,"|") ^= "");
     transcript_id2=compress(scan(transcript_id,i,"|"));
     output;
     end;
   keep miso_id event_type gene_id exc_inc transcript_id2;
run;






proc sort data=miso_xs_se_inc_frag;
   by miso_id event_type gene_id exc_inc transcript_id2;
proc freq data=miso_xs_se_inc_frag noprint;
   by miso_id event_type gene_id exc_inc;
   tables transcript_id2 / out=xs_count_per_frag;
run;

proc sort data=xs_count_per_frag;
   by miso_id gene_id exc_inc descending count;
run;

data xs_count_per_frag2;
  set xs_count_per_frag;
   by miso_id gene_id exc_inc;
   if first.exc_inc then output;
   drop PERCENT transcript_id2;
run;

proc sort data=xs_count_per_frag2;
  by miso_id gene_id exc_inc count;
proc sort data=xs_count_per_frag;
  by miso_id gene_id exc_inc count;
run;

data miso_xs_se_inc_frag2;
  merge xs_count_per_frag (in=in1) xs_count_per_frag2 (in=in2);
  by miso_id gene_id exc_inc count;
  if in1 and in2;
run;


data miso_xs_ri_exc;
   set miso_junc2xs2;
   where event_type in ("RI") and exc_inc="E";
   length transcript_id2 $100.;
   do i = 1 by 1 while(scan(transcript_id,i,"|") ^= "");
     transcript_id2=compress(scan(transcript_id,i,"|"));
     output;
     end;
   keep miso_id event_type gene_id exc_inc transcript_id2;
run;




data miso_xs_ri_inc_frag;
   set miso_frag2xs2;
   where event_type in ("RI") and exc_inc="I";
   length transcript_id2 $100.;
   do i = 1 by 1 while(scan(transcript_id,i,"|") ^= "");
     transcript_id2=compress(scan(transcript_id,i,"|"));
     output;
     end;
   keep miso_id event_type gene_id exc_inc transcript_id2;
run;


data miso_xs_se_inc_frag3;
   set miso_xs_se_inc_frag2;
   drop percent;
run;

proc sort data=miso_xs_se_inc_frag3;
   by miso_id event_type gene_id exc_inc transcript_id2;
proc sort data=miso_xs_se_inc2;
   by miso_id event_type gene_id exc_inc transcript_id2;
run;

data miso_xs_se_inc3;
  merge miso_xs_se_inc2 (in=in1) miso_xs_se_inc_frag3 (in=in2);
  by miso_id event_type gene_id exc_inc transcript_id2;
  if in1 and in2;
run;


data miso_xs_se_inc3_drop;
  set miso_xs_se_inc3;
  keep miso_id exc_inc;
run;

proc sort data=miso_xs_se_inc3_drop nodup;
  by miso_id exc_inc;
proc sort data=miso_xs_se_inc_frag2;
  by miso_id exc_inc;
run;

data miso_xs_se_inc_frag_only;
  merge miso_xs_se_inc_frag2 (in=in1) miso_xs_se_inc3_drop (in=in2);
  by miso_id exc_inc; 
  if in2 then delete;
run;





proc sort data=miso_xs_ri_inc_frag;
   by miso_id event_type gene_id exc_inc transcript_id2;
proc freq data=miso_xs_ri_inc_frag noprint;
   by miso_id event_type gene_id exc_inc;
   tables transcript_id2 / out=xs_count_per_frag;
run;

proc sort data=xs_count_per_frag;
   by miso_id gene_id exc_inc descending count;
run;

data xs_count_per_frag2;
  set xs_count_per_frag;
   by miso_id gene_id exc_inc;
   if first.exc_inc then output;
   drop PERCENT transcript_id2;
run;

proc sort data=xs_count_per_frag2;
  by miso_id gene_id exc_inc count;
proc sort data=xs_count_per_frag;
  by miso_id gene_id exc_inc count;
run;

data miso_xs_ri_inc_frag2;
  merge xs_count_per_frag (in=in1) xs_count_per_frag2 (in=in2);
  by miso_id gene_id exc_inc count;
  if in1 and in2;
run;

data miso_xs_stack;
    set miso_xs_a5ss_a3ss2 miso_xs_mxe2 miso_xs_ri_exc miso_xs_ri_inc_frag2 miso_xs_se_exc
    miso_xs_se_inc_frag_only miso_xs_se_inc3;
    keep miso_id event_type gene_id exc_inc transcript_id2;
    rename transcript_id2=transcript_id;
run;




/* Stack!! */


data xs2gene;
  set hg19.hg19_xscript2gene;
run;

proc sort data=miso_xs_stack;
  by gene_id transcript_id;
proc sort data=xs2gene;
  by gene_id transcript_id;
run;

data miso_xs_stack2;
  merge miso_xs_stack (in=in1) xs2gene (in=in2);
  by gene_id transcript_id;
  if in1 and in2;
run;


proc sort data=miso_xs_stack2 nodup;
   by transcript_id miso_id;
proc freq data=miso_xs_stack2 noprint;
   tables transcript_id / out=miso_count;
proc sort data=miso_count;
  by descending count;
run; * 11 max miso events per fragment;


data miso_xs_stack2cat;
   array miso[11] $20.;
   retain miso1-miso11;
   set miso_xs_stack2 ;
   by transcript_id;
   if first.transcript_id then do;
       call missing(of miso1-miso6);
       records=0;
       end;
   records + 1;
   miso[records]=miso_id;
   if last.transcript_id then output;
run;

data miso_xs_stack2cat2;
  set miso_xs_stack2cat;
  length miso_id_cat $200.;
  miso_id_cat = catx("|", OF miso1-miso6);
  keep transcript_id miso_id_cat;
run;

data xs_counts;
  set cc.hg19_rsem_100p_apn2_xs_isoPerc;
  where cell_type=425;
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
  keep sample_id transcript_id isopct;
run;

proc sort data=miso_xs_stack2cat2;
   by transcript_id;
proc sort data=xs_counts;
  by transcript_id;
run;

data miso_xs_stack2cat_cnt;
  merge miso_xs_stack2cat2 (in=in1) xs_counts (in=in2);
  by transcript_id;
  if in1 and in2;
run;

data miso_xs_stack_cnt;
   set miso_xs_stack2cat_cnt;
   length miso_id $20.;
   do i=1 by 1 while(scan(miso_id_cat,i,"|") ^= "");
     miso_id=compress(scan(miso_id_cat,i,"|"));
     output; end;
   drop miso_id_cat;
run;

proc sort data=miso_xs_stack_cnt;
   by miso_id transcript_id;
proc sort data=miso_xs_stack2;
   by miso_id transcript_id;
run;

data miso_xs_stack_cnt2;
  merge miso_xs_stack_cnt (in=in1) miso_xs_Stack2 (in=in2);
  by miso_id transcript_id;
  if in1 and in2;
run;


proc sort data=miso_xs_stack_cnt2;
  by miso_id gene_id exc_inc sample_id;
proc means data=miso_xs_stack_cnt2 noprint;
  by miso_id gene_id exc_inc sample_id;
  var isopct;
  output out=miso_sum_isopct sum=;
run;

proc sort data=miso_sum_isopct;
  by sample_id miso_id gene_id;
proc transpose data=miso_sum_isopct out=miso_sum_isopct_sbys;
   by sample_id miso_id gene_id;
  id exc_inc;
  var isopct;
run;

data miso_sum_isopct_sbys2;
  set miso_sum_isopct_sbys;
  iso_diff = I / (E + I);
run;


proc sort data=miso_psi_vs_ea2;
   by miso_id sample_id;
proc sort data=miso_sum_isopct_sbys2;
   by miso_id sample_id;
run;

data miso_sum_isopct_sbys3;
  merge miso_sum_isopct_sbys2 (in=in1) miso_psi_vs_ea2 (in=in2);
  by miso_id sample_id;
  if in1 and in2;
run;


proc corr data=miso_sum_isopct_sbys3 pearson;
   var iso_diff psi_q3_depth_sum psi_apn_avg psi_q3_apn_avg miso_PSI;
run; 

/* psi_q3_apn_avg agress best with iso_diff */

/* For each GENE and MISO EVENT, I need:
GeneID
Splicing type (A3SS, A5SS, MXE, RI, SE)
Exonic sequence inclusion (A|B)
PSI control (q3 apn)
PSI case (q3 apn)
*/


data miso_annot;
  set miso_index_w_gene;
  if gene_id="" then delete;
  length spliced_in_coord $300.;
  if event_type="A3SS" then
       spliced_in_coord=catt("chr",scan(annotation,1,":"),":",scan(scan(annotation,5,":"),1,"|"),"-",scan(scan(annotation,5,":"),2,"|"),":",scan(annotation,7,":"));



  if event_type="A5SS" then
       spliced_in_coord=catt("chr",scan(annotation,1,":"),":",scan(scan(annotation,3,":"),1,"|"),"-",scan(scan(annotation,3,":"),2,"|"),":",scan(annotation,7,":"));



  if event_type="MXE" then
       spliced_in_coord=catx("/", catt("chr",scan(annotation,1,":"),":",scan(annotation,5,":"),"-",scan(annotation,6,":"),":",scan(annotation,13,":")),
                                  catt("chr",scan(annotation,1,":"),":",scan(annotation,8,":"),"-",scan(annotation,9,":"),":",scan(annotation,13,":")));


  if event_type="RI" then
       spliced_in_coord=catt("chr",scan(annotation,1,":"),":",scan(scan(annotation,2,":"),2,"-"),"-",scan(scan(annotation,4,":"),1,"-"),":",scan(annotation,5,":"));


  if event_type="SE" then
       spliced_in_coord=catt("chr",scan(annotation,1,":"),":",scan(annotation,5,":"),"-",scan(annotation,6,":"),":",scan(annotation,10,":"));



  keep annotation event_type spliced_in_coord gene_id;
run;

proc sort data=miso_annot nodup;
  by annotation gene_id event_type spliced_in_coord;
run;

/*************************************************************/
data miso_mean_psi;
   set mean_psi_by_event;
   keep miso_id status psi_q3_apn_avg;
run;

proc sort data=miso_mean_psi;
  by miso_id status;
proc transpose data=miso_mean_psi out=miso_mean_psi_sbys;
  by miso_id;
  id status;
  var psi_q3_apn_avg;
run;

data miso_mean_psi_sbys2;
  set miso_mean_psi_Sbys;
  percent_spliced_in_CTL= ctl *100;
  percent_spliced_in_t1d= t1d *100;
  diff_percent_spliced_in=(T1D -CTL) * 100;
  drop _NAME_ CTL T1D;
run;

data sum_inc_exc_depth;
  set sum_inc_exc_depth_by_miso;
  keep miso_id status depth_frag_I depth_junc_I depth_junc_E;
  rename depth_frag_I=sum_depth_frag_I depth_junc_I=sum_depth_junc_I depth_junc_E=sum_depth_junc_E;
run;


proc sort data=sum_inc_exc_depth;
  by miso_id status;
proc transpose data=sum_inc_exc_depth out=sum_inc_frag_sbys;
  by miso_id;
  id status;
  var sum_depth_frag_I;
run;

proc transpose data=sum_inc_exc_depth out=sum_inc_junc_sbys;
  by miso_id;
  id status;
  var sum_depth_junc_I;
run;

proc transpose data=sum_inc_exc_depth out=sum_exc_junc_sbys;
  by miso_id;
  id status;
  var sum_depth_junc_E;
run;



data sum_inc_frag_sbys2;
  set sum_inc_frag_sbys;
  rename ctl=ctl_sum_depth_frag_I t1d=t1d_sum_depth_frag_I;
  drop _NAME_ ;
run;


data sum_inc_junc_sbys2;
  set sum_inc_junc_sbys;
  rename ctl=ctl_sum_depth_junc_I t1d=t1d_sum_depth_junc_I;
  drop _NAME_ ;
run;

data sum_exc_junc_sbys2;
  set sum_exc_junc_sbys;
  rename ctl=ctl_sum_depth_junc_E t1d=t1d_sum_depth_junc_E;
  drop _NAME_ ;
run;




data mean_inc_exc_depth;
  set mean_inc_exc_depth_by_miso;
  keep miso_id status depth_frag_I depth_junc_I depth_junc_E;
  rename depth_frag_I=mean_depth_frag_I depth_junc_I=mean_depth_junc_I depth_junc_E=mean_depth_junc_E;
run;


proc sort data=mean_inc_exc_depth;
  by miso_id status;
proc transpose data=mean_inc_exc_depth out=mean_inc_frag_sbys;
  by miso_id;
  id status;
  var mean_depth_frag_I;
run;

proc transpose data=mean_inc_exc_depth out=mean_inc_junc_sbys;
  by miso_id;
  id status;
  var mean_depth_junc_I;
run;

proc transpose data=mean_inc_exc_depth out=mean_exc_junc_sbys;
  by miso_id;
  id status;
  var mean_depth_junc_E;
run;



data mean_inc_frag_sbys2;
  set mean_inc_frag_sbys;
  rename ctl=ctl_mean_depth_frag_I t1d=t1d_mean_depth_frag_I;
  drop _NAME_ ;
run;


data mean_inc_junc_sbys2;
  set mean_inc_junc_sbys;
  rename ctl=ctl_mean_depth_junc_I t1d=t1d_mean_depth_junc_I;
  drop _NAME_ ;
run;

data mean_exc_junc_sbys2;
  set mean_exc_junc_sbys;
  rename ctl=ctl_mean_depth_junc_E t1d=t1d_mean_depth_junc_E;
  drop _NAME_ ;
run;


proc sort data=miso_mean_psi_sbys2;
  by miso_id;

proc sort data=sum_inc_frag_sbys2;
  by miso_id;
proc sort data=sum_inc_junc_sbys2;
  by miso_id;
proc sort data=sum_exc_junc_sbys2;
  by miso_id;

proc sort data=mean_inc_frag_sbys2;
  by miso_id;
proc sort data=mean_inc_junc_sbys2;
  by miso_id;
proc sort data=mean_exc_junc_sbys2;
  by miso_id;


proc sort data=miso_shortid;
  by miso_id;
run;

data miso_Test_mean;
  merge miso_mean_psi_sbys2 (in=in1)  miso_shortid 
  sum_inc_frag_sbys2 sum_inc_junc_sbys2 sum_exc_junc_sbys2
  mean_inc_frag_sbys2 mean_inc_junc_sbys2 mean_exc_junc_sbys2 ;
  by miso_id;
  if in1 ;
run;

data miso_keep;
  set flag_event_len_ok_sbys2;
  where flag_psi_estimable=1;
  keep miso_id;
run;

proc sort data=miso_test_mean ;
  by miso_id;
proc sort data=miso_keep;
  by miso_id;
run;


data miso_test_mean2;
  merge miso_test_mean (in=in1) miso_keep (in=in2);
  by miso_id;
  if in1 and in2;
run;

proc sort data=miso_test_mean2;
  by annotation;
proc sort data=miso_annot;
  by annotation;
run;

data miso_annot_w_psi;
  merge miso_annot (in=in1) miso_test_mean2 (in=in2);
  by annotation;
  if in1 and in2;
  drop annotation;
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

proc sort data=miso_annot_w_psi;
  by gene_id;
proc sort data=ds_test;
  by gene_id;
proc sort data=orank_test;
  by gene_id;
proc sort data=trank_test;
  by gene_id;
run;

data DS_table;
  merge ds_test (in=in1) orank_test (in=in2) trank_test (in=in3) miso_annot_w_psi (in=in4);
  by gene_id;
  if not in4 then event_type="n/a";
  if in1 then output;
run;


data cc.miso_ds425_genes_table;
  set ds_table;
run;


data cc.miso_annot_w_psi;
set miso_annot_w_psi;
run;


