libname cc '!PATCON/case_control/sas_data/fixed';

ods listing;

proc import datafile="!PATCON/useful_human_data/aceview_hg19/SpliceAid_factors_aceview_id.csv"
   out=splicing_factor_list dbms=csv replace;
      guessingrows=max;
      run;



      proc import datafile="!PATCON/useful_human_data/aceview_hg19/genes_with_rrm1_domains.txt"
         out=rrm1_genes dbms=csv replace;
	    guessingrows=max;
	    run;


	    data ds_genes;
	      set cc.m37_anova_fdr_ds_case4;
	        run;
		  
		  proc sort data=ds_genes;
		      by gene_id;
		        proc sort data=rrm1_genes nodup;
			    by gene_id;
			      proc sort data=splicing_factor_list nodup;
			          by gene_id;
				      run;
				          
				       data ds_Genes2;
				          merge ds_genes (in=in1) splicing_factor_list (in=in2) rrm1_genes (in=in3);
					    by gene_id;
					      if in2 then flag_sf=1; else flag_sf=0;
					         if in3 then flag_rrm1=1; else flag_rrm1=0;
						   if in1 then output;
						     run;
						       
						       proc sort data=ds_Genes2;
						           by cell_Type;
							       run;
							           
							         proc freq data=ds_genes2;
								      by cell_type;
								           tables flag_fdr05*flag_rrm1 / chisq;
									        tables flag_fdr05*flag_sf / chisq;
										run;


