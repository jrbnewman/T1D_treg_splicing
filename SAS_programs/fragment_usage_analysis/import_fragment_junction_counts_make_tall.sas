ods html close; ods listing;

libname cc '!PATCON/case_control/sas_data/fixed';
libname cclocal '/TB14/TB14/immigrans/store/cc_sandbox/sas_data';


/* Import fragment and junction counts in wide format and convert to tall and save */

/* Fragment counts */

     data WORK.FRAG_TALL    ;
     %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
     infile '/TB14/TB14/immigrans/store/cc_sandbox/event_analysis_output/hpc/bri_cc_fragment_counts_tall.tsv'
 delimiter='09'x MISSOVER DSD lrecl=32767 firstobs=2 ;
        informat event_id $13. ;
        informat sample_id $22. ;
        informat apn best32. ;
        format event_id $13. ;
        format sample_id $22. ;
      	format apn best12. ;
     input
                 event_id $
     Sample_id $
                 apn
     ;
     if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
     run;

/* Junction counts - 425 CTL */

    data WORK.JUNC_CTL425_WIDE    ;
    %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
    infile
'/TB14/TB14/immigrans/SAS_WRK1/cc_sandbox/event_analysis_output/bri_cc_junction_counts_wide_425_control.tsv'
delimiter='09'x MISSOVER DSD lrecl=32767 firstobs=2 ;
       informat event_id $16. ;       informat Sample_P2_B1_50066425 best32. ;       
informat Sample_P2_G4_26089425 best32. ;       informat Sample_P4_F10_51558425 best32. ;
       informat Sample_P4_A12_31588425 best32. ;       informat Sample_P4_H11_17451425 best32. ;
       informat Sample_P2_H2_11333425 best32. ;       informat Sample_P4_G13_28547425 best32. ;
       informat Sample_P2_F1_60935425 best32. ;       informat Sample_P4_H7_76979425 best32. ;
       informat Sample_P4_C12_24886425 best32. ;       informat Sample_P4_A9_30214425 best32. ;
      informat Sample_P2_C5_23619425 best32. ;       informat Sample_P2_C2_49925425 best32. ;   
    informat Sample_P2_F6_63435425 best32. ;       informat Sample_P6_D3_24509425 best32. ;     
  informat Sample_P2_D1_43495425 best32. ;       informat Sample_P2_F4_84251425 best32. ;       
informat Sample_P2_B5_16970425 best32. ;       informat Sample_P2_E5_67530425 best32. ;     
  informat Sample_P2_A1_97721425 best32. ;       informat Sample_P4_D8_39397425 best32. ;   
    informat Sample_P6_E3_20952425 best32. ;       informat Sample_P4_B12_48429425 best32. ;  
     informat Sample_P6_E2_48902425 best32. ;       informat Sample_P4_B9_32639425 best32. ;  
     informat Sample_P4_H10_48429425 best32. ;        informat Sample_P4_E9_87886425 best32. ; 
       informat Sample_P6_D4_85413425 best32. ;        informat Sample_P2_A6_29941425 best32. ;
 informat Sample_P2_F2_95777425 best32. ; informat Sample_P4_G9_37099425 best32. ; 
informat Sample_P2_E3_20773425 best32. ; informat Sample_P2_B2_70252425 best32. ; 
informat Sample_P6_D1_65486425 best32. ; informat Sample_P2_C4_25560425 best32. ; 
informat Sample_P4_A8_23429425 best32. ; informat Sample_P2_D2_29052425 best32. ; 
informat Sample_P2_C1_77517425 best32. ; informat Sample_P2_D3_25527425 best32. ; 
informat Sample_P4_C11_46780425 best32. ; informat Sample_P5_C8_72730425 best32. ;
 informat Sample_P2_E2_59869425 best32. ; informat Sample_P6_C3_53681425 best32. ; 
informat Sample_P5_A7_33881425 best32. ; informat Sample_P2_A4_21854425 best32. ; 
informat Sample_P4_A10_66496425 best32. ; informat Sample_P2_F3_20952425 best32. ; 
format event_id $16. ; format Sample_P2_B1_50066425 best12. ; format Sample_P2_G4_26089425 best12. ;
 format Sample_P4_F10_51558425 best12. ; format Sample_P4_A12_31588425 best12. ;
 format Sample_P4_H11_17451425 best12. ; format Sample_P2_H2_11333425 best12. ; 
format Sample_P4_G13_28547425 best12. ; format Sample_P2_F1_60935425 best12. ; 
format Sample_P4_H7_76979425 best12. ; format Sample_P4_C12_24886425 best12. ;
 format Sample_P4_A9_30214425 best12. ; format Sample_P2_C5_23619425 best12. ;
 format Sample_P2_C2_49925425 best12. ; format Sample_P2_F6_63435425 best12. ; 
format Sample_P6_D3_24509425 best12. ; format Sample_P2_D1_43495425 best12. ; 
format Sample_P2_F4_84251425 best12. ; format Sample_P2_B5_16970425 best12. ; 
format Sample_P2_E5_67530425 best12. ;   format Sample_P2_A1_97721425 best12. ;
   format Sample_P4_D8_39397425 best12. ;   format Sample_P6_E3_20952425 best12. ;  
 format Sample_P4_B12_48429425 best12. ;   format Sample_P6_E2_48902425 best12. ; 
  format Sample_P4_B9_32639425 best12. ;   format Sample_P4_H10_48429425 best12. ; 
  format Sample_P4_E9_87886425 best12. ;   format Sample_P6_D4_85413425 best12. ;   
format Sample_P2_A6_29941425 best12. ;   format Sample_P2_F2_95777425 best12. ;   
format Sample_P4_G9_37099425 best12. ;   format Sample_P2_E3_20773425 best12. ;  
 format Sample_P2_B2_70252425 best12. ;   format Sample_P6_D1_65486425 best12. ;   
format Sample_P2_C4_25560425 best12. ;   format Sample_P4_A8_23429425 best12. ;   
format Sample_P2_D2_29052425 best12. ;   format Sample_P2_C1_77517425 best12. ;   
format Sample_P2_D3_25527425 best12. ;   format Sample_P4_C11_46780425 best12. ; 
  format Sample_P5_C8_72730425 best12. ;   format Sample_P2_E2_59869425 best12. ; 
  format Sample_P6_C3_53681425 best12. ;   format Sample_P5_A7_33881425 best12. ;  
 format Sample_P2_A4_21854425 best12. ;   format Sample_P4_A10_66496425 best12. ;  
 format Sample_P2_F3_20952425 best12. ;
input event_id $ Sample_P2_B1_50066425 Sample_P2_G4_26089425 Sample_P4_F10_51558425 Sample_P4_A12_31588425
 Sample_P4_H11_17451425 Sample_P2_H2_11333425 Sample_P4_G13_28547425 Sample_P2_F1_60935425 Sample_P4_H7_76979425
 Sample_P4_C12_24886425 Sample_P4_A9_30214425 Sample_P2_C5_23619425 Sample_P2_C2_49925425 Sample_P2_F6_63435425
 Sample_P6_D3_24509425 Sample_P2_D1_43495425 Sample_P2_F4_84251425 Sample_P2_B5_16970425 Sample_P2_E5_67530425
 Sample_P2_A1_97721425 Sample_P4_D8_39397425 Sample_P6_E3_20952425 Sample_P4_B12_48429425 Sample_P6_E2_48902425
 Sample_P4_B9_32639425 Sample_P4_H10_48429425 Sample_P4_E9_87886425 Sample_P6_D4_85413425 Sample_P2_A6_29941425
 Sample_P2_F2_95777425 Sample_P4_G9_37099425 Sample_P2_E3_20773425 Sample_P2_B2_70252425 Sample_P6_D1_65486425
 Sample_P2_C4_25560425 Sample_P4_A8_23429425 Sample_P2_D2_29052425 Sample_P2_C1_77517425 Sample_P2_D3_25527425
 Sample_P4_C11_46780425 Sample_P5_C8_72730425 Sample_P2_E2_59869425 Sample_P6_C3_53681425 Sample_P5_A7_33881425
 Sample_P2_A4_21854425 Sample_P4_A10_66496425 Sample_P2_F3_20952425 ;
if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
run;


/* Junction counts - 425 T1D */

     data WORK.JUNC_T1D425_WIDE    ;
     %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
     infile
 '/TB14/TB14/immigrans/SAS_WRK1/cc_sandbox/event_analysis_output/bri_cc_junction_counts_wide_425_cases.tsv'
 delimiter='09'x MISSOVER DSD lrecl=32767 firstobs=2 ;
        informat event_id $16. ;        informat Sample_P4_D10_10109425 best32. ;        informat Sample_P4_E7_1641425 best32. ;        informat Sample_P6_A3_63105425 best32. ;        informat Sample_P2_B4_86451425 best32. ;        informat Sample_P2_A5_95508425 best32. ;        informat Sample_P2_C3_77594425 best32. ;        informat Sample_P2_A2_33161425 best32. ;        informat Sample_P6_A1_12286425 best32. ;        informat Sample_P4_D11_74450425 best32. ;        informat Sample_P2_D5_1018425 best32. ;        informat Sample_P2_H1_35307425 best32. ;        informat Sample_P4_C8_92302425 best32. ;        informat Sample_P4_B11_80342425 best32. ;        informat Sample_P2_B3_27902425 best32. ;        informat Sample_P2_A3_84593425 best32. ;        informat Sample_P6_F2_88536425 best32. ;        informat Sample_P2_H4_86185425 best32. ;        informat Sample_P4_H9_69075425 best32. ;        informat Sample_P2_D6_77109425 best32. ;        informat Sample_P4_F12_9648425 best32. ;        informat Sample_P2_B6_32193425 best32. ;        informat Sample_P2_D4_27959425 best32. ;        informat Sample_P6_E1_70009425 best32. ;        informat Sample_P4_G11_30502425 best32. ;        informat Sample_P4_F9_16654425 best32. ;        informat Sample_P4_G7_44881425 best32. ;        informat Sample_P4_D12_92702425 best32. ;        informat Sample_P4_G10_72894425 best32. ;        informat Sample_P2_G1_89130425 best32. ;        informat Sample_P4_D7_92179425 best32. ;        informat Sample_P6_B2_87313425 best32. ;        informat Sample_P4_E10_51768425 best32. ;        informat Sample_P4_E11_75147425 best32. ;        informat Sample_P4_B10_56485425 best32. ;        informat Sample_P4_E12_41581425 best32. ;        informat Sample_P2_E1_11083425 best32. ;        informat Sample_P4_F8_2326425 best32. ;        informat Sample_P6_C2_41409425 best32. ;        informat Sample_P2_G6_72468425 best32. ;        informat Sample_P4_F11_87313425 best32. ;        informat Sample_P2_H5_84738425 best32. ;        informat Sample_P6_A2_79550425 best32. ;        informat Sample_P5_B8_62765425 best32. ;        informat Sample_P4_B8_18894425 best32. ;        informat Sample_P4_D9_93515425 best32. ;        informat Sample_P2_H6_77571425 best32. ;        informat Sample_P4_H12_96401425 best32. ;        informat Sample_P4_C9_57963425 best32. ;        informat Sample_P6_B3_63925425 best32. ;        informat Sample_P2_E4_21239425 best32. ;        informat Sample_P2_E6_31241425 best32. ;        informat Sample_P2_F5_98838425 best32. ;        informat Sample_P4_G8_16359425 best32. ;        informat Sample_P4_H8_4251425 best32. ;        informat Sample_P4_A11_16020425 best32. ;        informat Sample_P4_E8_60525425 best32. ;        informat Sample_P6_B1_28588425 best32. ;

format event_id $16. ;        format Sample_P4_D10_10109425 best12. ;        format Sample_P4_E7_1641425 best12. ;        format Sample_P6_A3_63105425 best12. ;        format Sample_P2_B4_86451425 best12. ;        format Sample_P2_A5_95508425 best12. ;        format Sample_P2_C3_77594425 best12. ;        format Sample_P2_A2_33161425 best12. ;        format Sample_P6_A1_12286425 best12. ;        format Sample_P4_D11_74450425 best12. ;        format Sample_P2_D5_1018425 best12. ;        format Sample_P2_H1_35307425 best12. ;        format Sample_P4_C8_92302425 best12. ;        format Sample_P4_B11_80342425 best12. ;        format Sample_P2_B3_27902425 best12. ;        format Sample_P2_A3_84593425 best12. ;        format Sample_P6_F2_88536425 best12. ;        format Sample_P2_H4_86185425 best12. ;        format Sample_P4_H9_69075425 best12. ;        format Sample_P2_D6_77109425 best12. ;        format Sample_P4_F12_9648425 best12. ;         format Sample_P2_B6_32193425 best12. ;         format Sample_P2_D4_27959425 best12. ;         format Sample_P6_E1_70009425 best12. ;         format Sample_P4_G11_30502425 best12. ;         format Sample_P4_F9_16654425 best12. ;         format Sample_P4_G7_44881425 best12. ;         format Sample_P4_D12_92702425 best12. ;         format Sample_P4_G10_72894425 best12. ;         format Sample_P2_G1_89130425 best12. ;         format Sample_P4_D7_92179425 best12. ;         format Sample_P6_B2_87313425 best12. ;         format Sample_P4_E10_51768425 best12. ;         format Sample_P4_E11_75147425 best12. ;         format Sample_P4_B10_56485425 best12. ;         format Sample_P4_E12_41581425 best12. ;         format Sample_P2_E1_11083425 best12. ;         format Sample_P4_F8_2326425 best12. ;         format Sample_P6_C2_41409425 best12. ;         format Sample_P2_G6_72468425 best12. ;         format Sample_P4_F11_87313425 best12. ;         format Sample_P2_H5_84738425 best12. ;         format Sample_P6_A2_79550425 best12. ;         format Sample_P5_B8_62765425 best12. ;         format Sample_P4_B8_18894425 best12. ;         format Sample_P4_D9_93515425 best12. ;         format Sample_P2_H6_77571425 best12. ;         format Sample_P4_H12_96401425 best12. ;         format Sample_P4_C9_57963425 best12. ;         format Sample_P6_B3_63925425 best12. ;         format Sample_P2_E4_21239425 best12. ;         format Sample_P2_E6_31241425 best12. ;         format Sample_P2_F5_98838425 best12. ;         format Sample_P4_G8_16359425 best12. ;         format Sample_P4_H8_4251425 best12. ;         format Sample_P4_A11_16020425 best12. ;         format Sample_P4_E8_60525425 best12. ;         format Sample_P6_B1_28588425 best12. ;

input                 event_id $                  Sample_P4_D10_1010945                  Sample_P4_E7_1641425                  Sample_P6_A3_63105425              Sample_P2_B4_86451425              Sample_P2_A5_95508425              Sample_P2_C3_77594425              Sample_P2_A2_33161425              Sample_P6_A1_12286425              Sample_P4_D11_74450425              Sample_P2_D5_1018425              Sample_P2_H1_35307425              Sample_P4_C8_92302425              Sample_P4_B11_80342425              Sample_P2_B3_27902425              Sample_P2_A3_84593425              Sample_P6_F2_88536425              Sample_P2_H4_86185425              Sample_P4_H9_69075425              Sample_P2_D6_77109425              Sample_P4_F12_9648425              Sample_P2_B6_32193425              Sample_P2_D4_27959425              Sample_P6_E1_70009425              Sample_P4_G11_30502425              Sample_P4_F9_16654425              Sample_P4_G7_44881425              Sample_P4_D12_92702425              Sample_P4_G10_72894425              Sample_P2_G1_89130425              Sample_P4_D7_92179425              Sample_P6_B2_87313425              Sample_P4_E10_51768425              Sample_P4_E11_75147425              Sample_P4_B10_56485425              Sample_P4_E12_41581425              Sample_P2_E1_11083425              Sample_P4_F8_2326425              Sample_P6_C2_41409425              Sample_P2_G6_72468425              Sample_P4_F11_87313425              Sample_P2_H5_84738425              Sample_P6_A2_79550425              Sample_P5_B8_62765425              Sample_P4_B8_18894425              Sample_P4_D9_93515425               Sample_P2_H6_77571425               Sample_P4_H12_96401425               Sample_P4_C9_57963425               Sample_P6_B3_63925425               Sample_P2_E4_21239425               Sample_P2_E6_31241425               Sample_P2_F5_98838425               Sample_P4_G8_16359425               Sample_P4_H8_4251425               Sample_P4_A11_16020425               Sample_P4_E8_60525425               Sample_P6_B1_28588425   ;
   if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
   run;

/* Junction counts - 426 CTL */

     data WORK.JUNC_CTL426_WIDE    ;
     %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
     infile
 '/TB14/TB14/immigrans/SAS_WRK1/cc_sandbox/event_analysis_output/bri_cc_junction_counts_wide_426_control.tsv'
 delimiter='09'x MISSOVER DSD lrecl=32767 firstobs=2 ;
        informat event_id $16. ;
        informat Sample_P3_D5_26349426 best32. ;
        informat Sample_P4_H4_10788426 best32. ;
        informat Sample_P3_C9_29941426 best32. ;
        informat Sample_P4_H5_53681426 best32. ;
        informat Sample_P4_A3_12290426 best32. ;
        informat Sample_P4_E4_11661426 best32. ;
        informat Sample_P3_G1_47601426 best32. ;
        informat Sample_P4_E5_81333426 best32. ;
        informat Sample_P3_E6_84251426 best32. ;
        informat Sample_P3_B10_4385426 best32. ;
        informat Sample_P3_D6_33881426 best32. ;
        informat Sample_P3_F6_26089426 best32. ;
        informat Sample_P3_E9_22378426 best32. ;
        informat Sample_P3_F7_57113426 best32. ;
        informat Sample_P3_C1_50006426 best32. ;
        informat Sample_P3_E12_87886426 best32. ;
        informat Sample_P4_G6_78636426 best32. ;
        informat Sample_P3_A3_59869426 best32. ;
        informat Sample_P4_G3_25634426 best32. ;
        informat Sample_P3_C12_17573426 best32. ;
        informat Sample_P3_H6_27268426 best32. ;
        informat Sample_P4_G4_96988426 best32. ;
        informat Sample_P3_A9_12735426 best32. ;
        informat Sample_P4_B5_81816426 best32. ;
        informat Sample_P4_H6_97955426 best32. ;
        informat Sample_P4_B4_34096426 best32. ;
        informat Sample_P3_B3_95777426 best32. ;
        informat Sample_P3_E4_20773426 best32. ;
        informat Sample_P3_G9_17615426 best32. ;
        informat Sample_P6_C5_78636426 best32. ;
        informat Sample_P3_F4_20952426 best32. ;
        informat Sample_P3_G3_11333426 best32. ;
        informat Sample_P3_H2_49925426 best32. ;
        informat Sample_P3_G6_85413426 best32. ;
        informat Sample_P4_C4_12528426 best32. ;
        informat Sample_P4_D5_26395426 best32. ;
       informat Sample_P3_A7_75705426 best32. ;
       informat Sample_P3_D7_23619426 best32. ;
       informat Sample_P4_E3_51558426 best32. ;
       informat Sample_P3_G2_48902426 best32. ;
       informat Sample_P3_D10_40585426 best32. ;
       informat Sample_P3_E3_58652426 best32. ;
       informat Sample_P3_F2_70252426 best32. ;
       informat Sample_P4_C6_4017426 best32. ;
       informat Sample_P3_A2_28547426 best32. ;
       informat Sample_P3_H10_23429426 best32. ;
       informat Sample_P4_A6_22121426 best32. ;
       informat Sample_P3_A10_37683426 best32. ;
       informat Sample_P4_F3_17451426 best32. ;
       informat Sample_P4_F5_24886426 best32. ;
       informat Sample_P3_H4_72730426 best32. ;
       informat Sample_P3_A6_25560426 best32. ;
       informat Sample_P4_F7_37099426 best32. ;
       informat Sample_P3_G10_25755426 best32. ;
       informat Sample_P3_G7_42742426 best32. ;
       informat Sample_P3_D4_25527426 best32. ;
       informat Sample_P4_D6_72636426 best32. ;
       informat Sample_P4_C7_57113426 best32. ;
       informat Sample_P3_B5_76432426 best32. ;
       informat Sample_P6_H2_25560426 best32. ;
       informat Sample_P4_H3_64359426 best32. ;
       informat Sample_P3_H8_22378426 best32. ;
       informat Sample_P3_D1_48966426 best32. ;
       informat Sample_P3_B1_29014426 best32. ;
       informat Sample_P3_E5_21854426 best32. ;
       informat Sample_P3_B4_37683426 best32. ;
       informat Sample_P4_D4_76022426 best32. ;
       informat Sample_P3_C2_60935426 best32. ;
       informat Sample_P3_E7_16970426 best32. ;
       informat Sample_P3_E10_63435426 best32. ;
       informat Sample_P3_H7_76104426 best32. ;
       informat Sample_P6_G1_75705426 best32. ;
       informat Sample_P5_F10_48429426 best32. ;
       informat Sample_P4_A5_63145426 best32. ;
       informat Sample_P3_D3_25527426 best32. ;
       format event_id $16. ;
       format Sample_P3_D5_26349426 best12. ;
       format Sample_P4_H4_10788426 best12. ;
        format Sample_P3_C9_29941426 best12. ;
        format Sample_P4_H5_53681426 best12. ;
        format Sample_P4_A3_12290426 best12. ;
        format Sample_P4_E4_11661426 best12. ;
        format Sample_P3_G1_47601426 best12. ;
        format Sample_P4_E5_81333426 best12. ;
        format Sample_P3_E6_84251426 best12. ;
        format Sample_P3_B10_4385426 best12. ;
        format Sample_P3_D6_33881426 best12. ;
        format Sample_P3_F6_26089426 best12. ;
        format Sample_P3_E9_22378426 best12. ;
        format Sample_P3_F7_57113426 best12. ;
        format Sample_P3_C1_50006426 best12. ;
        format Sample_P3_E12_87886426 best12. ;
        format Sample_P4_G6_78636426 best12. ;
        format Sample_P3_A3_59869426 best12. ;
        format Sample_P4_G3_25634426 best12. ;
        format Sample_P3_C12_17573426 best12. ;
        format Sample_P3_H6_27268426 best12. ;
        format Sample_P4_G4_96988426 best12. ;
        format Sample_P3_A9_12735426 best12. ;
        format Sample_P4_B5_81816426 best12. ;
        format Sample_P4_H6_97955426 best12. ;
        format Sample_P4_B4_34096426 best12. ;
        format Sample_P3_B3_95777426 best12. ;
        format Sample_P3_E4_20773426 best12. ;
        format Sample_P3_G9_17615426 best12. ;
        format Sample_P6_C5_78636426 best12. ;
        format Sample_P3_F4_20952426 best12. ;
        format Sample_P3_G3_11333426 best12. ;
        format Sample_P3_H2_49925426 best12. ;
        format Sample_P3_G6_85413426 best12. ;
        format Sample_P4_C4_12528426 best12. ;
        format Sample_P4_D5_26395426 best12. ;
        format Sample_P3_A7_75705426 best12. ;
        format Sample_P3_D7_23619426 best12. ;
        format Sample_P4_E3_51558426 best12. ;
        format Sample_P3_G2_48902426 best12. ;
        format Sample_P3_D10_40585426 best12. ;
        format Sample_P3_E3_58652426 best12. ;
        format Sample_P3_F2_70252426 best12. ;
        format Sample_P4_C6_4017426 best12. ;
       format Sample_P3_A2_28547426 best12. ;
       format Sample_P3_H10_23429426 best12. ;
       format Sample_P4_A6_22121426 best12. ;
       format Sample_P3_A10_37683426 best12. ;
       format Sample_P4_F3_17451426 best12. ;
       format Sample_P4_F5_24886426 best12. ;
       format Sample_P3_H4_72730426 best12. ;
       format Sample_P3_A6_25560426 best12. ;
       format Sample_P4_F7_37099426 best12. ;
       format Sample_P3_G10_25755426 best12. ;
       format Sample_P3_G7_42742426 best12. ;
       format Sample_P3_D4_25527426 best12. ;
       format Sample_P4_D6_72636426 best12. ;
       format Sample_P4_C7_57113426 best12. ;
       format Sample_P3_B5_76432426 best12. ;
       format Sample_P6_H2_25560426 best12. ;
       format Sample_P4_H3_64359426 best12. ;
       format Sample_P3_H8_22378426 best12. ;
       format Sample_P3_D1_48966426 best12. ;
       format Sample_P3_B1_29014426 best12. ;
       format Sample_P3_E5_21854426 best12. ;
       format Sample_P3_B4_37683426 best12. ;
       format Sample_P4_D4_76022426 best12. ;
       format Sample_P3_C2_60935426 best12. ;
       format Sample_P3_E7_16970426 best12. ;
       format Sample_P3_E10_63435426 best12. ;
       format Sample_P3_H7_76104426 best12. ;
       format Sample_P6_G1_75705426 best12. ;
       format Sample_P5_F10_48429426 best12. ;
       format Sample_P4_A5_63145426 best12. ;
       format Sample_P3_D3_25527426 best12. ;
    input
                event_id $
                Sample_P3_D5_26349426
                Sample_P4_H4_10788426
                Sample_P3_C9_29941426
                Sample_P4_H5_53681426
                Sample_P4_A3_12290426
                Sample_P4_E4_11661426
                Sample_P3_G1_47601426
                Sample_P4_E5_81333426
                Sample_P3_E6_84251426
              Sample_P3_B10_4385426
              Sample_P3_D6_33881426
              Sample_P3_F6_26089426
              Sample_P3_E9_22378426
              Sample_P3_F7_57113426
              Sample_P3_C1_50006426
              Sample_P3_E12_87886426
              Sample_P4_G6_78636426
              Sample_P3_A3_59869426
              Sample_P4_G3_25634426
              Sample_P3_C12_17573426
              Sample_P3_H6_27268426
              Sample_P4_G4_96988426
              Sample_P3_A9_12735426
              Sample_P4_B5_81816426
              Sample_P4_H6_97955426
              Sample_P4_B4_34096426
              Sample_P3_B3_95777426
              Sample_P3_E4_20773426
              Sample_P3_G9_17615426
              Sample_P6_C5_78636426
              Sample_P3_F4_20952426
              Sample_P3_G3_11333426
              Sample_P3_H2_49925426
              Sample_P3_G6_85413426
              Sample_P4_C4_12528426
              Sample_P4_D5_26395426
              Sample_P3_A7_75705426
              Sample_P3_D7_23619426
              Sample_P4_E3_51558426
              Sample_P3_G2_48902426
              Sample_P3_D10_40585426
              Sample_P3_E3_58652426
              Sample_P3_F2_70252426
              Sample_P4_C6_4017426
              Sample_P3_A2_28547426
              Sample_P3_H10_23429426
              Sample_P4_A6_22121426
              Sample_P3_A10_37683426
              Sample_P4_F3_17451426
              Sample_P4_F5_24886426
              Sample_P3_H4_72730426
              Sample_P3_A6_25560426
              Sample_P4_F7_37099426
              Sample_P3_G10_25755426
              Sample_P3_G7_42742426
              Sample_P3_D4_25527426
              Sample_P4_D6_72636426
              Sample_P4_C7_57113426
              Sample_P3_B5_76432426
              Sample_P6_H2_25560426
              Sample_P4_H3_64359426
              Sample_P3_H8_22378426
              Sample_P3_D1_48966426
              Sample_P3_B1_29014426
              Sample_P3_E5_21854426
              Sample_P3_B4_37683426
              Sample_P4_D4_76022426
              Sample_P3_C2_60935426
              Sample_P3_E7_16970426
              Sample_P3_E10_63435426
              Sample_P3_H7_76104426
              Sample_P6_G1_75705426
              Sample_P5_F10_48429426
              Sample_P4_A5_63145426
              Sample_P3_D3_25527426  ;

  if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
  run;



/* Junction counts - 426 T1D */

    data WORK.JUNC_T1D426_WIDE    ;
    %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
    infile
'/TB14/TB14/immigrans/SAS_WRK1/cc_sandbox/event_analysis_output/bri_cc_junction_counts_wide_426_cases.tsv'
delimiter='09'x MISSOVER DSD lrecl=32767 firstobs=2 ;
       informat event_id $16. ;
       informat Sample_P6_F3_90726426 best32. ;
       informat Sample_P4_B6_57951426 best32. ;
       informat Sample_P3_G5_86451426 best32. ;
       informat Sample_P4_C10_2326426 best32. ;
       informat Sample_P4_A7_1694426 best32. ;
       informat Sample_P3_D2_35307426 best32. ;
       informat Sample_P4_D3_87313426 best32. ;
       informat Sample_P4_B7_16654426 best32. ;
       informat Sample_P3_C7_95508426 best32. ;
       informat Sample_P4_D1_99706426 best32. ;
       informat Sample_P6_C4_83970426 best32. ;
       informat Sample_P3_F9_32193426 best32. ;
       informat Sample_P4_F4_73605426 best32. ;
       informat Sample_P3_A5_62765426 best32. ;
       informat Sample_P3_C4_77594426 best32. ;
       informat Sample_P3_B9_84738426 best32. ;
       informat Sample_P4_F6_56485426 best32. ;
       informat Sample_P3_F1_65680426 best32. ;
       informat Sample_P3_D9_39396426 best32. ;
       informat Sample_P3_F12_41409426 best32. ;
       informat Sample_P4_C1_72894426 best32. ;
       informat Sample_P4_E2_78111426 best32. ;
       informat Sample_P3_F10_72468426 best32. ;
       informat Sample_P4_H1_70009426 best32. ;
       informat Sample_P3_H3_84593426 best32. ;
       informat Sample_P3_B2_11083426 best32. ;
       informat Sample_P6_A4_4251426 best32. ;
       informat Sample_P4_D2_97288426 best32. ;
       informat Sample_P3_A4_27902426 best32. ;
       informat Sample_P4_A1_43138426 best32. ;
       informat Sample_P3_C5_83970426 best32. ;
       informat Sample_P3_H1_45311426 best32. ;
      informat Sample_P3_E2_33161426 best32. ;
      informat Sample_P3_H9_77109426 best32. ;
      informat Sample_P6_H1_92179426 best32. ;
      informat Sample_P3_A1_92179426 best32. ;
      informat Sample_P4_C2_16020426 best32. ;
      informat Sample_P3_G13_43594426 best32. ;
      informat Sample_P4_F1_84189426 best32. ;
      informat Sample_P4_E1_12286426 best32. ;
      informat Sample_P4_B2_70658426 best32. ;
      informat Sample_P6_C6_62579426 best32. ;
      informat Sample_P4_G1_79550426 best32. ;
      informat Sample_P3_C6_21239426 best32. ;
      informat Sample_P4_B1_30502426 best32. ;
      informat Sample_P3_C3_44884426 best32. ;
      informat Sample_P3_H12_18894426 best32. ;
      informat Sample_P3_B6_27959426 best32. ;
      informat Sample_P4_G5_92702426 best32. ;
      informat Sample_P4_A2_47020426 best32. ;
      informat Sample_P3_G4_28992426 best32. ;
      informat Sample_P6_B4_88927426 best32. ;
      informat Sample_P6_H3_88536426 best32. ;
      informat Sample_P3_F5_41724426 best32. ;
      informat Sample_P4_F2_80342426 best32. ;
      informat Sample_P3_H5_82865426 best32. ;
      informat Sample_P3_B7_86185426 best32. ;
      informat Sample_P6_B5_84189426 best32. ;
      informat Sample_P4_E6_69075426 best32. ;
      informat Sample_P4_G2_63105426 best32. ;
      informat Sample_P3_D12_93515426 best32. ;
      informat Sample_P4_H2_5255426 best32. ;
      informat Sample_P3_F3_11699426 best32. ;
      informat Sample_P4_C5_1641426 best32. ;
      informat Sample_P4_B3_74450426 best32. ;
      informat Sample_P4_A4_50316426 best32. ;
      informat Sample_P4_C3_75147426 best32. ;
      informat Sample_P3_C10_4700426 best32. ;
      format event_id $16. ;
      format Sample_P6_F3_90726426 best12. ;
      format Sample_P4_B6_57951426 best12. ;
      format Sample_P3_G5_86451426 best12. ;
      format Sample_P4_C10_2326426 best12. ;
      format Sample_P4_A7_1694426 best12. ;
    format Sample_P3_D2_35307426 best12. ;
    format Sample_P4_D3_87313426 best12. ;
    format Sample_P4_B7_16654426 best12. ;
    format Sample_P3_C7_95508426 best12. ;
    format Sample_P4_D1_99706426 best12. ;
    format Sample_P6_C4_83970426 best12. ;
    format Sample_P3_F9_32193426 best12. ;
    format Sample_P4_F4_73605426 best12. ;
    format Sample_P3_A5_62765426 best12. ;
    format Sample_P3_C4_77594426 best12. ;
    format Sample_P3_B9_84738426 best12. ;
    format Sample_P4_F6_56485426 best12. ;
    format Sample_P3_F1_65680426 best12. ;
    format Sample_P3_D9_39396426 best12. ;
    format Sample_P3_F12_41409426 best12. ;
    format Sample_P4_C1_72894426 best12. ;
    format Sample_P4_E2_78111426 best12. ;
    format Sample_P3_F10_72468426 best12. ;
    format Sample_P4_H1_70009426 best12. ;
    format Sample_P3_H3_84593426 best12. ;
    format Sample_P3_B2_11083426 best12. ;
    format Sample_P6_A4_4251426 best12. ;
    format Sample_P4_D2_97288426 best12. ;
    format Sample_P3_A4_27902426 best12. ;
    format Sample_P4_A1_43138426 best12. ;
    format Sample_P3_C5_83970426 best12. ;
    format Sample_P3_H1_45311426 best12. ;
    format Sample_P3_E2_33161426 best12. ;
    format Sample_P3_H9_77109426 best12. ;
    format Sample_P6_H1_92179426 best12. ;
    format Sample_P3_A1_92179426 best12. ;
    format Sample_P4_C2_16020426 best12. ;
    format Sample_P3_G13_43594426 best12. ;
    format Sample_P4_F1_84189426 best12. ;
    format Sample_P4_E1_12286426 best12. ;
    format Sample_P4_B2_70658426 best12. ;
    format Sample_P6_C6_62579426 best12. ;
    format Sample_P4_G1_79550426 best12. ;
    format Sample_P3_C6_21239426 best12. ;
    format Sample_P4_B1_30502426 best12. ;
    format Sample_P3_C3_44884426 best12. ;
    format Sample_P3_H12_18894426 best12. ;
       format Sample_P3_B6_27959426 best12. ;
       format Sample_P4_G5_92702426 best12. ;
       format Sample_P4_A2_47020426 best12. ;
       format Sample_P3_G4_28992426 best12. ;
       format Sample_P6_B4_88927426 best12. ;
       format Sample_P6_H3_88536426 best12. ;
       format Sample_P3_F5_41724426 best12. ;
       format Sample_P4_F2_80342426 best12. ;
       format Sample_P3_H5_82865426 best12. ;
       format Sample_P3_B7_86185426 best12. ;
       format Sample_P6_B5_84189426 best12. ;
       format Sample_P4_E6_69075426 best12. ;
       format Sample_P4_G2_63105426 best12. ;
       format Sample_P3_D12_93515426 best12. ;
       format Sample_P4_H2_5255426 best12. ;
       format Sample_P3_F3_11699426 best12. ;
       format Sample_P4_C5_1641426 best12. ;
       format Sample_P4_B3_74450426 best12. ;
       format Sample_P4_A4_50316426 best12. ;
       format Sample_P4_C3_75147426 best12. ;
       format Sample_P3_C10_4700426 best12. ;
    input
                event_id $
                Sample_P6_F3_90726426
                Sample_P4_B6_57951426
                Sample_P3_G5_86451426
                Sample_P4_C10_2326426
                Sample_P4_A7_1694426
                Sample_P3_D2_35307426
                Sample_P4_D3_87313426
                Sample_P4_B7_16654426
                Sample_P3_C7_95508426
                Sample_P4_D1_99706426
                Sample_P6_C4_83970426
                Sample_P3_F9_32193426
                Sample_P4_F4_73605426
                Sample_P3_A5_62765426
                Sample_P3_C4_77594426
                Sample_P3_B9_84738426
                Sample_P4_F6_56485426
                Sample_P3_F1_65680426
                Sample_P3_D9_39396426
            Sample_P3_F12_41409426
            Sample_P4_C1_72894426
            Sample_P4_E2_78111426
            Sample_P3_F10_72468426
            Sample_P4_H1_70009426
            Sample_P3_H3_84593426
            Sample_P3_B2_11083426
            Sample_P6_A4_4251426
            Sample_P4_D2_97288426
            Sample_P3_A4_27902426
            Sample_P4_A1_43138426
            Sample_P3_C5_83970426
            Sample_P3_H1_45311426
            Sample_P3_E2_33161426
            Sample_P3_H9_77109426
            Sample_P6_H1_92179426
            Sample_P3_A1_92179426
            Sample_P4_C2_16020426
            Sample_P3_G13_43594426
            Sample_P4_F1_84189426
            Sample_P4_E1_12286426
            Sample_P4_B2_70658426
            Sample_P6_C6_62579426
            Sample_P4_G1_79550426
            Sample_P3_C6_21239426
            Sample_P4_B1_30502426
            Sample_P3_C3_44884426
            Sample_P3_H12_18894426
            Sample_P3_B6_27959426
            Sample_P4_G5_92702426
            Sample_P4_A2_47020426
            Sample_P3_G4_28992426
            Sample_P6_B4_88927426
            Sample_P6_H3_88536426
            Sample_P3_F5_41724426
            Sample_P4_F2_80342426
            Sample_P3_H5_82865426
            Sample_P3_B7_86185426
            Sample_P6_B5_84189426
            Sample_P4_E6_69075426
            Sample_P4_G2_63105426
            Sample_P3_D12_93515426
             Sample_P4_H2_5255426
             Sample_P3_F3_11699426
             Sample_P4_C5_1641426
             Sample_P4_B3_74450426
             Sample_P4_A4_50316426
             Sample_P4_C3_75147426
             Sample_P3_C10_4700426
 ;
 if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
 run;


    data WORK.JUNC_T1D127_WIDE    ;
    %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
    infile
'/TB14/TB14/immigrans/SAS_WRK1/cc_sandbox/event_analysis_output/bri_cc_junction_counts_wide_127_cases.tsv'
delimiter='09'x MISSOVER DSD lrecl=32767 firstobs=2 ;
informat event_id $16. ;
informat Sample_P5_A10_78111127 best32.;
informat Sample_P5_G4_89130127 best32.;
informat Sample_P5_A9_92702127 best32.;
informat Sample_P5_F3_28992127 best32.;
informat Sample_P5_C4_60525127 best32.;
informat Sample_P5_H1_33161127 best32.;
informat Sample_P1_C1_97288 best32.;
informat Sample_P5_F12_27902127 best32.;
informat Sample_P5_G13_77109127 best32.;
informat Sample_P5_E3_37602127 best32.;
informat Sample_P5_D10_92302127 best32.;
informat Sample_P6_F4_70658127 best32.;
informat Sample_P6_E5_64368127 best32.;
informat Sample_P5_B1_97706127 best32.;
informat Sample_P5_H6_62579127 best32.;
informat Sample_P5_C3_4700127 best32.;
informat Sample_P1_G3_11083 best32.;
informat Sample_P6_F6_99706127 best32.;
informat Sample_P5_A1_5255127 best32.;
informat Sample_P5_D3_4251127 best32.;
informat Sample_P5_G2_45311127 best32.;
informat Sample_P1_F3_16359 best32.;
informat Sample_P5_F8_90726127 best32.;
informat Sample_P5_E6_77571127 best32.;
informat Sample_P1_B1_47020 best32.;
informat Sample_P5_B5_75147127 best32.;
informat Sample_P5_C11_31241111 best32.;
informat Sample_P5_E5_92302111 best32.;
informat Sample_P1_H2_63928 best32.;
informat Sample_P1_E1_63333 best32.;
informat Sample_P5_B10_73605127 best32.;
informat Sample_P5_E10_7419127 best32.;
informat Sample_P1_G1_74450 best32.;
informat Sample_P5_D6_69075127 best32.;
informat Sample_P5_G1_35307127 best32.;
informat Sample_P5_H7_83970127 best32.;
informat Sample_P1_A1_57963 best32.;
informat Sample_P5_E11_32193127 best32.;
informat Sample_P5_H12_82865127 best32.;
informat Sample_P5_A6_88536127 best32.;
informat Sample_P1_D1_80342 best32.;
informat Sample_P6_H4_16654127 best32.;
informat Sample_P5_E9_44884127 best32.;
informat Sample_P5_B7_1018127 best32.;
informat Sample_P6_G6_28588127 best32.;
informat Sample_P5_G6_21239127 best32.;
informat Sample_P1_G2_84738 best32.;
informat Sample_P5_D5_43155127 best32.;
informat Sample_P5_A4_99706127 best32.;
informat Sample_P5_E1_1641127 best32.;
informat Sample_P5_H11_86451127 best32.;
informat Sample_P5_F4_88927127 best32.;
informat Sample_P6_D6_86451127 best32.;
informat Sample_P5_A2_43138127 best32.;
informat Sample_P5_C9_31241127 best32.;
informat Sample_P5_F5_65680127 best32.;
informat Sample_P5_E2_43594127 best32.;

format event_id $16. ;
format Sample_P5_A10_78111127 best12.;
format Sample_P5_G4_89130127 best12.;
format Sample_P5_A9_92702127 best12.;
format Sample_P5_F3_28992127 best12.;
format Sample_P5_C4_60525127 best12.;
format Sample_P5_H1_33161127 best12.;
format Sample_P1_C1_97288 best12.;
format Sample_P5_F12_27902127 best12.;
format Sample_P5_G13_77109127 best12.;
format Sample_P5_E3_37602127 best12.;
format Sample_P5_D10_92302127 best12.;
format Sample_P6_F4_70658127 best12.;
format Sample_P6_E5_64368127 best12.;
format Sample_P5_B1_97706127 best12.;
format Sample_P5_H6_62579127 best12.;
format Sample_P5_C3_4700127 best12.;
format Sample_P1_G3_11083 best12.;
format Sample_P6_F6_99706127 best12.;
format Sample_P5_A1_5255127 best12.;
format Sample_P5_D3_4251127 best12.;
format Sample_P5_G2_45311127 best12.;
format Sample_P1_F3_16359 best12.;
format Sample_P5_F8_90726127 best12.;
format Sample_P5_E6_77571127 best12.;
format Sample_P1_B1_47020 best12.;
format Sample_P5_B5_75147127 best12.;
format Sample_P5_C11_31241111 best12.;
format Sample_P5_E5_92302111 best12.;
format Sample_P1_H2_63928 best12.;
format Sample_P1_E1_63333 best12.;
format Sample_P5_B10_73605127 best12.;
format Sample_P5_E10_7419127 best12.;
format Sample_P1_G1_74450 best12.;
format Sample_P5_D6_69075127 best12.;
format Sample_P5_G1_35307127 best12.;
format Sample_P5_H7_83970127 best12.;
format Sample_P1_A1_57963 best12.;
format Sample_P5_E11_32193127 best12.;
format Sample_P5_H12_82865127 best12.;
format Sample_P5_A6_88536127 best12.;
format Sample_P1_D1_80342 best12.;
format Sample_P6_H4_16654127 best12.;
format Sample_P5_E9_44884127 best12.;
format Sample_P5_B7_1018127 best12.;
format Sample_P6_G6_28588127 best12.;
format Sample_P5_G6_21239127 best12.;
format Sample_P1_G2_84738 best12.;
format Sample_P5_D5_43155127 best12.;
format Sample_P5_A4_99706127 best12.;
format Sample_P5_E1_1641127 best12.;
format Sample_P5_H11_86451127 best12.;
format Sample_P5_F4_88927127 best12.;
format Sample_P6_D6_86451127 best12.;
format Sample_P5_A2_43138127 best12.;
format Sample_P5_C9_31241127 best12.;
format Sample_P5_F5_65680127 best12.;
format Sample_P5_E2_43594127 best12.;
input 
event_id $ Sample_P5_A10_78111127 Sample_P5_G4_89130127 Sample_P5_A9_92702127 
Sample_P5_F3_28992127 Sample_P5_C4_60525127 Sample_P5_H1_33161127 Sample_P1_C1_97288
Sample_P5_F12_27902127 Sample_P5_G13_77109127 Sample_P5_E3_37602127 Sample_P5_D10_92302127
Sample_P6_F4_70658127 Sample_P6_E5_64368127 Sample_P5_B1_97706127 Sample_P5_H6_62579127 
Sample_P5_C3_4700127 Sample_P1_G3_11083 Sample_P6_F6_99706127 Sample_P5_A1_5255127 Sample_P5_D3_4251127
 Sample_P5_G2_45311127 Sample_P1_F3_16359 Sample_P5_F8_90726127 Sample_P5_E6_77571127 Sample_P1_B1_47020
 Sample_P5_B5_75147127 Sample_P5_C11_31241111 Sample_P5_E5_92302111 Sample_P1_H2_63928 Sample_P1_E1_63333
 Sample_P5_B10_73605127 Sample_P5_E10_7419127 Sample_P1_G1_74450 Sample_P5_D6_69075127 Sample_P5_G1_35307127
 Sample_P5_H7_83970127 Sample_P1_A1_57963 Sample_P5_E11_32193127 Sample_P5_H12_82865127 Sample_P5_A6_88536127
 Sample_P1_D1_80342 Sample_P6_H4_16654127 Sample_P5_E9_44884127 Sample_P5_B7_1018127 Sample_P6_G6_28588127
 Sample_P5_G6_21239127 Sample_P1_G2_84738 Sample_P5_D5_43155127 Sample_P5_A4_99706127 Sample_P5_E1_1641127
 Sample_P5_H11_86451127 Sample_P5_F4_88927127 Sample_P6_D6_86451127 Sample_P5_A2_43138127 Sample_P5_C9_31241127
 Sample_P5_F5_65680127 Sample_P5_E2_43594127
;
 if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
 run;




    data WORK.JUNC_CTL127_WIDE    ;
    %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
    infile
'/TB14/TB14/immigrans/SAS_WRK1/cc_sandbox/event_analysis_output/bri_cc_junction_counts_wide_127_control.tsv'
delimiter='09'x MISSOVER DSD lrecl=32767 firstobs=2 ;
informat event_id $16. ;
informat Sample_P5_B11_22121127 best32. ;
informat Sample_P5_C12_28425127 best32. ;
informat Sample_P1_A2_31588 best32. ;
informat Sample_P5_F6_57113127 best32. ;
informat Sample_P5_D2_64359127 best32. ;
informat Sample_P5_A8_76022127 best32. ;
informat Sample_P1_F1_46780 best32. ;
informat Sample_P5_C10_91359127 best32. ;
informat Sample_P5_H2_65486127 best32. ;
informat Sample_P5_A3_37250127 best32. ;
informat Sample_P5_G7_97721127 best32. ;
informat Sample_P5_B4_26871127 best32. ;
informat Sample_P1_H1_17451 best32. ;
informat Sample_P5_C1_48902127 best32. ;
informat Sample_P5_C7_24886127 best32. ;
informat Sample_P6_G4_25560127 best32. ;
informat Sample_P6_F5_4017127 best32. ;
informat Sample_P6_H12_11333127 best32. ;
informat Sample_P5_G10_43495127 best32. ;
informat Sample_P5_H9_21854127 best32. ;
informat Sample_P5_D12_76104127 best32. ;
informat Sample_P5_H3_70252127 best32. ;
informat Sample_P5_G8_22378127 best32. ;
informat Sample_P5_F9_81417127 best32. ;
informat Sample_P5_G3_60935127 best32. ;
informat Sample_P5_G11_29014127 best32. ;
informat Sample_P5_C5_26349127 best32. ;
informat Sample_P5_A5_25527127 best32. ;
informat Sample_P6_G3_59869127 best32. ;
informat Sample_P5_B3_81816127 best32. ;
informat Sample_P1_D2_31384 best32. ;
informat Sample_P6_H6_84251127 best32. ;
informat Sample_P6_H5_25560127 best32. ;
informat Sample_P5_B2_12290127 best32. ;
informat Sample_P5_F11_29894127 best32. ;
informat Sample_P5_E4_40585127 best32. ;
informat Sample_P1_E2_33881 best32. ;
informat Sample_P5_E12_17615127 best32. ;
informat Sample_P5_F1_49925127 best32. ;
informat Sample_P5_B6_81900127 best32. ;
informat Sample_P5_H10_33881127 best32. ;
informat Sample_P6_E6_37099127 best32. ;
informat Sample_P6_D5_53681127 best32. ;
informat Sample_P5_D1_25755127 best32. ;
informat Sample_P5_H5_47601127 best32. ;
informat Sample_P5_G9_50006127 best32. ;
informat Sample_P1_B2_78685 best32. ;
informat Sample_P1_C3_23619 best32. ;
informat Sample_P6_A5_95777127 best32. ;
informat Sample_P5_F2_20952127 best32. ;
informat Sample_P5_E7_87886127 best32. ;
informat Sample_P5_D4_23619127 best32. ;
informat Sample_P1_C2_48429 best32. ;
informat Sample_P5_F7_37260127 best32. ;
informat Sample_P5_D9_42742127 best32. ;
informat Sample_P6_E4_29941127 best32. ;
informat Sample_P1_B3_16970 best32. ;
informat Sample_P5_C6_34096127 best32. ;
informat Sample_P5_H4_48966127 best32. ;
informat Sample_P5_C2_20773127 best32. ;
informat Sample_P5_H8_76432127 best32. ;
informat Sample_P5_A11_67530127 best32. ;
format event_id $16. ;
format Sample_P5_B11_22121127 best12. ;
format Sample_P5_C12_28425127 best12. ;
format Sample_P1_A2_31588 best12. ;
format Sample_P5_F6_57113127 best12. ;
format Sample_P5_D2_64359127 best12. ;
format Sample_P5_A8_76022127 best12. ;
format Sample_P1_F1_46780 best12. ;
format Sample_P5_C10_91359127 best12. ;
format Sample_P5_H2_65486127 best12. ;
format Sample_P5_A3_37250127 best12. ;
format Sample_P5_G7_97721127 best12. ;
format Sample_P5_B4_26871127 best12. ;
format Sample_P1_H1_17451 best12. ;
format Sample_P5_C1_48902127 best12. ;
format Sample_P5_C7_24886127 best12. ;
format Sample_P6_G4_25560127 best12. ;
format Sample_P6_F5_4017127 best12. ;
format Sample_P6_H12_11333127 best12. ;
format Sample_P5_G10_43495127 best12. ;
format Sample_P5_H9_21854127 best12. ;
format Sample_P5_D12_76104127 best12. ;
format Sample_P5_H3_70252127 best12. ;
format Sample_P5_G8_22378127 best12. ;
format Sample_P5_F9_81417127 best12. ;
format Sample_P5_G3_60935127 best12. ;
format Sample_P5_G11_29014127 best12. ;
format Sample_P5_C5_26349127 best12. ;
format Sample_P5_A5_25527127 best12. ;
format Sample_P6_G3_59869127 best12. ;
format Sample_P5_B3_81816127 best12. ;
format Sample_P1_D2_31384 best12. ;
format Sample_P6_H6_84251127 best12. ;
format Sample_P6_H5_25560127 best12. ;
format Sample_P5_B2_12290127 best12. ;
format Sample_P5_F11_29894127 best12. ;
format Sample_P5_E4_40585127 best12. ;
format Sample_P1_E2_33881 best12. ;
format Sample_P5_E12_17615127 best12. ;
format Sample_P5_F1_49925127 best12. ;
format Sample_P5_B6_81900127 best12. ;
format Sample_P5_H10_33881127 best12. ;
format Sample_P6_E6_37099127 best12. ;
format Sample_P6_D5_53681127 best12. ;
format Sample_P5_D1_25755127 best12. ;
format Sample_P5_H5_47601127 best12. ;
format Sample_P5_G9_50006127 best12. ;
format Sample_P1_B2_78685 best12. ;
format Sample_P1_C3_23619 best12. ;
format Sample_P6_A5_95777127 best12. ;
format Sample_P5_F2_20952127 best12. ;
format Sample_P5_E7_87886127 best12. ;
format Sample_P5_D4_23619127 best12. ;
format Sample_P1_C2_48429 best12. ;
format Sample_P5_F7_37260127 best12. ;
format Sample_P5_D9_42742127 best12. ;
format Sample_P6_E4_29941127 best12. ;
format Sample_P1_B3_16970 best12. ;
format Sample_P5_C6_34096127 best12. ;
format Sample_P5_H4_48966127 best12. ;
format Sample_P5_C2_20773127 best12. ;
format Sample_P5_H8_76432127 best12. ;
format Sample_P5_A11_67530127 best12. ;
input event_id $  Sample_P5_B11_22121127 Sample_P5_C12_28425127 Sample_P1_A2_31588 Sample_P5_F6_57113127
Sample_P5_D2_64359127 Sample_P5_A8_76022127 Sample_P1_F1_46780 Sample_P5_C10_91359127 Sample_P5_H2_65486127
Sample_P5_A3_37250127 Sample_P5_G7_97721127 Sample_P5_B4_26871127 Sample_P1_H1_17451 Sample_P5_C1_48902127
Sample_P5_C7_24886127 Sample_P6_G4_25560127 Sample_P6_F5_4017127 Sample_P6_H12_11333127 Sample_P5_G10_43495127
Sample_P5_H9_21854127 Sample_P5_D12_76104127 Sample_P5_H3_70252127 Sample_P5_G8_22378127 Sample_P5_F9_81417127
Sample_P5_G3_60935127 Sample_P5_G11_29014127 Sample_P5_C5_26349127 Sample_P5_A5_25527127 Sample_P6_G3_59869127
Sample_P5_B3_81816127 Sample_P1_D2_31384 Sample_P6_H6_84251127 Sample_P6_H5_25560127 Sample_P5_B2_12290127
Sample_P5_F11_29894127 Sample_P5_E4_40585127 Sample_P1_E2_33881 Sample_P5_E12_17615127 Sample_P5_F1_49925127
Sample_P5_B6_81900127 Sample_P5_H10_33881127 Sample_P6_E6_37099127 Sample_P6_D5_53681127 Sample_P5_D1_25755127
Sample_P5_H5_47601127 Sample_P5_G9_50006127 Sample_P1_B2_78685 Sample_P1_C3_23619 Sample_P6_A5_95777127
Sample_P5_F2_20952127 Sample_P5_E7_87886127 Sample_P5_D4_23619127 Sample_P1_C2_48429 Sample_P5_F7_37260127
Sample_P5_D9_42742127 Sample_P6_E4_29941127 Sample_P1_B3_16970 Sample_P5_C6_34096127 Sample_P5_H4_48966127
Sample_P5_C2_20773127 Sample_P5_H8_76432127 Sample_P5_A11_67530127
;
 if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
 run;






    data WORK.JUNC_400_WIDE    ;
    %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
    infile
'/TB14/TB14/immigrans/SAS_WRK1/cc_sandbox/event_analysis_output/bri_cc_junction_counts_wide_400.tsv'
delimiter='09'x MISSOVER DSD lrecl=32767 firstobs=2 ;
informat event_id $16.;
informat Sample_P6_G10_50316400 best32.;
informat Sample_P6_B10_76002400 best32.;
informat Sample_P6_G9_79484400 best32.;
informat Sample_P6_F12_63928400 best32.;
informat Sample_P6_E12_78636400 best32.;
informat Sample_P6_C9_11333400 best32.;
informat Sample_P6_F7_76432400 best32.;
informat Sample_P6_H7_23619400 best32.;
informat Sample_P6_E10_22121400 best32.;
informat Sample_P6_A11_81417400 best32.;
informat Sample_P6_D10_34096400 best32.;
informat Sample_P6_C11_47601400 best32.;
informat Sample_P6_D8_67530400 best32.;
informat Sample_P6_B9_62765400 best32.;
informat Sample_P6_C7_32193400 best32.;
informat Sample_P6_H8_37099400 best32.;
informat Sample_P6_E8_64368400 best32.;
informat Sample_P6_A7_75147400 best32.;
informat Sample_P6_F9_77694400 best32.;
informat Sample_P6_B7_29941400 best32.;
informat Sample_P6_E9_27902400 best32.;
informat Sample_P6_G11_35307400 best32.;
informat Sample_P6_D11_28547400 best32.;
informat Sample_P6_H11_48902400 best32.;
informat Sample_P6_E7_17615400 best32.;
informat Sample_P6_A9_28992400 best32.;
informat Sample_P6_D7_29894400 best32.;
informat Sample_P6_B12_29014400 best32.;
informat Sample_P6_F11_89130400 best32.;
informat Sample_P6_G7_95508400 best32.;
informat Sample_P6_F8_86185400 best32.;
informat Sample_P6_C8_61784400 best32.;
informat Sample_P6_H9_5255400 best32.;
informat Sample_P6_A12_49925400 best32.;
informat Sample_P6_A10_81333400 best32.;
informat Sample_P6_E11_11083400 best32.;
informat Sample_P6_B8_98626400 best32.;
informat Sample_P6_A8_16970400 best32.;
informat Sample_P6_H10_50006400 best32.;
informat Sample_P6_D12_31241400 best32.;
informat Sample_P6_G8_57185400 best32.;
informat Sample_P6_D9_59869400 best32.;
informat Sample_P6_F10_45750400 best32.;
informat Sample_P6_C12_97721400 best32.;
informat Sample_P6_B11_65680400 best32.;
informat Sample_P6_C10_18977400 best32.;
format event_id $16.;
format Sample_P6_G10_50316400 best12.;
format Sample_P6_B10_76002400 best12.;
format Sample_P6_G9_79484400 best12.;
format Sample_P6_F12_63928400 best12.;
format Sample_P6_E12_78636400 best12.;
format Sample_P6_C9_11333400 best12.;
format Sample_P6_F7_76432400 best12.;
format Sample_P6_H7_23619400 best12.;
format Sample_P6_E10_22121400 best12.;
format Sample_P6_A11_81417400 best12.;
format Sample_P6_D10_34096400 best12.;
format Sample_P6_C11_47601400 best12.;
format Sample_P6_D8_67530400 best12.;
format Sample_P6_B9_62765400 best12.;
format Sample_P6_C7_32193400 best12.;
format Sample_P6_H8_37099400 best12.;
format Sample_P6_E8_64368400 best12.;
format Sample_P6_A7_75147400 best12.;
format Sample_P6_F9_77694400 best12.;
format Sample_P6_B7_29941400 best12.;
format Sample_P6_E9_27902400 best12.;
format Sample_P6_G11_35307400 best12.;
format Sample_P6_D11_28547400 best12.;
format Sample_P6_H11_48902400 best12.;
format Sample_P6_E7_17615400 best12.;
format Sample_P6_A9_28992400 best12.;
format Sample_P6_D7_29894400 best12.;
format Sample_P6_B12_29014400 best12.;
format Sample_P6_F11_89130400 best12.;
format Sample_P6_G7_95508400 best12.;
format Sample_P6_F8_86185400 best12.;
format Sample_P6_C8_61784400 best12.;
format Sample_P6_H9_5255400 best12.;
format Sample_P6_A12_49925400 best12.;
format Sample_P6_A10_81333400 best12.;
format Sample_P6_E11_11083400 best12.;
format Sample_P6_B8_98626400 best12.;
format Sample_P6_A8_16970400 best12.;
format Sample_P6_H10_50006400 best12.;
format Sample_P6_D12_31241400 best12.;
format Sample_P6_G8_57185400 best12.;
format Sample_P6_D9_59869400 best12.;
format Sample_P6_F10_45750400 best12.;
format Sample_P6_C12_97721400 best12.;
format Sample_P6_B11_65680400 best12.;
format Sample_P6_C10_18977400 best12.;
input event_id $ Sample_P6_G10_50316400 Sample_P6_B10_76002400 Sample_P6_G9_79484400 Sample_P6_F12_63928400
Sample_P6_E12_78636400 Sample_P6_C9_11333400 Sample_P6_F7_76432400 Sample_P6_H7_23619400 Sample_P6_E10_22121400
Sample_P6_A11_81417400 Sample_P6_D10_34096400 Sample_P6_C11_47601400 Sample_P6_D8_67530400 Sample_P6_B9_62765400
Sample_P6_C7_32193400 Sample_P6_H8_37099400 Sample_P6_E8_64368400 Sample_P6_A7_75147400 Sample_P6_F9_77694400
Sample_P6_B7_29941400 Sample_P6_E9_27902400 Sample_P6_G11_35307400 Sample_P6_D11_28547400 Sample_P6_H11_48902400
Sample_P6_E7_17615400 Sample_P6_A9_28992400 Sample_P6_D7_29894400 Sample_P6_B12_29014400 Sample_P6_F11_89130400
Sample_P6_G7_95508400 Sample_P6_F8_86185400 Sample_P6_C8_61784400 Sample_P6_H9_5255400 Sample_P6_A12_49925400
Sample_P6_A10_81333400 Sample_P6_E11_11083400 Sample_P6_B8_98626400 Sample_P6_A8_16970400 Sample_P6_H10_50006400
Sample_P6_D12_31241400 Sample_P6_G8_57185400 Sample_P6_D9_59869400 Sample_P6_F10_45750400 Sample_P6_C12_97721400
Sample_P6_B11_65680400 Sample_P6_C10_18977400 ;

 if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
 run;


/* Transpose and stack junction counts */

proc transpose data=junc_ctl425_wide out=junc_ctl425_tall(rename=(COL1=apn _NAME_=sample_id));
by event_id;
run;

proc transpose data=junc_t1d425_wide out=junc_t1d425_tall(rename=(COL1=apn _NAME_=sample_id));
by event_id;
run;

proc transpose data=junc_ctl426_wide out=junc_ctl426_tall(rename=(COL1=apn _NAME_=sample_id));
by event_id;
run;

proc transpose data=junc_t1d426_wide out=junc_t1d426_tall(rename=(COL1=apn _NAME_=sample_id));
by event_id;
run;

proc transpose data=junc_ctl127_wide out=junc_ctl127_tall(rename=(COL1=apn _NAME_=sample_id));
by event_id;
run;

proc transpose data=junc_t1d127_wide out=junc_t1d127_tall(rename=(COL1=apn _NAME_=sample_id));
by event_id;
run;

proc transpose data=junc_400_wide out=junc_400_tall(rename=(COL1=apn _NAME_=sample_id));
by event_id;
run;


data junc_counts;
  set junc_ctl425_tall junc_t1d425_tall junc_ctl426_tall junc_t1d426_tall
      junc_ctl127_tall junc_t1d127_tall junc_400_tall;
run;

/* Make permenant */

data cclocal.fragment_counts_by_sample;
  set frag_tall;
run;

data cclocal.junction_counts_by_sample;
   set junc_counts;
run;

