/****************************************************************************************
BIOSTAT 713 - Examples from Cox models Lecture 3
"Modeling covariate effects - part 2"
****************************************************************************************/

options formdlim='-' ps=120;

* Create library that contains example data;
libname dukecath "/informatics/BIOS713_Fall2020/data/DukeCathR/data" access=readonly;

* Create library for formats;
libname library "/informatics/BIOS713_Fall2020/data/DukeCathR/fmtlib" access=readonly;


/***************************************************************************************
Look at some variables from Framingham Risk model + NUMDZV, and how these relate 
to survival time in patients with non-interventional index cath in 2003-2006 presenting with no ACS;
****************************************************************************************/
data sub;
   set dukecath.dukecathr;
   if YRCATH_G=5 and ACS=0 and INTVCATH^=1 and RSEQCATHNUM=1;
   
   * collapse AGE_G to classify patients as >=65 or younger;
   agege65 = (age_g>=10) + 0*age_g;
   label agege65 = 'Age >= 65 y'
         gender = 'Female sex';
   
   * for now subset on patients with non-missing data;
   if GENDER>. and AGE_G>. and AGEGE65>. and HXDIAB>. and HXSMOKE>. and SYSBP_R>. and 
   TOTCHOL_R>. and HDL_R>. and NUMDZV>.;
   
   keep RSUBJID YRCATH_G ACS INTVCATH GENDER AGE_G AGEGE65 HXDIAB HXSMOKE SYSBP_R TOTCHOL_R HDL_R NUMDZV DEATH DAYS2LKA;
   run;

   
/***************************************************************************************
Distribution of class variables;
****************************************************************************************/
proc freq data=sub;
   tables GENDER AGE_G AGEGE65 HXDIAB HXSMOKE NUMDZV;
   run;
   
/***************************************************************************************
Distribution of continuous variables;
****************************************************************************************/
proc means data=sub n nmiss min q1 median q3 max mean stddev;
   var SYSBP_R TOTCHOL_R HDL_R days2lka;
   run;
   
/***************************************************************************************
Center and scale continuous variables
****************************************************************************************/
data sub;
   set sub;
   sysbp0 = (sysbp_r - 147)/10;
   totchol0 = (totchol_r - 179)/10;
   hdl0 = (hdl_r - 47)/10;
   label sysbp0 = 'Systolic BP, per 10 mmHg'
         totchol0 = 'Total cholesterol, per 10 mg/dL'
         hdl0 = 'HDL cholesterol, per 10 mg/dL'
         NUMDZV='NUMDZV'
         ;
   run;

/***************************************************************************************
Contrast and Estimate statements
Cox model including interaction between GENDER and NUMDZV
****************************************************************************************/     
title 'Model including GENDER, NUMDZV, and interaction';
proc phreg data=sub;
   class GENDER (param=ref ref='0') NUMDZV (param=ref ref='1');
   model days2lka * death(0) = GENDER NUMDZV GENDER*NUMDZV / ties=efron ;
            format GENDER NUMDZV;
   hazardratio NUMDZV / diff=ref;
   hazardratio GENDER / diff=ref;
   
   estimate 'Female vs. Male, NUMDZV=1' GENDER 1;
   estimate 'Female vs. Male, NUMDZV=2' GENDER 1 GENDER*NUMDZV 1 0;
   estimate 'Female vs. Male, NUMDZV=3' GENDER 1 GENDER*NUMDZV 0 1;
   
   contrast 'interaction' GENDER*NUMDZV 1 0, 
                          GENDER*NUMDZV 0 1;                          
   run;  
 
/***************************************************************************************
Example of comparing models;
****************************************************************************************/   
title 'Model including GENDER, TOTCHOL0, and interaction';
proc phreg data=sub;
   class GENDER (param=ref ref='0');
   model days2lka * death(0) = GENDER totchol0 GENDER*totchol0 / ties=efron ;
            format GENDER;
   hazardratio totchol0 / unit=1;
   run;
   
* comparison with model just including GENDER and totchol0;   
title 'Model including just GENDER, and TOTCHOL0';
proc phreg data=sub;
   class GENDER (param=ref ref='0');
   model days2lka * death(0) = GENDER totchol0  / ties=efron ;
            format GENDER;
   run;

* get p-value;
* In this code I am hardcoding the -2 log partial likelihood values obtained from above model fits;
* To apply this to a different problem you would need to change the values under datalines;
data LRT;
   input test $ parms1 fit1 parms2 fit2;
   LRT = fit1 - fit2;
   df = parms2 - parms1; 
   pval = 1 - cdf('CHISQUARE',LRT,df);
   datalines;
   interaction 2 3578.663 3 3577.853
   ;
   run;
proc print data=LRT;
run; 
 
/***************************************************************************************
Check linearity of relationship with TOTCHOL0
****************************************************************************************/   
%include "/informatics/BIOS713_Fall2020/lecture/RCS_Reg.sas" / nosource;
title 'Check linearity of relationship with TOTCHOL0';
* first get rid of some extreme values';
data sub;
   set sub;
   if TOTCHOL0 > 25 then TOTCHOL0=.;
   run;
%RCS_Reg(INFILE = sub,
         MAIN_SPLINE_VAR = TOTCHOL0,
         AVK_MSV = 0,
         KNOTS_MSV = 5 35 65 95,
         TYP_REG = cox,
         DEP_VAR = DEATH,
         SURV_TIME_VAR = DAYS2LKA,
         EXP_BETA = 0,
         PRINT_OR_HR = 0,
         NO_GRAPH = 0,
         Y_REF_LINE = 1
);

* try broken stick model with knot at 2;
data sub;
   set sub;
   
   if totchol0>. then do;
      x1 = min(totchol0,2);
      x2 = max(totchol0,2)-2;
      end;
   label x1 = 'Total cholesterol, per 10 mg/dL, below ~200 mg/dL'
         x2 = 'Total cholesterol, per 10 mg/dL, above ~200 mg/dL'
         ;
   run;
title 'Relationship with total cholesterol, accounting for non-linearity';   
proc phreg data=sub;
   model days2lka * death(0) = 
            x1 x2/ ties=efron ;
   Total_Cholesterol: test x1=0, x2=0; *overall test for cholesterol;
   run;  