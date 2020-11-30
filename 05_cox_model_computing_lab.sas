/****************************************************************************************
BIOSTAT 713 - Computing Lab
Cox Model
****************************************************************************************/

options formdlim='-' ps=120;

* Create library that contains example data;
libname dukecath "/informatics/BIOS713_Fall2020/data/DukeCathR/data" access=readonly;

* Create library for formats;
libname library "/informatics/BIOS713_Fall2020/data/DukeCathR/fmtlib" access=readonly;


/***************************************************************************************
Get index cath between 1999-2002 (YRCATH_G=4) for patients presenting without ACS;
****************************************************************************************/
data sub;
   set dukecath.dukecathr;
   if ACS=0 and YRCATH_G=4 and RSEQCATHNUM=1;
   
   * delete records where stenosis measures are missing;
   if RCAST=. or LMST=. or LADST=. or LCXST=. or PRXLADST=. then delete;  
  
   * collapse AGE_G to classify patients as >=65 or younger;
   agege65 = (age_g>=10) + 0*age_g;
   label agege65 = 'Age >= 65 y'
         gender = 'Female sex';
         
   * rescale survival time to years (and add one day);
   yr2dth = (DAYS2LKA + 1) / 365.25;
   label yr2dth = "Years from Index Cath to Death or End of Follow Up";
           
   * patient revascularized within 30 days (treat this as intended treatment strategy);
   revasc30=0;
   if . < dspci <=30 or . < dscabg<= 30 then revasc30=1;
   label revasc30 = "Initial treatment strategy = revasc";
   
   keep RSUBJID GENDER AGE_G AGEGE65 HXDIAB HXSMOKE 
        RCAST LMST LADST LCXST PRXLADST 
        NUMDZV DEATH DAYS2LKA yr2dth revasc30;
   run;

* number of records and number with missing values;   
proc means data=sub n nmiss;
   run;

title 'KM curves and logrank test by REVASC30';
proc lifetest data=sub plots=(survival (nocensor atrisk) lls) method=KM notable;
   strata REVASC30 / test=logrank;
   time yr2dth * death(0);
   run;

/***************************************************************************************
Basic syntax of Cox model for binary variable (REVASC30);
Two types of syntax for specifying event time:
   1) Time_variable * event_variable (censor value) = covariates
   2) (t1, t2) * event_variable (censor value) = covariates ("counting process syntax")
****************************************************************************************/title 'Standard model syntax';
title2 'Breslow method for tied event times (=SAS default)';
proc phreg data=sub;
   model yr2dth * death(0) = REVASC30 / ties=breslow rl=wald type3(wald);
   run;
title2 'Discrete method for tied event times';
proc phreg data=sub;
   model yr2dth * death(0) = REVASC30 / ties=discrete rl=wald type3(wald);
   run;
title2 'Efron method for tied event times';
proc phreg data=sub;
   model yr2dth * death(0) = REVASC30 / ties=efron rl=wald type3(wald);
   run;

title2 'Discrete method for tied event times, type3 score test';
proc phreg data=sub;
   model yr2dth * death(0) = REVASC30 / ties=discrete rl=wald type3(score);
   run;

title2 'Discrete method for tied event times, type3 likelihood ratio test';
proc phreg data=sub;
   model yr2dth * death(0) = REVASC30 / ties=discrete rl=wald type3(lr);
   run;

  
title 'Counting process syntax';
* Counting process syntax;
* For this example, subject is at risk from time of cath (time zero). Define this as tZero;
data sub;
   set sub;
   tZero = 0;
   run;
proc phreg data=sub;
   model (tZero,yr2dth) * death(0) = REVASC30 / ties=efron rl=wald type3(wald);
   run;
  
  
/***************************************************************************************
Set REVASC30 as class variable and try different parameterizations
****************************************************************************************/   
title 'Specify REVASC30 as class variable';
title2 'Default settings';
proc phreg data=sub;
   class REVASC30;
   model yr2dth * death(0) = REVASC30 / ties=efron rl=wald type3(wald);
   run;
title2 'Explicit settings for class variable, e.g., specify 0 as reference level';
proc phreg data=sub;
   class REVASC30 (param=ref order=internal ref='0');
   model yr2dth * death(0) = REVASC30 / ties=efron rl=wald type3(wald);
   hazardratio REVASC30 / diff=ref;
   run;
title2 'GLM coding';
proc phreg data=sub;
   class REVASC30 / param=GLM order=internal ref=FIRST;
   model yr2dth * death(0) = REVASC30 / ties=efron rl=wald type3(wald);
   hazardratio REVASC30 / diff=all;
   run;  
title2 'Effect coding';
proc phreg data=sub;
   class REVASC30 (param=effect order=internal ref='0');
   model yr2dth * death(0) = REVASC30 / ties=efron rl=wald type3(wald);
   hazardratio REVASC30 / diff=ref;
   run;     
 
 
/***************************************************************************************
Continuous variable: LCXST
****************************************************************************************/   
title "Univariable analysis of LCXST";
* first check distribution in case anything weird;
* output percentiles of distribution;
title2 "Check distribution of LCXST";
proc univariate data=sub;
   var LCXST ;
   output out=lcxperc pctlpre=P_ pctlpts= 5 to 95 by 5;
   run;
title2 "Percentiles of distribution of LCXST";   
proc print data= lcxperc;
run;
  
title2 "Fit model that is linear in LCXST ";   
proc phreg data=sub;
   model yr2dth * death(0) = LCXST / ties=efron rl=wald type3(wald);
   hazardratio LCXST / units=10;
   run;

* check for nonlinearity;
title2 "Check for nonlinearity in LCXST using RCS_Reg macro";
%include "/informatics/BIOS713_Fall2020/lecture/RCS_Reg.sas" / nosource;
* run macro with 4 knots placed at 10th 40th 60th and 80th percentiles;
* cannot place knots at recommended locations because LCXST is discrete;
%RCS_Reg(INFILE = sub,
         MAIN_SPLINE_VAR = LCXST ,
         AVK_MSV = 0,
         KNOTS_MSV = 10 40 60 80,
         TYP_REG = cox,
         DEP_VAR = DEATH,
         SURV_TIME_VAR = yr2dth,
         EXP_BETA = 0,
         PRINT_OR_HR = 0,
         NO_GRAPH = 0,
         Y_REF_LINE = 1
);   

* get estimates of HR from model including RCS for LCXST;
title2 "Get estimates of HR for LCXST when modeled with RCS";
%RCS_Reg(INFILE = sub,
         MAIN_SPLINE_VAR = LCXST ,
         AVK_MSV = 0,
         KNOTS_MSV = 10 40 60 80,
         TYP_REG = cox,
         DEP_VAR = DEATH,
         SURV_TIME_VAR = yr2dth,
         EXP_BETA = 1,
         PRINT_OR_HR = 0,
         NO_GRAPH = 0,
         Y_REF_LINE = 1,
         REF_VAL=75,
         SPECIF_VAL=10 20 30 40 50 60 70 75 80 90 
);   
   


/***************************************************************************************
Fit similar model using EFFECT statement in PROC PHREG;
****************************************************************************************/   
title "Use EFFECT statement to fit spline transformation of LCXST";   
proc phreg data=sub;
   effect lcx_spl = spline(LCXST / naturalcubic basis=TPF(noint) knotmethod=list(5 50 75 95) );
   model yr2dth * death(0) = lcx_spl / ties=efron rl=wald type3(wald);

   * get tests for overall association and non-linearity;
   contrast 'Overall association with LCX' lcx_spl 1 0 0, lcx_spl 0 1 0, lcx_spl 0 0 1 / E;
   contrast 'Test for non-linearity' lcx_spl 0 1 0, lcx_spl 0 0 1 / E;
   * the warnings in the SAS log about this contrast statement can be ignored *;
   
   * try to get estimates of hazard ratios compared to reference of 75;
   hazardratio '75 vs 10' LCXST / units=65 at (LCXST=10);
   hazardratio '75 vs 70' LCXST / units=5 at (LCXST=70);
   hazardratio '80 vs 75' LCXST / units=5 at (LCXST=75);
   hazardratio '90 vs 75' LCXST / units=15 at (LCXST=75);   
   
   run;  
   
/***************************************************************************************
Include interaction with NUMDZV;
****************************************************************************************/   
title "Interaction between LCXST and NUMDZV";
title2 "First check distribution of LCXST within NUMDZV";
proc means data=sub n q1 median q3 mean min max ;
   class NUMDZV;
   var LCXST;
   run;
proc means data=sub n q1 median q3 mean min max ;
   var LCXST;
run;   
title2 "Fit Cox model with interaction between NUMDZV and LCXST - reference cell coding";   
proc phreg data=sub;
   strata AGEGE65 REVASC30;
   class NUMDZV (param=ref ref='1' order=internal);
   format NUMDZV;
   model yr2dth * death(0) = NUMDZV|LCXST / ties=efron rl=wald type3(wald);
   hazardratio LCXST / units=10;
   contrast 'Effect of LCXST at NUMDZV=1' LCXST 1 / E;
   contrast 'Effect of LCXST, averaged across NUMDZV' LCXST 1  LCXST*NUMDZV 0.33333333333 0.33333333333 / E;
   contrast 'Effect of NUMDZV at LCXST=0' NUMDZV 1 0, NUMDZV 0 1 / E;
   run;  
title2 "Fit Cox model with interaction between NUMDZV and LCXST - GLM coding";   
proc phreg data=sub;
   strata AGEGE65 REVASC30;
   class NUMDZV /param=GLM order=internal ref=FIRST;
   format NUMDZV;
   model yr2dth * death(0) = NUMDZV|LCXST / ties=efron rl=wald type3(wald);
   hazardratio LCXST / units=10;
   contrast 'Main effect of LCXST' LCXST 1 / E;
   contrast 'Effect of NUMDZV at LCXST=0' NUMDZV 1 -1 0, NUMDZV 0 1 -1 / E;   
   contrast 'Effect of NUMDZV at LCXST=60.8' NUMDZV 1 -1 0 LCXST*NUMDZV 60.8 -60.8 0, NUMDZV 0 1 -1 LCXST*NUMDZV 0 60.8 -60.8 / E;      
   run;     
   
/***************************************************************************************
Multivariable analysis for several stenosis variables;
****************************************************************************************/   
title "Multivariable analyis for several stenosis variables, with stratification for AGEGE65 and REVASC30";   
title2 "First check correlations between variables";
proc corr data=sub;
   var LMST PRXLADST LADST LCXST RCAST;
   run;
title2 "Multivariable Cox regression - include main effects for all variables";   
proc phreg data=sub;
   strata AGEGE65 REVASC30;
   model yr2dth * death(0) = LMST PRXLADST LADST LCXST RCAST / ties=efron rl=wald type3(wald);
   hazardratio LMST / units=10;
   hazardratio PRXLADST / units=10;
   hazardratio LADST / units=10;
   hazardratio LCXST / units=10;
   hazardratio RCAST / units=10;
   run;  
   
title2 "Multivariable Cox regression - forward selection of main effects and interactions";   
proc phreg data=sub;
   strata AGEGE65 REVASC30;
   model yr2dth * death(0) = LMST|PRXLADST|LADST|LCXST|RCAST / ties=efron 
      selection=forward slentry=0.10;
   run;  

title2 "Multivariable Cox regression - backward selection of main effects and interactions";   
proc phreg data=sub;
   strata AGEGE65 REVASC30;
   model yr2dth * death(0) = LMST|PRXLADST|LADST|LCXST|RCAST / ties=efron 
      selection=backward slstay=0.10;
   run;  
   
title2 "Multivariable Cox regression - stepwise selection of main effects and interactions";   
proc phreg data=sub;
   strata AGEGE65 REVASC30;
   model yr2dth * death(0) = LMST|PRXLADST|LADST|LCXST|RCAST / ties=efron 
      selection=stepwise slentry=0.10 slstay=0.10;
   run;  
   
 
 
/***************************************************************************************
Compare model fit statistics and c-index for the three models;
Unfortunately, c-index not calculated when STRATA is used;
Switch stratification factors to class variables in model;
****************************************************************************************/   
title "Multivariable model - from forward selection";
proc phreg data=sub concordance=harrell ;
   class AGEGE65 REVASC30 /param=ref order=INTERNAL ref=FIRST;
   model yr2dth * death(0) = AGEGE65 REVASC30 LCXST LMST / ties=efron ;
   run;  
   
title "Multivariable model - from backward selection";
proc phreg data=sub concordance=harrell ;
   class AGEGE65 REVASC30 /param=ref order=INTERNAL ref=FIRST;
   model yr2dth * death(0) = AGEGE65 REVASC30
          LMST PRXLADST LMST*PRXLADST 
          LADST LMST*LADST PRXLADST*LADST LMST*PRXLADST*LADST 
          LCXST LMST*LCXST PRXLADST*LCXST LMST*PRXLADST*LCXST LADST*LCXST LMST*LADST*LCXST PRXLADST*LADST*LCXST 
          LMST*PRXLADST*LADST*LCXST 
          RCAST LADST*RCAST LCXST*RCAST LADST*LCXST*RCAST
   / ties=efron ;
   run;  

title "Multivariable model - from stepwise selection";
proc phreg data=sub concordance=harrell ;
   class AGEGE65 REVASC30 /param=ref order=INTERNAL ref=FIRST;
   model yr2dth * death(0) = AGEGE65 REVASC30 LMST LCXST / ties=efron ;
   run;  
   
/***************************************************************************************
Fit multivariable model with predictors for GENDER, AGEGE65, HXSMOKE, HXDIAB, REVASC30, ;
Get predicted survival for patients of interest;
****************************************************************************************/     
title "Multivariable model - GENDER, AGEGE65, HXSMOKE, HXDIAB, REVASC30";
proc phreg data=sub concordance=harrell ;
   class GENDER AGEGE65 HXSMOKE HXDIAB REVASC30 /param=ref order=INTERNAL ref=FIRST;
   model yr2dth * death(0) = GENDER AGEGE65 HXSMOKE HXDIAB REVASC30 / ties=efron ;
   run;  
   

/***************************************************************************************
Get predicted survival for patients of interest, stratifying on REVASC30;
****************************************************************************************/     
title "Predicted Survival for patients of interest";
* make covariates data that contains the unique distinct covariates of interest;
proc sort data=sub out=covs (keep=GENDER AGEGE65 HXSMOKE HXDIAB ) nodupkey;
   by GENDER AGEGE65 HXSMOKE HXDIAB;
   run;

proc phreg data=sub;
   strata REVASC30;
   class GENDER AGEGE65 HXSMOKE HXDIAB  /param=ref order=INTERNAL ref=FIRST;
   model yr2dth * death(0) =  GENDER AGEGE65 HXSMOKE HXDIAB  / ties=efron ;
   baseline out=outp covariates=covs survival=surv lower=low upper=upp / method=PL;
   run;
   
* make description for type of patient;
data outp;
   length description $60;
   set outp;
   description = "GENDER="||put(GENDER,1.)||", AGEGE65="||put(AGEGE65,1.)||", HXSMOKE="||put(HXSMOKE,1.)||", HXDIAB="||put(HXDIAB, 1.);
   run;
   
* graph survival estimates for patients of interest;
title 'Predicted survival curves for patients of interest';  
proc sgpanel data=outp noautolegend; 
   panelby description / novarname columns=2;
   step x=yr2dth y=surv / group=REVASC30 lineattrs=(thickness=2)  justify=LEFT name='pred_surv'  ; 
   colaxis min=0 max=15 values=(0 to 15 by 1 ) labelattrs=(Weight=Bold) label="Years from Cath"; 
   rowaxis min=0 max=0.1 values=(0 to 1 by 0.1 ) labelattrs=(Weight=Bold) label="Survival Probability";    
   keylegend 'pred_surv' /title="REVASC30" position=bottom across=2  NOBORDER; 
run;   


/***************************************************************************************
Get direct adjusted survival for REVASC30;
Adjust for GENDER, AGEGE65, HXDIAB, HXSMOKE;
Calculates this by averaging up survival predictions made for suitable reference pop.
****************************************************************************************/     
title "Adjusted survival for REVASC30";
* create reference population based on adjustment covariates in most recent year of cath;
data covsa;
   set sub;
   keep GENDER AGEGE65 HXDIAB HXSMOKE;
   run;

proc phreg data=sub;
   strata REVASC30;
   class GENDER AGEGE65 HXDIAB HXSMOKE /param=ref order=INTERNAL ref=FIRST;
   model yr2dth * death(0) =  GENDER AGEGE65 HXDIAB HXSMOKE / ties=efron ;
   baseline out=outa covariates=covsa survival=adjusted_survival / diradj method=BRESLOW;
   run;
   
  
* graph adjusted survival estimates for patients of interest;
title 'Survival curves for REVASC30, adjusted for GENDER, AGEGE65, HXDIAB, HXSMOKE';  
proc sgplot data=outa noautolegend; 
   step x=yr2dth y=adjusted_survival / group=REVASC30 lineattrs=(thickness=2)  justify=LEFT name='adj_surv'; 
   xaxis min=0 max=15 values=(0 to 15 by 1 ) labelattrs=(Weight=Bold) label="Years from Cath"; 
   yaxis min=0 max=0.1 values=(0 to 1 by 0.1 ) labelattrs=(Weight=Bold) label="Survival Probability";  
   keylegend 'adj_surv' /title="REVASC30" position=bottom across=2  NOBORDER; 
run;   
   
   
* compare this to unadjusted curves;
title 'Unadjusted survival curves for REVASC30';  
proc lifetest data=sub plots(only)=survival (atrisk nocensor ) notable;
   strata REVASC30;
   time yR2DTH * DEATH (0);
   run;
   
