/***************************************************************************** 
Analysis of competing risks (regression models)
Collett, Example 12.7 - Survival of laboratory mice 
*****************************************************************************/

/***************************************************************************** 
Read data set
******************************************************************************/

filename mice 
"/informatics/BIOS713_Fall2020/data/Collett/Survival of laboratory mice.dat"
     ;
proc format;
   value cause
      0 = 'Censored'
      1 = 'Thymic lymphoma'
      2 = 'Reticulum cell sarcoma'
      3 = 'Other causes'
      ;
   value environ
      1 = 'Standard'
      2 = 'Germ-free'
      ;
run;      
      
data mice;
   infile mice firstobs=3;
   input environment causeofdeath  time ;
   format causeofdeath cause. environment environ.;
   label time = 'Survival time in days';
   
   * calculate 0/1 indicator for germ-free environment;
   if environment=2 then germfree=1;
   else if environment=1 then germfree=0;
   label germfree = 'Germ-free environment';
   run;

title 'Summarize survival times by group and cause';   
title2 'These summaries are OK because no censoring of death times in this dataset';
proc freq data=mice;
   tables causeofdeath * germfree / nopercent norow;
   run;
proc means data=mice n median q1 q3;
   class causeofdeath germfree;
   var time;
   run;  
* analysis of time to death (all-cause) by environment;   
ods graphics on;
proc lifetest data=mice plots=survival(failure) notable;
   strata germfree / test=logrank;
   time time * causeofdeath (0);
   run;
ods graphics off; 
  
/***************************************************************************** 
Estimate and plot the cumulative incidence functions for Thymic Lymphoma in
standard and germ-free environments
Calculate Gray's test to test whether the CIF differ by group
******************************************************************************/  
title 'CIF estimates';
ods graphics on;
ods select cifplot;
proc lifetest data=mice plots=CIF(test);
   strata germfree;
   * specify the FAILCODE option to get CIF and Grays test;
   time time * causeofdeath (0) / failcode=1;
   run;
ods graphics off;

/***************************************************************************** 
Compare CIF estimates to naive KM event rate estimates
******************************************************************************/  
ods graphics on;
proc lifetest data=mice plots=survival(failure);
   strata germfree;
   * for KM estimates, censor events coded as 0, 2, 3;
   time time * causeofdeath (0 2 3);  * <-- censoring other causes of death;
   run;
ods graphics off;  
   
/***************************************************************************** 
Repeat the CIF analyses for other causes of death
Repeat the KM analyses too if you want to see how they compare
******************************************************************************/ 

/***************************************************************************** 
Cause-specific hazards models
Fit one model for each cause
******************************************************************************/ 
title 'Cause-specific hazard ratio (germ-free vs. standard)';
title2 'For Thymic Lympohoma';
proc phreg data=mice;
   * to get the cause-specific estimates for Thymic Lymphoma;
   * we censor the other causes of death (2 3); 
   model time * causeofdeath(0 2 3) = germfree / rl;
   run;

* newer version of SAS has a specific failcode option for requesting cause-specific model;
* with this option we can get estimated CIF functions using the baseline statement;

* first make a dataset containing the values of the covariates;
* for which we want to predict CIF estimates;
data covs;
   input germfree;
   datalines;
   0
   1
   ;
   run;
proc phreg data=mice;
   * use failcode(cox) option to get the cause-specific estimates for Thymic Lymphoma;
   model time * causeofdeath(0) = germfree / 
      failcode(cox)=1 
      rl;
   baseline out=lymph_CIF1 cif=cif covariates=covs; * estimated CIF;
   run;
* use STEP function in PROC SGPLOT to plot the model-based CIF estimates;   
title 'Model-based CIF estimates (cause-specific hazards model)';
title2 'For Thymic Lympohoma';
proc sgplot data=lymph_CIF1;
   step x=time y=cif / group=germfree justify=left;
   xaxis min=0 max=1000;
   yaxis min=0 max=1;
   run;   

/***************************************************************************** 
Repeat the cause-specific analysis for event types 2 and 3
******************************************************************************/ 
   
/***************************************************************************** 
Subdistribution hazards models (Fine-Gray models)
******************************************************************************/ 
title 'Subdistribution hazards ratio (germ-free vs. standard)';
title2 'For Thymic Lympohoma';
proc phreg data=mice;
   * to get subdistribution hazards model use the FAILCODE option;
   * and specify FAILCODE=1 for Thymic Lymphoma;
   model time * causeofdeath(0) = germfree / 
      failcode(fg(breslow)) = 1 
      rl;
   run;
   
/***************************************************************************** 
Re-run the Fine-Gray model and output the CIF estimates from the model.
Plot these model-based estimates and notice how these estimated curves
differ from the non-parametric CIF estimates;
******************************************************************************/ 
title 'Model-based CIF estimates (subdistribution hazards model)';
title2 'For Thymic Lympohoma';
* first make a dataset containing the values of the covariates;
* for which we want to predict CIF estimates;
data covs;
   input germfree;
   datalines;
   0
   1
   ;
   run;
proc phreg data=mice;
   model time * causeofdeath(0) = germfree / 
      failcode(fg(breslow)) = 1 
      rl;
   * use BASELINE statement to output model-based CIF estimates;
   * for Thymic Lymphoma for the two environments;   
   baseline covariates=covs out=lymph_CIF CIF;
   run;  
* use STEP function in PROC SGPLOT to plot the model-based CIF estimates;   
proc sgplot data=lymph_CIF;
   step x=time y=cif / group=germfree justify=left;
   xaxis min=0 max=1000;
   yaxis min=0 max=1;
   run;

/***************************************************************************** 
Repeat the subdistribution analysis for event types 2 and 3
******************************************************************************/  


/***************************************************************************** 
Model checking diagnostics for cause-specific model
Because we only have one binary covariate, the main assumption to check
is PROPORTIONAL HAZARDS.
******************************************************************************/ 
title 'Model checking: cause-specific hazards analysis';
title2 'For Thymic Lympohoma';
title3 'Proportional hazards';
proc phreg data=mice ZPH(transform=IDENTITY);
   model time * causeofdeath(0 2 3) = germfree / rl;
   run;    

/***************************************************************************** 
For the subdistribution hazards model, the only residuals currently
available in PHREG are the score and the (unweighted) Schoenfeld residuals.
In this example, the main assumption to check is proportional hazards 
which is not easily checked with these.
How else could you check the assumption that subdistribution hazard functions
 are proportional?
******************************************************************************/ 

