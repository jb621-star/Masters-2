/****************************************************************************************
BIOSTAT 713 - SAS Computing Lab
Parametric Models
****************************************************************************************/

/****************************************************************************************
Using a subset of the DukeCath dataset, examine the relationship between smoking
and death or myocardial infarction within 1 year after cath
****************************************************************************************/

options formdlim='-' ps=120;

* Create library that contains example data;
libname dukecath "/informatics/BIOS713_Fall2020/data/DukeCathR/data" access=readonly;

* Create library for formats;
libname library "/informatics/BIOS713_Fall2020/data/DukeCathR/fmtlib" access=readonly;


/***************************************************************************************
Input data
Subset on patients whose index cath was 1985-1990 and who presented without ACS; 
****************************************************************************************/
data sub;
       
   set dukecath.dukecathr;
   if YRCATH_G=1 and acs=0 and RSEQCATHNUM=1;

   * simplify some of the covariates to be binary;
   white = (race_g=1) + 0*race_g;
   agege65 = (age_g>=10) + 0*age_g;
   female = (gender=1) + 0*gender;
   cad3 = (numdzv>=3) + 0*numdzv;
   label white = 'White race'
         agege65 = 'Age >= 65'
         female = 'Female sex'
         cad3 = '3-vessel disease'
;   

   * define a new endpoint for death or myocardial infarction by 1 year of cath;
   * rescale time so that first day of study is day 1 not day 0;
   dthmi1=0;
   d2dthmi1=1;
   label dthmi1 = 'Death or MI within 1 year'
         d2dthmi1 = 'Days from cath to death or MI within 1 year'
         ;
   if . < dsmi <= 366 or (death=1 and . < days2lka <= 366) then do;
      dthmi1 = 1;
      d2dthmi1 = min(dsmi, days2lka) + 1;
      end;
   else d2dthmi1 = min(days2lka, 366) + 1;
   
   keep rsubjid rseqcathnum yrcath_g acs 
        hxsmoke white agege65 female cad3 dthmi1 d2dthmi1 death days2lka dsmi;
   
   run;
   
title 'Baseline characteristics in smokers vs. non-smokers';
proc freq data=sub;
   tables (agege65 female white cad3 ) * hxsmoke / nopercent norow;
   run;
title;

title 'Check derivation of composite endpoint';
proc means data=sub n nmiss min max;
   class dthmi1 death;
   var days2lka dsmi d2dthmi1;
   run;
   

/***************************************************************************************
* PART 1: Get to know your data
* KM curve and cumulative hazard plot
****************************************************************************************/
title 'Nonparametric estimates of survival, cumulative hazard, and hazard functions';
PROC LIFETEST DATA=sub METHOD=KM PLOTS=(S(atrisk) logsurv loglogs h ) notable;
  STRATA hxsmoke; 
  TIME d2dthmi1*dthmi1(0);
  ODS OUTPUT SurvivalPlot=S_KM;
RUN;
title;   
   
*** Based on this, how might you model these data?

/***************************************************************************************
* PART 2:  Fit exponential model to estimate average hazard rates and hazar ratio
*          Show that we can estimate the same thing by fitting Poisson model
****************************************************************************************/

/***************************************************************************************
Although exponential does not look appropriate (hazard high initially then decreasing)
we might choose to fit the exponential to summarize the average incidence rate
through 1 year of follow up
****************************************************************************************/

*** Do this by fitting exponential model in proc lifereg;
*** and estimating the hazard rate and 95% confidence interval in each group;
ods output Estimates=expests;
proc lifereg data=sub;
   model d2dthmi1*dthmi1(0) = hxsmoke / dist=exponential;
   estimate 'Non-smokers' intercept 1 hxsmoke 0 / cl;
   estimate 'Smokers' intercept 1 hxsmoke 1 / cl;
   estimate 'Difference' intercept 0 hxsmoke 1 / cl;
   run;
ods output close;

*** estimate the hazard rates = exp(-alpha_i) where alpha_i are estimates
*** from the AFT model (PROC LIFEREG);
*** also back transform the confidence limits;
*** scale so that units are per 100 patient years;
*** also estimate the hazard ratio for smokers vs. non-smokers;
data hazards;
   set expests;
   if label ^= 'Difference';
   hazard = round( exp(-estimate) * 100 * 366, .1);
   hazard_low = round(exp(-upper) * 100 * 366, .1);
   hazard_upp = round(exp(-lower) * 100 * 366, .1);
   label hazard = 'Hazard rate, per 100 patient years'
         hazard_low = 'Lower limit of 95% CI'
         hazard_upp = 'Upper limit of 95% CI'
         ;
   run;
title 'Hazard estimates and 95% CI, per 100 patient years';
proc print data=hazards noobs;
   var label hazard hazard_low hazard_upp;
   run;  
title;

title 'Compare to raw %';
proc freq data=sub;
   tables hxsmoke * dthmi1 / nopercent nocol;
   run;
title;

data hrs;
   set expests;
   if label = 'Difference';
   label = 'HR';
   hrest = round( exp(-estimate), 0.01);
   hr_low = round(exp(-upper), 0.01);
   hr_upp = round(exp(-lower), 0.01);
   label hrest = 'Hazard ratio (Smokers vs. Non-smokers)'
         hr_low = 'Lower limit of 95% CI'
         hr_upp = 'Upper limit of 95% CI'
         ;
   run;
title 'Hazard ratio estimate and 95% CI';
proc print data=hrs noobs;
   var label hrest hr_low hr_upp;
   run;  
title;
    
/***************************************************************************************
An alternative approach is to fit a Poisson model, adjusting for follow up time
This works because of the relationship between Exponential and Poisson
****************************************************************************************/
data sub1;
   set sub;
   * need to calculate log(d2dthmi1) to use as an offset;
   * this takes care of different follow up times for different patients; 
   * also rescale time so that units will be in years;
   logt = log(d2dthmi1/366);
   run;
proc genmod data=sub1;
   model dthmi1 = hxsmoke / dist=poisson link=log offset=logt;
   estimate 'Non-smokers' intercept 1 hxsmoke 0 / exp;
   estimate 'Smokers' intercept 1 hxsmoke 1 / exp;
   estimate 'Difference' intercept 0 hxsmoke 1 / exp;
   run;

/***************************************************************************************
* PART 3: Fit AFT model to quantify effect of smoking on death/MI in 1 year 
*         Start by checking if AFT assumption is reasonable using QQ plot
****************************************************************************************/

/***************************************************************************************
* First let's check the quantiles in smokers and non-smokers to make sure AFT reasonable;
Pick off reasonable set of percentiles, and plot for one group versus the other;
****************************************************************************************/
proc sort data=s_km (where=(survival>.));
   by survival;
   run;

* define percentile as the earliest time when survival drops below 1 - (perc/100);
* after sorting data set by increasing survival, this will be the last time point where S(t)<1-p/100;
%macro getperc(indat, perc);

   data s0_&perc;
   set &indat (where=(stratum='No' and survival < (1-(&perc/100)) )) end=last;
   if last; * grab the first time that survival drops below threshold;
   mergev=1;
   keep time mergev;
   run;
   data s1_&perc;
   set &indat (where=(stratum='Yes' and survival < (1-(&perc/100)) )) end=last;
   if last;
   mergev=1;
   keep time mergev;
   run;
   data p_&perc;
      merge s0_&perc (rename=(time=t0))
             s1_&perc (rename=(time=t1))
       ;
      by mergev;
      percentile=&perc;
      drop mergev;
      run;

%mend getperc;

%getperc(indat=s_km, perc=01);
%getperc(indat=s_km, perc=02);
%getperc(indat=s_km, perc=03);
%getperc(indat=s_km, perc=04);
%getperc(indat=s_km, perc=05);
%getperc(indat=s_km, perc=06);
%getperc(indat=s_km, perc=07);
%getperc(indat=s_km, perc=08);
%getperc(indat=s_km, perc=09);
%getperc(indat=s_km, perc=10);

data percentiles;
   set p_:;
   run;

proc print data=percentiles;
run;

proc sgplot data=percentiles ;
title "Q-Q plot for smoking groups (1st-90th percentiles of survival time)";
scatter x=t0 y=t1 ;
yaxis label="Percentile (days) for Smokers";
xaxis label="Percentile (days) for Non-smokers";
run;
title;

** AFT model looks OK;

/***************************************************************************************
* PART 4: Fit AFT model using LOG-LOGISTIC distribution
*         Calculate and examine Cox-Snell residuals to check model fit
****************************************************************************************/
   
/***************************************************************************************
Now lets fit an accelerated failure time model
including covariates for age>=65, female sex, and cad3
The Weibull looks reasonable, but compare results fitting log-logistic
****************************************************************************************/
** Fit log-logistic AFT and check fit;
title 'Fit Log-logistic AFT - all covariates';
proc lifereg data=sub;
   model d2dthmi1*dthmi1(0) = hxsmoke agege65 female cad3 / dist=llogistic;
   output out=outll cres=cs; *output cox-snell residuals;
   run;

title2 'Check fit of log-logistic by checking if Cox-Snell residuals are ~exponential(1)'; 
proc lifetest data=outll plots(only)=(LS) notable; 
time cs*dthmi1(0); 
run; 


** Hmmm... model fit looks OK based on this plot;
  
/***************************************************************************************
* PART 5: Fit AFT model using WEIBULL distribution
*         Calculate and examine Cox-Snell residuals to check model fit
*         Save model fit statistics so we can use them later for Likelihood Ratio Test
****************************************************************************************/

** Fit Weibull AFT and check fit;
title 'Fit Weibull AFT - smoking + covariates (agege65, female, cad3)';
ods output fitstatistics=fit1 ModelInfo=model1;
proc lifereg data=sub;
   model d2dthmi1*dthmi1(0) = hxsmoke agege65 female cad3 / dist=weibull;
   output out=outw cres=cs; *output cox-snell residuals;
   run;
title;   

** save -2 log likelihood statistic and number of parameters;
data fit1;
   merge fit1 (where=(Criterion='-2 Log Likelihood'))
         model1 (where=(Label1='Number of Parameters') rename=(nvalue1=nparm1))
         ;
   mergev=1;
   keep mergev nparm1 value;
   run;
   
title2 'Check fit of Weibull by checking if Cox-Snell residuals are ~exponential(1)'; 
proc lifetest data=outw plots(only)=(LS) notable;
time cs*dthmi1(0);
run;

** Weibull seems to fit OK;

/***************************************************************************************
* PART 6: Try a simpler model and test if it fits as well as model from PART 5
****************************************************************************************/
** After adjusting for covariates, smoking is no longer associated with increased risk;
** We could stop here...  but lets try to simplify our model a bit;

** Female sex doesnt seem to be strongly associated so lets drop it from model
** and test if smaller model fits substantially worse;
title 'Fit Weibull AFT - drop female sex';
ods output FitStatistics=fit0 ModelInfo=model0;
proc lifereg data=sub;
   model d2dthmi1*dthmi1(0) = hxsmoke agege65 cad3 / dist=weibull;
   run;
title;   

** save -2 log likelihood test;
data fit0;
   merge fit0 (where=(Criterion='-2 Log Likelihood'))
         model0 (where=(Label1='Number of Parameters') rename=(nvalue1=nparm0))
         ;
   mergev=1;
   keep mergev nparm0 value;
   run;
   
title 'Calculate -2 log likelihood ratio test for dropping female from model';
data fits;
   length comparison $25;
   merge fit1 (in=in1 rename=(value=stat1)) 
         fit0 (in=in0 rename=(value=stat0))
         ;
   by mergev;
   
   comparison = 'Drop female from model';
   testLL = stat0 - stat1;  * difference in -2 log L;
   dfLL = nparm1 - nparm0;  * difference in number of parameters from the two models;
   pvalLL = 1 - cdf('CHISQUARE', testLL, dfLL);
   
   label testLL = 'Likelihood Ratio Test Statistic'
         dfLL = 'Degrees of freedom'
         pvalLL = 'P-value'
         ;
   format pvalLL pvalue6.3;
   run;
proc print data=fits;
   var comparison testLL dfll pvalLL;
run;   
   
** The smaller model without female sex and white race does not fit significantly worse;

/***************************************************************************************
* PART 7: Using model from PART 6, calculate the hazard ratio for smoking, adjusting for 
*         other covariates in the model.
****************************************************************************************/
** Lets estimate the hazard ratio for smoking adjusting for age>=65 and CAD3;
** To estimate the hazard ratio we recall that a Weibull AFT model is also a PH model;
** with betaj = - alphaj / sigma;
** where alphaj and sigma are estimated from the AFT;
** But we will have to use the formulas for variance given by Collett (5.6.3);
** in order to get a 95% CI;
title 'Fit Weibull AFT - hxsmoke agege65 cad3';
ods output covb = covm ParameterEstimates=ests;
proc lifereg data=sub;
   model d2dthmi1*dthmi1(0) = hxsmoke agege65 cad3 / dist=weibull covb;
   run;
title;   
ods output close;

* pull the needed results into one dataset to facilitate calculations;
data all;
   merge ests (where=(lowcase(Parameter)='hxsmoke') rename=(Estimate=alpha) keep=Parameter Estimate )
         ests (where=(lowcase(Parameter)='scale') rename=(Estimate=sigma) keep=Parameter Estimate )
         covm (where=(lowcase(Variable)='hxsmoke') rename=(HXSMOKE=v_alpha Scale=cov_as) keep=Variable HXSMOKE Scale)
         covm (where=(lowcase(Variable)='scale') rename=(Scale=v_sigma) keep=Variable Scale)
         ;
   
   drop Variable Parameter;
   
   ** calculate beta, the log hazard ratio for HXSMOKE;
   beta = - (alpha / sigma);
   
   ** calculate the approximate variance for beta (Collett, formula 5.49);
   v_beta = ( (sigma**2) * v_alpha + (alpha**2) * v_sigma - 2 * alpha * sigma * cov_as ) / (sigma**4);      
   
   ** calculate approximate 95% CI for beta;
   lower = beta - 1.96 * sqrt(v_beta);
   upper = beta + 1.96 * sqrt(v_beta);
   
   ** exponentiate to get back to HR scale;
   smoke_hr = exp(beta);
   smoke_low = exp(lower);
   smoke_upp = exp(upper);
   label smoke_hr = 'Hazard Ratio (Smokers vs. Non-smokers)'
         smoke_low = 'Lower CI Limit'
         smoke_upp = 'Upper CI Limit'
         ;
   
   run;

title2 'Hazard ratio for Smoking, adjusted for age>=65 and 3-vessel CAD';
proc print data=all;
   var smoke:;
run;




