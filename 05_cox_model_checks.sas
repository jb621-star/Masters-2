/****************************************************************************************
BIOSTAT 713
Cox Model Checks
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

/*
 * Define a new patient id to use for plots, counting from 1 to 1666;
 */
data sub;
   set sub;
   patient = _n_;
   run;

/*
Fit Cox PH Model including AGEGE65 HXDIAB LMST LCXST
*/
title "Fit Cox PH Model including AGEGE65 HXDIAB LMST LCXST";
title2 "Use OUTPUT statement to save residuals and predictions";
proc phreg data=sub;
model yr2dth*death(0) = AGEGE65 HXDIAB LMST LCXST;
output out=RES logsurv=rescoxs2 survival=surv resdev=resdev resmart=resmart xbeta=xbeta;
run;

/*
Output the following and store in dataset RES
logsurv: logS(t)=log(S_0(t)^exp(xB))
survival: S(t)=S_0(t)^exp(xB)
resdev:r_Di
resmart:r_Mi
xbeta: XB-->linear predictor or risk score 
*/

*Manipulate the dataset created by the 
OUTPUT statement to manually add Cox Snell residuals r_Ci;
data RES;
set RES;
*Create Cox-Snell Residuals: -logS(t);
cxsnres = -rescoxs2;
run;


/************************************************************
Diagnostic 1: Assess overall model adequacy using Cox Snell residuals.
The plot H(r_Ci) vs r_Ci can be manually obtained using 
the KM estimator where "time" is now the Cox-Snell Residuals.
************************************************************/
title "Check distribution of Cox-Snell residuals from model fit";
title2 "If cumulative hazard plot is linear with slope 1 then residuals are OK (exponential(1))";
proc lifetest data=RES method=KM plots=(logsurv lls) notable;
time cxsnres*death(0);
run;
/*Conclusion: In our example, there's no reason to think that the overall fit of our model is poor.*/


/************************************************************ 
Diagnostic 2: Assess model adequacy for each individual.
The plot r_Mi vs subject index can be manually obtained using 
PROC SGPLOT.
************************************************************/
title "Check Martingale residuals by subject (e.g., to identify outlier patients)";
proc sgplot data=RES;
scatter Y=resmart X=patient;
loess Y=resmart X=patient;
run;
/*Conclusion: The model appears to be a good fit for all subjects. 
However, if I wanted to better identify those few subjects with slightly longer survival than expected
from the model, I could print them out to examine */
proc print data=RES;
where resmart < -2;
run;


/************************************************************ 
Diagnostic 3: Assess model adequacy for each individual.
Alternatively, we can plot r_Di vs subject index using PROC SGPLOT to accomplish this goal.
************************************************************/
title "Check Deviance residuals by subject (e.g., to identify outlier patients)";
proc sgplot data=RES;
scatter Y=resdev X=patient;
loess Y=resdev X=patient;
run;
/*Conclusion: Our conclusion here is the same as that from Diagnostic 2.*/


/************************************************************ 
Diagnostic 4: Assess model adequacy across the range of predictions.
Alternatively, we can plot r_Di vs risk score using PROC SGPLOT to accomplish this goal.
************************************************************/
title "Check Deviance residuals versus X*beta to check model fit across range of predictions";
proc sgplot data=RES;
scatter Y=resdev X=xbeta;
loess Y=resdev X=xbeta;
run;
/*Conclusion: On average the model does not seem to be over- or under-predicting 
across the range of xbeta */




/************************************************************ 
Diagnostic 5: Assess proportional hazards assumption for each variable.
Plot scaled Schoenfeld resdiual for covariates by time; 
*Option ZPH provides plots of of a measure that is based on scaled Schoenfeld residuals;
*Use IDENTITY transform to test for slope with respect to time
***B_j(t) = r*_Sij + \hat{B}_j is a time-varying coefficient
***Source: Grambsch, P. M. and Therneau, T. M. (1994), “Proportional Hazards Tests and Diagnostics Based on Weighted Residuals,” Biometrika, 81, 515–526.
************************************************************/
title "Check proportional hazards assumption for each variable using scaled Schoenfeld residuals";
ods graphics on;
proc phreg data=sub ZPH(fit=LOESS transform=IDENTITY); 
model yr2dth*death(0) = AGEGE65 HXDIAB LMST LCXST;
output out=RES2 wtressch=_ALL_;*request weighted schoenfeld residuals for all variables;
run;
ods graphics off;

/*Conclusion: See some evidence of non-PH (e.g., for AGEGE65), but I am not concerned
 *because beta(t) does not change that much over follow up time
 */


/************************************************************ 
Diagnostic 6: Assess proportional hazards assumption for each variable.
Add interaction with time to model
These are defined using programming statements in proc phreg
Use linear time scale to be consistent
Center time at 5 years
************************************************************/
title "Check proportional hazards assumption for each variable by checking interaction with study time";
proc phreg data=sub; 
model yr2dth*death(0) = AGEGE65 HXDIAB LMST LCXST
              age_t diab_t lm_t lcx_t;
* programming statements for defining time interactions;              
age_t = AGEGE65 * (yr2dth - 5);
diab_t = HXDIAB * (yr2dth - 5);
lm_t = LMST * (yr2dth - 5);
lcx_t = LCXST * (yr2dth - 5);
* end of programming statements;              
run;

/*Conclusion: Same conclusion as from Diagnostic 5
 */


/************************************************************ 
Diagnostic 7: Assess functional form of continuous covariates.
1. Fit the null Cox PH model and obtain martingale residuals.
2. Plot r_Mi vs continuous predictor 
************************************************************/
title "Check Martingale residuals from null model versus X to determine how to transform X";
proc phreg data=sub;
model yr2dth*death(0) = ;
output out=RES3 resmart=resmart ;
run;

proc sgplot data=RES3;
scatter Y=resmart X=LMST;
loess Y=resmart X=LMST;
run;
proc sgplot data=RES3;
scatter Y=resmart X=LCXST;
loess Y=resmart X=LCXST;
run;

/*Conclusion: linear association OK for LMST, some nonlinearity for LCXST. */
