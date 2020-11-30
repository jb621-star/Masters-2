/****************************************************************************************
Data example from Collett - Table 3.6

Survival times for 36 patients with a malignant tumor in the kidney (hypernephroma).
Patients were all treated with a combination of chemotherapy and immunotherapy.
Surgical removal of the kidney (nephrectomy) was performed on some patients.
Age at time of diagnosis was also recorded.

months: months from diagnosis to death or end of follow up
death: 1 = death
       0 = end of follow-up
agegp: age (grouped) at time of diagnosis
surg:  1 = nephrectomy (kidney removed)
       0 = no nephrectomy
****************************************************************************************/

options formdlim='-' ps=120;

/***************************************************************************************
Input data
****************************************************************************************/
proc format;
   value agefmt 
     1 = "<60"
     2 = "60-70"
     3 = ">70"
;
   value yesno
     0 = "No"
     1 = "Yes"
;    
run;
       
data neph;
   input months death agegp surg;
   format agegp agefmt. death yesno. surg yesno.;
   label months = "Months from diagnosis to death or last follow-up"
         death = "Death"
         agegp = "Age at diagnosis (years)"
         surg = "Nephrectomy performed"
;
   datalines;
9   1  1  0
6   1  1  0
21  1  1  0
15  1  2  0
8   1  2  0
17  1  2  0
12  1  3  0
104 0  1  1  
9   1  1  1  
56  1  1  1  
35  1  1  1  
52  1  1  1  
68  1  1  1  
77  0  1  1  
84  1  1  1  
8   1  1  1  
38  1  1  1  
72  1  1  1  
36  1  1  1  
48  1  1  1  
26  1  1  1  
108 1  1  1  
5   1  1  1  
108 0  2  1  
26  1  2  1  
14  1  2  1  
115 1  2  1  
52  1  2  1  
5   0  2  1  
18  1  2  1  
36  1  2  1  
9   1  2  1  
10  1  3  1  
9   1  3  1   
18  1  3  1  
6   1  3  1  
;
run;

/*
proc print data=neph;
   by surg agegp;
   run;
*/
title 'Number of patients by age group and surgery';
proc freq data=neph;
   tables surg * agegp / nopercent norow nocol;
run;
title;

/***************************************************************************************
KM curve and cumulative hazard plot
****************************************************************************************/
ods graphics on;
PROC LIFETEST DATA=neph METHOD=KM PLOTS=(S(atrisk) logsurv loglogs ) notable;
  STRATA surg;
  TIME months*death(0);
  ODS OUTPUT SurvivalPlot=S_KM;
RUN;
ods graphics off;


/***************************************************************************************
Fit exponential model
****************************************************************************************/
proc means data=neph n sum;
   class surg;
   var death months;
   run;

PROC LIFEREG DATA=neph ;
  MODEL months*death(0)= surg / DISTRIBUTION=exponential ;
RUN;


/***************************************************************************************
Fit Weibull model
****************************************************************************************/
title 'Weibull model with covariate for surgery';
PROC LIFEREG DATA=neph ;
  MODEL months*death(0)= surg / DISTRIBUTION=Weibull ;
RUN;


/***************************************************************************************
Fit Weibull model with additional covariates
****************************************************************************************/
title 'Weibull model with covariates for surgery and age';
PROC LIFEREG DATA=neph ;
  class agegp;
  MODEL months*death(0)= surg agegp / DISTRIBUTION=Weibull ;
RUN;

title 'Weibull model with covariates for surgery, age, and their interaction';
PROC LIFEREG DATA=neph ;
  class agegp;
  MODEL months*death(0)= surg|agegp / DISTRIBUTION=Weibull ;
RUN;


** manually compute the likelihood ratio test results;
title 'Calculate likelihood ratio test statistics for comparing nested models';
data test;
   length comp $25.;
   input comp $ base new pb pn;
   lrchi = base - new;
   diffp = pn - pb;
   pval = 1 - cdf('CHISQUARE', lrchi, diffp);
   label comp = 'Comparison'
         base = '-2 Log L for base model'
         new = '-2 Log L for new model'
         pb = 'Number of parameters for base model'
         pn = 'Number of parameters for new model'
         lrchi = 'Chi-squared test statistic'
         diffp = 'Degrees of freedom'
         pval = 'P-value';
   format pval pvalue6.3;
   datalines;
   Mod1_vs_Mod2 94.384 87.758 2 4
   Mod2_vs_Mod3 87.758 83.064 4 6
;
run;

proc print data=test;
run;
