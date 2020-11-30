/****************************************************************************************
BIOSTAT 713 - Examples from Cox models Lecture 4
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
Fitting Cox model including AGEGE65 as strata variable
****************************************************************************************/  
/*
proc lifetest data=sub plots=(S LLS) notable;
strata AGEGE65;
time days2lka * death(0);
run;
*/
proc phreg data=sub;
   model days2lka * death(0) = AGEGE65 GENDER HXDIAB HXSMOKE 
             / ties=efron ;
   run;
proc phreg data=sub;
   strata AGEGE65 GENDER;
   model days2lka * death(0) = HXDIAB HXSMOKE
             / ties=efron ;
   run;


/***************************************************************************************
Estimating the survival function at 1, 2, ..., 5 years for a particular subject
****************************************************************************************/   
* make dataset with covariates of patient that we want estimated survival for;
* for example, arbitarily selected record (n=100) from input data;
title 'Estimated survival probability for patient 100';
data subj100;
   set sub;
   if _n_ = 100;
   keep  RSUBJID GENDER AGEGE65 HXDIAB HXSMOKE ;
   run;
proc print data=subj100;
run;   

proc phreg data=sub;
   model days2lka * death(0) = GENDER AGEGE65 HXDIAB HXSMOKE
             / ties=efron ;
   baseline out=s100 covariates=subj100 timelist=365 730 1095 1460 1825 survival=surv lower=low upper=upp /method=breslow; 
   run;
proc print data=  s100;
   var rsubjid GENDER AGEGE65 HXDIAB HXSMOKE
            days2lka surv ;*low upp;
run; 

/***************************************************************************************
Adjusted survival
****************************************************************************************/   
* raw survival curves by diabetes;
proc lifetest data=sub notable plots(only)=survival;
   strata hxdiab;
   time days2lka * death(0);
   run;

* make dataset with covariates to represent the whole population;
data ref;
   set sub;
   keep  RSUBJID GENDER AGEGE65 HXSMOKE ;
   run;

* fit cox model with HXDIAB as strata, and adjustment for other factors;
* get directed adjusted estimate of survival by HXDIAB, using whole population as reference;
proc phreg data=sub;
   strata HXDIAB; * include HXDIAB as a strata varible;
   model days2lka * death(0) = GENDER AGEGE65 HXSMOKE
             / ties=efron ;
   baseline out=adjs covariates=ref survival=surv lower=low upper=upp / diradj method=breslow; 
   run;

* graph the adjusted curves;
title 'Survival curves for HXDIAB, adjusted for differences in GENDER, AGEGE65, HXSMOKE';  
title2 'Direct Adjusted Method';
proc sgplot data=adjs noautolegend; 
   step x=days2lka y=surv /group=hxdiab lineattrs=(thickness=2)  justify=LEFT name='adj_surv'; 

    xaxis min=0 max=4000 values=(0 to 4000 by 180 ) labelattrs=(Weight=Bold) label="Days from Cath"; 
    yaxis min=0 max=0.1 values=(0 to 1 by 0.1 ) labelattrs=(Weight=Bold) label="Siurvival Probability"; 
    keylegend 'adj_surv' /title="HXDIAB" location=outside position=bottom across=2  NOBORDER; 
run;

* fit cox model with HXDIAB as strata, and adjustment for other factors;
* get adjusted estimate of survival by HXDIAB, for the average patient;
* make dataset with covariates to represent the whole population;
proc means data=sub;
   var GENDER AGEGE65 HXSMOKE;
   output out=avg mean(GENDER AGEGE65 HXSMOKE)=GENDER AGEGE65 HXSMOKE;
   run;

proc phreg data=sub;
   strata HXDIAB; * include HXDIAB as a strata varible;
   model days2lka * death(0) = GENDER AGEGE65 HXSMOKE
             / ties=efron ;
   baseline out=adjs2 covariates=avg survival=surv lower=low upper=upp / method=breslow; 
   run;

* graph the adjusted curves;
title 'Survival curves for HXDIAB, adjusted for differences in GENDER, AGEGE65, HXSMOKE';  
title2 'Average patient method';
proc sgplot data=adjs2 noautolegend; 
   step x=days2lka y=surv /group=hxdiab lineattrs=(thickness=2)  justify=LEFT name='adj_surv'; 

    xaxis min=0 max=4000 values=(0 to 4000 by 180 ) labelattrs=(Weight=Bold) label="Days from Cath"; 
    yaxis min=0 max=0.1 values=(0 to 1 by 0.1 ) labelattrs=(Weight=Bold) label="Siurvival Probability"; 
    keylegend 'adj_surv' /title="HXDIAB" location=outside position=bottom across=2  NOBORDER; 
run;

* graph the two methods of adjustment on the same figure;
data adjc;
   set adjs (in=in1)
   adjs2 (in=in2);
   if in1 then method = 'Approach 1: Average Survival';
   else if in2 then method = 'Approach 2: Average Patient';
   run;

title 'Survival curves for HXDIAB, adjusted for differences in GENDER, AGEGE65, HXSMOKE';  
proc sgpanel data=adjc noautolegend; 
   panelby Method / columns=2 novarname;
   step x=days2lka y=surv /group=HXDIAB lineattrs=(thickness=2)  justify=LEFT name='adj_surv'; 

    colaxis min=0 max=4000 values=(0 to 4000 by 180 ) labelattrs=(Weight=Bold) label="Days from Cath"; 
    rowaxis min=0 max=0.1 values=(0 to 1 by 0.1 ) labelattrs=(Weight=Bold) label="Siurvival Probability"; 
    keylegend 'adj_surv' /title="HXDIAB" position=bottom across=2  NOBORDER; 
run;
   
/***************************************************************************************
C-index for Cox model including GENDER, AGEGE65, HXDIAB, HXSMOKE, sysbpo, totchol0 (x1, x2), hdl0;
****************************************************************************************/   
title 'Harrells c-index';
proc phreg data=sub concordance=harrell;
   model days2lka * death(0) = GENDER AGEGE65 HXDIAB HXSMOKE
            sysbp0 x1 x2 hdl0 / ties=efron ;
   run;

title 'UNOs c-index';
proc phreg data=sub concordance=uno;
   model days2lka * death(0) = GENDER AGEGE65 HXDIAB HXSMOKE
            sysbp0 x1 x2 hdl0 / ties=efron ;
   run;
