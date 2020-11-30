/***********************************************************************************************************/
/* BIOSTAT 713                                                                                             */
/* Nonparametric methods - estimating S(t), h(t), H(t)                                                     */
/* Computing Lab                                                                                           */
/***********************************************************************************************************/

/*******************************************Prepare Data****************************************************/
*Create SAS library that contains example data;
libname DUKECATH "/informatics/BIOS713_Fall2020/data/DukeCathR/data" access=readonly;
libname library "/informatics/BIOS713_Fall2020/data/DukeCathR/fmtlib" access=readonly;

*Create survival dataset in WORK library;
data SURV;
  set DUKECATH.dukecathr;
  where RSEQCATHNUM = 1 AND YRCATH_G = 7; 
  keep RSUBJID YRCATH_G RSEQCATHNUM GENDER NUMDZV DAYS2LKA DEATH;
run;

proc print data=SURV(obs=5) ;
   var RSUBJID YRCATH_G RSEQCATHNUM GENDER NUMDZV DAYS2LKA DEATH;
run;


/************************Questions to explore today************************/
/*
1. What is the probability that a patient in this cohort survives to 2 years or beyond?
2. What is the cumulative hazard of death by 2 years in this cohort?
3. What is the instantaneous risk of death at 2 years for a patient in this cohort?
*/

*************************** Kaplan-Meier Method ****************************;
*First answer Q1, Q2, Q3 using Kaplan-Meier Method;
PROC LIFETEST DATA=SURV METHOD=KM PLOTS=(S(atrisk cl nocensor) LS H(cl));
  TIME DAYS2LKA*DEATH(0);
  ODS OUTPUT SurvivalPlot=S_KM
             NegLogSurvivalPlot=CH_KM
             HazardPlot=H_KM; 
RUN;
*Manipulate the S_KM dataset to view the data at each
t_j in the same way it was presented in class;
data S_KM(keep=time atrisk event pj qj survival sdf_ucl sdf_lcl stratum);
keep time atrisk event pj qj survival sdf_ucl sdf_lcl stratum;
set S_KM(where=(Survival ^=.));
pj = event/atrisk;
qj = 1 - pj;
run;

proc print data=S_KM(obs=10);
run; 

*1. What is the probability that a patient in this cohort survives to 2 years or beyond?;
data Q1_KM;
set S_KM(where=(Survival ^= . and .< Time <=730)) end=last;
if last;
run;
proc print data=Q1_KM;
run;

*2. What is the cumulative hazard of death by 2 years in this cohort?;
data Q2_KM;
set CH_KM(where=(Survival ^= . and .< Time <=730)) end=last;
if last;
run;
proc print data=Q2_KM;
run;

*3. What is the instantaneous risk of death at 2 years for a patient in this cohort?;
*We'll just select the closet time prior to the time that we are interested in just due to the 
behavior of the kernel smoother.;
data Q3_KM;
set H_KM(where=(Hazard ^= . and .< Time <=730)) end=last;
if last;
run;
proc print data=Q3_KM;
run;

*************************** Life-Table Method ****************************;
*Answer Q1, Q2, Q3 using Life-Table Method;
PROC LIFETEST DATA=SURV METHOD=LT INTERVALS=(0 to 1500 by 180) PLOTS=(S(atrisk cl) LS H(cl));
  TIME DAYS2LKA*DEATH(0);
  ODS OUTPUT SurvivalPlot=S_LT
             NegLogSurvivalPlot=CH_LT
             HazardPlot=H_LT; 
RUN;
*1. What is the probability that a patient in this cohort survives to 2 years or beyond?;
data Q1_LT(where=(ind=1));
set S_LT;
ind = (Time <= 730 < Time + 180);
run;
proc print data=Q1_LT;
run;

*2. What is the cumulative hazard of death by 2 years in this cohort?;
data Q2_LT(where=(ind=1));
set CH_LT;
ind = (Time <= 730 < Time + 180);
run;
proc print data=Q2_LT;
run;

*3. What is the instantaneous risk of death at 2 years for a patient in this cohort?;
data Q3_LT(where=(ind=1));
set H_LT;
ind = (Midpoint - 90 <= 730 < Midpoint + 90);
run;
proc print data=Q3_LT;
run;

*************************** Nelson-Aalen Method ****************************;
*Now answer Q1, Q2 using Nelson-Aalen Method;
PROC LIFETEST DATA=SURV METHOD=BRESLOW NELSON PLOTS=(S(atrisk cl nocensor) LS);
  TIME DAYS2LKA*DEATH(0);
  ODS OUTPUT SurvivalPlot=S_NA
             NegLogSurvivalPlot=CH_NA; 
RUN;
*1. What is the probability that a patient in this cohort survives to 2 years or beyond?;
data Q1_NA;
set S_NA(where=(Survival ^= . and .< Time <=730)) end=last;
if last;
run;
proc print data=Q1_NA;
run;

*2. What is the cumulative hazard of death by 2 years in this cohort?;
data Q2_NA;
set CH_NA(where=(Survival ^= . and .< Time <=730)) end=last;
if last;
run;
proc print data=Q2_NA;
run;


/***********************************************************************************************************/
/* In-class practice                                                                                       */
/* What is the probability of surviving to 6 months or beyond (6 months = 182.5 days)                      */
/* Report this estimate using the Kaplan-Meier, Life-Table, and Nelson-Aalen methods                       */
/***********************************************************************************************************/

