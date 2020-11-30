/***********************************************************************************************************/
/* BIOSTAT 713                                                                                             */
/* Nonparametric methods for comparing survival curves                                                     */
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

*Number of women vs. men, distribution of NUMDZV;
proc freq data=SURV;
   tables GENDER NUMDZV ;
   run;


/*************************************** Questions to explore today****************************************/
/*
1. Does the distribution of survival time differ between women and men in this cohort? 
2. Does the distribution of survival time differ by number of diseased vessels determined at index cath?   
3. Does the distribution of survival time differ between women and men in this cohort after controlling 
   for number of diseased vessels?
*/   

/***********************************************************************************************************/
*1. Does the distribution of survival time differ between women and men in this cohort? ;
* Using K-M method, calculate survival curves, cumulative hazards, and hazard functions by gender;
* Compare survival distributions using Log-rank and Wilcoxon tests;
PROC LIFETEST DATA=SURV METHOD=KM PLOTS=(S(atrisk cl nocensor test) LS H(cl)) NOTABLE;
  STRATA GENDER /TEST=(LOGRANK WILCOXON);
  TIME DAYS2LKA*DEATH(0);
RUN;


/***********************************************************************************************************/
*2. Does the distribution of survival time differ by number of diseased vessels determined at index cath?   ;
* IN-CLASS PRACTICE:  
*   Repeat analysis, but this time compare the NUMDZV groups instead of GENDER;
* 


/***********************************************************************************************************/
*3. Does the distribution of survival time differ between women and men in this cohort after controlling 
   for number of diseased vessels?;
* Compare survival distribution in women vs. men, stratifying by NUMDZV;
* Calculate stratified Log-rank and Wilcoxon tests;

* first look at distribution of NUMDZV by GENDER to see if there might be confounding;
proc freq data=SURV;
   tables NUMDZV * GENDER / nopercent norow chisq;
   run;

* now compare survival curves by GENDER, controlling for NUMDZV;   
PROC LIFETEST DATA=SURV METHOD=KM NOTABLE;
  STRATA NUMDZV / GROUP=GENDER TEST=(LOGRANK WILCOXON);
  TIME DAYS2LKA*DEATH(0);
  format NUMDZV;  * strip character formatting to get more sensible ordering of curves;
RUN;

* to investigate further - plot survival curves for women and men, by NUMDZV;   
proc sort data=SURV;
   by NUMDZV GENDER;
   run;
PROC LIFETEST DATA=SURV METHOD=KM NOTABLE PLOTS=(S(atrisk cl nocensor)) NOTABLE;
  by NUMDZV;
  STRATA GENDER;
  TIME DAYS2LKA*DEATH(0);
RUN;