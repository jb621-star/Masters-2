/****************************************************************************************
Using a subset of the DukeCath dataset, examine the relationship between diabetes and survival.
****************************************************************************************/

options formdlim='-' ps=120;

libname dukecath "/informatics/BIOS713_Fall2020/data/DukeCathR/data" access=readonly;
libname library "/informatics/BIOS713_Fall2020/data/DukeCathR/fmtlib" access=readonly;

/***************************************************************************************
Input data
Subset on patients whose index cath was 2003-2006 and who presented without ACS; 
****************************************************************************************/
data sub;
   set dukecath.dukecathr;
   if YRCATH_G=5 and acs=0 and RSEQCATHNUM=1;
   run;


title 'Baseline characteristics in diabetics vs. non-diabetics';
proc freq data=sub;
   tables (age_g gender race_g ) * hxdiab / nopercent norow;
   run;
title;

/***************************************************************************************
KM curve and cumulative hazard plot
****************************************************************************************/
title 'Nonparametric estimates of survival and cumulative hazards function';
PROC LIFETEST DATA=sub METHOD=KM PLOTS=(S(atrisk) logsurv loglogs ) notable;
  STRATA hxdiab; 
  TIME days2lka*death(0);
  ODS OUTPUT SurvivalPlot=S_KM
             LogNegLogSurvivalPlot=LLS_KM;
RUN;
title;

* make prettier version of LLS plot (to deal with logtime<0);
title 'Plot of nonparametric estimate of log(H(t)) vs. log(t)';
PROC SGPLOT DATA=LLS_KM;
   SERIES X=LOG_TIME_ Y=LOG___LOG_SURVIVAL___ / GROUP=Stratum MARKERS;
   WHERE log_time_>=0;
   XAXIS LABEL="Log (time)";
   YAXIS LABEL="Log(-log(Survival))";
   run;
title;   

/***************************************************************************************
Pick off reasonable set of percentiles, and plot for one group versus the other;
****************************************************************************************/
proc sort data=s_km (where=(survival>.));
   by survival;
   run;

* define percentile as the earliest time when survival drops below 100-perc;
%macro getperc(indat, perc);

   data s0_&perc;
   set &indat (where=(stratum='No' and survival < ((100 - &perc)/100))) end=last;
   if last;
   mergev=1;
   keep time mergev;
   run;
   data s1_&perc;
   set &indat (where=(stratum='Yes' and survival < ((100 - &perc)/100))) end=last;
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
%getperc(indat=s_km, perc=05);
%getperc(indat=s_km, perc=10);
%getperc(indat=s_km, perc=20);
%getperc(indat=s_km, perc=25);
%getperc(indat=s_km, perc=30);
%getperc(indat=s_km, perc=40);

data percentiles;
   set p_:;
   run;

proc print data=percentiles;
run;

title "Q-Q plot for diabetes groups (1st-40th percentiles of survival time)";
proc sgplot data=percentiles ;
scatter x=t0 y=t1 ;
yaxis label="Percentile (days) for Diabetes group";
xaxis label="Percentile (days) for No Diabetes group";
run;
title;

/***************************************************************************************
Fit Weibull model to explore associations between various factors and survival
****************************************************************************************/
data sub;
   set sub;
   white = (race_g=1) + 0*race_g;
   agege80 = (age_g=13) + age_g*0;
   female = (gender=1) + 0*gender;
   label white = 'White race'
         agege80 = 'Age >= 80'
         female = 'Female sex'
;
   run;

* check derivations;
proc freq data=sub;
   tables white * race_g / list missing nocum nopercent;
   tables agege80 * age_g / list missing nocum nopercent;
   tables female * gender / list missing nocum nopercent;
   run;

title 'Fit Weibull AFT';
proc lifereg data=sub;
   model days2lka*death(0) = hxdiab female white agege80;
   run;
