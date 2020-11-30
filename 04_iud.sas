/****************************************************************************************
Lecture: Fitting a parametric model to one group
Data example from Collett - Table 1.1

Data from 18 women aged between 18-35 years, who had experienced two previous pregnancies. 
Measurements were made on the number of weeks from commencement of use of a particular 
intrauterine device (IUD), until discontinuation due to menstrual bleeding problems. 

weeks: weeks from starting use to discontinuation
delta: 1 = discontinuation due to mentrual bleeding problems
       0 = censoring (due to other reasons for discontinuation of lost to follow-up) 
****************************************************************************************/

options formdlim='-';

/***************************************************************************************
Input data
****************************************************************************************/
data iud;
   input weeks delta;
   datalines;
   10 1
   13 0
   18 0
   19 1
   23 0
   30 1
   36 1
   38 0
   54 0 
   56 0
   59 1
   75 1
   93 1
   97 1
   104 0
   107 1
   107 0
   107 0
;
run;


/***************************************************************************************
KM curve and cumulative hazard plot
****************************************************************************************/
ods graphics on;
PROC LIFETEST DATA=iud METHOD=KM PLOTS=(S(atrisk) logsurv loglogs H(cl));
  TIME weeks*delta(0);
  ODS OUTPUT SurvivalPlot=S_KM;
RUN;
ods graphics off;

/***************************************************************************************
Simple summary stats
****************************************************************************************/
proc means data=iud n nmiss sum mean;
   class delta;
   var weeks;
   run; 


/***************************************************************************************
Fit exponential model
****************************************************************************************/
ods graphics on;
PROC LIFEREG DATA=iud ;
  MODEL weeks*delta(0)= / DISTRIBUTION=exponential ;
RUN;
ods graphics off;


/***************************************************************************************
Fit Weibull model
****************************************************************************************/
ods graphics on;
PROC LIFEREG DATA=iud ;
  MODEL weeks*delta(0)= / DISTRIBUTION=Weibull ;
RUN;
ods graphics off;


/***************************************************************************************
Plot the estimated survival curves on top of the KM curves
****************************************************************************************/
* generate curve 'data' for fitted exponential and weibull models;
%let lambda_e = 0.0086;
%let lambda_g = 0.000454;
%let gamma_g = 1.676;
data exp;
   do weeks = 1 to 110 ;
      survival = exp(-&lambda_e * weeks);
      Estimate = 'Exponential';
      output;
   end;
   run;
data weib;
   do weeks = 1 to 110 ;
      survival = exp(-&lambda_g * ( weeks )**&gamma_g );
      Estimate = 'Weibull';
      output;
   end;
   run;

data all;
   set S_KM (keep=time survival rename=(time=weeks))
       exp (keep=weeks survival estimate)
       weib (keep=weeks survival estimate)
;
   if Estimate = "" then Estimate="KM";
   
   * calculate some transformations of survival;
   * log(-log(s));
   LLS = log(-log(survival));
   label LLS = "Log(-Log(S(t)))";
   
   * log(t);
   logt = log(weeks);
   label logt = "Log(weeks)";
run;


proc print data=all (obs=100);
run;

*Plotting the Survival Functions in PROC SGPLOT;
ods graphics on;
proc sgplot data=all ;
title "Survival Curve Estimates";
step x=weeks y=Survival /group=Estimate lineattrs=(thickness=2) legendlabel="Survival Curve Estimates";
yaxis label="Survival Probability" labelattrs=(weight=bold size=14) values=(0 to 1 by 0.1);
xaxis label="Time to IUD discontinuation (weeks)" labelattrs=(weight=bold size=14) values=(0 to 120 by 20);
keylegend / titleattrs=(weight=bold) valueattrs=(weight=bold);
run;
ods graphics off;

*Plotting the Log(-Log(S)) function versus log(weeks) in PROC SGPLOT;
ods graphics on;
proc sgplot data=all (where=(weeks>10));
title "Log(-Log(S)) Curve Estimates versus Log(weeks)";
series x=logt y=LLS /group=Estimate lineattrs=(thickness=2) legendlabel="Log(-Log(S)) Curve Estimates";
yaxis label="Log(-Log(S))" labelattrs=(weight=bold size=14) ;
xaxis label="Log weeks" labelattrs=(weight=bold size=14) ;
keylegend / titleattrs=(weight=bold) valueattrs=(weight=bold);
run;
ods graphics off;
