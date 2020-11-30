/***************************************************************************** 
BIOS 713
Cox models - time-dependent covariates
*****************************************************************************/

/***************************************************************************** 
Collett, Example 8.1
Bone marrow transplantation in treatment of leukemia
*****************************************************************************/

/***************************************************************************** 
Read data
*****************************************************************************/
filename bone 
"/informatics/BIOS713_Fall2020/data/Collett/Bone marrow transplantation in the treatment of leukaemia.dat"
     ;

data bm;
infile bone firstobs=3;
input patient time status group page dage precovery ptime
;
label patient = 'Patient'
      time = 'Survival time in days'
      status = 'Event indicator (0=censored, 1=event)'
      group = 'Disease group (1=ALL, 2=low risk AML, 3=high risk AML)'
      page = 'Age of patient'
      dage = 'Age of donor'
      precovery = 'Platelet recovery indicator (0=no, 1=yes)'
      ptime = 'Time in days to return of platelets to normal level (if precovery=1)'
      ;
run;

proc print data=bm; run;

/***************************************************************************** 
Example: Cox model with time-varying covariate for platelet recovery
*****************************************************************************/

/***************************************************************************** 
First method: use programming statements within proc phreg
*****************************************************************************/
proc phreg data=bm;
   model time*status(0) = p1 / rl ties=efron;
   * programming statements to turn on indicator at time of platelet recovery;
   if (ptime=. or time<=ptime)then p1=0;
   else if (ptime>. and time>ptime) then p1=1;
run;


/***************************************************************************** 
QUESTIONS - SOLUTIONS;
****************************************************************************/
* What happens if you run a model including the variable ‘precovery’ without accounting for the fact that it is time-dependent? ;
* ----;
* We get a different estimate of the hazard ratio;
* This estimate will tend to be biased too low, because patients with recovery are put in the risk set for
* the recovery group from time zero, not from time of recovery;
* This leads to over-estimating the hazard of death in the non=recovery group (risk set is too small);
proc phreg data=bm;
   model time*status(0) = precovery / rl ties=efron;
run;

* try adding the ZPH option to check if effect of p1 is constant over time;
* ----;
* SAS ignores the ZPH option;
proc phreg data=bm ZPH(fit=LOESS transform=IDENTITY);
   model time*status(0) = p1 / rl ties=efron;
   * programming statements to turn on indicator at time of platelet recovery;
   if (ptime=. or time<=ptime)then p1=0;
   else if (ptime>. and time>ptime) then p1=1;
run;


/***************************************************************************** 
Second method: restructure data in counting process format
*****************************************************************************/
data bm_count; 
   set bm;
   * if platelets dont recover then output a single record (time interval);
   * and set p1=0 on this interval;
   if ptime=. then do;
   t1=0; t2=time; event=status; p1=0;
   output; end;
   * if platelets do recover then output two records (time intervals);
   * set p1=0 on the first interval and p1=1 on the second;	
   * event indicator is 0 on first interval and equals status on the second;
    else if ptime>. then do;
	t1=0; t2=ptime; event=0; p1=0; output; * first interval;
	t1=ptime; t2=time; event=status; p1=1; output; * second interval;
	end;
run; 

proc print data=bm_count; run;

proc phreg data=bm_count;
model (t1,t2)*event(0) = p1 / rl ties=efron;
run;

/***************************************************************************** 
QUESTIONS
****************************************************************************/
* Try adding the ZPH option to check if the effect of variable P1 is constant over time;
* ----;
proc phreg data=bm_count ZPH(fit=LOESS transform=IDENTITY);
model (t1,t2)*event(0) = p1 / rl ties=efron;
run;

* This time it works!;
* So we can check for non-PH for variable p1 using ZPH, but only if we use the counting process syntax;

* CAUTION:  just because an option runs in SAS, doesn't mean it always makes sense to use it;
* In the above example it does make sense to run the ZPH test for time-dependent covariate p1;
* However, it is NOT sensible to make survival predictions using the BASELINE statement with ;
* time dependent covariates.  If will run but it will give you the incorrect answer.;


/***************************************************************************** 
Example - Change time scale to age and account for left truncation
****************************************************************************/
data bma;
   set bm;
   
   * recalculate event time so that it is in years of age;
   timea = page + (time/365);
   
   * define 'start' (t1) and 'stop' (t2) variables;
   * t1=age when patient entered the study;
   * t2=age when patient had event or was censored;
   t1_a = page;
   t2_a = timea;
   
 
run;

proc print data=bma;
run;


/***************************************************************************** 
Fit Cox model using the age as time scale ;
Need (start, stop) method because follow up is left truncated;
Examine effect of donor age;
*****************************************************************************/
proc phreg data=bma;
   model (t1_a, t2_a)*status(0) = dage / rl ties=efron;
run;


/***************************************************************************** 
QUESTIONS
****************************************************************************/
* Estimate the HR for 1 year increase in donor age when the hazard is a function of days since transplant;
* ----;
proc phreg data=bma;
   model time*status(0) = dage / rl ties=efron;
run;

* This gives a slightly different estimate for the hazard ratio when;
* we are comparing hazards as a function of days since transplant;
* The results from the above two analyses may be similar or different, depending on how different the ;
* hazard functions look as a function of age, versus as a function of days since transplant;


* Try fitting the model for the effect of donor age on hazard as a function of patient age using the following syntax.;
proc phreg data=bma;
   model timea*status(0) = dage / rl ties=efron;
run;

* What do you get?;
* ----;
* This gives a slightly different estimate of the hazard ratio compared to the model accounting for left truncation;

* What is wrong with this approach?;
* ----;
* This approach is incorrect because we incorrectly assume that patients are at risk for events since birth;
* However, we are only able to observe events for patients who surive long enough to enter our study;
* Not accounting for left truncation will cause the hazards to be under-estimated (in all groups);
* This may cause bias (in either direction) in the hazard ratio estimate;
