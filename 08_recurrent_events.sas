/***************************************************************************** 
Example - analyzing recurrent events using several different models
*****************************************************************************/

/***************************************************************************** 
From Therneau & Grambsch, 8.5.4: Recurrence of bladder cancer

Data originally listed in Wei, Lin, and Wasserfield (1989), Regression
analysis of multivariate incomplete failure time data by modeling the 
marginal distributions. JASA, 84: 1065-1073.
*****************************************************************************/

/***************************************************************************** 
Read data
*****************************************************************************/
filename blad 
"/informatics/BIOS713_Fall2020/data/Therneau/bladder.dat";

*** Below is the SAS code provided by Dr. Terry Therneau;
*
* Read in the bladder data set
*   drop the one useless subject who has no follow-up time
* ;
* In the dataset there are a max of 4 recurrences per patient
* ;

data temp;
    infile blad missover;
    retain id 0;
    input rx  futime number size  r1-r4;
    if (futime =0) then delete; * patient with no follow-up;
    id = id +1;  * generate an ID number for patients;
	rx = rx-1; * convert rx to 0 (placebo) / 1 (thiotepa);
	label 
       id = 'Subject id (1-85)'
       futime = 'Followup time (months)'
       number = 'Initial number of tumors'
       size = 'Initial size of tumors'
       rx = 'Treat code, 0=placebo, 1=thiotepa'
       r1 = 'Months to 1st recurrence'
       r2 = 'Months to 2nd recurrence'
       r3 = 'Months to 3rd recurrence'
       r4 = 'Months to 4th recurrence'
       ;

run;

proc print data=temp; run;

/***************************************************************************** 
Poisson model
This assumes constant event rate over time, and variance equal to mean;
*****************************************************************************/
title 'Poisson regression';

* Prepare the data set;
*    Calculate the total number of events per subject;
*    Get the follow-up time for the last event;

data bladder0;
   set temp;
   by id;
   lfup = log(futime);
   nevents = n(r1, r2, r3, r4); * count non-missing event times;
run;

proc print data=bladder0; run;

title2 'Estimate the relative event rate (# events/month) for Thiotepa vs. Placebo';
* Fit the Poisson model with log link, using log(futime) as offset;
proc genmod data=bladder0;
model nevents = rx size number/ dist=poisson link=log offset=lfup ;
estimate 'effect of thiotepa' rx 1 / exp; 
run;

* Poisson model assumes that mean and variance are equal;
* This may not be true;
* Negative binomial provides a more flexible variance;

/***************************************************************************** 
Negative binomial model
This assumes constant event rate over time;
*****************************************************************************/
title 'Negative binomial regression';

title2 'Estimate the relative event rate (# events/month) for Thiotepa vs. Placebo';
* Fit the Negative binomial model with log link, using log(futime) as offset;

proc genmod data=bladder0;
model nevents = rx size number/ dist=negbin link=log offset=lfup ;
estimate 'effect of thiotepa' rx 1 / exp;
run;


/***************************************************************************** 
Andersen-Gill Analysis
*****************************************************************************/
title 'Andersen-Gill (AG) Analysis';

* Step 1: make Andersen-Gill style data;
data bladder1;
    set temp;

    drop futime r1-r4;

    time1 =0;
    enum  =0;

    * create interval for recurrence 1;
    if (r1 ne .) then do;
	time2 = r1;
	status= 1;
	enum  = 1;
	output;
	time1 = r1;
	end;

    * create interval for recurrence 2;
    if (r2 ne .) then do;
	time2 = r2;
	status= 1;
	enum  = 2;
	output;
	time1 = r2;
	end;

    * create interval for recurrence 3;
    if (r3 ne .) then do;
	time2 = r3;
	status= 1;
	enum  = 3;
	output;
	time1 = r3;
	end;

    * create interval for recurrence 4;
    if (r4 ne .) then do;
	time2 = r4;
	status= 1;
	enum  = 4;
	output;
	time1 = r4;
	end;

    * create a final interval until end of follow-up;
    * if subject has no recurrences then this will be only interval;
    if (futime > time1) then do;
	time2 = futime;
	status= 0;
	enum  = enum +1;
	output;
	end;
	
	label time1 = 'Start of interval'
	      time2 = 'End of interval'
	      enum = 'Event number'
	      status = 'Event indicator (1=event, 0=censored)'
	      ;
run;

proc print data=bladder1;
   var id rx number size enum time1 time2 status;
run;

proc freq data=bladder1;
   tables enum * status / nopercent norow nocol;
   run;
   
* Step 2: Fit Andersen-Gill model using PROC PHREG;

* Note the use of the COVS(AGGREGATE) option and ID statement;
* This is necessary to account for correlation within a subject;
title2 'Andersen-Gill model with common treatment effect for each recurrence';
proc phreg data=bladder1 covs(aggregate);
model (time1,time2)*status(0) = rx size number /rl;
id id;
run;

/*
title2 'Compare with Cox regression for first event only';
proc phreg data=bladder1;
where enum=1;
model (time1,time2)*status(0) = rx size number /rl;
run;

title2 'Incorrect AG analysis: without the robust sandwich variance estimate';
proc phreg data=bladder1;
model (time1,time2)*status(0) = rx size number /rl;
id id;
run;
*/

/***************************************************************************** 
Use PROC PHREG to generate nonparametric cumulative mean function estimates by RX
Uses the AG-style data as input;
*****************************************************************************/  
* relabel time2 to improve axis labeling;
data bladder1g;
   set bladder1;
   label time2 = "Months since study entry";
   run;
ods graphics on;
proc phreg data=bladder1g covs(aggregate) 
      plots(overlay=byrow)=mcf;
   id id;
   strata rx;
   model (time1,time2)*status(0)=;
   run; 
ods graphics off;   
   

/***************************************************************************** 
Prentice, Williams, and Peterson (PWP) Analysis
*****************************************************************************/
title 'PWP Analysis';

* Step 1: make PWP style data;
* We can use the same data structure created above for the AG analysis;

* Step 2: Fit the PWP model using PROC PHREG;
* This is the same as the AG model except ENUM (event number) 
* is used as a STRATA variable;
* This is to give the 1st, 2nd, ...  events their own baseline hazard function;

* Note the use of the COVS(AGGREGATE) option and ID statement;
* This is necessary to account for correlation within a subject;

title2 'PWP model with common treatment effect for each recurrence';
proc phreg data=bladder1 covs(aggregate);
strata enum;
model (time1,time2)*status(0) = rx size number / rl;
id id;
run;
* In the above model we are assuming the coefficients for covariates;
* are the same for all sequential events;


title2 'PWP model with different treatment effect for each recurrence';

* Add to the dataset a treatment variable specific to each event recurrence;
data bladder1;
   set bladder1;
   * define a treatment variable separately for each event;
   rx1 = rx*(enum=1);
   rx2 = rx*(enum=2);
   rx3 = rx*(enum=3);
   rx4 = rx*(enum=4);
   label rx1 = 'Treatment, for 1st recurrence'
         rx2 = 'Treatment, for 2nd recurrence'
         rx3 = 'Treatment, for 3rd recurrence'
         rx4 = 'Treatment, for 4th recurrence'
         ;
   run;
         
* >>> NOW WRITE THE PHREG CODE TO FIT THIS MORE FLEXIBLE MODEL;


/***************************************************************************** 
Wei, Lin, and Weissfeld (WLW) analysis;
*****************************************************************************/
title 'WLW Analysis';

* Step 1: make WLW style data;

* There are a maximum of 4 recurrences per subject;
* Make a dataset where EACH subject has 4 records (time intervals);
* The event time for each records is the time of the jth event (or end of followup);

data bladder2;
    set bladder1 (keep=id rx number size time1 time2 enum status);
    by id;

    futime = time2;
    if (1 <= enum <= 4) then output;
    
    * if subject doesnt have 4 intervals already then add intervals;
    * each of the new intervals will end at the max event time or end of follow up;
    if (last.id =1) then do;
	temp = enum +1;
	do enum = temp to 4;
	    status =0;
	    output;
	    end;
	end;
	drop temp time1 time2 ;
	label futime = 'Time (months) of jth event or end of follow up'
        ;
run;


proc print data=bladder2; run;

* Step 2: Fit the WLW model using PROC PHREG;

* This model does not use the (start, stop) syntax because the interval
* for each event begins at 0;
* This is the key difference between this analysis and the PWP analysis;

* The STRATA statement gives the 1st, 2nd, ...  events;
* their own baseline hazard function;

* Note the use of the COVS(AGGREGATE) option and ID statement;
* This is necessary to account for correlation within a subject;

title2 'WLW model with common treatment effect for each recurrence';
proc phreg data=bladder2 covs(aggregate);
strata enum;
model futime*status(0) = rx size number / rl;
id id;
run;

title2 'WLW model with separate treatment effect for each recurrence';

data bladder2; 
set bladder2; 
   * define treatment indicators for each event;
   rx1 = rx*(enum=1);
   rx2 = rx*(enum=2);
   rx3 = rx*(enum=3);
   rx4 = rx*(enum=4);
   label rx1 = 'Treatment, for 1st recurrence'
         rx2 = 'Treatment, for 2nd recurrence'
         rx3 = 'Treatment, for 3rd recurrence'
         rx4 = 'Treatment, for 4th recurrence'
         ;
run;

* >>> NOW WRITE THE PHREG CODE TO FIT THIS MORE FLEXIBLE MODEL;