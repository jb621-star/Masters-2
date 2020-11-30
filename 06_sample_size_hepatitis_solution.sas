/***************************************************************************** 
Collett, Example 15.1 - Sample size calculation for Active Chronic Hepatitis 
*****************************************************************************/

/***************************************************************************** 
Read data set to get estimate of survival in standard therapy group 
******************************************************************************/

filename hept "/informatics/BIOS713_Fall2020/data/Collett/Chronic active hepatitis.dat"
     ;

data hept0;
   infile hept firstobs=3;
   input treatment time status;
   
   * for this example we only need the standard group (treatment=2);
   if treatment=2;
   run;
   
/***************************************************************************** 
Get estimate of Survival in standard treatment group 
Note that these estimates are slightly different from those given by Collett,
Example 15.2, which were read off the graph.
******************************************************************************/
title 'Survival in standard treatment group';
proc lifetest data=hept0 method=pl outsurv=surv0 timelist=24 33 42 48 72;
   time time*status(0);
   run;   
  
proc print data=surv0;
run;   


/***************************************************************************** 
Determine required sample size using the same settings as hand calculation
in Collett, Example 15.1 - 15.2
Provide as many points from Standard survival curve as possible
Assume no loss to follow up
******************************************************************************/
title 'Sample size required to detect HR=0.57 with 0.90 power and alpha=0.05';
proc power;
   twosamplesurvival test=logrank
      curve("Standard") = 
      0:1
      2:0.95
      3:0.91
      4:0.86
      7:0.82
      10:0.77
      22:0.73
      28:0.68
      29:0.64
      32:0.59
      37:0.55
      40:0.50
      41:0.45
      54:0.41
      61:0.36
      63:0.32
      71:0.27
      refsurvival = "Standard"
      hazardratio = 0.57
      accrualtime = 18
      followuptime = 24
      grouplossexphazards = (0 0)
      power = 0.9
      alpha = 0.05
      sides = 2
      groupweights = (1 1)
      ntotal = .;
run;

/***************************************************************************** 
A sample size analysis relies heavily on the assumptions.
Try some of the following modifications to see how they affect the sample size/power.
Modifications:
1. Output the expected total number of deaths instead of the total sample size
2. How does the required sample change for a range of hazard ratios between 0.6 and 0.8
3. How does the sample size change if randomization to the standard:new group was 1:2
4. How would power change if you extended follow up to 48 months with 408 patients?
******************************************************************************/
title1 "Modifications";
title2 '#1 - Total number events required to detect HR=0.57 with 0.90 power and alpha=0.05';
proc power;
   twosamplesurvival test=logrank
      curve("Standard") = 
      0:1
      2:0.95
      3:0.91
      4:0.86
      7:0.82
      10:0.77
      22:0.73
      28:0.68
      29:0.64
      32:0.59
      37:0.55
      40:0.50
      41:0.45
      54:0.41
      61:0.36
      63:0.32
      71:0.27
      refsurvival = "Standard"
      hazardratio = 0.57
      accrualtime = 18
      followuptime = 24
      grouplossexphazards = (0 0)
      power = 0.9
      alpha = 0.05
      sides = 2
      groupweights = (1 1)
      eventstotal = .;
run;


title2 '#2 - Sample size required to detect HR=0.6-0.8 with 0.90 power and alpha=0.05';
proc power;
   twosamplesurvival test=logrank
      curve("Standard") = 
      0:1
      2:0.95
      3:0.91
      4:0.86
      7:0.82
      10:0.77
      22:0.73
      28:0.68
      29:0.64
      32:0.59
      37:0.55
      40:0.50
      41:0.45
      54:0.41
      61:0.36
      63:0.32
      71:0.27
      refsurvival = "Standard"
      hazardratio = 0.6 0.7 0.8
      accrualtime = 18
      followuptime = 24 
      grouplossexphazards = (0 0)
      power = 0.9
      alpha = 0.05
      sides = 2
      groupweights = (1 1)
      ntotal = .;
run;

title2 '#3 - Sample size when randomization is 1:2';
proc power;
   twosamplesurvival test=logrank
      curve("Standard") = 
      0:1
      2:0.95
      3:0.91
      4:0.86
      7:0.82
      10:0.77
      22:0.73
      28:0.68
      29:0.64
      32:0.59
      37:0.55
      40:0.50
      41:0.45
      54:0.41
      61:0.36
      63:0.32
      71:0.27
      refsurvival = "Standard"
      hazardratio =0.57
      accrualtime = 18
      followuptime = 24 
      grouplossexphazards = (0 0)
      power = 0.9
      alpha = 0.05
      sides = 2
      groupweights = (1 2)
      ntotal = .;
run;

title2 '#4 - Power to detect HR=0.57 with N=408 patients, 0.90 power, and alpha=0.05 and 48 months follow-up';
proc power;
   twosamplesurvival test=logrank
      curve("Standard") = 
      0:1
      2:0.95
      3:0.91
      4:0.86
      7:0.82
      10:0.77
      22:0.73
      28:0.68
      29:0.64
      32:0.59
      37:0.55
      40:0.50
      41:0.45
      54:0.41
      61:0.36
      63:0.32
      71:0.27
      refsurvival = "Standard"
      hazardratio = 0.57
      accrualtime = 18
      followuptime = 48
      grouplossexphazards = (0 0)
      power = .
      alpha = 0.05
      sides = 2
      groupweights = (1 1)
      ntotal = 408;
run;