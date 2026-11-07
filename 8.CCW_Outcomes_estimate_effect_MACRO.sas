libname source "C:\Users\Iftekhar Ahmed\Box\Foster-Ahmed-bup_depression\Analytic files";
libname oud "C:\Users\Iftekhar Ahmed\Box\Foster-Ahmed-bup_depression\Analytic files\dataset_output";
libname icd "C:\Users\Iftekhar Ahmed\Box\Foster-Ahmed-bup_depression\Analytic files\ICD_9_10_combined_final";



/*ESTIMATING EFFECTS*/

/*We now have two weighted data sets representing potential outcome from day 0 to 365 under two distinct hypothetical interventions:
1) No treatment initiation 
2) Treatment initiation within day 0 to 30

We can contrast the outcomes within each of these data sets using our effect measure of interest, the risk difference and risk ratio.*/

/*
*dataset= 
hosp
overdose
opi_use
bup_use
composite
mortality


*outcome= 
hospitalization
opioid_overdose
opioid_use
bup_use
composite_outcome
all_cause_mortality

*/

/*
%let dataset= overdose;
%let outcome= opioid_overdose;

*/




/*untreated dataset*/

DATA no_treatment;
	SET oud.&dataset._trunc_notrt;
	arm='untreated';
RUN;



/*We can do the same thing to the treatment initiation within day 0 to 30 group*/

DATA treatment_30;
	SET oud.&dataset._trunc_trt30;
	arm='treated';
RUN;


/*

proc freq data=treatment_30;
tables long_&outcome;

run;

proc freq data=no_treatment;
tables long_&outcome;

run;

*/

/*combining the treated and untreated satasets*/

Data &dataset._combined; 
length arm $10;
set  TREATMENT_30 no_TREATMENT;

run;




proc freq data= &dataset._combined;
tables arm*long_&outcome;
run;



/*creating estimates of survival*/

proc phreg data=&dataset._combined noprint;
    by  arm;                /* Key: PHREG runs separately within each arm and replicate */
    model end_interval*long_&outcome(0)= / entrytime=start_interval;
    weight Cum_IPCW_unstabl_trunc;
    baseline out=&dataset._surv survival=surv timelist=0 to 365 by 5;
run;



/*
proc freq data=&dataset._combined;
    tables arm*end_interval / missing;
run;
*/



/*calculating risks from survival data*/

data &dataset._risk;
    set &dataset._surv;
    by arm;
    if end_interval=365 then do;
        risk = 1 - surv;
        output;
    end;
    keep arm risk;
run;





PROC PRINT DATA=&dataset._risk;
	
RUN;



proc transpose data=&dataset._risk out=&dataset._risk_wide prefix=risk_;
    id arm;
    var risk;
run;

data &dataset._riskdiff;
    set &dataset._risk_wide;
    RiskDiff_30_vs_none = risk_treated - risk_untreated;
	Risk_ratio_30_vs_none= risk_treated/risk_untreated;
run;




*This gives the risks in ech group, risk difference, and risk ratio. For confidence interval, we need to bootstrap;

PROC PRINT DATA=&dataset._riskdiff;
	
RUN;




/* Step 3: Plot both survival curves */
PROC SGPLOT DATA=&dataset._surv;
    SERIES X=end_interval Y=surv / GROUP=arm LINEATTRS=(THICKNESS=2);
    XAXIS LABEL="Follow-up Time (Days)";
    YAXIS LABEL="Survival Probability"  VALUES=(0.6 TO 1 BY 0.1);
    TITLE "Survival Curves: Treatment vs. No Treatment";
    KEYLEGEND / LOCATION=inside POSITION=topright;
RUN;




*Risk difference treated-untreated= 0.076926
Risk ratio treated/untreated= 1.29158



*|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||;


*BOOTSRAPPING CONFIDENCE INTERVAL;



proc sort data=&dataset._combined out=subjects(keep=studyid) nodupkey;
    by studyid;
run;


*we are creating replicates at the subject level;


proc surveyselect data=subjects
    out=boot_ids
    seed=30459584
    method=urs
    samprate=1
    outhits
    rep=1000;
run;




proc sql;
    create table boot_all as
    select a.Replicate,
           b.*,
           a.NumberHits as boot_wt
    from boot_ids as a
    inner join &dataset._combined as b
        on a.studyid = b.studyid;
quit;


proc sort data=boot_all;
    by Replicate arm;
run;



/* bootstrapping at the interval level is not appropriate*/

/* STEP 2: Bootstrap the merged dataset 
proc surveyselect data=&dataset._combined out=boot_all
     seed=30459584
     method=urs
     samprate=1
     outhits
     rep=1000;
run;


*/





/*: Fit model by arm for each replicate */
proc phreg data=boot_all noprint;
    by Replicate arm;                /* Key: PHREG runs separately within each arm and replicate */
    model end_interval*long_&outcome(0)= / entrytime=start_interval;
    weight Cum_IPCW_unstabl_trunc;
    baseline out=boot_surv survival=surv timelist=0 to 365 by 5;
run;

/*  Extract final survival and compute risk per arm and replicate */
data boot_risk;
    set boot_surv;
    by Replicate arm;
    if end_interval=365 then do;
        risk = 1 - surv;
        output;
    end;
    keep Replicate arm risk;
run;

/*  Reshape to wide and compute risk difference */
proc sort data=boot_risk; by Replicate arm; run;

proc transpose data=boot_risk out=risk_wide prefix=risk_;
    by Replicate;
    id arm;
    var risk;
run;

data riskdiff;
    set risk_wide;
    RiskDiff_30_vs_none = risk_treated - risk_untreated;
	Risk_ratio_30_vs_none= risk_treated/risk_untreated;
run;




/*THIS IS WHERE YOU WILL FIND THE CI FOR TREATED, UNTREATED, RISK DIFFERENCE, RISK RATIO*/


/*  Compute bootstrap percentile CI */
proc univariate data=riskdiff noprint;
    var risk_treated risk_untreated RiskDiff_30_vs_none Risk_ratio_30_vs_none;
    output out=ci_final
        pctlpts=2.5, 97.5
        pctlpre=CI_treat_ CI_untreat_ CI_RD_ CI_RR_;
run;


/* Display results */
proc print data=ci_final;
run;











*------------------------------------------------------------------------------------------------------------------

Create cumulatuive incidence plot;

*===================================================================================================================;



/* STEP 0: Convert survival to cumulative risk */
data boot_cumrisk;
    set boot_surv;
    cumrisk = 1 - surv; /* 1 - survival = cumulative risk */
run;

/* STEP 1: Compute bootstrap percentile CI at each time point for each arm */
proc sort data=boot_cumrisk;
    by arm end_interval;
run;

proc univariate data=boot_cumrisk noprint;
    by arm end_interval;
    var cumrisk;
    output out=cumrisk_ci
        pctlpts=2.5 97.5
        pctlpre=CI_;
run;


/*
get the risks from the original analysis before bootstrapping
*/




Data risk_original; set &dataset._surv;
risk=1-surv;
drop surv;
run;


proc sort data=risk_original;
    by arm end_interval;
run;

proc sort data=cumrisk_ci;
    by arm end_interval;
run;



*combine the original risks data and the CI data from bootstrapping;

Data   cumrisk_ci_final; 
merge risk_original cumrisk_ci;
by arm end_interval;

run;





/* STEP 2: Plot cumulative risk curve with bootstrap CI */


/*
ods listing gpath="C:\Users\Iftekhar Ahmed\OneDrive - purdue.edu\PhD Research\Disseration Research Project\Aim 1";
ods graphics / reset imagename="cumrisk_curve" imagefmt=png 
               width=6in height=4in dpi=300;

*/

ods pdf file="C:\Users\Iftekhar Ahmed\OneDrive - purdue.edu\PhD Research\Disseration Research Project\Aim 1\cumrisk_curve.pdf" notoc;
ods graphics / reset width=6in height=4in dpi=300;


proc sgplot data=cumrisk_ci_final;
    /* shaded CI band */
    band x=end_interval lower=CI_2_5 upper=CI_97_5 / group=arm transparency=0.3;

    /* median cumulative risk curve */
    series x=end_interval y=risk / group=arm lineattrs=(thickness=2);

    xaxis label="Follow-up time (Days)";
    yaxis label="Absolute risk (95% CI)" values=(0 to .1 by 0.01)
    valuesformat=percent8.1  ;

	title "All-cause mortality";
    *keylegend / location=inside position=topright across=1;
run;


ods pdf close;

/*
ods graphics / reset;
ods listing close;




