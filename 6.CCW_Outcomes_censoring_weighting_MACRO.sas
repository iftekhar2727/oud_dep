libname source "C:\Users\Iftekhar Ahmed\Box\Foster-Ahmed-bup_depression\Analytic files";
libname oud "C:\Users\Iftekhar Ahmed\Box\Foster-Ahmed-bup_depression\Analytic files\dataset_output";
libname icd "C:\Users\Iftekhar Ahmed\Box\Foster-Ahmed-bup_depression\Analytic files\ICD_9_10_combined_final";



/*MACRO values.

dataset= 
overdose
opi_use
bup_use


outcome= 
opioid_overdose
opioid_use
bup_use


outcome_date=
overdose_date
opioid_date
bup_date 

*/


%let dataset= overdose;
%let outcome= opioid_overdose;
%let outcome_date= overdose_date;




Data &dataset._ccw; 
  retain studyid depression_date order_date_antidep rx_start;
  set oud.outcome_&dataset ;
  rx_start= order_date_antidep-depression_date;
  drop treatment_grp date_evidence_oud baseline_start baseline_end  ;
run;




*need to investigate the rx_start dates. when do patients get rx?;


PROC univariate data= &dataset._ccw; 
var rx_start;
histogram rx_start;
run;



Proc freq data=&dataset._ccw;
tables rx_start;
run;








Data &dataset._ccw_2;
    merge &dataset._ccw (In=A)    oud.outcome_mortality (keep= studyid date_of_death all_cause_mortality); 
	 by studyid;
     if A;
	 If missing (all_cause_mortality) then all_cause_mortality=0;

run;


Data &dataset._ccw_3;
  set &dataset._ccw_2;
  if all_cause_mortality=0 AND &outcome=0 then followup=(followup_end-followup_start)+1;
  else do;
    if all_cause_mortality=1 AND &outcome=1 AND date_of_death<&outcome_date then delete;
	else followup_end= min (&outcome_date, date_of_death);
	followup= (followup_end-followup_start)+1;
  end;
run;



PROC FREQ DATA= &dataset._ccw_3;
tables &outcome all_cause_mortality followup;
run;




proc print data=&dataset._ccw_3;
where missing (&outcome);
run;


data test; set &dataset._ccw_2;
where not missing(date_of_death) AND date_of_death<&outcome_date;
run;






*|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
*-------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------

STEP 1: CLONING AND CENSORING

*-------------------------------------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------------------------------------
|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||;



*CLONING;
*==================================================================================================================;
*###################################################################################################################;


DATA &dataset._trt0_30  &dataset._notrt; 
/*Each of these data sets represents a potential clone. We will be cloning everyone into all groups, even those they are incompatible
with at baseline, though we will be setting a "t0-censoring" flag before cloning. An alternative approach is to clone people only into
groups they are compatible with at baseline, though you will need to deal with potential confounding separately during your final comparative
analysis if you choose to do this*/
	SET &dataset._ccw_3;
	IF RX_start = 0 THEN DO;
		t0censnotrt = 1;	
	END;
	ELSE DO;
		t0censnotrt=0;
	END;
	OUTPUT &dataset._trt0_30;
	OUTPUT &dataset._notrt;
RUN;





*CENSORING 0-30 day treatment group;
*==================================================================================================================;
*###################################################################################################################;


/*Now we need to create censoring within the cohorts.*/
		
/*Let's do the trt0_30 cohort first. We need to censor everyone who doesn't start by day 30 at day 30*/

DATA oud.&dataset._trt_0_30_with_cens;
	SET &dataset._trt0_30;
	/*First, let's sort out the people with less than 31 days of followup. They should use their original values*/
	IF followup < 31 THEN DO;
		Cens_followup = followup;
		if all_cause_mortality=0 then do;
            Cens_&outcome = &outcome;
		    Cens_mortality=0;
		end;
        else if all_cause_mortality=1 AND  &outcome=0  then do;
		    Cens_&outcome = 0;
		    Cens_mortality=1;
		end;
        else if all_cause_mortality=1 AND &outcome=1 AND date_of_death>=&outcome_date then do;
		    Cens_&outcome = &outcome;
		    Cens_mortality=0;
		end;
		Cens_nostart = 0; /*We also want to make it clear this is a separate censoring mechanism if we decided to weight for it*/
	END;
	ELSE IF RX_start = . OR RX_start > 30 THEN DO; /*This will get everyone who doesn't start or starts after day 30 but survives to that day*/
		Cens_followup = 30;
		Cens_&outcome = 0;
		Cens_mortality=0;
		Cens_nostart = 1; /*This flag makes it clear they were censored for not starting*/
	END;
	ELSE IF RX_start < 31 THEN DO; /*Now we assign the people who started within the first 30 days their original values*/
		Cens_followup = followup;
		if all_cause_mortality=0 then do;
            Cens_&outcome = &outcome;
		    Cens_mortality=0;
		end;
        else if all_cause_mortality=1 AND  &outcome=0  then do;
		    Cens_&outcome = 0;
		    Cens_mortality=1;
		end;
        else if all_cause_mortality=1 AND  &outcome=1 AND date_of_death>=&outcome_date then do;
		    Cens_&outcome = &outcome;
		    Cens_mortality=0;
		end;
		Cens_nostart = 0;
	END;
RUN;





*CENSORING no treatment group;
*==================================================================================================================;
*###################################################################################################################;


/*Additionally, the notrt cohort for those interested. These patients need to have their follow-up censored
on the date they start an RX.*/

DATA oud.&dataset._notrt_with_cens;
	SET &dataset._notrt;
	IF RX_start = . THEN DO; /*Everyone who is missing an RX start date uses their old outcome and follow-up information*/
		Cens_followup = followup;
		if all_cause_mortality=0 then do;
            Cens_&outcome = &outcome;
		    Cens_mortality=0;
		end;
        else if all_cause_mortality=1 AND  &outcome=0  then do;
		    Cens_&outcome = 0;
		    Cens_mortality=1;
		end;
        else if all_cause_mortality=1 AND &outcome=1 AND date_of_death>=&outcome_date then do;
		    Cens_&outcome = &outcome;
		    Cens_mortality=0;
		end;;
		Cens_start = 0; /*We want to make sure it's clear that anyone who was censored here was not censored for starting an RX*/
	END;
	ELSE IF RX_start ^= . AND followup > RX_start THEN DO; /*Everyone with an RX_start time on or before their followup date should be censored*/
		Cens_followup = RX_start;
		Cens_&outcome = 0; /*We set their outcome to 0, since they were censored*/
        Cens_mortality= 0;
        Cens_start = 1; /*We want this flag to make sure it's clear they were censored for starting*/
	END;
	ELSE IF RX_start ^= . AND followup <= RX_start THEN DO; /*Everyone with an earlier or same follow-up date than RX start should get their original followup*/
		Cens_followup = followup;
		if all_cause_mortality=0 then do;
            Cens_&outcome = &outcome;
		    Cens_mortality=0;
		end;
        else if all_cause_mortality=1 AND  &outcome=0  then do;
		    Cens_&outcome = 0;
		    Cens_mortality=1;
		end;
        else if all_cause_mortality=1 AND &outcome=1 AND date_of_death>=&outcome_date then do;
		      Cens_&outcome = &outcome;
		      Cens_mortality=0;
		end;
	    Cens_start = 0;
	END;
RUN;






*|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
*-------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------

STEP 2: WEIGHTING

*-------------------------------------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------------------------------------
|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||;




*WEIGHTING 0-30 day treatment group;
*==================================================================================================================;
*###################################################################################################################;




/*CENSORING WEIGHTS FOR THE 0-30 DAY treatment strategy*/

/*This strategy involves two forms of censoring: general loss to follow-up and censoring because patients did not start taking any medications by day 30. We will
assume that general loss to follow-up is random, meaning we only have potential selection bias from "cens_nostart=1" at day 30. Handling differential loss to
follow-up requires a separate set of weights.*/

/*First, let's make the data set long. We only need two observations: one covering day 0-30 (where everyone contributes) and one covering day 30-365 (where
only those who initiated on or prior to day 30 contribute).*/


DATA &dataset._long_trt_0_30;
	SET oud.&dataset._trt_0_30_with_cens;
	DO start_interval = 0 to 30 by 30; /*This will generate up to 2 observations per individuals with our desired start timings*/

		date = depression_date + start_interval; /*Find the date for the start of the interval*/
        format date mmddyy10.;

		IF start_interval = 0 THEN DO; /*The first interval to handle*/

			IF Cens_followup > 30 THEN DO; /*If someone's censored follow-up is greater than 30 days...*/
				long_&outcome = 0; /*They didn't have the outcome during the interval*/
				long_mortality=0;
				long_cens_nostart = 0; /*They didn't get censored during the interval*/
				interval_fu = 30; /*They should get credit for all 30 days during the interval*/
				end_interval = 30; /*The interval will end after 30 days*/
			END;

			ELSE IF Cens_followup <= 30 THEN DO; /*Otherwise, if people were followed less than 30 days...*/
				long_&outcome = cens_&outcome; /*They should get the value of the outcome from the "short" data set*/
                long_mortality= cens_mortality;
				long_cens_nostart = 0; /*They definitely weren't censored due to not starting treatment*/
				interval_fu = Cens_followup; /*They should get credit for however long we followed them*/
				end_interval = Cens_followup; /*Their interval will end after their followup duration*/
			END;

			OUTPUT;
			IF Cens_followup <= 30 AND cens_nostart = 0 THEN GOTO thend; /*Making sure nobody who has less than 30 days of follow-up gets a second observation*/

		END; /*Ending the loop for start_interval = 0*/

		ELSE IF start_interval = 30 THEN DO; /*The second interval is a bit trickier*/
			
			IF cens_nostart = 1 THEN DO; /*First, we should create observations for the people censored due to not starting at day 30*/
				long_&outcome = 0; /*They shouldn't be counted as having the outcome during this interval*/
				long_mortality= 0;
				long_cens_nostart = 1; /*They were censored due to not starting treatment*/
				interval_fu = 0; /*They should get 0 follow-up*/
				end_interval = 30; /*They shouldn't be followed past day 30*/
				recentstart = 1; /*This is for convenience and to signify that these people's hypothetical treatment regimen would involve starting at day 30, rather than earlier*/
			END;

			ELSE IF cens_nostart = 0 THEN DO; /*We also need observations for those who did start by day 30*/
				
				IF RX_start >= 0 THEN recentstart = 1; /*We need to create a "recent start" variable to ensure we only upweight those who started near the end of the interval. since antidepressants takes nearly 4-6 weeks to show meaningful clinical effects, we will include everyone receiving treatment between day 0 and 30*/
				ELSE recentstart = 0;


				IF Cens_followup <= 365 THEN DO; /*What if their follow-up ended in this period?*/
					long_&outcome = cens_&outcome; /*They should receive their true outcome*/;
					long_mortality= cens_mortality;
					long_cens_nostart = 0; /*They were not censored due to not starting treatment*/
					interval_fu = Cens_followup - 30; /*This interval covers as many days as they had past 30*/
					end_interval = Cens_followup; /*Their interval ends at their follow-up time*/
				END;

			END; /*Ending the loop for cens_nostart = 0*/

			OUTPUT;

		END; /*Ending the loop for start_interval = 30*/

	END;

	thend: ; /*This is just a label to help get the people who should only have 1 observation out of the DO loop*/
RUN;




* MERGE THE TIME FIXED COVARIATES;


Data &dataset._long_trt_0_30_cov_f;
merge &dataset._long_trt_0_30 (in=A) oud.cov_time_fixed;
by studyid;
if A;
run;


* MERGE THE TIME VARYING COVARIATES;

Data &dataset._long_trt_0_30_cov;
merge &dataset._long_trt_0_30_cov_f (in=A) oud.cov_time_vary;
by studyid start_interval;
if A;
run;





/*First, let's split this long data set into 3 different data sets.*/

DATA &dataset._first_30_days &dataset._other_days &dataset._other_days_not_recent;
	SET &dataset._long_trt_0_30_cov;
	IF start_interval = 0 THEN OUTPUT &dataset._first_30_days; /*The first data set will have all of the observations covering time 0 to 30.*/
	ELSE IF start_interval = 30 AND recentstart = 1 THEN OUTPUT &dataset._other_days; /*The second has all the observations we will be using to fit our censoring weights*/
	ELSE OUTPUT &dataset._other_days_not_recent;
RUN;

/*For the first data set, the weights are just 1 for every observation, because all censoring is random.*/

DATA &dataset._first_30_weights;
	SET &dataset._first_30_days;
	Cum_IPCW_unstabilized = 1;
	
RUN;



/*For the second data set, we need to fit a logistic regression model predicting the probability of being uncensored at the start of the interval among those with recentstart=1*/



/*
PROC FREQ DATA= &dataset._other_days;
 table long_cens_nostart;
 run;

*/


*Denominator for IPCW;

PROC GENMOD DATA=&dataset._other_days;
	
	MODEL long_cens_nostart =  gender_2 age race_2  homelessness ADHD bipolar_disorder personality_disorder
	                           schizophrenia GAD OCD PTSD any_mental_illness selfharm AUD DUD TUD pain 
							   cov_opi_overdose cov_hosp cov_opioid cov_bup cov_benzos cov_antipsychotics cov_gabapentinoids
							   cov_outp_enc cov_emergency combinedscore year_index  intv_outp_enc intv_emerg_enc intv_benzos intv_antipsychotics intv_gabapentinoids 

							   /*interaction terms*/ 
                               
                               cov_hosp*intv_gabapentinoids TUD*cov_hosp cov_hosp*intv_antipsychotics cov_hosp*cov_benzos  cov_hosp*cov_bup
                               cov_benzos*year_index combinedscore*intv_emerg_enc 

                               / link=logit dist=binomial; 

	OUTPUT OUT=denominator prob=FU_uncens;
RUN;







*IPCW modelling for competing risk MORTALITY;



PROC GENMOD DATA=denominator;
	
	MODEL long_mortality =     gender_2 age race_2  homelessness ADHD bipolar_disorder personality_disorder
	                           schizophrenia GAD OCD PTSD any_mental_illness selfharm AUD DUD TUD pain 
							   cov_opi_overdose cov_hosp cov_opioid cov_bup cov_benzos cov_antipsychotics cov_gabapentinoids
							   cov_outp_enc cov_emergency combinedscore year_index  intv_outp_enc intv_emerg_enc intv_benzos intv_antipsychotics intv_gabapentinoids 

							   /*interaction terms*/ 
                               
                               
                               / link=logit dist=binomial; 

	OUTPUT OUT=denominator_mort prob=FU_uncens_mort;
RUN;









*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

*test for interactions;

/*

%let var_inter= cov_hosp ;



PROC GENMOD DATA=&dataset._other_days;
	CLASS gender race_final insurance_type phq_cat/ param=ref;
	MODEL long_cens_nostart =  gender age race_final insurance_type homelessness ADHD bipolar_disorder personality_disorder
	                           schizophrenia GAD OCD PTSD any_mental_illness selfharm AUD DUD TUD pain phq_cat
							   cov_opi_overdose cov_hosp cov_opioid cov_bup cov_benzos cov_antipsychotics cov_gabapentinoids
							   cov_outp_enc cov_emergency combinedscore year_index  intv_outp_enc intv_emerg_enc intv_benzos intv_antipsychotics intv_gabapentinoids 

    /* interaction terms with &var_inter *
    &var_inter*gender
    &var_inter*age
    &var_inter*race_final
    &var_inter*insurance_type
    &var_inter*homelessness
    &var_inter*ADHD
    &var_inter*bipolar_disorder
    &var_inter*personality_disorder
    &var_inter*schizophrenia
    &var_inter*OCD
    &var_inter*PTSD
    &var_inter*any_mental_illness
    &var_inter*selfharm
    &var_inter*AUD
    &var_inter*DUD
    &var_inter*intv_gabapentinoids
    &var_inter*GAD
    &var_inter*phq_cat
    &var_inter*cov_opi_overdose
    &var_inter*pain
    &var_inter*cov_opioid
    &var_inter*cov_bup
    &var_inter*combinedscore
    &var_inter*cov_antipsychotics
    &var_inter*cov_gabapentinoids
    &var_inter*cov_outp_enc
    &var_inter*cov_emergency
    &var_inter* cov_benzos
    &var_inter*year_index
    &var_inter*intv_outp_enc
    &var_inter*intv_emerg_enc
    &var_inter*TUD
    &var_inter*intv_antipsychotics
    &var_inter* intv_benzos

    / dist=binomial link=logit;

run;


*/



*/

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;








/* Merge */
PROC SORT DATA=denominator_mort; BY studyid start_interval; RUN;


DATA combined_weight;
    set denominator_mort;
	combined_weight= FU_uncens_mort*FU_uncens;
	label= "";
RUN;



/*Compute unstabilized IPCW*/;

DATA &dataset._other_days_weights;
    MERGE combined_weight &dataset._other_days_not_recent;
    BY studyid start_interval;

    IF recentstart=1 AND long_cens_nostart=0 THEN  Cum_IPCW_unstabilized = 1 / combined_weight;
        
    ELSE IF recentstart=0 AND long_cens_nostart=0 THEN Cum_IPCW_unstabilized = 1;
    
    ELSE Cum_IPCW_unstabilized = 0;
        
 attrib _all_ label='';

RUN;





/*And now we can save our final analytic dataset, as well as a data set that includes observations we can use to examine the performance of IPCW.*/

DATA oud.&dataset._wted_trt_0_30 oud.&dataset._zeros_0_30;
	SET &dataset._other_days_weights &dataset._first_30_weights;
	BY studyid start_interval;
	IF cum_IPCW_unstabilized ^= 0 THEN OUTPUT oud.&dataset._wted_trt_0_30;
	ELSE OUTPUT oud.&dataset._zeros_0_30;
RUN;





*88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888;

/*

proc sgplot data=oud.&dataset._wted_trt_0_30;
    vbox Cum_IPCW_stabilized;
	where start_interval = 30 AND recentstart = 1;
run;

proc sgplot data=oud.&dataset._wted_trt_0_30;
    histogram Cum_IPCW_stabilized;
    density Cum_IPCW_stabilized;
    where start_interval = 30 AND recentstart = 1;
run;


*/




proc univariate data=oud.&dataset._wted_trt_0_30 ;
    var Cum_IPCW_unstabilized;
    output out=percentiles 
        pctlpts = 0 .1 .5 1 5 10 50 90 95 99 99.5 99.9 100
        pctlpre = P_;
	*where recentstart=1;
run;






*8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888;

*TWO SIDED TRUNCATION.  USE UNSTABILIZED WEIGHTS;

/*

proc univariate data=oud.&dataset._wted_trt_0_30 noprint;
    var Cum_IPCW_unstabilized;
    output out=wq 
        pctlpts = 1 99
        pctlpre = wt_;
	where recentstart=1;
run;

* Save percentiles into macro variables ;

data _null_;
    set wq;
    call symputx('wt_low',  wt_1);
    call symputx('wt_high', wt_99);
	
run;

%put &=wt_low &=wt_high;


* Apply truncation using the new cut points. Although the cut points have been set using only recentstart data, the truncation is being done for the entire dataset; 


data oud.&dataset._trunc;
    set oud.&dataset._wted_trt_0_30;

    * Two-sided truncation ;
    Cum_IPCW_unstabl_trunc = Cum_IPCW_unstabilized;

    if Cum_IPCW_unstabl_trunc <  &wt_low  then Cum_IPCW_unstabl_trunc = &wt_low;
    if Cum_IPCW_unstabl_trunc >  &wt_high then Cum_IPCW_unstabl_trunc = &wt_high;
	
run;


proc univariate data=oud.&dataset._trunc;
    var Cum_IPCW_unstabl_trunc;   
run;


*/


*8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888;



*UPPER TRUNCATION ONLY;



proc univariate data=oud.&dataset._wted_trt_0_30 noprint;
    var Cum_IPCW_unstabilized;
    output out=wq 
        pctlpts = 99.9
        pctlpre = wt_;
	*where recentstart=1;
run;

* Save percentiles into macro variables ;

data _null_;
    set wq;
    call symputx('wt_high', wt_99_9);
	
run;

%put &=wt_high;



* Apply truncation using the new cut points. the truncation is being done for the entire dataset; 



data oud.&dataset._trunc_trt30;
    set oud.&dataset._wted_trt_0_30;

    * one-sided truncation ;
    Cum_IPCW_unstabl_trunc = Cum_IPCW_unstabilized;

    if Cum_IPCW_unstabl_trunc >  &wt_high then Cum_IPCW_unstabl_trunc = &wt_high;
	
run;



proc univariate data=oud.&dataset._trunc_trt30;
    var Cum_IPCW_unstabl_trunc;   
run;





*/








*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$;
*00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000;

*WEIGHTING No treatment group;

*==================================================================================================================;
*###################################################################################################################;




/*CENSORING WEIGHTS FOR THE "UNTREATED" ARM*/

/*This set has two forms of censoring: general loss to follow-up and censoring because of starting a prescription. While we can assume the former 
is basically random, we should adjust for censoring due to the latter by fitting inverse probability of remaining uncensored weights. We also
need to make sure we account for the "t0" censoring that occurs. 

Let's take care of the time 0 censoring first in its own step.*/

/*First, we fit a model in the notrt data set predicting the probability of being uncensored at time 0 (i.e., t0censnotrt = 0) based on baseline 
covariates. We'll use GENMOD for multivariable logistic regression. Because we want to predict the probability of being UNCENSORED, we should not use a DESCENDING statement.*/





* MERGE THE TIME FIXED AND TIME VARYING COVARIATES. only values for interval 0 should be used;


Data cov_time_vary_zero; set oud.cov_time_vary;
if start_interval=0;
run;


Data &dataset._notrt_with_cens_cov;
merge oud.&dataset._notrt_with_cens (in=A) oud.cov_time_fixed cov_time_vary_zero;
by studyid ;
if A;
run;



/*
Data t0_mortality; set oud.&dataset._notrt_with_cens;
where all_cause_mortality=1;
run;

proc freq Data= t0_mortality;
tables followup;
run;


**** There are no deaths on t0; no weighting is needed for death;

*/



*Denominator for IPCW;

PROC GENMOD DATA=&dataset._notrt_with_cens_cov;
	
	MODEL t0censnotrt = gender_2 age race_2 homelessness  ADHD bipolar_disorder personality_disorder
	                    schizophrenia GAD PTSD OCD any_mental_illness selfharm AUD DUD TUD pain 
					    cov_opi_overdose cov_hosp cov_bup cov_opioid cov_benzos cov_antipsychotics cov_gabapentinoids cov_outp_enc
					    cov_emergency combinedscore year_index intv_outp_enc intv_emerg_enc intv_benzos intv_antipsychotics intv_gabapentinoids 

                               /* interaction effects*/

                               / link=logit dist=binomial; 

	OUTPUT OUT=denominator prob=t0_uncens;
RUN;





*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;






/*Now we can generate a new data set with those who were not censored at t0, assigning them weights based on these conditional probabilities to counter
the t0 censoring. We can call those weights "t0IPCW" to represent their action at t0.*/


/*Compute unstabilized IPCW*/;

DATA &dataset._notrt_with_cens_no_0;
    SET denominator;
	WHERE t0censnotrt ^= 1;
    
        t0IPCW_unstabilized = 1 / t0_uncens;
        
RUN;







/*
DATA &dataset._notrt_with_cens_no_0;
	SET &dataset._notrt_t0_uncens;
	WHERE t0censnotrt ^= 1;
	t0IPCW=1/t0_uncens;
RUN;
*/




/*This new data set still has potentially informative censoring, however, from those with cens_start = 1. We will need to fix this using IPCW applied
to a "long" data set of observations, basing the IPCW on a pooled multivariable logistic model predicting the probability of censoring.*/

/*First, let's make the data set long. Since we're looking at 365 days, we can create 25 intervals of 15 days each*/


DATA &dataset._long_no_trt;
	SET &dataset._notrt_with_cens_no_0;
	DO start_interval = 0 to 360 by 15; /*This will generate up to 25 observations per individual, each covering 15 days*/

        date = depression_date + start_interval; /*Find the date for the start of the interval*/
		format date mmddyy10.;

		
		IF Cens_followup > (start_interval + 15) THEN DO; /*If someone's censored follow-up is greater than the start of the interval plus 15...*/
			long_&outcome = 0; /*They didn't have the outcome during the interval*/
			long_mortality= 0;
			long_cens_start = 0; /*They didn't get censored because of starting the treatment in that interval*/
			interval_fu = 15; /*They should get credit for all 15 days of person-time in the interval*/
			end_interval = start_interval + 15; /*The interval will end after 15 days*/
		END;

		ELSE IF Cens_followup <= (start_interval + 15) THEN DO; /*On the other hand, if someone's censored follow-up falls within the interval...*/
			long_&outcome = cens_&outcome; /*They should get the value of the outcome from the "short" data set*/
			long_mortality= cens_mortality;
			long_cens_start = cens_start; /*They should get the value of cens_start from the short data set. This prevents people who experienced the outcome or were lost to follow-up people from being miscategorized*/
			interval_fu = Cens_followup - start_interval; /*They should get the days of follow-up based on their actual censored followup*/
			end_interval = Cens_followup; /*Their interval should end at the actual censored followup*/
		END;

		OUTPUT;
		IF Cens_followup <= (start_interval + 15) THEN LEAVE; /*Now those who met the second "IF" criteria should stop generating observations*/
	END;


RUN;





* MERGE THE TIME FIXED COVARIATES;


Data &dataset._long_no_trt_cov_f;

merge &dataset._long_no_trt (in=A) oud.cov_time_fixed;

by studyid ;
if A;
run;






* MERGE THE TIME VARYING COVARIATES;

Data &dataset._long_no_trt_cov;
merge &dataset._long_no_trt_cov_f (in=A) oud.cov_time_vary;
by studyid start_interval;
if A;
run;




/*
PROC FREQ DATA= &dataset._long_no_trt_cov;
 table long_cens_start;
 run;


*/

/*Now, we can fit a pooled multivariable logistic model predicting the probability of long_cens_start being equal to 0 (i.e., remaining uncensored)
across the entire grace period. We probably want to include a term for start_interval with both linear and quadratic terms to capture time trends somewhat flexibly.*/

PROC GENMOD DATA=&dataset._long_no_trt_cov;
	
	MODEL long_cens_start =  gender_2 age race_2 homelessness  ADHD bipolar_disorder personality_disorder
	                           schizophrenia GAD PTSD OCD any_mental_illness selfharm AUD DUD TUD pain 
							   cov_opi_overdose cov_hosp cov_bup cov_opioid cov_benzos cov_antipsychotics cov_gabapentinoids cov_outp_enc 
							   cov_emergency combinedscore year_index  intv_outp_enc intv_emerg_enc intv_benzos intv_antipsychotics intv_gabapentinoids 
                               start_interval start_interval*start_interval

							   /* interaction effects*/
                               

                               / link=logit dist=binomial; 

	OUTPUT OUT=denominator prob=FU_uncens;
RUN;
			





*IPCW modelling for competing risk MORTALITY;



PROC GENMOD DATA=denominator;
	
	MODEL long_mortality =  gender_2 age race_2  homelessness ADHD bipolar_disorder personality_disorder
	                           schizophrenia GAD OCD PTSD any_mental_illness selfharm AUD DUD TUD pain 
							   cov_opi_overdose cov_hosp cov_opioid cov_bup cov_benzos cov_antipsychotics cov_gabapentinoids
							   cov_outp_enc cov_emergency combinedscore year_index  intv_outp_enc intv_emerg_enc intv_benzos intv_antipsychotics intv_gabapentinoids 
                               start_interval start_interval*start_interval

							   /*interaction terms*/ 
                               
                              
                               / link=logit dist=binomial; 

	OUTPUT OUT=denominator_mort prob=FU_uncens_mort;
RUN;







/* Merge */
PROC SORT DATA=denominator_mort; BY studyid start_interval; RUN;


DATA combined_weight;
    set denominator_mort;
	combined_weight= FU_uncens_mort*FU_uncens;
	label= "";
RUN;







/*We can now generate our interval-specific and cumulative IPCWs*/

DATA &dataset._long_no_trt_weights;
	SET combined_weight;
	RETAIN Cum_IPCW_unstabilized ; /*Allow this variable to cross observations*/
	BY studyid;

	Interval_IPCW_unstabilized = 1/combined_weight;
   

	IF first.studyid THEN Cum_IPCW_unstabilized = t0IPCW_unstabilized; /*For the first observation of each patient, initiate their cumulative IPCW as their t0IPCW*/
      
	ELSE Cum_IPCW_unstabilized = Interval_IPCW_unstabilized * Cum_IPCW_unstabilized; /*For every subsequent observation, their cumulative IPCW is equal to the product of their last cumulative IPCW and the interval IPCW*/  
     
    IF long_cens_start =1 THEN Cum_IPCW_unstabilized = 0;
	   
RUN;






/*And now we can save our final analytic dataset, as well as a data set that includes observations we can use to examine the performance of IPCW.*/

DATA oud.&dataset._wted_no_trt oud.&dataset._zeros_no_trt;
	SET &dataset._long_no_trt_weights;
	BY studyid start_interval;
	IF cum_IPCW_unstabilized ^= 0 THEN OUTPUT oud.&dataset._wted_no_trt;
	ELSE OUTPUT oud.&dataset._zeros_no_trt;
RUN;







proc univariate data=oud.&dataset._wted_no_trt ;
    var Cum_IPCW_unstabilized;
    output out=percentiles 
        pctlpts = 0 .1 .5 1 5 10 50 90 95 99 99.5 99.9 100
        pctlpre = P_;
	*where long_cens_start =0;
run;







*UPPER TRUNCATION ONLY;



proc univariate data=oud.&dataset._wted_no_trt noprint;
    var Cum_IPCW_unstabilized;
    output out=wq 
        pctlpts = 99.9
        pctlpre = wt_;
	*where long_cens_start =0;
run;

* Save percentiles into macro variables ;

data _null_;
    set wq;
    call symputx('wt_high', wt_99_9);
	
run;

%put &=wt_high;



* Apply truncation using the new cut points. Although the cut points have been set using only recentstart data, the truncation is being done for the entire dataset; 



data oud.&dataset._trunc_notrt;
    set oud.&dataset._wted_no_trt;

    * one-sided truncation ;
    Cum_IPCW_unstabl_trunc = Cum_IPCW_unstabilized;

    if Cum_IPCW_unstabl_trunc >  &wt_high then Cum_IPCW_unstabl_trunc = &wt_high;
	
run;



proc univariate data=oud.&dataset._trunc_notrt;
    var Cum_IPCW_unstabl_trunc;   
run;
