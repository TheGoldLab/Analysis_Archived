README file for code associated with:

A bias-variance trade-off governs individual differences in on-line learning in an unpredictable environment

Christopher M. Glaze, Alexandre L.S. Filipowicz, Joseph W. Kable, Vijay Balasubramanian, and Joshua I. Gold

Nature Human Behavior
2018

Model-fitting routines. These routines work but are not organized or commented as clearly as the figure-generating code. Please contact the authors with questions.

Notes (from C.Glaze):

Routines for fitting data and saved parameter estimates.

Creates *.mat files beginning with ‘allparam’ that contain fit parameters, and have names with basic syntax
‘allparam_fixedvarfits_[learning model]_[fixed parameters].mat’. (‘fixedvarfits’ refers to the fact that subjective generative variance is written as a constant linear function of the task generative variance vs. previous learning models CG analyzed in which variance is learned online within each session).

*.m files beginning with ‘analyze’ contain routines for fitting parameters with syntax ‘analyze_everyone_2AFC_[learning model].mat’

All files ending in ‘crossval’ are dedicated to estimating particular learning model parameters in one of two randomly generated partitions of experimental sessions for each subject. Session assignments to partitions were random but fixed for all analysis so that results could be most accurately compared across models. Partition information is in ‘random_session_assignments.mat’.

All fits included two basic kinds of structures:
fixedparamstrct: 
	-params0: initial parameter estimates
	-lb, ub: lower, upper boundaries on parameter estimates
	-free: binary vector indicating which parameters are free. those that are not are fixed at the 	params0 values.

paramstrct: structure with fixedparamstrct information and additional fields:
	-params: fit parameters
	-LL: model LL 
	-LLrl (sampling model fits only): model LL re-calculated after fitting was complete

NOTE: I’ve only included final models for manuscript, in which the task generative variance is assumed to be known and used correctly. This is based on (1) the fact that subjects were shown the true Gaussian distributions before each session and directly instructed on each variance, (2) allowed 120 practice trials before each experimental session and (3) previous model analysis showing that allowing subjective generative variance to freely vary resulted in lower model evidence.

MODELS:

***cmex files***
-All c-mex files contain the basic learning algorithms investigated for the research, pre-compiled for Mac with the default gcc compiler in the Xcode associated with El Captain. Windows users will have to re-compile: 

-Bayesian_mixedXX.c has the ideal observer algorithm, with each version storing different kinds of information (e.g. Bayesian_6.c only stores and returns posterior over state space for speed of calculation in fitting algorithm).

-fixedLR_learnHXX.c has a reduced learning algorithm that makes point-estimates of hazard rate trial by trial with a delta rule that operates on change point probability (2nd moment of posterior in state space).

-particle_filter_learnHXX.c: novel sampling method for learning hazard rates as reported in manuscript.

-dsprt_mixed_Hvc2.c: non-learning algorithm implementing HMM solution with a hazard rate that can vary trial-by-trial, used for fitting the initial regression model and blocked estimates of subjective hazard rate.

***geffit_XXX.m***
m-files that get the data likelihood for a given model. Each of these take in arrays of task and choice data so that entire subjects can be fit all at once. Each uses Matlab’s ‘parfor’ command to loop through arrays for speed. Those without the parallel processing toolbox can simply change this to ‘for’.

***fit_XXX.m***
m-files that do the actual fits of choice data to learning models. Every function has the same syntax:
-first 3 arguments are arrays of task and choice data, 
-2nd two indicate task generative variances and the average data mean (always 30 in this task)
-last argument is a structure fixedparamstrct with the following parallel fields:
	-params0: initial parameter estimates
	-lb, ub: lower, upper boundaries on parameter estimates
	-free: binary vector indicating which parameters are free. those that are not are fixed at the 	params0 values.
