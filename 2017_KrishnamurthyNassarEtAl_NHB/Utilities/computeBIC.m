function [BIC, AIC]=computeBIC(logLikelihood, numParams, n)

%% computes BIC given penalizes extra parameters (numParams) a

BIC=-2.*logLikelihood+numParams.*log(n);
AIC=-2.*logLikelihood+2.*numParams;

