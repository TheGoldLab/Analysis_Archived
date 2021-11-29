function val_ = ddOUnb_val(fits, data, param)
% function vals_ = ddOUn_val(fits, data)
%
% Compute the DV given A and lambda (the drift rate and the leak term)
% based on an Ornstein-Uhlenback model with a bound.
% The DV is given by
%
%       mu = A*COH./LAMBDA.*(exp(LAMBDA*T-TAU)-1)   if mu<bound
%       mu = bound                                  if mu>bound
%
% , where TAU is the neural delay and bound is the 99th-percentile of the
% max. 
%
%
% INPUTS:
%   FITS - [A, LAMBDA, TAU]
%   DATA - DATA(:,1) ... coh [0 ... 1]
%          DATA(:,2) ... time (sec)
%
%

tau   = param(1);
bound = param(2);
acm   = fits(1).*data(:,1);
val_  = acm./fits(2).*(exp(fits(2).*(data(:,2)-tau))-1)+fits(3).*data(:,2);

val_(data(:,2)<tau) = 0;
val_(val_>bound)    = bound;



