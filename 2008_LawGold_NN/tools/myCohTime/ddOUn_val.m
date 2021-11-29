function val_ = ddOUn_val(fits, data, param)
% function vals_ = ddOUn_val(fits, data)
%
% Compute the DV given A and lambda (the drift rate and the leak term)
% based on an Ornstein-Uhlenback model. The DV is given by
%
%       mu = A*COH./LAMBDA.*(exp(LAMBDA*T-TAU)-1),
%
% , where TAU is the neural delay.
%
%
% INPUTS:
%   FITS - [A, LAMBDA, PARAM]
%   DATA - DATA(:,1) ... coh [0 ... 1]
%          DATA(:,2) ... time (sec)
%
%
   
tau  = param(1);
acm  = fits(1).*data(:,1);
val_ = acm./fits(2).*(exp(fits(2).*(data(:,2)-tau))-1);
val_(data(:,2)<tau) = 0;


