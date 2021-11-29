function mu = ddOUn(fits, data)
% function vals_ = ddOUn(fits, data)
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
%   FITS - [A, LAMBDA, TAU]
%   DATA - DATA(:,1) ... coh [0 ... 1]
%          DATA(:,2) ... time (sec)
%
%
   
acm = fits(1).*data(:,1);
mu  = acm./fits(2).*(exp(fits(2).*(data(:,2)-fits(3)))-1);   
