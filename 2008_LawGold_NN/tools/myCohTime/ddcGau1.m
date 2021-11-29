function vals_ = ddcGau1(fits, data)
% function vals_ = ddcGau1(fits, data)
%
% 1 parameters: A
%
% Computes performance predicted by DDM with perfect accumulation.
%   Assumes DD to a fixed TIME (given as col 2 of data).
%   Thus, pct correct is fraction of gaussian above
%       0 at TIME.
%
%   at values in "data":
%       data(1)   ... coh [0 ... 1]
%       data(2)   ... time (sec)
%
%   given parameters in "fits":
%       fits(1) ... A      (coh scale)
% 

% return initial values (init, min, max)
if nargin < 1 || isempty(fits)

    vals_ = [10 0.0001  1000];
else
    
    PHI   = 1.5;
    acm   = fits(1).*data(:,1);
    mu    = acm.*data(:,2);
    nu    = sqrt(PHI*acm.*data(:,2));
    vals_ = 1-normcdf(0,mu,nu);    
    

end
