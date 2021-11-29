function vals_ = ddPow3FixLapse(fits, data, lapse)
% function vals_ = ddPow3FixLapse(fits, data, lapse)
%
% 2 parameters: A, alpha
%
% Computes TIME-based DD function with LAPSE and BIAS.
%   Assumes DD to a fixed TIME (given as col 2 of data).
%   Thus, pct correct is fraction of gaussian above
%       0 at TIME.
%   Max Drift rate depends on coherence as a power law:
%       A_max(coh) = A*coh^1.25
%   And on time as a power function
%       T.^N
%
%   at values in "data":
%       data(1)   ... coh [0 ... 1]
%       data(2)   ... time (sec)
%
%   given parameters in "fits":
%       fits(1) ... A      (coh scale)
%       fits(2) ... N      (time exponent)
%       fits(3) ... lambda ("lapse")
%
% lapse uses abbott's law
% P = 0.5 + (0.5 - L)P*
% here  P* = erf

% return initial values (init, min, max)
if nargin < 2 || isempty(fits)

    vals_ = [20  0.0001 1000];

else

    PHI   = 0.3;
    R0    = 10;
    M     = 1;

    acm   = fits(1).*data(:,1).^M;
    tn    = data(:,2);
    mu    = acm.*tn;
    nu    = sqrt(PHI.*tn.*(2*R0 + acm));
    
    if lapse<0.05
        vals_ = 0.5 + (0.5 - lapse).* ...
                erf(mu./(nu.*sqrt(2)));
    else
        vals_ = 0.5 + (0.5 - 0.05).* ...
                erf(mu./(nu.*sqrt(2)));
    end    
    
end
