function vals_ = ddOU2(fits, data)
% function vals_ = ddOU2(fits, data)
%
% 2 parameters: A, alpha
%
% Computes Ornstein-Uhlenback DD function with constant drift rate.
%   Assumes DD to a fixed TIME (given as col 2 of data).
%   Thus, pct correct is fraction of gaussian above
%       0 at TIME.
%   Max Drift rate depends on coherence as a power law:
%       A_max(coh) = A*coh^M, where M is fixed, below
%   And on time as a decaying exponential
%       A(coh, t) = A_max * exp(-alpha*t)
%
%   at values in "data":
%       data(1)   ... coh [0 ... 1]
%       data(2)   ... time (sec)
%       data(3)   ... dot dir (-1/1)
%
%   given parameters in "fits":
%       fits(1) ... A      (coh scale)
%       fits(2) ... alpha  (time exponent)
% 

% return initial values (init, min, max)
if nargin < 1 || isempty(fits)
    
    % data matrix includes pcor (last column)
    % guess lapse from high coherence trials    
    if nargin < 2 || isempty(data)
        lapse = 0;
    else
        Lmax = data(:,1) == max(data(:,1));
        lapse = min([0.49, 1.0 - sum(data(Lmax,end))./sum(Lmax)]);
    end

%     vals_ = [ ...
%         30        30  30; ...
%          0.001 -100       100; ...
%      lapse        0.0001    0.18];

    vals_ = [ ...
          2        0.0001  150; ...
        -20      -20        -0.00000001];

else
    
    R0    = 10;
    M     = 1;
    PHI   = 0.3;
    
    acm   = fits(1).*data(:,1).^M;
    mu    = acm./fits(2).*(exp(fits(2).*data(:,2))-1);
    nu    = sqrt((2*R0 + acm)*PHI./(2*fits(2)).*(exp(2*fits(2).*data(:,2))-1));
    vals_ = 0.5 + 0.5.*erf(mu./(nu.*sqrt(2)));    
    
    
%     R0    = 10;
%     M     = 1;
%     PHI   = 0.3;
%     
%     acm   = fits(1).*data(:,1).^M;
%     mu    = acm./fits(2).*(exp(fits(2).*data(:,2))-1);
%     nu    = sqrt((2*R0 + acm)*PHI./(2*fits(2)).*(exp(2*fits(2).*data(:,2))-1));
%     vals_ = 0.5 + 0.5.*erf(mu./(nu.*sqrt(2)));    

end
