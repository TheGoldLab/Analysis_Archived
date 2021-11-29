function vals_ = ddOU3(fits, data)
% function vals_ = ddOU3(fits, data)
%
% 3 parameters: A, alpha, lapse
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
%       fits(3) ... lambda ("lapse")
% 
% lapse uses abbott's law
% P = G + (1 - G - L)P*
% where here P* = erf and G (guess) = 0.5

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
        -20      -20        -0.00000001; ...
      lapse        0         0];

else
    
    R0    = 10;
    M     = 1.3;    % use 1.3 instead of 1, because if I set M as a free parameter,
                    % the fitted value for both monkeys is about 1.3. This
                    % means that coherence dependence of drift rate more for high
                    % coh than low coh.
    PHI   = 0.3;
    
    acm   = fits(1).*data(:,1).^M;
    mu    = acm./fits(2).*(exp(fits(2).*data(:,2))-1);
    nu    = sqrt((2*R0 + acm)*PHI./(2*fits(2)).*(exp(2*fits(2).*data(:,2))-1));
    vals_ = 0.5 + (0.5 - fits(3)).*erf(mu./(nu.*sqrt(2)));    
end
