function [wk_, wkr_, lag_, r4_] = getBIAS_wkFromFIRA(tin, lag, r4)
% function [wk_, wkr_, lag_, r4_] = getBIAS_wkFromFIRA(tin, lag, r4)
%
% tin:      1=right, -1=left
%
% returns:
%   wk_     ... wk-filtered choices wrt tin
%   wkr_    ... correlation coefficient (& CIs) between wk_ & choices

if nargin < 1 || isempty(tin)
    tin = 0;
end

if nargin < 2
    lag = [];
end

% recompute wk-filtered choices
% cohs are 0..1
% view time in seconds
% dir=-1/0/1, chc=-1/1 signed by TOUT(-1)/TIN(1)
% 0% coh -> dir = 0
Lgood = getFIRA_ecodesByName('task')==6 & getFIRA_ecodesByName('correct')>=0;
ddat  = getFIRA_ecodesByName({'dot_coh', 'dot_dur', 'dot_dir', 'choice', 'correct'}, [], Lgood);
ddat(:,1)   = ddat(:,1)./100;
ddat(:,2)   = ddat(:,2)./1000;
ddat(:,3:4) = sign(cos(ddat(:,3:4).*pi./180)).*sign(cos(tin*pi/180));
L0          = ddat(:,1) == 0;
ddat(L0,3)  = 0;

% fit to unbiased model
if nargin < 3 || isempty(r4)
    [f,s,t,p,r] = ctPsych_fit(@ddExp3, ddat(:,1:4), ddat(:,5), [], []);
    r4 = r(:,4);
    % ctPsych_plot(@ddExp3, f, ddat(:,1:4), ddat(:,5), 3, []);
end

% get wk data
Lrew          = getFIRA_ecodesByName('num_rewards', [], Lgood) > 0;
vins          = ddat(:,[4 4]);
vins(~Lrew,1) = 0;
vins(Lrew,2)  = 0;

[a,b,c] = getBIAS_wkF(vins(1:end-1,:), r4(2:end), lag);

% get wk
wk_        = nans(size(Lgood));
wk_(Lgood) = [0; c]; % filtered by cor/err "tosum"

% correlation coef between local bias and resids
wkr_       = b(1:3)';

% return lag
if nargout > 2
    lag_ = b(4);
end

% return resids
if nargout > 3
    r4_ = [find(Lgood) r4];
end

return

if 0
    r1 = R(2,1);

    Lg       = getFIRA_ecodesByName('task')==6 & getFIRA_ecodesByName('correct')>0;
    chc      = sign(cos(getFIRA_ecodesByName('choice')*pi/180));
    chc(~Lg) = 0;
    filt     = 1./exp((1:20)./6);
    wk_      = filter(filt, sum(filt), chc);
    wk_      = [0; wk_(1:end-1)];

    [R, P, RLO, RUP] = corrcoef(wk_(Lgood), r(:,3));
    wkr_             = [R(2,1) RLO(2,1) RUP(2,1)];

    disp([r1 R(2,1)])
end
