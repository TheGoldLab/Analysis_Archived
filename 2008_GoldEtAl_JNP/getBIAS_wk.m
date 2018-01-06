function [wk_, r_, fout_, fits_, stds_] = getBIAS_wk(vin, vout, lag_max, lag_step, tosum, tofit, finit)
% function [wk_, r_, fout_, fits_, stds_] = getBIAS_wk(vin, vout, lag_max, lag_step, tosum, tofit, finit)
%
% wk_ is matrix of kernels
%   D1 (rows)    are lags
%   D2 (columns) are vins
%   D3 are "tofit" ... first is raw
%
% r_ is matrix of correlation coefficients
%   D1 rows are lags
%   D2 columns are vins; last is "tosum"
%   D3 are "tofit" ... first is raw
%
% fout_ is choices filtered by given wk
%   D1 rows are trials
%   D2 columns are wks associated with given vin; last is "tosum"
%   D3 are wks from "tofit" ... first is raw
% 
% stds_ is matrix of stds of raw wks
%   D1 (rows)    are lags
%   D2 (columns) are vins
%
% fits_ is cell array of best-fitting parameters from "tofit"
%   
% default "to fit" are:
%   kPow
%   kExp1
%   kExp2
%   kExp3
if nargin < 3 || isempty(lag_max)
    lag_max = 800;
    r_max   = lag_max;
elseif isscalar(lag_max)
    r_max   = lag_max;
elseif length(lag_max == 2);
    r_max   = lag_max(2);
    lag_max = lag_max(1);
end

if nargin < 4 || isempty(lag_step)
    lag_step = 10;
end

% get kernels
nv  = size(vin, 2);
wk_ = nans(lag_max, nv);
%st  = nans(lag_max, nv);
for vv = 1:nv
    wk_(:,vv) = getWienerK(vin(:,vv), vout, lag_max);
%    [wk_(:,vv),r,st(:,vv)] = getWienerK(vin(:,vv), vout, lag_max);
end

if nargin < 5 || isempty(tosum)
    if size(vin, 2) > 1
        tosum = 1:size(vin,2);
    else
        tosum = [];
    end
end

% if necessary, compute correlation coefficients
if nargout > 1
    [r_,f] = get_rs(vin, vout, wk_, r_max, lag_step, tosum);
end

% return filtered values
if nargout > 2
    fout_ = f;
end

if nargout > 3
    fits_ = {};
end

if nargout > 4
    stds_ = st;
end

% get fit
if nargin >= 6 && ~isempty(tofit)
    if ischar(tofit) && strcmp(tofit, 'all')
        tofit = {@kPow, @kExp1, @kExp2, @kExp3, @kExp2a};
    elseif ~iscell(tofit)
        tofit = {tofit};
    end
    if nargin < 7
        finit = {};
    end
    nfits  = length(tofit);
    wk_    = cat(3, wk_,   nans(lag_max, nv, nfits));
    r_     = cat(3, r_,    nans(size(r_,1), size(r_,2), nfits));
    fout_  = cat(3, fout_, nans(size(fout_,1), size(fout_,2), nfits));
    for fi = 1:nfits
        [inits,A,b] = feval(tofit{fi});
        is          = (1:lag_max)';
        fits        = nans(size(inits,1), nv);        
        for vi = 1:nv
            if ~isempty(finit)
                inits(:,1) = finit{fi}(:,vi);
            end
            fits(:,vi) = fmincon(@werr, inits(:,1), A, b, [], [], ...
                inits(:,2), inits(:,3), [], ...
                optimset('LargeScale', 'off', 'Display', 'off', 'Diagnostics', 'off'), ...
                wk_(:,vi), tofit{fi}, is);
            wk_(:, vi, fi+1) = feval(tofit{fi}, fits(:,vi), is);
        end
        [r_(:,:,fi+1),fout_(:,:,fi+1)] = ...
            get_rs(vin, vout, wk_(:,:,fi+1), r_max, lag_step, tosum);
        fits_ = cat(1, fits_, fits);
    end
end

%% SUBFUNCTION get_rs
%
% returns:
%   rs_   ... matrix of correlation coefficients
%   fout_ ... matrix of filtered choices -- uses the one
%               with the highest correlation coefficient
function [rs_, fout_] = get_rs(vin, vout, wk, r_max, lag_step, tosum)

vr  = size(vin, 1);
vc  = size(vin, 2) + double(~isempty(tosum));
if lag_step == 0
    steps = r_max;
    rs_   = nans(3, vc);
else
    steps = 1:lag_step:r_max;
    rs_   = nans(r_max, vc);
end
fout_ = nans(vr, vc);
fout  = nans(vr, vc);
X     = [nans(vr, 1) vout-mean(vout)];
rmax  = nans(vc,2);
for ss = steps
    % compute relative weights
    % if ~isempty(tosum)
    %     wts = sum(abs(wk(1:ss,tosum)),1);
    %     wts = wts./sum(wts);
    % end
    for vv = 1:vc
        if ~isempty(tosum) && vv == vc
            % fout(:,vv) = sum(fout(:,tosum).*repmat(wts,size(fout,1),1),2);
            fout(:,vv) = sum(fout(:,tosum),2);
        else
            fout(:,vv) = filter(wk(1:ss,vv), 1, vin(:,vv));
        end
        if lag_step == 0
            [R,P,RLO,RUP] = corrcoef(fout(:,vv), vout);
            rs_(:,vv)     = [R(2,1), RLO(2,1), RUP(2,1)];
            fout_(:,vv)   = fout(:,vv);
        else
            X(:,1)     = fout(:,vv) - sum(fout(:,vv))/vr;
            c          = (X' * X) / (vr - 1);
            d          = sqrt(diag(c));
            R          = c./(d*d');
            rs_(ss,vv) = R(2,1);
            if isnan(rmax(vv)) || rs_(ss,vv) > rmax(vv,1)
                rmax(vv,:)  = [rs_(ss,vv) ss];
                fout_(:,vv) = fout(:,vv);
            end
        end
    end
end
%disp(rmax)

% err func for fitting
function err_ = werr(fits, wk, fun, is)

%    err_ = -sum(normpdf(feval(fun,fits,is),wk,st));
err_ = sum(abs(wk-feval(fun, fits, is)));
%    err_ = sum((wk-feval(fun, fits, is)).^2);

%% SUBFUNCTION: kExp1
%
% func is:
%   a./n.*exp(-i/tau)
%
%   where n = sum(exp(-i/tau))
%
%   fits(1) = a
%   fits(2) = tau
function [ks_, A_, b_] = kExp1(fits, is)

if nargin < 2
    ks_ = [ ...
        1 -100    100; ...
        1    0.1 1000];
    A_  = [];
    b_  = [];

else
    e1  = exp(-is./fits(2));
    ks_ = fits(1)./sum(e1).*e1;
end

%% SUBFUNCTION: kExp2
%
% func is:
%   a1./n1.*exp(-i/tau1) + a2./n2.*exp(-i/tau2)
%
%   fits(1) = a1
%   fits(2) = a2
%   fits(3) = tau1
%   fits(4) = tau2
function [ks_, A_, b_] = kExp2(fits, is)

if nargin < 2
    ks_ = [ ...
        1 -100    100; ...   % 1
        1 -100    100; ...   % 1
        1   0.1  100; ...   % 1
        20   0.1 1000];      % 20
    A_  = [0 0 1 -1];
    b_  = 0;

else
    e1  = exp(-is./fits(3));
    e2  = exp(-is./fits(4));
    ks_ = fits(1)./sum(e1).*e1 + fits(2)./sum(e2).*e2;
end


%% SUBFUNCTION: kExp2a
%
% func is:
%   f5+ a1./n1.*exp(-i/tau1) + a2./n2.*exp(-i/tau2)
%
%   fits(1) = a1
%   fits(2) = a2
%   fits(3) = tau1
%   fits(4) = tau2
%   fits(5) = lower asymptote 1
function [ks_, A_, b_] = kExp2a(fits, is)

if nargin < 2
    ks_ = [ ...
        1 -100    100; ...
        1 -100    100; ...
        1    0.1  500; ...
        20   0.1  500; ...
        0 -100    100];
    A_  = [0 0 1 -1 0];
    b_  = 0;

else
    e1  = exp(-is./fits(3));
    e2  = exp(-is./fits(4));
    ks_ = fits(5) + fits(1)./sum(e1).*e1 + fits(2)./sum(e2).*e2;
end

%% SUBFUNCTION: kExp3
%
% func is:
%   a1./n1.*exp(-i/tau1) + ...
%   a2./n2.*exp(-i/tau2) + ...
%   a3./n3.*exp(-i/tau3)
%
%   fits(1) = a1
%   fits(2) = a2
%   fits(3) = a3
%   fits(4) = tau1
%   fits(5) = tau2
%   fits(6) = tau3
function [ks_, A_, b_] = kExp3(fits, is)

if nargin < 2
    ks_ = [ ...
        1 -100    100; ...
        1 -100    100; ...
        1 -100    100; ...
        1    0.1 1000; ...
        10   0.1 1000; ...
        10   0.1 1000];
    A_  = [];
    b_  = [];
else
    e1  = exp(-is./fits(4));
    e2  = exp(-is./fits(5));
    e3  = exp(-is./fits(6));
    ks_ = ...
        fits(1)./sum(e1).*e1 + ...
        fits(2)./sum(e2).*e2 + ...
        fits(3)./sum(e3).*e3;
end

%% SUBFUNCTION: kPow
%
% func is:
%   a.*i.^n
%
%   fits(1) = a
%   fits(2) = n
function [ks_, A_, b_] = kPow(fits, is)

if nargin < 2
    ks_ = [ ...
        1 -100 100; ...
        1 -100 100];
    A_  = [];
    b_  = [];

else
    e1  = exp(-is./fits(2));
    ks_ = fits(1).*is.^fits(2);
end
