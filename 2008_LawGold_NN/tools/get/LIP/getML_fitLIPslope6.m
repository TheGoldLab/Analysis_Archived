function [B, Bsem, ym] = getML_fitLIPslope6(x, coh, y, ys, B0, Bcon)
% Fit LIP responses to a piecewise-linear model. The model assumes that LIP responses (spike rate or roc areas) depends on
% both coh and time, and is given by:
%
%
%          / = b0                     , if T<t1
%   R(C,T) - = b0 + (b1*C+b2)(T-t1)   , if t1<=T<t2
%          \ = b0 + (b1*C+b2)(t2-t1)  , if T>=t2
%
%
%  INPUTS:
%   x    - time in second       
%   coh  - fraction coherence (0...1)
%   y    - LIP responses
%   ys   - se of LIP responses
%   B0   - initial conditions
%   Bcon - upper and lower bound for model parameters B
%
%
%  OUTPUTS:
%   B    - model parameters [b0; b1; b2; t1; t2]
%   Bsem - standard error for B
%   ym   - predicted values from model

% check inputs
if nargin < 3
    error('getML_fitLIPslope6 needs at least three inputs')
    return;
elseif size(x)~=size(y)
    error('x and y should have the same size')
    return;
end

if nargin < 4 | isempty(ys)
    % same as normal least sq fit this case
    ys = ones(size(y));
end

if nargin < 6 | isempty(Bcon)
    % set boundary conditions
    Bcon = [0 -3 -2 0; 100 3 2 100];
end

if nargin < 5 | isempty(B0) | any(isnan(B0))
    %%%  guess initial parameters  %%%    
    % base line firing (B0)
    L  = x==min(x);
    bl = nanmean(y(L));
    
    % coh dependence (B1) and slope at 0% coh (B2)
    c     = unique(coh);
    slope = nans(2,length(c));
    for i = 1:length(c)
        L  = coh==c(i)&~isnan(y);
        slope(:,i) = [ones(size(x(L))) x(L)]\y(L);
    end
    L = ~isnan(slope(2,:));
    a = [ones(sum(L),1), c(L)]\(slope(2,L)');
    
    % begin and end time for ramp-like response
    L99 = coh==0.999;
    L51 = coh==0.512;
    if sum(L99)>=sum(L51)
        Lgd = L99;
        bl  = nanmean(y(Lgd&x<=200&x>=0));
        bHi = sramp_fitW(x(Lgd), y(Lgd), ys(Lgd), bl, [],...
            [min(x(Lgd)) 0 min(x(Lgd)); max(x(Lgd)) 30 max(x(Lgd))]);
    else
        Lgd = L51;
        bl  = nanmean(y(Lgd&x<=200&x>=0));
        bHi = sramp_fitW(x(Lgd), y(Lgd), ys(Lgd), bl, [],...
            [min(x(Lgd)) 0 min(x(Lgd)); max(x(Lgd)) 30 max(x(Lgd))]);
    end
    
    % initial conditions
    B0in = B0;
    B0   = [bl a(2) a(1) bHi(1) nanmax(y)];
    if any(isnan(B0in))
        B0(~isnan(B0in)) = B0in(~isnan(B0in));
    end
end



warning off
[B, f, eflag, output, l, g, h] = ...
        fmincon(@ctslope_err, B0, [], [], [], [],...
                Bcon(1,:), Bcon(2,:), [],...
                optimset('LargeScale', 'off', 'Display', 'off', 'Diagnostics', 'off'),...
                {x coh y ys});

[B; B0]            
% return ym
m              = coh*B(2)+B(3);
ym             = m.*(x-B(4))+B(1);
if any(x<B(4))
    ym(x<B(4)) = B(1);
end

if any(x>=B(5))
    ym(x>=B(5)) = m(x>=B(5)).*(B(5)-B(4))+B(1);
end


%%    
% the following computation assume that errors are gaussian        

% return standard error of b using inverse of hessian matrix
if nargout > 1
    Bsem = diag(h^-0.5);
end

warning on

return;


function NLL = ctslope_err(b_, data_)
    x_   = data_{1};
    coh_ = data_{2};
    y_   = data_{3};
    ys_  = data_{4};    
    Lgd  = ~isnan(data_{1}) & ~isnan(data_{2}) & ~isnan(data_{3});
    
    % parts for computing err and g
    m              = coh_*b_(2)+b_(3);
    ym             = m.*(x_-b_(4))+b_(1);
    if any(x_<b_(4))
        ym(x_<b_(4))  = b_(1);
    end
    if any(x_>=b_(5))
        ym(x_>=b_(5)) = m(x_>=b_(5)).*(b_(5)-b_(4))+b_(1);
    end
    
    % return negative of log likelihood
    NLL = -sum(log(normpdf(y_(Lgd),ym(Lgd),ys_(Lgd))));
return;



