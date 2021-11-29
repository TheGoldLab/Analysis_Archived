function [b, bsem, gof, ym] = getML_fitLIPslope2(x, coh, y, ys, b0, bcon)
% fit data to slope model without plateau, coh is a power function


% check inputs
if nargin < 3
    error('sramp_fitW needs at least three inputs')
    return;
elseif size(x)~=size(y)
    error('x and y has to be of the same size')
    return;
end

if nargin < 4 | isempty(ys)
    % same as normal least sq fit this case
    ys = ones(size(y));
end


if nargin < 5 | isempty(b0)
    % guess initial parameters
    L  = x==min(x);
    bl = nanmean(y(L));
    L  = coh==(min(coh));
    m0 = [ones(size(x(L))) x(L)]\y(L);  % slope at lowest coh
    L  = x==median(x);
    a  = [ones(size(coh(L))) log(coh(L))]\y(L);  % pick the median time to estimate coh dependence
    b0 = [bl m0(2) a(2) 1]';
end

if nargin < 6 | isempty(bcon)
    % set boundary conditions
    bcon = [0 -3 -1000 -10; 100 3 1000 10];
end

coh(coh==0) = 1;

warning off
% declare data_ to compute err
[b, f, eflag, output, l, g, h] = ...
        fmincon(@ctslope_err, b0, [], [], [], [],...
                bcon(1,:), bcon(2,:), [],...
                optimset('LargeScale', 'off', 'Display', 'off', 'Diagnostics', 'off'),...
                {x coh y ys});

m  = b(2)+b(3)*(coh.^b(4));
ym = m.*(x-x(1))+b(1);
     

%%    
% the following computation assume that errors are gaussian        

% return standard error of b using inverse of hessian matrix
if nargout > 1
    bsem = diag(h^-0.5);
end
% return goodness of fit
% gof = probability of obtaining observations given the model is true
if nargout > 2
    gof = 1-chi2cdf(ctslope_err(b, {x coh y ys}), length(x)-length(b));
end
warning on

return;


function chi = ctslope_err(b_, data_)
    x_   = data_{1};
    coh_ = data_{2};
    y_   = data_{3};
    ys_  = data_{4};    
    Lgd  = ~isnan(data_{1}) & ~isnan(data_{2}) & ~isnan(data_{3});
    
    % parts for computing err and g
    m              = b_(2)+b_(3)*(coh_.^b_(4));
    ym             = m.*(x_-x_(1))+b_(1);
    
    sqy  = (y_(Lgd)-ym(Lgd)).^2;
    w    = 1./(ys_(Lgd).^2);
    
    chi  = sum(sqy.*w);    
return;

