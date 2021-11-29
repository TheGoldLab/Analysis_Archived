function [b, bsem, gof, ym] = getML_fitLIPslope4(x, coh, y, ys, b0, bcon)
% fit data to slope model, log dependence on coh


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

% fix coh for log scale
if max(coh)==0.999
    if any(coh==0)
        coh(coh==0)=0.01;
    end
elseif max(coh)==99.9
    if any(coh==0)
        coh(coh==0)=1;
    end
end

if nargin < 5 | isempty(b0)
    %%%  guess initial parameters  %%%    
    % base line firing
    L  = x==min(x);
    bl = nanmean(y(L));
    
    % estimate slope at all coh
    c     = unique(coh);
    slope = nans(2,length(c));
    for i = 1:length(c)
        L  = coh==c(i)&~isnan(y);
        slope(:,i) = [ones(size(x(L))) x(L)]\y(L);
    end
    
    % estimate a, the coh dependence term
    L = ~isnan(slope(2,:));
    a = [ones(sum(L),1), log(c(L))]\(slope(2,L)');
    
    % max firing (plateau)
    L  = coh==(max(coh));
    pk = myprctile(y(L),95);
    
    % initial conditions
    b0 = [bl slope(2,1) a(2) pk]';
end

if nargin < 6 | isempty(bcon)
    % set boundary conditions
    bcon = [0 -3 -2 0; 100 3 2 100];
end

%%% NEW %%%
% limit plateau values
bcon(:,4) = [0.8*pk; 1.2*pk];

warning off
% declare data_ to compute err
[b, f, eflag, output, l, g, h] = ...
        fmincon(@ctslope_err, b0, [], [], [], [],...
                bcon(1,:), bcon(2,:), [],...
                optimset('LargeScale', 'off', 'Display', 'off', 'Diagnostics', 'off'),...
                {x coh y ys});

m  = b(2)+log(coh)*b(3);
ym = m.*(x-x(1))+b(1);
if any(ym>b(4))
    ym(ym>b(4)) = b(4);
end
         

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
    m              = b_(2)+log(coh_)*b_(3);
    ym             = m.*(x_-x_(1))+b_(1);
    if any(ym>b_(4))
        ym(ym>b_(4)) = b_(4);
    end
    
    sqy  = (y_(Lgd)-ym(Lgd)).^2;
    w    = 1./(ys_(Lgd).^2);
    
    chi  = sum(sqy.*w);    
return;

