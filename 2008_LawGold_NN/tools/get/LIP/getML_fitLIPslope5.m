function [b, bsem, p, ym] = getML_fitLIPslope5(x, coh, y, ys, b0, bcon, bl)
% fit data to slope model, nonlinear dependence on coh


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
    %%%  guess initial parameters  %%%    
    % estimate slope at all coh
    c     = unique(coh);
    slope = nans(2,length(c));
    for i = 1:length(c)
        L  = coh==c(i)&~isnan(y);
        slope(:,i) = [ones(size(x(L))) x(L)]\y(L);
    end
    
    % estimate a, the coh dependence term
    L = ~isnan(slope(2,:));
    a = [ones(sum(L),1), c(L)]\(slope(2,L)');
    
    % max firing (plateau)
    L  = coh==(max(coh));
    pk = myprctile(y(L),95);
    
    % initial conditions
    b0 = [0 a(2) 1 pk]';
end

if nargin < 6 | isempty(bcon)
    % set boundary conditions
    bcon = [0 -3 -2 0; 100 3 2 100];
end

if nargin < 7 | isempty(bl)
    L  = x==min(x);
    bl = nanmean(y(L));
end

warning off
% declare data_ to compute err
[b, f, eflag, output, l, g, h] = ...
        fmincon(@ctslope_err, b0, [], [], [], [],...
                bcon(1,:), bcon(2,:), [],...
                optimset('LargeScale', 'off', 'Display', 'off', 'Diagnostics', 'off'),...
                {x coh y ys bl});

m  = b(1)+b(2).*(coh.^b(3));
ym = m.*(x-x(1))+bl;
if any(ym>b(4))
    ym(ym>b(4)) = b(4);
end
         

%%    
% the following computation assume that errors are gaussian        

% return standard error of b using inverse of hessian matrix
if nargout > 1
    bsem = diag(h^-0.5);
end

% return significance of the coh term
if nargout > 2
    [c, ss2] = ctslope_err(b, {x coh y ys bl});
    
    % fit data to null model w/o coh term
    [b_null, f, eflag, output, l, g, h] = ...
        fmincon(@ctslopenull_err, b0([1,2,4]), [], [], [], [],...
        bcon(1,[1,2,4]), bcon(2,[1,2,4]), [],...
        optimset('LargeScale', 'off', 'Display', 'off', 'Diagnostics', 'off'),...
        {x y ys});
    [c, ss2null] = ctslopenull_err(b_null, {x y ys});
    
    % get statistics
    f = (sum(ss2null)-sum(ss2))/(sum(ss2)/(length(ss2)-length(b)));
    p = 1-fcdf(f,1,length(ss2)-length(b)-1);
end
warning on

return;


function [chi, sqy] = ctslope_err(b_, data_)
    x_   = data_{1};
    coh_ = data_{2};
    y_   = data_{3};
    ys_  = data_{4};    
    bl_  = data_{5};
    
    Lgd  = ~isnan(data_{1}) & ~isnan(data_{2}) & ~isnan(data_{3});
    
    % parts for computing err and g
    m              = b_(1)+b_(2).*(coh_.^b_(3));
    ym             = m.*(x_-x_(1))+bl_;
    if any(ym>b_(4))
        ym(ym>b_(4)) = b_(4);
    end
    
    sqy  = (y_(Lgd)-ym(Lgd)).^2;
    w    = 1./(ys_(Lgd).^2);
    
    chi  = sum(sqy.*w);
return;


function [chi sqy] = ctslopenull_err(b_, data_)
    x_   = data_{1};
    y_   = data_{2};
    ys_  = data_{3};    
    Lgd  = ~isnan(data_{1}) & ~isnan(data_{2});
    
    % parts for computing err and g
    m              = b_(2);
    ym             = m.*(x_-x_(1))+b_(1);
    if any(ym>b_(3))
        ym(ym>b_(3)) = b_(3);
    end
    
    sqy  = (y_(Lgd)-ym(Lgd)).^2;
    w    = 1./(ys_(Lgd).^2);
    
    chi  = sum(sqy.*w);    
return;
