function [b, bsem, p, ym] = getML_fitLIPslope3(x, coh, y, ys, b0, bcon, ci_args)
% fit data to slope model, linear dependence on coh


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


if nargin < 5 | isempty(b0) | any(isnan(b0))
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
    a = [ones(sum(L),1), c(L)]\(slope(2,L)');
    
    % max firing (plateau)
    L  = coh==(max(coh));
    pk = myprctile(y(L),95);
    
    % initial conditions
    %b0 = [bl slope(2,1) a(2) pk]';
    b0 = [bl 0 a(2) pk]';
    
    % initial conditions
    b0in = b0;
    if any(isnan(b0in))
        b0(~isnan(b0in)) = b0in(~isnan(b0in));
    end

end

if nargin < 6 | isempty(bcon)
    % set boundary conditions
    bcon = [0 -3 -2 0; 100 3 2 100];
end


warning off
% declare data_ to compute err
[b, f, eflag, output, l, g, h] = ...
        fmincon(@ctslope_err, b0, [], [], [], [],...
                bcon(1,:), bcon(2,:), [],...
                optimset('LargeScale', 'off', 'Display', 'off', 'Diagnostics', 'off'),...
                {x coh y ys});

m  = b(2)+coh*b(3);
ym = m.*(x-x(1))+b(1);
if any(ym>b(4))
    ym(ym>b(4)) = b(4);
end
         

%%    
% the following computation assume that errors are gaussian        
if nargout > 1

    % confidence intervals
    if nargin < 7 || ~iscell(ci_args)

        % Compute Standard errors
        %   The covariance matrix is the negative of the inverse of the
        %   hessian of the natural logarithm of the probability of observing
        %   the data set given the optimal parameters.
        %   For now, use the numerically-estimated HESSIAN returned by fmincon
        %   (which remember is computed from -log likelihood)
        % -H because we used -logL in quick_err
        bsem = sqrt(diag(-((-h)^(-1))));

    else

        % citype is cell array with:
        %   number of simulated data sets to create
        %   confidence interval
        if isempty(ci_args)
            ci_args = {[], []};
        end

        % Check for input etype, which indicates that we will use
        %    Monte Carlo resampling (parametric bootstrap).
        %   See Wichmann & Hill (2001)
        %       The psychometric function: II. Bootstrap-
        %           based confidence intervals and sampling
        % citype{1} is number of simulated data sets to use
        if isempty(ci_args{1})
            mcn = 100;
        else
            mcn = ci_args{1};
        end
        if mcn > 0
            mcfits = zeros(mcn, length(b));
            for ii = 1:mcn
                if ~mod(ii,10)
                    disp(sprintf('Bootstrap CIs, set %d', ii))
                end
                % compute fit on simulated data set
                ysim = normrnd(y,ys);
                mcb  = fmincon(@ctslope_err, b0, [], [], [], [],...
                            bcon(1,:), bcon(2,:), [],...
                            optimset('LargeScale', 'off', 'Display', 'off', 'Diagnostics', 'off'),...
                            {x coh ysim ys});
   
                mcfits(ii,:) = mcb;
            end
            % ci_args{2} is confidence interval
            %   'sem' is 68
            %   'iqr' is 50
            %   90 is default
            if length(ci_args) < 2 || isempty(ci_args{2})
                CI = 90; % confidence interval
            elseif ischar(ci_args{2})
                switch ci_args{2}
                    case 'sem'
                        CI = 68;
                    case 'iqr'
                        CI = 50;
                    otherwise
                        CI = 90;
                end
            elseif isscalar(ci_args{2})
                CI = ci_args{2};
            end
            
%            bsem = [myprctile(mcfits,50-CI/2)' myprctile(mcfits,50+CI/2)']
            mcse = nanmean([b'-myprctile(mcfits,50-CI/2); myprctile(mcfits,50+CI/2)-b']);
            bsem = [b-mcse' b+mcse'];
  
        else
            bsem = [];
        end
    end
end

% return significance of the coh term
if nargout > 2
    [c, ss2] = ctslope_err(b, {x coh y ys});
    
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
    Lgd  = ~isnan(data_{1}) & ~isnan(data_{2}) & ~isnan(data_{3});
    
    % parts for computing err and g
    m              = b_(2)+coh_*b_(3);
    ym             = m.*(x_-x_(1))+b_(1);
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
