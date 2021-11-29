function [fits, sems, p] = getML_neuroACT(data)
% fit neural responses to A*Cm*Tn model using maximum likelihood regression
% data inputs are:
%   1 ... fraction coherence [0...1]
%   2 ... time in sec        [0...1]
%   3 ... mean response      [spikes/sec]
%   4 ... standard error     [spikes/sec]


% estimate initial parameters
% estimate baseline using mean of the first bin
Lt = data(:,2)==nanmin(data(:,2));
bl = nanmean(data(Lt,3));

% estimate A
Uc = nonanunique(data(:,1));
Lt = data(:,2)==nanmax(data(:,2));
r  = nans(length(Uc),1);
for i = 1:length(Uc)
    Lc   = data(:,1)==Uc(i);
    r(i) = nanmean(data(Lt&Lc,3));
end
A  = Uc\r;
    
bcon = [bl  0 100;  ...
        A   0 1000; ...
        1   1    1; ...
        1   1    1];


warning off

% declare data_ to compute err
[fits, f, eflag, output, l, g, h] = ...
        fmincon(@ACT_err, bcon(:,1), [], [], [], [],...
                bcon(:,2), bcon(:,3), [],...
                optimset('LargeScale', 'off', 'Display', 'off', 'Diagnostics', 'off'),...
                data);


% return standard error of b using inverse of hessian matrix
if nargout > 1
    sems = diag(h^-0.5);
end
 

% return p-value for fitted parameters, the null hypothesis is that the
% fitted parameters is equal to 0, except for A. H0 for A is A=1;
if nargout > 2
    p = nans(3,1);
    for i = 1:3
        LLF  = -ACT_err(fits, data);
        LLR  = -ACTReduced_err(fits, data, i);
        
        if LLR < LLF
            p(i) = 1-chi2cdf(2*(LLF-LLR),1);
            
        else
            p(i) = 1;
            
		end
    end
end

return;


function NLLF = ACT_err(b, d)
    % get good trials
    Lgd  = ~isnan(d(:,1)) & ~isnan(d(:,2)) & ~isnan(d(:,3)) & ~isnan(d(:,4));   % good trials
    tbe  = nanmean(d(:,2));
    
    % computer negative of the log likelihood value
    pred = b(1) + b(2) .* (d(Lgd,1).^b(3)) .* ((d(Lgd,2)-tbe).^b(4));
    NLLF = -sum(log(normpdf(d(Lgd,3), pred, d(Lgd,4))));
return;


function NLLF = ACTReduced_err(b, d, i) % get log likelihood for reduced model
    % get good trials
    Lgd  = ~isnan(d(:,1)) & ~isnan(d(:,2)) & ~isnan(d(:,3)) & ~isnan(d(:,4));   % good trials
    tbe  = nanmean(d(:,2));
    
    switch i
        case 1
            pred = b(1) + (d(Lgd,1).^b(3)).*((d(Lgd,2)-tbe).^b(4));
          
        case 2
            pred = b(1) + b(2).*((d(Lgd,2)-tbe).^b(4));
            
        case 3
            pred = b(1) + b(2).*(d(Lgd,1).^b(3));
    end
    
    % computer negative of the log likelihood value
    NLLF = -sum(log(normpdf(d(Lgd,3), pred, d(Lgd,4))));
return;


