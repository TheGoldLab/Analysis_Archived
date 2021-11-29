function rateN = getML_normalizedRate(mbins, rate)
% USAGE: Subtract all rate with mean baseline firing before dots on,
%        and divide all rate with 95-percentile of the firing rate (across
%        all bins and coh)

% normalize
% subtract everything with baseline (average of all bins smaller than 50ms)
% divide everything with 95 percentile of the max 
bl   = nanmean(rate(:,mbins<=50));
bl(abs(bl)==inf) = nan;
while length(bl)>1
    bl = nanmean(bl);
end

xx   = rate(5:7,mbins>=0 & mbins<=1200);
rmax = myprctile(xx(:), 97.5);
rateN = (rate-bl)/(rmax-bl);

