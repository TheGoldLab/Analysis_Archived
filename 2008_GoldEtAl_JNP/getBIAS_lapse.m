function [lapse_,lse_] = getBIAS_lapse(ldat)
% function [lapse_,lse_] = getBIAS_lapse(ldat)
%
% ldat is
%   correct
%   coh
%   time

% get high coh, long time
Llapse = ldat(:,2) == max([ldat(:,2); 0.512]) & ...
    ldat(:,3) > max(0.3, prctile(ldat(:,3),40));
if sum(Llapse) > 2
    lapse_ = 1-sum(Llapse&ldat(:,1)==1)./sum(Llapse);
    lse_   = sqrt(lapse_.*(1-lapse_)./sum(Llapse));
else
    lapse_ = 0;
    lse_   = nan;
end
