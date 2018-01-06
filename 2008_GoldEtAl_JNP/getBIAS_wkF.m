function [wk_, r_, fout_] = getBIAS_wkF(vin, vout, lag_m2)
% function [wk_, r_, fout_] = getBIAS_wkF(vin, vout, lag_m2)
%
% see getBIAS_wk for details.
%   returns the kernel, r [r ci_lo ci_hi len], filtered outputs
%   for the kernel that provides the highest r

lag_max    = 600;
lag_step   = 10;
tosum      = [1 2];
tofit      = {'kExp2'};

% find max r
if nargin < 3 || isempty(lag_m2)
    [wk,r,c,d] = getBIAS_wk(vin, vout, lag_max, lag_step, tosum, tofit);
    lag_m2     = find(r(:,end,end)==max(r(:,end,end)),1)+lag_step;
end

% recompute wks
[wk, r, ff, d] = getBIAS_wk(vin, vout, lag_m2, 0, tosum, tofit);

% return just fit kernels, rstuff, fout
wk_                 = zeros(lag_max, 2);
wk_(1:size(wk,1),1) = wk(:,2,1);
wk_(1:size(wk,1),2) = wk(:,2,2);
r_                  = [r(:,end,end); lag_m2];
fout_               = ff(:,end,end);

return

subplot(2,1,1);cla reset; hold on;
plot(wk(:,1,1), 'k-');
plot(wk(:,1,2), 'r-');

subplot(2,1,2);cla reset; hold on;
plot(wk(:,2,1), 'k-');
plot(wk(:,2,2), 'r-');

r = input('next');