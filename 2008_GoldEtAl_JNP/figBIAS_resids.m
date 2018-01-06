function fig_ = figBIAS_resids(num, ex_monk, ex_ss)
% function fig_ = figBIAS_resids(num, ex_monk, ex_ss)
%
% Plots:
%   1 % correct pmf vs coh
%   2 % rt choice pmf vs signed coh
%   3 resids

if nargin < 1 || isempty(num)
    num = 2;
end

if nargin < 2 || isempty(ex_monk)
    %ex_monk = 'Ava';
    ex_monk = 'Atticus';
    %ex_monk = 'Cyrus';
    %ex_monk = 'ZZ';
end

if nargin < 3 || isempty(ex_ss)
    ex_ss = 185;
end
fun = @ddExp3;

%%
% set up the figure
%%
% units should be in inches, from wysifig
wid  = 4.0; % total width
hts  = 1.3;
cols = {[.4 .6]};
[axs,fig_] = getBIAS_axes(num, wid, hts, cols);

%%
% Get the data
%%
% data rows are trials, columns are:
%   1 coherence (0...1)
%   2 time (seconds)
%   3 dot dir (-1/0/1)
%   4 choice (-1/1)
%   5 correct (0/1)
%
data      = FS_getDotsTrainingData(ex_monk);
sessions  = unique(data(:,1));
Lgood     = data(:,1) == sessions(ex_ss) & data(:,2) <= 2 & data(:,3)>=0;
data      = data(Lgood, [5 6 8 9 3]);

% do the fit
[f,s,t,p,r] = ctPsych_fit(fun, data(:,1:4), data(:,5));
os = ones(1,3);
sp = {'line_style', '-', 'colors', {0.8.*os 0.4.*os 0.*os}};
%% PANEL xx: % Correct PMF
%
% axes(axs(1)); cla reset; hold on;
% ctPsych_plot(fun, f, data(:,1:2), data(:,5), 3, [], axs(1), sp);

%% PANEL 1: % rt choices PMF
%
axes(axs(1)); cla reset; hold on;
ctPsych_plot(fun, f, data(:,1:4), data(:,5), 3, [], axs(1), sp);
ylabel('% right choices')

%% PANEL 2: sample residuals
%
axes(axs(2)); cla reset; hold on;
set(axs(2), 'FontSize', 14);
exi   = 140:154;
srdat = data(exi,:);
srdat(srdat(:,4)==-1,4)=0;
Lrc   = srdat(:,4)==1;
Lld   = srdat(:,3)==-1;
Lcor  = srdat(:,5)==1;
inds  = (1:size(srdat,1))';

plot([0 inds(end)+1], [0.5 0.5], 'k:')

preds      = p(exi);
preds(Lld) = 1-preds(Lld);
plot(inds, preds, 'k.')

plot([inds(Lrc) inds(Lrc)]', [preds(Lrc) srdat(Lrc,4)]', 'k-')
plot(inds(Lrc&Lcor),   srdat(Lrc&Lcor,4 ),   'ko', 'MarkerFaceColor', 'w');
plot(inds(Lrc&~Lcor),  srdat(Lrc&~Lcor, 4),  'kx');

plot([inds(~Lrc) inds(~Lrc)]', [preds(~Lrc) srdat(~Lrc,4)]', 'k--')
plot(inds(~Lrc&Lcor),  srdat(~Lrc&Lcor, 4),  'ko', 'MarkerFaceColor', 'w');
plot(inds(~Lrc&~Lcor), srdat(~Lrc&~Lcor, 4), 'kx');

xlim([0 inds(end)+1])
xlabel('Trial number')

return
%% PANEL 3: histogram of residuals
%
axes(axs(3)); cla reset; hold on;
res = r(:,4);
[na, xa] = hist(res, -1:.1:1, 'k');
hmax = ceil(max(na)*1.1);
plot([0 hmax], repmat(median(res), 1, 2), 'r-', 'LineWidth', 2);
hold on;
barh(xa, na, 'FaceColor', 'k');
set(axs(3), 'XLim', [0 hmax], 'YLim', [-1 1], 'XTickLabel', [], ...
    'XTick', [], 'YTickLabel', [], 'Box', 'off');
disp(median(res))
xlabel('Count')

