%% plotML_fig7.m
% Matlab script for plotting Fig 7  (Specificity of learning effects in behavioral, MT and LIP data)
% Final version for publication to Nature Neuroscience
%
%   a) Example data showing the specificity of learning
%   b) Specificity of learning of behavioral, MT and LIP data
%

%% load data

[rzp_cy, dthp_cy, sp_cy]       = getML_psySpecificity('Cy', 0);
[rzp_zz, dthp_zz, sp_zz, pdat] = getML_psySpecificity('ZZ', 0);

[rzm_cy, dthm_cy, sm_cy] = getML_neuroSpecificity('Cy', 'MT', 0);
[rzm_zz, dthm_zz, sm_zz] = getML_neuroSpecificity('ZZ', 'MT', 0);

[rzl_cy, dthl_cy, sl_cy] = getML_neuroSpecificity('Cy', 'LIP', 0);
[rzl_zz, dthl_zz, sl_zz] = getML_neuroSpecificity('ZZ', 'LIP', 0);



%% plot data
%%%%%%%%%%%%%%%%%%%%%%%%%% plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create figure
fh = figure; 
mm2inch = 25.4;
pw = 85/mm2inch;    % paper width, make it nature size of one column, 85mm x 85mm
ph = 85/mm2inch;    % paper height
wysifig(fh, pw, ph) % set it to US letter size

% define figure parameters
alw = 0.5;  % axis line width
afs = 5;    % axis font size
lfs = 5;    % label font size
ffs = 12;   % figure font size
tfs = 6;    % text font size
lw  = 1.5; % line width
elw = 0.5;
ms  = 2;    % marker size

% define line color
lc1  = [0 0 0];
lc1e = [0.5 0.5 0.5];
lc2  = [1 0.3 0.3];
lc2e = [1 0.6 0.6];
lc3  = [0 1 1];
lc3e = [0.5 1 1];


asize = 0.3;
%% Fig 7a
% plot DV color coded by motion directions
th = annotation('textbox', [0.07-0.07 0.58+asize+0.02 0.05 0.05]);
set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'String', 'a', 'FontSize', ffs, 'FontWeight', 'bold', 'FontName', 'Helvetica') 

set(axes, 'Units', 'Normalized', 'Position', [0.11 0.565 asize asize*pw/ph], 'FontName', 'Helvetica')
td  = pdat{1};
ses = pdat{2};
DV  = pdat{3};
DVd = pdat{4};
minses = 30;
maxses = 50;

hold on
rmw = 10;
% line(repmat(ses(td>30)',2,1), [DV(td>30)-DVd(td>30); DV(td>30)+DVd(td>30)], 'Color', [0.6 0.6 0.6], 'LineWidth', 1)
% line(repmat(ses(td>-30 & td<=30)',2,1), [DV(td>-30 & td<=30)-DVd(td>-30 & td<=30); DV(td>-30 & td<=30)+DVd(td>-30 & td<=30)], 'Color', [0.6 0.6 0.6], 'LineWidth', 1)
% line(repmat(ses(td<=-30)',2,1), [DV(td<=-30)-DVd(td<=-30); DV(td<=-30)+DVd(td<=-30)], 'Color', [0.6 0.6 0.6], 'LineWidth', 1)

h1 = plot(ses(td>30), DV(td>30), '+', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', ms);
h2 = plot(ses(td>-30 & td<=30), DV(td>-30 & td<=30), 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', ms);
h3 = plot(ses(td<=-30), DV(td<=-30), 'o', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'k', 'MarkerSize', ms);

plot(ses, nanrunmean(DV,rmw), '-', 'Color', 'k', 'LineWidth', lw);
hold off
set(gca, 'LineWidth', alw, 'FontSize', afs, ...
    'XLim', [minses maxses], 'XTick', [30 40 50], ...
    'YLim', [0 15], 'YTick', [0 5 10 15])  % general
XLabel('Session (d)', 'FontSize', lfs)
YLabel('\ita\rm_{be} (spikes per s^2 per coh)', 'FontSize', lfs)
[lx ly] = axisXY2figXY(36.8,12);
h = legend([h1 h2 h3], '\fontsize{5}\theta > 30\circ',...
    '\fontsize{5}-30\circ < \theta \leq 30\circ',...
    '\fontsize{5}\theta \leq -30\circ', ...
    'Location', [lx ly 0.1*asize 0.1*asize]);
legend('boxoff')


%% Fig 7b
% plot specificity of learning data for cyrus

asize=0.3;
th = annotation('textbox', [0.5-0.1 0.58+asize+0.02 0.05 0.05]);
set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'String', 'b', 'FontSize', ffs, 'FontWeight', 'bold', 'FontName', 'Helvetica') 

set(axes, 'Units', 'Normalized', 'Position', [0.56 0.565 asize asize], 'FontName', 'Helvetica')
mx = 1.6;
% MT
hold on
plot(abs(rzm_cy), dthm_cy, 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', lc3, 'MarkerSize', ms)
plot(abs(rzm_cy), dthm_cy, 'o', 'MarkerFaceColor', lc3, 'MarkerEdgeColor', 'none', 'MarkerSize', ms)
L = ~isnan(dthm_cy);
b = regress(dthm_cy(L),[ones(size(rzm_cy(L))) abs(rzm_cy(L))]);
plot([0:0.1:mx], b(1)+b(2).*[0:0.1:mx], '-', 'LineWidth', lw, 'Color', lc3);

% LIP
plot(abs(rzl_cy), 1000*dthl_cy, 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', lc2, 'MarkerSize', ms)
plot(abs(rzl_cy), 1000*dthl_cy, 'o', 'MarkerFaceColor', lc2, 'MarkerEdgeColor', 'none', 'MarkerSize', ms)
L = ~isnan(dthl_cy);
b = regress(1000*dthl_cy(L),[ones(size(rzl_cy(L))) abs(rzl_cy(L))]);
plot([0:0.1:mx], b(1)+b(2).*[0:0.1:mx], '-', 'LineWidth', lw, 'Color', lc2);

set(gca, 'XLim', [0 mx], 'YAxisLocation', 'right') 
hold off
set(gca, 'LineWidth', alw, 'FontSize', afs, ...
         'XLim', [0 mx], 'XTick', [0 0.5 1 1.5 2.0 2.5], 'XTickLabel', {'0', '0.5', '1.0', '1.5', '2.0', '2.5'}, ...
         'YLim', [-250 250], 'YTick', [-250 -125 0 125 250])
hold off

% behavior
set(axes, 'Units', 'Normalized', 'Position', [0.56 0.565 asize asize])
hold on
plot(abs(rzp_cy), dthp_cy, 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', lc1, 'MarkerSize', ms)
plot(abs(rzp_cy), dthp_cy, 'o', 'MarkerFaceColor', lc1, 'MarkerEdgeColor', 'none', 'MarkerSize', ms)
L = ~isnan(dthp_cy);
b = regress(dthp_cy(L),[ones(size(rzp_cy(L))) abs(rzp_cy(L))]);
plot([0:0.1:mx], b(1)+b(2).*[0:0.1:mx], '-', 'LineWidth', lw, 'Color', lc1);
hold off

set(gca, 'Color', 'none', 'LineWidth', alw, 'FontSize', afs, ...
          'XLim', [0 mx], 'XTick', [0 0.5 1 1.5], 'XTickLabel', {'0', '0.5', '1.0', '1.5'}, ...
          'YLim', [-10 8], 'YTick', [-10 -5 0 5 10])
text(0.05, -10, 'Monkey C', 'FontSize', tfs, ...
     'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom')
text(-0.5, -14.3, ['Difference in \ita\rm_{be} with respect to running average (spikes per s^2 per coh)'], ...
      'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
      'FontSize', lfs, 'Rotation', 90)
text(2.1, -14.3, ['Difference in \color[rgb]{0 1 1}\ita\rm_{MT}\color[rgb]{0 0 0} & \color[rgb]{1 0.3 0.3}\ita\rm_{LIP}\color[rgb]{0 0 0} with respect to running average (spikes per s^2 per coh)'], ...
      'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
      'FontSize', lfs, 'Rotation', 90)
 
 

% Fig 5c
% plot specificity of learning data for ZsaZsa
set(axes, 'Units', 'Normalized', 'Position', [0.56 0.125 asize asize])
mx = 2.5;
% MT
hold on
plot(abs(rzm_zz), dthm_zz, 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', lc3, 'MarkerSize', ms)
plot(abs(rzm_zz), dthm_zz, 'o', 'MarkerFaceColor', lc3, 'MarkerEdgeColor', 'none', 'MarkerSize', ms)
L = ~isnan(dthm_zz);
b = regress(dthm_zz(L),[ones(size(rzm_zz(L))) abs(rzm_zz(L))]);
plot([0:0.1:mx], b(1)+b(2).*[0:0.1:mx], '-', 'LineWidth', lw, 'Color', lc3);

% LIP
plot(abs(rzl_zz), 1000*dthl_zz, 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', lc2, 'MarkerSize', ms)
plot(abs(rzl_zz), 1000*dthl_zz, 'o', 'MarkerFaceColor', lc2, 'MarkerEdgeColor', 'none', 'MarkerSize', ms)
L = ~isnan(dthl_zz);
b = regress(1000*dthl_zz(L),[ones(size(rzl_zz(L))) abs(rzl_zz(L))]);
plot([0:0.1:mx], b(1)+b(2).*[0:0.1:mx], '-', 'LineWidth', lw, 'Color', lc2);

set(gca, 'XLim', [0 mx], 'YAxisLocation', 'right') 
hold off
set(gca, 'LineWidth', alw, 'FontSize', afs, ...
         'XLim', [0 mx], 'XTick', [0 0.5 1 1.5 2.0 2.5], 'XTickLabel', {'0', '0.5', '1.0', '1.5', '2.0', '2.5'}, ...
         'YLim', [-200 200], 'YTick', [-200 -100 0 100 200])
hold off

% behavior
set(axes, 'Units', 'Normalized', 'Position', [0.56 0.125 asize asize])
hold on
plot(abs(rzp_zz), dthp_zz, 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', lc1, 'MarkerSize', ms)
plot(abs(rzp_zz), dthp_zz, 'o', 'MarkerFaceColor', lc1, 'MarkerEdgeColor', 'none', 'MarkerSize', ms)
L = ~isnan(dthp_zz);
b = regress(dthp_zz(L),[ones(size(rzp_zz(L))) abs(rzp_zz(L))]);
plot([0:0.1:mx], b(1)+b(2).*[0:0.1:mx], '-', 'LineWidth', lw, 'Color', lc1);
hold off

set(gca, 'Color', 'none', 'LineWidth', alw, 'FontSize', afs, ...
         'XLim', [0 mx], 'XTick', [0 0.5 1 1.5 2.0 2.5], 'XTickLabel', {'0', '0.5', '1.0', '1.5', '2.0', '2.5'}, ...
         'YLim', [-8 8], 'YTick', [-8 -4 0 4 8])
% xlabel('Absolute z score', 'FontSize', lfs)
% ylabel('Change in behavioral DV', 'FontSize', lfs)
text(0.05, -8, 'Monkey Z', 'FontSize', tfs, ...
     'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom')
xlabel(sprintf('Absolute z score\nof motion direction'), 'FontSize', lfs)



% get statistics
L              = ~isnan(dthm_cy) & ~isnan(rzm_cy);
[rm_cy, pm_cy] = corr(abs(rzm_cy(L)), dthm_cy(L));
L              = ~isnan(dthl_cy) & ~isnan(rzl_cy);
[rl_cy, pl_cy] = corr(abs(rzl_cy(L)), dthl_cy(L));
L              = ~isnan(dthp_cy) & ~isnan(rzp_cy);
[rp_cy, pp_cy] = corr(abs(rzp_cy(L)), dthp_cy(L));

L              = ~isnan(dthm_zz) & ~isnan(rzm_zz);
[rm_zz, pm_zz] = corr(abs(rzm_zz(L)), dthm_zz(L));
L              = ~isnan(dthl_zz) & ~isnan(rzl_zz);
[rl_zz, pl_zz] = corr(abs(rzl_zz(L)), dthl_zz(L));
L              = ~isnan(dthp_zz) & ~isnan(rzp_zz);
[rp_zz, pp_zz] = corr(abs(rzp_zz(L)), dthp_zz(L));




% % plot data
% % create figure
% fh = figure;
% pw = 8.5;    % paper width
% ph = 8.5;      % paper height
% wysifig(fh, pw, ph) % set it to US letter size
% 
% % define figure parameters
% alw = 1;    % axis line width
% afs = 14;   % axis font size
% lfs = 14;   % label font size
% tfs = 14;   % text font size
% lw  = 3;    % line width
% ms  = 8;    % marker size
% asize = 0.4;
% 
% 
% %%% plot psy specificity
% set(axes, 'Units', 'Normalized', 'Position', [0.5-asize/2 0.15+asize asize asize*pw/ph])
% hold on
% %mx = 
% 
% % for cyrus
% plot(abs(rzp_cy), dthp_cy, 'o', 'MarkerFaceColor', [0 0 1], 'MarkerEdgeColor', 'none', 'MarkerSize', 5)
% b = regress(dthp_cy,[ones(size(rzp_cy)) abs(rzp_cy)]);
% plot([0:0.1:ceil(max(abs(rzp_cy)))], b(1)+b(2).*[0:0.1:ceil(max(abs(rzp_cy)))], ...
%         '-', 'LineWidth', lw, 'Color', [0 0 1]);
% 
% % for zsa zsa
% plot(abs(rzp_zz), dthp_zz, 'o', 'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', 'none', 'MarkerSize', 5)
% b = regress(dthp_zz,[ones(size(rzp_zz)) abs(rzp_zz)]);
% plot([0:0.1:ceil(max(abs(rzp_zz)))], b(1)+b(2).*[0:0.1:ceil(max(abs(rzp_zz)))], ...
%         '-', 'LineWidth', lw, 'Color', [1 0 0]);
% hold off
% 
% % set(gca, 'LineWidth', alw, 'FontSize', afs, ...%'XLim', [0,2], 'YLim', [-60 110], ...
% %           'YTick', [-40 0 40 80])
% %XLabel('Absolute z score', 'FontSize', lfs)
% YLabel('% change in threshold', 'FontSize', lfs)
% xx = get(gca, 'XLim');
% yy = get(gca, 'YLim');
% text('Position', [0.05 0.98*yy(2)], 'String', 'Behavior', ...
%         'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', lfs+2)
% 
% 
% 
% %%% plot MT specificity
% set(axes, 'Units', 'Normalized', 'Position', [0.075 0.075 asize asize*pw/ph])
% hold on
% mx = 2;
% 
% % for cyrus
% plot(abs(rzm_cy), dthm_cy, 'o', 'MarkerFaceColor', [0 0 1], 'MarkerEdgeColor', 'none', 'MarkerSize', 5)
% b = regress(dthm_cy,[ones(size(rzm_cy)) abs(rzm_cy)]);
% plot([0:0.1:mx], b(1)+b(2).*[0:0.1:mx], '-', 'LineWidth', lw, 'Color', [0 0 1]);
% 
% % for zsa zsa
% plot(abs(rzm_zz), dthm_zz, 'o', 'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', 'none', 'MarkerSize', 5)
% b = regress(dthm_zz,[ones(size(rzm_zz)) abs(rzm_zz)]);
% plot([0:0.1:mx], b(1)+b(2).*[0:0.1:mx], '-', 'LineWidth', lw, 'Color', [1 0 0]);
% hold off
% 
% set(gca, 'LineWidth', alw, 'FontSize', afs, 'XLim', [0 2], 'YLim', [-100 100], ...
%           'YTick', [-80 -40 0 40 80])
% XLabel('Absolute z score', 'FontSize', lfs)
% YLabel('% change in ROC area', 'FontSize', lfs)
% xx = get(gca, 'XLim');
% yy = get(gca, 'YLim');
% text('Position', [0.05 0.98*yy(2)], 'String', 'MT', ...
%         'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', lfs+2)
% 
% 
% 
% %%% plot LIP specificity
% set(axes, 'Units', 'Normalized', 'Position', [0.15+asize 0.075 asize asize*pw/ph])
% hold on
% mx = 2.5;
% 
% % for cyrus
% plot(abs(rzl_cy), 1000*dthl_cy, 'o', 'MarkerFaceColor', [0 0 1], 'MarkerEdgeColor', 'none', 'MarkerSize', 5)
% b = regress(1000*dthl_cy,[ones(size(rzl_cy)) abs(rzl_cy)]);
% plot([0:0.1:mx], b(1)+b(2).*[0:0.1:mx], '-', 'LineWidth', lw, 'Color', [0 0 1]);
% 
% % for zsa zsa
% plot(abs(rzl_zz), 1000*dthl_zz, 'o', 'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', 'none', 'MarkerSize', 5)
% b = regress(1000*dthl_zz,[ones(size(rzl_zz)) abs(rzl_zz)]);
% plot([0:0.1:mx], b(1)+b(2).*[0:0.1:mx], '-', 'LineWidth', lw, 'Color', [1 0 0]);
% 
% set(gca, 'LineWidth', alw, 'FontSize', afs, 'XLim', [0 2.5], 'YLim', [-250 280], 'YTick', [-200 -100 0 100 200])
% XLabel('Absolute z score', 'FontSize', lfs)
% YLabel('% change in ROC area', 'FontSize', lfs)
% xx = get(gca, 'XLim');
% yy = get(gca, 'YLim');
% text('Position', [0.05 0.98*yy(2)], 'String', 'LIP', ...
%         'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', lfs+2)
% 
