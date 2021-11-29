% plotML_RLfig7.m
% plot changes in weights with training for fine discriminatino task

%%
% create figure
fh = figure; 
mm2inch = 25.4;
pw = 150/mm2inch;    % paper width, make it nature size of one column, 85mm x 85mm
ph = 100/mm2inch;   % paper height
wysifig(fh, pw, ph) % set it to US letter size

% define figure parameters
alw = 0.5;  % axis line width
afs = 7;    % axis font size
lfs = 7;    % label font size
ffs = 12;   % figure font size
tfs = 9;    % text font size
lw  = 1.5;  % line width
elw = 0.75;
ms  = 2;    % marker size

% define line color
lc1  = [1 0.2 0.2];
lc1e = [1 0.75 0.75];
lc2  = [0.2 0.2 1];
lc2e = [0.75 0.75 1];
lc3  = [0.2 0.2 0.2];
lc3e = [0.75 0.75 0.75];


asize = 0.18;



%% simulation parameters
recompute = 0;

Nc = 2*115985;
Nz = 77447;

TRIALS = [0 16000 80000 Nc];
bb     = 2.32;
be     = 6.53;
bw     = (6.53-2.32)/50;
THBINS = [2.^[bb:bw:be-bw]' 2.^[bb+bw:bw:be]'];
   
Monk      = 'Cy';
NSEN      = 200;
NTUNE     = 36;
ELN       = [7]*1e-8;
DIRS      = [-10,10];
SIMNUM    = [71:75];
[wavg_c, wse_c] = getML_RLFig7AvgWSortedByTh(Monk, NSEN, NTUNE, ELN, DIRS, SIMNUM, TRIALS, THBINS, recompute);


% optimal weights
NORMSCALE = sqrt(0.02);
[woavg_c, wod, woc_c, wocd] = getML_RLFig7OptimalWSortedByTh(Monk, NSEN, NTUNE, ELN, DIRS, SIMNUM, THBINS, NORMSCALE, recompute);


% % shift weights
wavg_c  = circshift(wavg_c,[0,-1]);
wse_c   = circshift(wse_c,[0,-1]);
woavg_c = circshift(woavg_c,[0,-1]);
DTUNE   = 360/NTUNE;
pdir    = linspace(-170,180,NTUNE);

%% plot weights (binned by thresholds)
% label
th = annotation('textbox', [0 0.95 0.05 0.05]);
set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'String', 'a', 'FontSize', ffs, 'FontWeight', 'bold', 'FontName', 'Helvetica') 



clim_ = [-0.004 0.004];
thm   = nangeomean(THBINS,2);
L     = thm>6 & thm<80;


% plot weights
set(axes, 'Units', 'Normalized', 'Position', [0.07 0.6 asize 1.3*asize*pw/ph], 'FontName', 'Helvetica')
imagesc(pdir, 1:sum(L), wavg_c(L,:,1))
z = wavg_c(L,:,1);  % transparency map for nans
z(isnan(wavg_c(L,:,1))) = 0;
z(~isnan(wavg_c(L,:,1))) = 1;
%alpha(z)
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, 'Color', [0.5 0.5 0.5], ...
         'XTick', [-170 -135 -90 -45 0 45 90 135 180], 'TickDir', 'in', ...
         'XTickLabel', {'-170', '', '', '', '0', '', '', '', '180'}, ...
         'YTick', [6 17 29 41], 'YTickLabel', {'8', '16', '32', '64'}, 'clim', clim_)
%xlabel('Direction tuning (deg)', 'FontSize', lfs)
ylabel('Threshold (% coh)', 'FontSize', lfs)
title('Initial', 'FontSize', tfs)
% % plot annotation arrow 
% [xb,yb] = axisXY2figXY(10,1.5);
% [xe,ye] = axisXY2figXY(10,0.2);
% anh     = annotation('arrow', [xb xe], [yb ye]);
% set(anh, 'Color', 'r', 'LineWidth', 0, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
%           'HeadWidth', 1.5*ms, 'HeadLength', 1.5*ms)
% [xb,yb] = axisXY2figXY(-10,1.5);
% [xe,ye] = axisXY2figXY(-10,0.2);
% anh     = annotation('arrow', [xb xe], [yb ye]);
% set(anh, 'Color', 'b', 'LineWidth', 0, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
%           'HeadWidth', 1.5*ms, 'HeadLength', 1.5*ms)
% plot annotation arrow 
[xb,yb] = axisXY2figXY(215,45-7+1);
[xe,ye] = axisXY2figXY(190,45-7+1);
anh     = annotation('arrow', [xb xe], [yb ye]);
set(anh, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
          'HeadWidth', ms, 'HeadLength', ms)
[xb,yb] = axisXY2figXY(215,45-18+1);
[xe,ye] = axisXY2figXY(190,45-18+1);
anh     = annotation('arrow', [xb xe], [yb ye]);
set(anh, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
          'HeadWidth', ms, 'HeadLength', ms)


 

set(axes, 'Units', 'Normalized', 'Position', [0.07+1*(0.039+asize) 0.6 asize 1.3*asize*pw/ph], 'FontName', 'Helvetica')
imagesc(pdir, 1:sum(L), wavg_c(L,:,2))
z = wavg_c(L,:,2);  % transparency map for nans
z(isnan(wavg_c(L,:,2))) = 0;
z(~isnan(wavg_c(L,:,2))) = 1;
%alpha(z)
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, 'Color', [0.5 0.5 0.5], ...
         'XTick', [-170 -135 -90 -45 0 45 90 135 180], 'TickDir', 'in', ...
         'XTickLabel', {'-170', '', '', '', '0', '', '', '', '180'}, ...
         'YTick', [6 17 29 41], 'YTickLabel', {}, 'clim', clim_)
%xlabel('Direction tuning (deg)', 'FontSize', lfs)
title('Early', 'FontSize', tfs)
hla = xlabel('Direction tuning (deg)', 'FontSize', lfs)
lpos = get(hla, 'position');
lpos(1) = 220;
set(hla, 'position', lpos)

% % plot annotation arrow 
% [xb,yb] = axisXY2figXY(10,1.5);
% [xe,ye] = axisXY2figXY(10,0.2);
% anh     = annotation('arrow', [xb xe], [yb ye]);
% set(anh, 'Color', 'r', 'LineWidth', 0, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
%           'HeadWidth', 1.5*ms, 'HeadLength', 1.5*ms)
% [xb,yb] = axisXY2figXY(-10,1.5);
% [xe,ye] = axisXY2figXY(-10,0.2);
% anh     = annotation('arrow', [xb xe], [yb ye]);
% set(anh, 'Color', 'b', 'LineWidth', 0, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
%           'HeadWidth', 1.5*ms, 'HeadLength', 1.5*ms)
% plot annotation arrow 
[xb,yb] = axisXY2figXY(215,45-7+1);
[xe,ye] = axisXY2figXY(190,45-7+1);
anh     = annotation('arrow', [xb xe], [yb ye]);
set(anh, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
          'HeadWidth', ms, 'HeadLength', ms)
[xb,yb] = axisXY2figXY(215,45-18+1);
[xe,ye] = axisXY2figXY(190,45-18+1);
anh     = annotation('arrow', [xb xe], [yb ye]);
set(anh, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
          'HeadWidth', ms, 'HeadLength', ms)

 
 

set(axes, 'Units', 'Normalized', 'Position', [0.07+2*(0.039+asize) 0.6 asize 1.3*asize*pw/ph], 'FontName', 'Helvetica')
imagesc(pdir, 1:sum(L), wavg_c(L,:,3))
z = wavg_c(L,:,3);  % transparency map for nans
z(isnan(wavg_c(L,:,3))) = 0;
z(~isnan(wavg_c(L,:,3))) = 1;
%alpha(z)
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, 'Color', [0.5 0.5 0.5], ...
         'XTick', [-170 -135 -90 -45 0 45 90 135 180], 'TickDir', 'in', ...
         'XTickLabel', {'-170', '', '', '', '0', '', '', '', '180'}, ...
         'YTick', [6 17 29 41], 'YTickLabel', {}, 'clim', clim_)
%xlabel('Direction tuning (deg)', 'FontSize', lfs)
title('Mid', 'FontSize', tfs)
% % plot annotation arrow 
% [xb,yb] = axisXY2figXY(10,1.5);
% [xe,ye] = axisXY2figXY(10,0.2);
% anh     = annotation('arrow', [xb xe], [yb ye]);
% set(anh, 'Color', 'r', 'LineWidth', 0, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
%           'HeadWidth', 1.5*ms, 'HeadLength', 1.5*ms)
% [xb,yb] = axisXY2figXY(-10,1.5);
% [xe,ye] = axisXY2figXY(-10,0.2);
% anh     = annotation('arrow', [xb xe], [yb ye]);
% set(anh, 'Color', 'b', 'LineWidth', 0, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
%           'HeadWidth', 1.5*ms, 'HeadLength', 1.5*ms)
% plot annotation arrow 
[xb,yb] = axisXY2figXY(215,45-7+1);
[xe,ye] = axisXY2figXY(190,45-7+1);
anh     = annotation('arrow', [xb xe], [yb ye]);
set(anh, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
          'HeadWidth', ms, 'HeadLength', ms)
[xb,yb] = axisXY2figXY(215,45-18+1);
[xe,ye] = axisXY2figXY(190,45-18+1);
anh     = annotation('arrow', [xb xe], [yb ye]);
set(anh, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
          'HeadWidth', ms, 'HeadLength', ms)

 

set(axes, 'Units', 'Normalized', 'Position', [0.07+3*(0.039+asize) 0.6 asize 1.3*asize*pw/ph], 'FontName', 'Helvetica')
imagesc(pdir, 1:sum(L), wavg_c(L,:,4))
z = wavg_c(L,:,4);  % transparency map for nans
z(isnan(wavg_c(L,:,4))) = 0;
z(~isnan(wavg_c(L,:,4))) = 1;
%alpha(z)
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, 'Color', [0.5 0.5 0.5], ...
         'XTick', [-170 -135 -90 -45 0 45 90 135 180], 'TickDir', 'in', ...
         'XTickLabel', {'-170', '', '', '', '0', '', '', '', '180'}, ...
         'YTick', [6 17 29 41], 'YTickLabel', {}, 'clim', clim_)
%%xlabel('Direction tuning (deg)', 'FontSize', lfs)
title('Late', 'FontSize', tfs)
% % plot annotation arrow 
% [xb,yb] = axisXY2figXY(10,1.5);
% [xe,ye] = axisXY2figXY(10,0.2);
% anh     = annotation('arrow', [xb xe], [yb ye]);
% set(anh, 'Color', 'r', 'LineWidth', 0, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
%           'HeadWidth', 1.5*ms, 'HeadLength', 1.5*ms)
% [xb,yb] = axisXY2figXY(-10,1.5);
% [xe,ye] = axisXY2figXY(-10,0.2);
% anh     = annotation('arrow', [xb xe], [yb ye]);
% set(anh, 'Color', 'b', 'LineWidth', 0, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
%           'HeadWidth', 1.5*ms, 'HeadLength', 1.5*ms)
% plot annotation arrow 
[xb,yb] = axisXY2figXY(215,45-7+1);
[xe,ye] = axisXY2figXY(190,45-7+1);
anh     = annotation('arrow', [xb xe], [yb ye]);
set(anh, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
          'HeadWidth', ms, 'HeadLength', ms)
[xb,yb] = axisXY2figXY(215,45-18+1);
[xe,ye] = axisXY2figXY(190,45-18+1);
anh     = annotation('arrow', [xb xe], [yb ye]);
set(anh, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
          'HeadWidth', ms, 'HeadLength', ms)

      
%% plot color scale
set(axes, 'Units', 'Normalized', 'Position', [0.07+4*(0.038+asize)+0.015 0.6+0.5*1.3*asize*pw/ph+0.02 0.2*asize 0.5*1.3*asize*pw/ph-0.02], 'FontName', 'Helvetica')
dat = clim_(2)*[1:-0.01:-1]';
imagesc(dat)
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, ...
         'XTick', [], 'TickDir', 'in', ...
         'XTickLabel', {}, ...
         'YTick', [1 101 201]-0.5, 'YTickLabel', {num2str(4) '0' num2str(-4) }, ...
         'clim', clim_)
text(0.5,-25, '\times10^{-3}', 'FontSize', afs)
 

%% plot optimal weight
% label
th = annotation('textbox', [0 0.45 0.05 0.05]);
set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'String', 'b', 'FontSize', ffs, 'FontWeight', 'bold', 'FontName', 'Helvetica') 


clim_ = [-0.004 0.004];
thm   = nangeomean(THBINS,2);


% plot weights
set(axes, 'Units', 'Normalized', 'Position', [0.07 0.1 asize 1.3*asize*pw/ph], 'FontName', 'Helvetica')
imagesc(pdir, 1:sum(L), woavg_c(L,:))
z = woavg_c(L,:);  % transparency map for nans
z(isnan(woavg_c(L,:))) = 0;
z(~isnan(woavg_c(L,:))) = 1;
%alpha(z)
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, 'Color', [0.5 0.5 0.5], ...
         'XTick', [-170 -135 -90 -45 0 45 90 135 180], 'TickDir', 'in', ...
         'XTickLabel', {'-170', '', '', '', '0', '', '', '', '180'}, ...
         'YTick', [6 17 29 41], 'YTickLabel', {'8', '16', '32', '64'}, 'clim', clim_)
xlabel('Direction tuning (deg)', 'FontSize', lfs)
ylabel('Threshold (% coh)', 'FontSize', lfs)
title('Optimal', 'FontSize', tfs)
% % plot annotation arrow 
% [xb,yb] = axisXY2figXY(10,1.5);
% [xe,ye] = axisXY2figXY(10,0.2);
% anh     = annotation('arrow', [xb xe], [yb ye]);
% set(anh, 'Color', 'r', 'LineWidth', 0, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
%           'HeadWidth', 1.5*ms, 'HeadLength', 1.5*ms)
% [xb,yb] = axisXY2figXY(-10,1.5);
% [xe,ye] = axisXY2figXY(-10,0.2);
% anh     = annotation('arrow', [xb xe], [yb ye]);
% set(anh, 'Color', 'b', 'LineWidth', 0, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
%           'HeadWidth', 1.5*ms, 'HeadLength', 1.5*ms)
% plot annotation arrow 
[xb,yb] = axisXY2figXY(215,45-7+1);
[xe,ye] = axisXY2figXY(190,45-7+1);
anh     = annotation('arrow', [xb xe], [yb ye]);
set(anh, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
          'HeadWidth', ms, 'HeadLength', ms)
[xb,yb] = axisXY2figXY(215,45-18+1);
[xe,ye] = axisXY2figXY(190,45-18+1);
anh     = annotation('arrow', [xb xe], [yb ye]);
set(anh, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
          'HeadWidth', ms, 'HeadLength', ms)

 
 

 
%% plot weights vs tuning for both simulations and optimal
% label
th = annotation('textbox', [0.28 0.45 0.05 0.05]);
set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'String', 'c', 'FontSize', ffs, 'FontWeight', 'bold', 'FontName', 'Helvetica') 

% sen1
set(axes, 'Units', 'Normalized', 'Position', [0.15+1*(0.038+asize) 0.1 asize 0.5*asize*pw/ph], 'FontName', 'Helvetica')
hold on
plot(pdir, zeros(size(pdir)), 'k:') 
plot(pdir, woavg_c(21,:,1), '-', 'Color', lc1, 'LineWidth', lw)
plot([pdir; pdir], [wavg_c(21,:,2)-wse_c(21,:,2); wavg_c(21,:,2)+wse_c(21,:,2)], '-', 'LineWidth', elw, 'Color', lc3e)
plot(pdir, wavg_c(21,:,2), 'o', 'MarkerFaceColor', lc3, 'MarkerEdgeColor', 'none', 'MarkerSize', ms)
hold off
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, ...
         'XLIM', [-170 180], 'YLIM', 1.3*clim_, ...
         'XTick', [-170 -135 -90 -45 0 45 90 135 180], 'TickDir', 'in', ...
         'XTickLabel', {'-170', '', '', '', '0', '', '', '', '180'}, ...
         'YTick', [-5 0 5]*1e-3)
%xlabel('Direction tuning (deg)', 'FontSize', lfs)
ylabel('Pooling weight', 'FontSize', lfs)
% % plot annotation arrow 
% [xb,yb] = axisXY2figXY(10,1.2*clim_(1));
% [xe,ye] = axisXY2figXY(10,1.3*clim_(1));
% anh     = annotation('arrow', [xb xe], [yb ye]);
% set(anh, 'Color', 'r', 'LineWidth', 0, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
%           'HeadWidth', 1.5*ms, 'HeadLength', 1.5*ms)
% [xb,yb] = axisXY2figXY(-10,1.2*clim_(1));
% [xe,ye] = axisXY2figXY(-10,1.3*clim_(1));
% anh     = annotation('arrow', [xb xe], [yb ye]);
% set(anh, 'Color', 'b', 'LineWidth', 0, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
%           'HeadWidth', 1.5*ms, 'HeadLength', 1.5*ms)
% text(0.95*-170, 0.98*1.3*clim_(2), '16.5% coh', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', tfs-3)

 
 
      

set(axes, 'Units', 'Normalized', 'Position', [0.15+2*(0.038+asize) 0.1 asize 0.5*asize*pw/ph], 'FontName', 'Helvetica')
hold on
plot(pdir, zeros(size(pdir)), 'k:') 
plot(pdir, woavg_c(21,:,1), '-', 'Color', lc1, 'LineWidth', lw)
plot([pdir; pdir], [wavg_c(21,:,3)-wse_c(21,:,3); wavg_c(21,:,3)+wse_c(21,:,3)], '-', 'LineWidth', elw, 'Color', lc3e)
plot(pdir, wavg_c(21,:,3), 'o', 'MarkerFaceColor', lc3, 'MarkerEdgeColor', 'none', 'MarkerSize', ms)
hold off
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, ...
         'XLIM', [-170 180], 'YLIM', 1.3*clim_, ...
         'XTick', [-170 -135 -90 -45 0 45 90 135 180], 'TickDir', 'in', ...
         'XTickLabel', {'-170', '', '', '', '0', '', '', '', '180'}, ...
         'YTick', [-5 0 5]*1e-3, 'YTickLabel', {})
xlabel('Direction tuning (deg)', 'FontSize', lfs)
% % plot annotation arrow 
% [xb,yb] = axisXY2figXY(10,1.2*clim_(1));
% [xe,ye] = axisXY2figXY(10,1.3*clim_(1));
% anh     = annotation('arrow', [xb xe], [yb ye]);
% set(anh, 'Color', 'r', 'LineWidth', 0, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
%           'HeadWidth', 1.5*ms, 'HeadLength', 1.5*ms)
% [xb,yb] = axisXY2figXY(-10,1.2*clim_(1));
% [xe,ye] = axisXY2figXY(-10,1.3*clim_(1));
% anh     = annotation('arrow', [xb xe], [yb ye]);
% set(anh, 'Color', 'b', 'LineWidth', 0, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
%           'HeadWidth', 1.5*ms, 'HeadLength', 1.5*ms)

 



set(axes, 'Units', 'Normalized', 'Position', [0.15+3*(0.038+asize) 0.1 asize 0.5*asize*pw/ph], 'FontName', 'Helvetica')
hold on
plot(pdir, zeros(size(pdir)), 'k:') 
plot(pdir, woavg_c(21,:,1), '-', 'Color', lc1, 'LineWidth', lw)
plot([pdir; pdir], [wavg_c(21,:,4)-wse_c(21,:,4); wavg_c(21,:,4)+wse_c(21,:,4)], '-', 'LineWidth', elw, 'Color', lc3e)
plot(pdir, wavg_c(21,:,4), 'o', 'MarkerFaceColor', lc3, 'MarkerEdgeColor', 'none', 'MarkerSize', ms)
hold off
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, ...
         'XLIM', [-170 180], 'YLIM', 1.3*clim_, ...
         'XTick', [-170 -135 -90 -45 0 45 90 135 180], 'TickDir', 'in', ...
         'XTickLabel', {'-170', '', '', '', '0', '', '', '', '180'}, ...
         'YTick', [-5 0 5]*1e-3, 'YTickLabel', {})
%xlabel('Direction tuning (deg)', 'FontSize', lfs)
% plot annotation arrow 
[xb,yb] = axisXY2figXY(10,1.4*clim_(1));
[xe,ye] = axisXY2figXY(10,1.3*clim_(1));
anh     = annotation('arrow', [xb xe], [yb ye]);
set(anh, 'Color', 'r', 'LineWidth', 0, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
          'HeadWidth', 1.5*ms, 'HeadLength', 1.5*ms)
[xb,yb] = axisXY2figXY(-10,1.4*clim_(1));
[xe,ye] = axisXY2figXY(-10,1.3*clim_(1));
anh     = annotation('arrow', [xb xe], [yb ye]);
set(anh, 'Color', 'b', 'LineWidth', 0, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
          'HeadWidth', 1.5*ms, 'HeadLength', 1.5*ms)

 




% sen2
set(axes, 'Units', 'Normalized', 'Position', [0.15+1*(0.038+asize) 0.1+0.5*asize*pw/ph+0.07 asize 0.5*asize*pw/ph], 'FontName', 'Helvetica')
hold on
plot(pdir, zeros(size(pdir)), 'k:') 
plot(pdir, woavg_c(10,:,1), '-', 'Color', lc1, 'LineWidth', lw)
plot([pdir; pdir], [wavg_c(10,:,2)-wse_c(10,:,2); wavg_c(10,:,2)+wse_c(10,:,2)], '-', 'LineWidth', elw, 'Color', lc3e)
plot(pdir, wavg_c(10,:,2), 'o', 'MarkerFaceColor', lc3, 'MarkerEdgeColor', 'none', 'MarkerSize', ms)
hold off
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, ...
         'XLIM', [-170 180], 'YLIM', 1.3*clim_, ...
         'XTick', [-170 -135 -90 -45 0 45 90 135 180], 'TickDir', 'in', ...
         'XTickLabel', {}, ...
         'YTick', [-5 0 5]*1e-3)
title('Early', 'FontSize', tfs)
ylabel('Pooling weight', 'FontSize', lfs)
text(0.95*-170, 0.98*1.3*clim_(2), '8.7% coh', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', tfs-3)



set(axes, 'Units', 'Normalized', 'Position', [0.15+2*(0.038+asize) 0.1+0.5*asize*pw/ph+0.07 asize 0.5*asize*pw/ph], 'FontName', 'Helvetica')
hold on
plot(pdir, zeros(size(pdir)), 'k:') 
plot(pdir, woavg_c(10,:,1), '-', 'Color', lc1, 'LineWidth', lw)
plot([pdir; pdir], [wavg_c(10,:,3)-wse_c(10,:,3); wavg_c(10,:,3)+wse_c(10,:,3)], '-', 'LineWidth', elw, 'Color', lc3e)
plot(pdir, wavg_c(10,:,3), 'o', 'MarkerFaceColor', lc3, 'MarkerEdgeColor', 'none', 'MarkerSize', ms)
hold off
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, ...
         'XLIM', [-170 180], 'YLIM', 1.3*clim_, ...
         'XTick', [-170 -135 -90 -45 0 45 90 135 180], 'TickDir', 'in', ...
         'XTickLabel', {}, ...
         'YTick', [-5 0 5]*1e-3, 'YTickLabel', {})
title('Mid', 'FontSize', tfs)


set(axes, 'Units', 'Normalized', 'Position', [0.15+3*(0.038+asize) 0.1+0.5*asize*pw/ph+0.07 asize 0.5*asize*pw/ph], 'FontName', 'Helvetica')
hold on
plot(pdir, zeros(size(pdir)), 'k:') 
plot(pdir, woavg_c(10,:,1), '-', 'Color', lc1, 'LineWidth', lw)
plot([pdir; pdir], [wavg_c(10,:,4)-wse_c(10,:,4); wavg_c(10,:,4)+wse_c(10,:,4)], '-', 'LineWidth', elw, 'Color', lc3e)
plot(pdir, wavg_c(10,:,4), 'o', 'MarkerFaceColor', lc3, 'MarkerEdgeColor', 'none', 'MarkerSize', ms)
hold off
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, ...
         'XLIM', [-170 180], 'YLIM', 1.3*clim_, ...
         'XTick', [-170 -135 -90 -45 0 45 90 135 180], 'TickDir', 'in', ...
         'XTickLabel', {}, ...
         'YTick', [-5 0 5]*1e-3, 'YTickLabel', {})
title('Late', 'FontSize', tfs)


     