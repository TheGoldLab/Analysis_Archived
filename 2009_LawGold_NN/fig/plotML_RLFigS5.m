% plotML_RLfigS9.m.m
% plot optimal weights for different correlation structures

%%
% create figure
fh = figure; 
mm2inch = 25.4;
pw = 150/mm2inch;    % paper width, make it nature size of one column, 85mm x 85mm
ph = 150/mm2inch;   % paper height
wysifig(fh, pw, ph) % set it to US letter size

% define figure parameters
alw = 0.5;  % axis line width
afs = 5;    % axis font size
lfs = 6;    % label font size
ffs = 12;   % figure font size
tfs = 10;    % text font size
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


%% get average weights
%% simulation parameters
% no corr
recompute = 0;

Nc = 115985;
Nz = 77447;

TRIALS = [0 2000 40000 Nc];
bb     = 2.32;
be     = 6.53;
bw     = (6.53-2.32)/50;
THBINS = [2.^[bb:bw:be-bw]' 2.^[bb+bw:bw:be]'];
   
Monk      = 'Cy';
NSEN      = 200;
NTUNE     = 36;
ELN       = [7]*1e-8;
DIRS      = [180 0];
SIMNUM    = [81];
TAG       = 'NoCorr';
[wnc, wse_c] = getML_RLFigS9AvgWSortedByTh(Monk, NSEN, NTUNE, ELN, DIRS, SIMNUM, TRIALS, THBINS, TAG, recompute);
wnc = circshift(wnc,[0,-10]);
DTUNE   = 360/NTUNE;
pdir    = linspace(-170,180,NTUNE);


% Shadlen corr
recompute = 0;

Nc = 115985;
Nz = 77447;

TRIALS = [0 2000 40000 Nc];
bb     = 2.32;
be     = 6.53;
bw     = (6.53-2.32)/50;
THBINS = [2.^[bb:bw:be-bw]' 2.^[bb+bw:bw:be]'];
   
Monk      = 'Cy';
NSEN      = 200;
NTUNE     = 36;
ELN       = [7]*1e-8;
DIRS      = [180 0];
SIMNUM    = [81];
TAG       = 'ShadCorr';
[wsc, wse_c] = getML_RLFigS9AvgWSortedByTh(Monk, NSEN, NTUNE, ELN, DIRS, SIMNUM, TRIALS, THBINS, TAG, recompute);
wsc = circshift(wsc,[0,-10]);


% my corr
recompute = 0;

Nc = 115985;
Nz = 77447;

TRIALS = [0 2000 40000 Nc];
bb     = 2.32;
be     = 6.53;
bw     = (6.53-2.32)/50;
THBINS = [2.^[bb:bw:be-bw]' 2.^[bb+bw:bw:be]'];
   
Monk      = 'Cy';
NSEN      = 200;
NTUNE     = 36;
ELN       = [7]*1e-8;
DIRS      = [180 0];
SIMNUM    = [51];
TAG       = '';
[wmc, wse_c] = getML_RLFigS9AvgWSortedByTh(Monk, NSEN, NTUNE, ELN, DIRS, SIMNUM, TRIALS, THBINS, TAG, recompute);
wmc = circshift(wmc,[0,-10]);


%% plot weights (binned by thresholds)
% label



clim_ = [-0.004 0.004];
thm   = nangeomean(THBINS,2);
L     = thm>6 & thm<80;


% plot weights
set(axes, 'Units', 'Normalized', 'Position', [0.07 0.6 asize 1.3*asize*pw/ph], 'FontName', 'Helvetica')
imagesc(pdir, 1:sum(L), wnc(L,:,4))
z = wnc(L,:,4);  % transparency map for nans
z(isnan(wnc(L,:,4))) = 0;
z(~isnan(wnc(L,:,4))) = 1;
alpha(z)
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, 'Color', [0.5 0.5 0.5], ...
         'XTick', [-170 -135 -90 -45 0 45 90 135 180], 'TickDir', 'in', ...
         'XTickLabel', {'-170', '', '', '', '0', '', '', '', '180'}, ...
         'YTick', [6 17 29 41], 'YTickLabel', {'8', '16', '32', '64'}, 'clim', clim_)
%xlabel('Direction tuning (degree)', 'FontSize', lfs)
ylabel('Threshold (% COH)', 'FontSize', lfs)
%title('Initial', 'FontSize', tfs)
% plot annotation arrow 
[xb,yb] = axisXY2figXY(0,1.5);
[xe,ye] = axisXY2figXY(0,0.2);
anh     = annotation('arrow', [xb xe], [yb ye]);
set(anh, 'Color', 'r', 'LineWidth', 0, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
          'HeadWidth', 1.5*ms, 'HeadLength', 1.5*ms)
[xb,yb] = axisXY2figXY(180,1.5);
[xe,ye] = axisXY2figXY(180,0.2);
anh     = annotation('arrow', [xb xe], [yb ye]);
set(anh, 'Color', 'b', 'LineWidth', 0, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
          'HeadWidth', 1.5*ms, 'HeadLength', 1.5*ms)



set(axes, 'Units', 'Normalized', 'Position', [0.07+1*(0.03+asize) 0.6 asize 1.3*asize*pw/ph], 'FontName', 'Helvetica')
imagesc(pdir, 1:sum(L), wsc(L,:,4))
z = wsc(L,:,4);  % transparency map for nans
z(isnan(wsc(L,:,4))) = 0;
z(~isnan(wsc(L,:,4))) = 1;
alpha(z)
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, 'Color', [0.5 0.5 0.5], ...
         'XTick', [-170 -135 -90 -45 0 45 90 135 180], 'TickDir', 'in', ...
         'XTickLabel', {'-170', '', '', '', '0', '', '', '', '180'}, ...
         'YTick', [6 17 29 41], 'YTickLabel', {}, ...
         'clim', clim_)
% title('Early', 'FontSize', tfs)
% hla = xlabel('Direction tuning (degree)', 'FontSize', lfs)
% lpos = get(hla, 'position');
% lpos(1) = 220;
%set(hla, 'position', lpos)

% plot annotation arrow 
[xb,yb] = axisXY2figXY(0,1.5);
[xe,ye] = axisXY2figXY(0,0.2);
anh     = annotation('arrow', [xb xe], [yb ye]);
set(anh, 'Color', 'r', 'LineWidth', 0, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
          'HeadWidth', 1.5*ms, 'HeadLength', 1.5*ms)
[xb,yb] = axisXY2figXY(180,1.5);
[xe,ye] = axisXY2figXY(180,0.2);
anh     = annotation('arrow', [xb xe], [yb ye]);
set(anh, 'Color', 'b', 'LineWidth', 0, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
          'HeadWidth', 1.5*ms, 'HeadLength', 1.5*ms)

      
set(axes, 'Units', 'Normalized', 'Position', [0.07+2*(0.03+asize) 0.6 asize 1.3*asize*pw/ph], 'FontName', 'Helvetica')
imagesc(pdir, 1:sum(L), wmc(L,:,4))
z = wmc(L,:,4);  % transparency map for nans
z(isnan(wmc(L,:,4))) = 0;
z(~isnan(wmc(L,:,4))) = 1;
alpha(z)
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, 'Color', [0.5 0.5 0.5], ...
         'XTick', [-170 -135 -90 -45 0 45 90 135 180], 'TickDir', 'in', ...
         'XTickLabel', {'-170', '', '', '', '0', '', '', '', '180'}, ...
         'YTick', [6 17 29 41], 'YTickLabel', {}, ...
         'clim', clim_)
%title('Mid', 'FontSize', tfs)
% plot annotation arrow 
[xb,yb] = axisXY2figXY(0,1.5);
[xe,ye] = axisXY2figXY(0,0.2);
anh     = annotation('arrow', [xb xe], [yb ye]);
set(anh, 'Color', 'r', 'LineWidth', 0, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
          'HeadWidth', 1.5*ms, 'HeadLength', 1.5*ms)
[xb,yb] = axisXY2figXY(180,1.5);
[xe,ye] = axisXY2figXY(180,0.2);
anh     = annotation('arrow', [xb xe], [yb ye]);
set(anh, 'Color', 'b', 'LineWidth', 0, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
          'HeadWidth', 1.5*ms, 'HeadLength', 1.5*ms)
     
% 
% set(axes, 'Units', 'Normalized', 'Position', [0.07+3*(0.03+asize) 0.6 asize 1.3*asize*pw/ph], 'FontName', 'Helvetica')
% imagesc(pdir, 1:sum(L), wavg_c(L,:,4))
% z = wavg_c(L,:,4);  % transparency map for nans
% z(isnan(wavg_c(L,:,4))) = 0;
% z(~isnan(wavg_c(L,:,4))) = 1;
% alpha(z)
% set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, 'Color', [0.5 0.5 0.5], ...
%          'XTick', [-170 -135 -90 -45 0 45 90 135 180], 'TickDir', 'in', ...
%          'XTickLabel', {'-170', '', '', '', '0', '', '', '', '180'}, ...
%          'YTick', [6 17 29 41], 'YTickLabel', {}, ...
%          'clim', clim_)
% title('Late', 'FontSize', tfs)
% % plot annotation arrow 
% [xb,yb] = axisXY2figXY(0,1.5);
% [xe,ye] = axisXY2figXY(0,0.2);
% anh     = annotation('arrow', [xb xe], [yb ye]);
% set(anh, 'Color', 'r', 'LineWidth', 0, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
%           'HeadWidth', 1.5*ms, 'HeadLength', 1.5*ms)
% [xb,yb] = axisXY2figXY(180,1.5);
% [xe,ye] = axisXY2figXY(180,0.2);
% anh     = annotation('arrow', [xb xe], [yb ye]);
% set(anh, 'Color', 'b', 'LineWidth', 0, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
%           'HeadWidth', 1.5*ms, 'HeadLength', 1.5*ms)
% % plot annotation arrow 
% [xb,yb] = axisXY2figXY(215,45-7+1);
% [xe,ye] = axisXY2figXY(190,45-7+1);
% anh     = annotation('arrow', [xb xe], [yb ye]);
% set(anh, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
%           'HeadWidth', ms, 'HeadLength', ms)
% [xb,yb] = axisXY2figXY(215,45-18+1);
% [xe,ye] = axisXY2figXY(190,45-18+1);
% anh     = annotation('arrow', [xb xe], [yb ye]);
% set(anh, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
%           'HeadWidth', ms, 'HeadLength', ms)
% 
%       
%       
%% plot color scale
set(axes, 'Units', 'Normalized', 'Position', [0.07+3*(0.03+asize)+0.03 0.6+0.5*1.3*asize*pw/ph+0.02 0.2*asize 0.5*1.3*asize*pw/ph-0.02], 'FontName', 'Helvetica')
dat = clim_(2)*[1:-0.01:-1]';
imagesc(dat)
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, ...
         'XTick', [], 'TickDir', 'in', ...
         'XTickLabel', {}, ...
         'YTick', [1 101 201]-0.5, 'YTickLabel', {num2str(4) '0' num2str(-4) }, ...
         'clim', clim_)
text(0.5,-25, '\times10^{-3}', 'FontSize', afs)

%%
% no corr
recompute = 0;

Nc = 115985;
Nz = 77447;

TRIALS = [Nc];
bb     = 2.32;
be     = 6.53;
bw     = (6.53-2.32)/50;
THBINS = [2.^[bb:bw:be-bw]' 2.^[bb+bw:bw:be]'];
   
Monk      = 'Cy';
NSEN      = 200;
NTUNE     = 36;
ELN       = [7]*1e-8;
DIRS      = [180 0];
SIMNUM    = [81];
NORMSCALE = sqrt(0.02);
TAG       = 'NoCorr';
[wnc0, wocd0] = getML_RLFigS9OptimalWSortedByTh0(Monk, NSEN, NTUNE, ELN, DIRS, SIMNUM, THBINS, NORMSCALE, TAG, recompute);

wnc0 = circshift(wnc0,[0,-10]);


% Shadlen corr
recompute = 0;

Nc = 115985;
Nz = 77447;

TRIALS = [Nc];
bb     = 2.32;
be     = 6.53;
bw     = (6.53-2.32)/50;
THBINS = [2.^[bb:bw:be-bw]' 2.^[bb+bw:bw:be]'];
   
Monk      = 'Cy';
NSEN      = 200;
NTUNE     = 36;
ELN       = [7]*1e-8;
DIRS      = [180 0];
SIMNUM    = [81];
NORMSCALE = sqrt(0.02);
TAG       = 'ShadCorr';
[wsc0, wocd0] = getML_RLFigS9OptimalWSortedByTh0(Monk, NSEN, NTUNE, ELN, DIRS, SIMNUM, THBINS, NORMSCALE, TAG, recompute);

wsc0 = circshift(wsc0,[0,-10]);


% My corr
recompute = 0;

Nc = 115985;
Nz = 77447;

TRIALS = [Nc];
bb     = 2.32;
be     = 6.53;
bw     = (6.53-2.32)/50;
THBINS = [2.^[bb:bw:be-bw]' 2.^[bb+bw:bw:be]'];
   
Monk      = 'Cy';
NSEN      = 200;
NTUNE     = 36;
ELN       = [7]*1e-8;
DIRS      = [180 0];
SIMNUM    = [51];
NORMSCALE = sqrt(0.02);
TAG       = '';
[wmc0, wocd0] = getML_RLFigS9OptimalWSortedByTh0(Monk, NSEN, NTUNE, ELN, DIRS, SIMNUM, THBINS, NORMSCALE, TAG, recompute);

wmc0 = circshift(wmc0,[0,-10]);



% %% simulation parameters
% recompute =0;
% 
% Nc = 115985;
% Nz = 77447;
% 
% TRIALS = [0 4000 40000 Nc];
% bb     = 2.32;
% be     = 6.53;
% bw     = (6.53-2.32)/50;
% THBINS = [2.^[bb:bw:be-bw]' 2.^[bb+bw:bw:be]'];
%    
% 
% % optimal weights
% Monk      = 'Cy';
% NSEN      = 200;
% NTUNE     = 36;
% ELN       = [7]*1e-8;
% DIRS      = [180 0];
% SIMNUM    = [81];
% NORMSCALE = sqrt(0.02);
% [woc0, wocd0] = getML_RLFigS9OptimalWSortedByTh0(Monk, NSEN, NTUNE, ELN, DIRS, SIMNUM, THBINS, NORMSCALE, recompute);
% woc0 = circshift(woc0, [0,-10]);
% wocd0 = circshift(wocd0, [0,-10]);
% 
% 
%% simulation parameters
recompute =0;

Nc = 115985;
Nz = 77447;

TRIALS = [0 4000 40000 Nc];
bb     = 2.32;
be     = 6.53;
bw     = (6.53-2.32)/50;
THBINS = [2.^[bb:bw:be-bw]' 2.^[bb+bw:bw:be]'];
   

% optimal weights
Monk      = 'Cy';
NSEN      = 200;
NTUNE     = 36;
ELN       = [7]*1e-8;
DIRS      = [180 0];
SIMNUM    = [51];
NORMSCALE = sqrt(0.02);
TAG       = '';
[woavg_c, wod, wmc2, wocd1] = getML_RLFigS9OptimalWSortedByTh1(Monk, NSEN, NTUNE, ELN, DIRS, SIMNUM, THBINS, NORMSCALE, TAG, recompute);
wmc2 = circshift(wmc2, [0,-10]);
wocd1 = circshift(wocd1, [0,-10]);
% 
% %%
recompute = 0;

Nc = 115985;
Nz = 77447;

TRIALS = [0 4000 40000 Nc];
bb     = 2.32;
be     = 6.53;
bw     = (6.53-2.32)/50;
THBINS = [2.^[bb:bw:be-bw]' 2.^[bb+bw:bw:be]'];
   

% optimal weights
Monk      = 'Cy';
NSEN      = 200;
NTUNE     = 36;
ELN       = [7]*1e-8;
DIRS      = [180 0];
SIMNUM    = [81];
NORMSCALE = sqrt(0.02);
[woavg_c, wod, wsc2, wocd2] = getML_RLFigS9OptimalWSortedByTh2(Monk, NSEN, NTUNE, ELN, DIRS, SIMNUM, THBINS, NORMSCALE, TAG, recompute);
wsc2 = circshift(wsc2, [0,-10]);
wocd2 = circshift(wocd2, [0,-10]);

% %%
% recompute = 0;
% 
% Nc = 115985;
% Nz = 77447;
% 
% TRIALS = [0 4000 40000 Nc];
% bb     = 2.32;
% be     = 6.53;
% bw     = (6.53-2.32)/50;
% THBINS = [2.^[bb:bw:be-bw]' 2.^[bb+bw:bw:be]'];
%    
% 
% % optimal weights
% Monk      = 'Cy';
% NSEN      = 200;
% NTUNE     = 36;
% ELN       = [7]*1e-8;
% DIRS      = [180 0];
% SIMNUM    = [81:82];
% NORMSCALE = sqrt(0.02);
% [woavg_c, wod, woc3, wocd3] = getML_RLFigS9OptimalWSortedByTh3(Monk, NSEN, NTUNE, ELN, DIRS, SIMNUM, THBINS, NORMSCALE, recompute);
% 
% woc3 = circshift(woc3, [0,-10]);
% wocd3 = circshift(wocd3, [0,-10]);









%% plot weights (binned by thresholds)
DTUNE   = 360/NTUNE;
pdir    = linspace(-170,180,NTUNE);

% label
% th = annotation('textbox', [0 0.95 0.05 0.05]);
% set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
%     'String', 'a', 'FontSize', ffs, 'FontWeight', 'bold', 'FontName', 'Helvetica') 



clim_ = [-0.004 0.004];
thm   = nangeomean(THBINS,2);
L     = thm>6 & thm<80;


% plot weights
set(axes, 'Units', 'Normalized', 'Position', [0.07 0.33 asize 1.3*asize*pw/ph], 'FontName', 'Helvetica')
imagesc(pdir, 1:sum(L), wnc0(L,:))
z = wnc0(L,:);  % transparency map for nans
z(isnan(wnc0(L,:))) = 0;
z(~isnan(wnc0(L,:))) = 1;
alpha(z)
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, 'Color', [0.5 0.5 0.5], ...
         'XTick', [-170 -135 -90 -45 0 45 90 135 180], 'TickDir', 'in', ...
         'XTickLabel', {'-170', '', '', '', '0', '', '', '', '180'}, ...
         'YTick', [6 17 29 41], 'YTickLabel', {'8', '16', '32', '64'}, 'clim', clim_)
%xlabel('Direction tuning (degree)', 'FontSize', lfs)
ylabel('Threshold (% COH)', 'FontSize', lfs)
%title('Initial', 'FontSize', tfs)
% plot annotation arrow 
[xb,yb] = axisXY2figXY(0,1.5);
[xe,ye] = axisXY2figXY(0,0.2);
anh     = annotation('arrow', [xb xe], [yb ye]);
set(anh, 'Color', 'r', 'LineWidth', 0, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
          'HeadWidth', 1.5*ms, 'HeadLength', 1.5*ms)
[xb,yb] = axisXY2figXY(180,1.5);
[xe,ye] = axisXY2figXY(180,0.2);
anh     = annotation('arrow', [xb xe], [yb ye]);
set(anh, 'Color', 'b', 'LineWidth', 0, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
          'HeadWidth', 1.5*ms, 'HeadLength', 1.5*ms)
% [xb,yb] = axisXY2figXY(215,45-18+1);
% [xe,ye] = axisXY2figXY(190,45-18+1);
% anh     = annotation('arrow', [xb xe], [yb ye]);
% set(anh, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
%           'HeadWidth', ms, 'HeadLength', ms)



set(axes, 'Units', 'Normalized', 'Position', [0.07+1*(0.03+asize) 0.33 asize 1.3*asize*pw/ph], 'FontName', 'Helvetica')
imagesc(pdir, 1:sum(L), wsc0(L,:))
z = wsc0(L,:);  % transparency map for nans
z(isnan(wsc0(L,:))) = 0;
z(~isnan(wsc0(L,:))) = 1;
alpha(z)
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, 'Color', [0.5 0.5 0.5], ...
         'XTick', [-170 -135 -90 -45 0 45 90 135 180], 'TickDir', 'in', ...
         'XTickLabel', {'-170', '', '', '', '0', '', '', '', '180'}, ...
         'YTick', [6 17 29 41], 'YTickLabel', {}, ...
         'clim', clim_)
%title('Early', 'FontSize', tfs)
% plot annotation arrow 
[xb,yb] = axisXY2figXY(0,1.5);
[xe,ye] = axisXY2figXY(0,0.2);
anh     = annotation('arrow', [xb xe], [yb ye]);
set(anh, 'Color', 'r', 'LineWidth', 0, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
          'HeadWidth', 1.5*ms, 'HeadLength', 1.5*ms)
[xb,yb] = axisXY2figXY(180,1.5);
[xe,ye] = axisXY2figXY(180,0.2);
anh     = annotation('arrow', [xb xe], [yb ye]);
set(anh, 'Color', 'b', 'LineWidth', 0, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
          'HeadWidth', 1.5*ms, 'HeadLength', 1.5*ms)
      
set(axes, 'Units', 'Normalized', 'Position', [0.07+2*(0.03+asize) 0.33 asize 1.3*asize*pw/ph], 'FontName', 'Helvetica')
imagesc(pdir, 1:sum(L), wmc0(L,:))
z = wmc0(L,:);  % transparency map for nans
z(isnan(wmc0(L,:))) = 0;
z(~isnan(wmc0(L,:))) = 1;
alpha(z)
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, 'Color', [0.5 0.5 0.5], ...
         'XTick', [-170 -135 -90 -45 0 45 90 135 180], 'TickDir', 'in', ...
         'XTickLabel', {'-170', '', '', '', '0', '', '', '', '180'}, ...
         'YTick', [6 17 29 41], 'YTickLabel', {}, ...
         'clim', clim_)
%title('Mid', 'FontSize', tfs)
% plot annotation arrow 
[xb,yb] = axisXY2figXY(0,1.5);
[xe,ye] = axisXY2figXY(0,0.2);
anh     = annotation('arrow', [xb xe], [yb ye]);
set(anh, 'Color', 'r', 'LineWidth', 0, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
          'HeadWidth', 1.5*ms, 'HeadLength', 1.5*ms)
[xb,yb] = axisXY2figXY(180,1.5);
[xe,ye] = axisXY2figXY(180,0.2);
anh     = annotation('arrow', [xb xe], [yb ye]);
set(anh, 'Color', 'b', 'LineWidth', 0, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
          'HeadWidth', 1.5*ms, 'HeadLength', 1.5*ms)
     

% set(axes, 'Units', 'Normalized', 'Position', [0.07+3*(0.03+asize) 0.58 asize 1.3*asize*pw/ph], 'FontName', 'Helvetica')
% imagesc(pdir, 1:sum(L), woc3)
% z = woc3;  % transparency map for nans
% z(isnan(woc3)) = 0;
% z(~isnan(woc3)) = 1;
% alpha(z)
% set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, 'Color', [0.5 0.5 0.5], ...
%          'XTick', [-170 -135 -90 -45 0 45 90 135 180], 'TickDir', 'in', ...
%          'XTickLabel', {'-170', '', '', '', '0', '', '', '', '180'}, ...
%          'YTick', [6 17 29 41], 'YTickLabel', {}, ...
%          'clim', clim_)
% %title('Late', 'FontSize', tfs)
% % plot annotation arrow 
% [xb,yb] = axisXY2figXY(0,1.5);
% [xe,ye] = axisXY2figXY(0,0.2);
% anh     = annotation('arrow', [xb xe], [yb ye]);
% set(anh, 'Color', 'r', 'LineWidth', 0, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
%           'HeadWidth', 1.5*ms, 'HeadLength', 1.5*ms)
% [xb,yb] = axisXY2figXY(180,1.5);
% [xe,ye] = axisXY2figXY(180,0.2);
% anh     = annotation('arrow', [xb xe], [yb ye]);
% set(anh, 'Color', 'b', 'LineWidth', 0, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
%           'HeadWidth', 1.5*ms, 'HeadLength', 1.5*ms)
% % plot annotation arrow 
% [xb,yb] = axisXY2figXY(215,45-7+1);
% [xe,ye] = axisXY2figXY(190,45-7+1);
% anh     = annotation('arrow', [xb xe], [yb ye]);
% set(anh, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
%           'HeadWidth', ms, 'HeadLength', ms)

      
      
%% plot weights (binned by thresholds)
DTUNE   = 360/NTUNE;
pdir    = linspace(-170,180,NTUNE);

% label
% th = annotation('textbox', [0 0.95 0.05 0.05]);
% set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
%     'String', 'a', 'FontSize', ffs, 'FontWeight', 'bold', 'FontName', 'Helvetica') 



clim_ = [-0.004 0.004];
thm   = nangeomean(THBINS,2);
L     = thm>6 & thm<80;


% plot weights
set(axes, 'Units', 'Normalized', 'Position', [0.07 0.06 asize 1.3*asize*pw/ph], 'FontName', 'Helvetica')
imagesc(pdir, 1:sum(L), wnc0(L,:))
z = wnc0(L,:);  % transparency map for nans
z(isnan(wnc0(L,:))) = 0;
z(~isnan(wnc0(L,:))) = 1;
alpha(z)
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, 'Color', [0.5 0.5 0.5], ...
         'XTick', [-170 -135 -90 -45 0 45 90 135 180], 'TickDir', 'in', ...
         'XTickLabel', {'-170', '', '', '', '0', '', '', '', '180'}, ...
         'YTick', [6 17 29 41], 'YTickLabel', {'8', '16', '32', '64'}, 'clim', clim_)
%xlabel('Direction tuning (degree)', 'FontSize', lfs)
ylabel('Threshold (% COH)', 'FontSize', lfs)
%title('Initial', 'FontSize', tfs)
% plot annotation arrow 
[xb,yb] = axisXY2figXY(0,1.5);
[xe,ye] = axisXY2figXY(0,0.2);
anh     = annotation('arrow', [xb xe], [yb ye]);
set(anh, 'Color', 'r', 'LineWidth', 0, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
          'HeadWidth', 1.5*ms, 'HeadLength', 1.5*ms)
[xb,yb] = axisXY2figXY(180,1.5);
[xe,ye] = axisXY2figXY(180,0.2);
anh     = annotation('arrow', [xb xe], [yb ye]);
set(anh, 'Color', 'b', 'LineWidth', 0, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
          'HeadWidth', 1.5*ms, 'HeadLength', 1.5*ms)



set(axes, 'Units', 'Normalized', 'Position', [0.07+1*(0.03+asize) 0.06 asize 1.3*asize*pw/ph], 'FontName', 'Helvetica')
imagesc(pdir, 1:sum(L), wsc2(L,:))
z = wsc2(L,:);  % transparency map for nans
z(isnan(wsc2(L,:))) = 0;
z(~isnan(wsc2(L,:))) = 1;
alpha(z)
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, 'Color', [0.5 0.5 0.5], ...
         'XTick', [-170 -135 -90 -45 0 45 90 135 180], 'TickDir', 'in', ...
         'XTickLabel', {'-170', '', '', '', '0', '', '', '', '180'}, ...
         'YTick', [6 17 29 41], 'YTickLabel', {}, ...
         'clim', clim_)
xlabel('Direction tuning (degree)', 'FontSize', lfs)
%title('Early', 'FontSize', tfs)
% plot annotation arrow 
[xb,yb] = axisXY2figXY(0,1.5);
[xe,ye] = axisXY2figXY(0,0.2);
anh     = annotation('arrow', [xb xe], [yb ye]);
set(anh, 'Color', 'r', 'LineWidth', 0, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
          'HeadWidth', 1.5*ms, 'HeadLength', 1.5*ms)
[xb,yb] = axisXY2figXY(180,1.5);
[xe,ye] = axisXY2figXY(180,0.2);
anh     = annotation('arrow', [xb xe], [yb ye]);
set(anh, 'Color', 'b', 'LineWidth', 0, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
          'HeadWidth', 1.5*ms, 'HeadLength', 1.5*ms)

      
      
set(axes, 'Units', 'Normalized', 'Position', [0.07+2*(0.03+asize) 0.06 asize 1.3*asize*pw/ph], 'FontName', 'Helvetica')
imagesc(pdir, 1:sum(L), wmc2(L,:))
z = wmc2(L,:);  % transparency map for nans
z(isnan(wmc2(L,:))) = 0;
z(~isnan(wmc2(L,:))) = 1;
alpha(z)
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, 'Color', [0.5 0.5 0.5], ...
         'XTick', [-170 -135 -90 -45 0 45 90 135 180], 'TickDir', 'in', ...
         'XTickLabel', {'-170', '', '', '', '0', '', '', '', '180'}, ...
         'YTick', [6 17 29 41], 'YTickLabel', {}, ...
         'clim', clim_)
%xlabel('Direction tuning (degree)', 'FontSize', lfs)
%title('Mid', 'FontSize', tfs)
% plot annotation arrow 
[xb,yb] = axisXY2figXY(0,1.5);
[xe,ye] = axisXY2figXY(0,0.2);
anh     = annotation('arrow', [xb xe], [yb ye]);
set(anh, 'Color', 'r', 'LineWidth', 0, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
          'HeadWidth', 1.5*ms, 'HeadLength', 1.5*ms)
[xb,yb] = axisXY2figXY(180,1.5);
[xe,ye] = axisXY2figXY(180,0.2);
anh     = annotation('arrow', [xb xe], [yb ye]);
set(anh, 'Color', 'b', 'LineWidth', 0, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
          'HeadWidth', 1.5*ms, 'HeadLength', 1.5*ms)

      














