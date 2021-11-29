%%
% create figure
mm2inch = 25.4;
pw = 150/mm2inch;    % paper width, make it nature size of one column, 85mm x 85mm
ph = 50/mm2inch;   % paper height

% define figure parameters
alw = 0.5;  % axis line width
afs = 5;    % axis font size
lfs = 6;    % label font size
ffs = 12;   % figure font size
tfs = 10;    % text font size
lw  = 1.5;  % line width
elw = 0.5;
ms  = 2;    % marker size

% define line color
lc1  = [0 0 0];
lc1e = [0 0 0];
lc2  = [1 0.2 0.2];
lc2e = [1 0.75 0.75];
lc3  = [0.5 0.5 0.5];
lc3e = [0.5 0.5 0.5];

fh = figure; 
wysifig(fh, pw, ph) % set it to US letter size



%% get fitted values for main models in figure 2
btmk_c = 21021;
blmk_c = 1184.1;


% get data
Monk      = 'Cy';
NSEN      = 200;
NTUNE     = 36;
ELN       = [2.5 3 5 7 9]*1e-8; %for correlated response
DIRS      = [-90 90];
SIMNUM    = [81:85];
bins      = [];

% get time constants
[bt_c, btd_c, bl_c, bld_c] = getML_RLFig2Tau(Monk, NSEN, NTUNE, ELN, DIRS, SIMNUM, bins, 0);

mlF2   = bl_c(3,:);
selF2  = nanmean([bl_c(3,:)'-shiftdim(bld_c(3,1,:)) shiftdim(bld_c(3,2,:))-bl_c(3,:)'],2)';
mtF2   = bt_c(3,:);
settF2 = nanmean([bt_c(3,:)'-shiftdim(btd_c(3,1,:)) shiftdim(btd_c(3,2,:))-bt_c(3,:)'],2)';
bF2    = regressXYW(mlF2',mtF2',ones(size(selF2')),ones(size(settF2')));


%% 2 pools of LIP neurons
% plot
asize = 0.72;

% 
% plot time constants
% 

% get data
SIMNAME   = 'TwoLIPPools';
Monk      = 'Cy';
NSEN      = 200;
NTUNE     = 36;
ELN       = [1 1.3 2 3.5 5]*1e-8; %for correlated response
DIRS      = [-90 90];
SIMNUM    = [81];

N      = 115985;


% get thresholds
bb     = 0;
be     = N;
bw     = 1000;
bs     = bw;
tbins_c   = [[bb:bs:be-bw]' [bb+bw:bs:be]'];
[fits_c, fitsd_c] = getML_RLSuppDocTh(Monk, NSEN, NTUNE, ELN, DIRS, SIMNUM, tbins_c, SIMNAME, 0);

% get lapse
bb     = 0;
be     = N;
bw     = 250;
bs     = bw;
lbins_c  = [[bb:bs:be-bw]' [bb+bw:bs:be]'];
[lp_c, lpd_c, lN_c] = getML_RLSuppDocLapse(Monk, NSEN, NTUNE, ELN, DIRS, SIMNUM, lbins_c, SIMNAME, 0);

% get time constants
[bt_c, btd_c, bl_c, bld_c] = getML_RLSuppDocTau(Monk, NSEN, NTUNE, ELN, DIRS, SIMNUM, lbins_c, SIMNAME, 0);

% plot
% label
th = annotation('textbox', [0.68 0.95 0.05 0.05]);
set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'String', 'c', 'FontSize', ffs, 'FontWeight', 'bold', 'FontName', 'Helvetica') 
set(axes, 'Units', 'Normalized', 'Position', [0.75 0.15 asize*ph/pw asize], 'FontName', 'Helvetica')
Monk = 'Cy';
ml   = bl_c(3,:);
sel  = nanmean([bl_c(3,:)'-shiftdim(bld_c(3,1,:)) shiftdim(bld_c(3,2,:))-bl_c(3,:)'],2)';
mt   = bt_c(3,:);
sett = nanmean([bt_c(3,:)'-shiftdim(btd_c(3,1,:)) shiftdim(btd_c(3,2,:))-bt_c(3,:)'],2)';

hold on
plot([mlF2; mlF2], [mtF2-settF2; mtF2+settF2], '-', 'LineWidth', 0.5, 'Color', lc3e)
plot([mlF2-selF2; mlF2+selF2], [mtF2; mtF2], '-', 'LineWidth', 0.5, 'Color', lc3e)
plot(mlF2,mtF2, 'ok', 'MarkerSize', 2*ms, 'MarkerFaceColor', lc3e, 'MarkerEdgeColor', 'none')
h1=plot([-500:100:7000], bF2(1)+bF2(2)*[-500:100:7000], '-', 'Color', lc3e', 'LineWidth', lw)

plot([ml; ml], [mt-sett; mt+sett], '-', 'LineWidth', 0.5, 'Color', lc1e)
plot([ml-sel; ml+sel], [mt; mt], '-', 'LineWidth', 0.5, 'Color', lc1e)
plot(ml,mt, 'ok', 'MarkerSize', 2*ms, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none')
b = regressXYW(ml',mt',ones(size(sel')),ones(size(sett')));
h2=plot([-500:100:7000], b(1)+b(2)*[-500:100:7000], '-k', 'LineWidth', lw)
[ycb p] = regresscb(ml',mt', b, [-500:100:7000], [blmk_c btmk_c]);
plot([-500:100:7000], b(1)+b(2)*[-500:100:7000], '-k', 'LineWidth', lw)
plot([-500:100:7000]', ycb(:,1), ':', 'LineWidth', 0.5, 'Color', 'k')
plot([-500:100:7000]', ycb(:,2), ':', 'LineWidth', 0.5, 'Color', 'k')
plot(blmk_c, btmk_c, 'r*', 'MarkerSize', 3*ms)
hold off
%legend([h1 h2], 'One pool', 'Two pools', 'Location', 'Best', 'FontSize', afs, 'Boxoff')
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, ...
         'xlim', [-500 7000], 'ylim', [0 60000], 'XTick', [0 3.5 7]*1e3, 'YTick', [0 3 6]*1e4, ...
         'XTickLabel', {'0' '3.5' '7'}, 'YTickLabel', {'0' '30' '60'}) 
%xlabel('\tau, lapse (\times 1000 trials)', 'FontSize', lfs)
xlabel('\tau_{la} (\times 1000 trials)', 'FontSize', lfs)
ylabel('\tau_{th} (\times 1000 trials)', 'FontSize', lfs)
hold on
plot(4000,5000, 'r*', 'MarkerSize', 3*ms)
text(6800, 5000, 'Monkey', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'Middle', 'FontSize', afs)
plot(4000,10000, 'o', 'MarkerSize', 2*ms, 'MarkerFaceColor', lc3e, 'MarkerEdgeColor', 'none')
text(6800,10000, 'One pool', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'Middle', 'FontSize', afs)
plot(4000,15000, 'o', 'MarkerSize', 2*ms, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none')
text(6800,15000, 'Two pools', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'Middle', 'FontSize', afs)
hold off



% 
% plot performance
% 
% one pool
recompute = 0;

Monk      = 'Cy';
NSEN      = 200;
NTUNE     = 36;
ELN       = [7]*1e-8; %for correlated response
DIRS      = [-90 90];
SIMNUM    = [81:90];

[fits_c, sems_c, tbins_c, lp_c, lpd_c, lN_c, lbins_c, blmd_c, blmdd_c] = getML_RLFig2ModelAvg(Monk, NSEN, NTUNE, ELN, DIRS, SIMNUM, recompute);

th = annotation('textbox', [0 0.95 0.05 0.05]);
set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'String', 'b', 'FontSize', ffs, 'FontWeight', 'bold', 'FontName', 'Helvetica') 

% plot Cy
set(axes, 'Units', 'Normalized', 'Position', [0.08 0.15 asize*ph/pw asize], 'FontName', 'Helvetica')
Monk = 'Cy';
load([Monk 'Combined.mat'])
N    = sum(data(:,3)>=0);
init = [0.07; 0.7; 20000];
bcon = [0.05 0.12; 0.45 1; 0 100000];
Lgd = fits_c(1,:)>0.03;
[bt_c, btd_c, gof, ymt_c] = exp_fitWd(tbins_c(Lgd), fits_c(1,Lgd)', sems_c(1,Lgd)', init, bcon, {100, 68});
hold on
plot([tbins_c(Lgd) tbins_c(Lgd)]',100*[fits_c(1,Lgd)-sems_c(1,Lgd); fits_c(1,Lgd)+sems_c(1,Lgd)], '-', 'Color', lc1e, 'LineWidth', elw)
plot(tbins_c(Lgd),100*fits_c(1,Lgd),'o', 'MarkerFaceColor', lc1, 'MarkerEdgeColor', 'none', 'MarkerSize', ms)
plot(tbins_c(Lgd), 100*ymt_c, 'k', 'LineWidth', lw)
hold off
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, ...
         'xlim', [0 N], 'ylim', [5 100], 'YScale', 'log', ...
         'XTick', [0 0.5 1]*10e4, 'XTickLabel', {'0' '50' '100'}, ...  
         'YTick', [5 10 20 40 80], 'YTickLabel', {'5' '10', '20', '40', '80'}) 
xlabel('Trial (\times 1000)', 'FontSize', lfs)
ylabel('Threshold (% coh)', 'FontSize', lfs)
title('One pool', 'FontSize', ffs-1)
text(0.98*N, 68, 'C', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'FontSize', tfs)


set(axes, 'Units', 'Normalized', 'Position', [0.08 0.15 asize*ph/pw asize], 'FontName', 'Helvetica')
Monk = 'Cy';
hold on
plot([lbins_c(lN_c>15) lbins_c(lN_c>15)]',1-[lp_c(lN_c>15)-lpd_c(lN_c>15) lp_c(lN_c>15)+lpd_c(lN_c>15)]', '-', 'Color', lc2e, 'LineWidth', elw)
plot(lbins_c(lN_c>15),1-lp_c(lN_c>15), 'v', 'MarkerFaceColor', lc2, 'MarkerEdgeColor', 'none', 'MarkerSize', ms)
plot(lbins_c, blmd_c(1)+(blmd_c(2)-blmd_c(1)).*exp(-1/blmd_c(3).*lbins_c), 'r', 'LineWidth', lw)
hold off
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, ...
         'Color', 'none', 'YColor', 'r', 'YAxisLocation', 'right', ...
         'xlim', [1 N], 'ylim', [0 0.5], 'XTick', [0 0.5 1]*10e4, 'YTick', [0 0.25 0.5], ...
         'XTickLabel', {}) 
%ylabel('Lapse', 'Color', 'r')




% Two pools
% get data
SIMNAME   = 'TwoLIPPools';
Monk      = 'Cy';
NSEN      = 200;
NTUNE     = 36;
ELN       = [1 1.3 2 3.5 5]*1e-8; %for correlated response
DIRS      = [-90 90];
SIMNUM    = [81];

N      = 115985;


% get thresholds
bb     = 0;
be     = N;
bw     = 1000;
bs     = bw;
tbins_c   = [[bb:bs:be-bw]' [bb+bw:bs:be]'];
[fits_c, fitsd_c] = getML_RLSuppDocTh(Monk, NSEN, NTUNE, ELN, DIRS, SIMNUM, tbins_c, SIMNAME, 0);

% get lapse
bb     = 0;
be     = N;
bw     = 250;
bs     = bw;
lbins_c  = [[bb:bs:be-bw]' [bb+bw:bs:be]'];
[lp_c, lpd_c, lN_c] = getML_RLSuppDocLapse(Monk, NSEN, NTUNE, ELN, DIRS, SIMNUM, lbins_c, SIMNAME, 0);

% get time constants
[bt_c, btd_c, bl_c, bld_c] = getML_RLSuppDocTau(Monk, NSEN, NTUNE, ELN, DIRS, SIMNUM, lbins_c, SIMNAME, 0);


SIMNAME   = 'TwoLIPPools';
Monk      = 'Cy';
NSEN      = 200;
NTUNE     = 36;
ELN       = [1 1.3 2 3.5 5]*1e-8; %for correlated response
DIRS      = [-90 90];
SIMNUM    = [81];
N         = 115985;
ind       = 3;

fits_c    = fits_c(:,:,ind);
sems_c    = fitsd_c(:,:,ind);
lp_c      = lp_c(:,ind);
lpd_c     = lpd_c(:,ind);
lN_c      = lN_c(:,ind);
tbins_c   = nanmean(tbins_c,2);
lbins_c   = nanmean(lbins_c,2);


% plot Cy
set(axes, 'Units', 'Normalized', 'Position', [0.39 0.15 asize*ph/pw asize], 'FontName', 'Helvetica')
Monk = 'Cy';
load([Monk 'Combined.mat'])
N    = sum(data(:,3)>=0);
hold on
plot([tbins_c tbins_c]',100*[fits_c(1,:)-sems_c(1,:); fits_c(1,:)+sems_c(1,:)], '-', 'Color', lc1e, 'LineWidth', elw)
plot(tbins_c,100*fits_c(1,:),'o', 'MarkerFaceColor', lc1, 'MarkerEdgeColor', 'none', 'MarkerSize', ms)
plot(tbins_c, 100*(bt_c(1,ind)+(bt_c(2,ind)-bt_c(1,ind)).*exp(-1/bt_c(3,ind).*tbins_c)), 'k', 'LineWidth', lw)
hold off
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, ...
         'xlim', [0 N], 'ylim', [5 100], 'YScale', 'log', ...
         'XTick', [0 0.5 1]*10e4, 'XTickLabel', {'0' '50' '100'}, ...  
         'YTick', [5 10 20 40 80], 'YTickLabel', {'5' '10', '20', '40', '80'}) 
xlabel('Trial (\times 1000)', 'FontSize', lfs)
%ylabel('Threshold (% coh)', 'FontSize', lfs)
title('Two pools', 'FontSize', ffs-1)
text(0.98*N, 68, 'C', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'FontSize', tfs)


set(axes, 'Units', 'Normalized', 'Position', [0.39 0.15 asize*ph/pw asize], 'FontName', 'Helvetica')
Monk = 'Cy';
hold on
plot([lbins_c(lN_c>15) lbins_c(lN_c>15)]',1-[lp_c(lN_c>15)-lpd_c(lN_c>15) lp_c(lN_c>15)+lpd_c(lN_c>15)]', '-', 'Color', lc2e, 'LineWidth', elw)
plot(lbins_c(lN_c>15),1-lp_c(lN_c>15), 'v', 'MarkerFaceColor', lc2, 'MarkerEdgeColor', 'none', 'MarkerSize', ms)
plot(lbins_c, bl_c(1,ind)+(bl_c(2,ind)-bl_c(1,ind)).*exp(-1/bl_c(3,ind).*lbins_c), 'r', 'LineWidth', lw)
hold off
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, ...
         'Color', 'none', 'YColor', 'r', 'YAxisLocation', 'right', ...
         'xlim', [1 N], 'ylim', [0 0.5], 'XTick', [0 0.5 1]*10e4, 'YTick', [0 0.25 0.5], ...
         'XTickLabel', {}) 
ylabel('Lapse', 'Color', 'r')


%% pooling weights
% create figure
mm2inch = 25.4;
pw = 150/mm2inch;    % paper width, make it nature size of one column, 85mm x 85mm
ph = 100/mm2inch;   % paper height

% define figure parameters
alw = 0.5;  % axis line width
afs = 5;    % axis font size
lfs = 6;    % label font size
ffs = 12;   % figure font size
tfs = 10;    % text font size
lw  = 1.5;  % line width
elw = 0.5;
ms  = 2;    % marker size

% define line color
lc1  = [0 0 0];
lc1e = [0 0 0];
lc2  = [1 0.2 0.2];
lc2e = [1 0.75 0.75];
lc3  = [0.5 0.5 0.5];
lc3e = [0.5 0.5 0.5];

fh = figure; 
wysifig(fh, pw, ph) % set it to US letter size

asize = 0.16

%% get data
recompute = 0;

SIMNAME   = 'TwoLIPPools';
Monk      = 'Cy';
NSEN      = 200;
NTUNE     = 36;
ELN       = [3.5]*1e-8; %for correlated response
DIRS      = [-90 90];
SIMNUM    = [81];

Nc        = 115985;
TRIALS = [Nc];
bb     = 2.32;
be     = 6.53;
bw     = (6.53-2.32)/50;
THBINS = [2.^[bb:bw:be-bw]' 2.^[bb+bw:bw:be]'];
   
[wavg2, wse2, wpavg2, wpse2, wnavg2, wnse2] = getML_RLFigS6AvgWSortedByTh(Monk, NSEN, NTUNE, ELN, DIRS, SIMNUM, TRIALS, THBINS, 'TWOPOOLS', recompute); 


Monk      = 'Cy';
NSEN      = 200;
NTUNE     = 36;
ELN       = [7]*1e-8; %for correlated response
DIRS      = [-90 90];
SIMNUM    = [81];

[wavg1, wse1, wpavg1, wpse1, wnavg1, wnse1] = getML_RLFigS6AvgWSortedByTh(Monk, NSEN, NTUNE, ELN, DIRS, SIMNUM, TRIALS, THBINS, 'ONEPOOL', recompute); 
    
% % shift weights
wavg2  = circshift(wavg2,[0,-1]);
wpavg2 = circshift(wpavg2,[0,-1]);
wnavg2 = circshift(wnavg2,[0,-1]);
wavg1  = circshift(wavg1,[0,-1]);

DTUNE   = 360/NTUNE;
pdir    = linspace(-170,180,NTUNE);

%% plot
% th = annotation('textbox', [0.42 0.95 0.05 0.05]);
% set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
%     'String', 'e', 'FontSize', ffs, 'FontWeight', 'bold', 'FontName', 'Helvetica') 

clim_ = [-0.008 0.008];
thm   = nangeomean(THBINS,2);
DTUNE = 360/NTUNE;
pdir  = linspace(-170,180,NTUNE);
L     = thm>6 & thm<80;

% plot weights     
set(axes, 'Units', 'Normalized', 'Position', [0.2+1*(0.03+asize) 0.58 asize 1.3*asize*pw/ph], 'FontName', 'Helvetica')

imagesc(pdir, 1:sum(L), wpavg2(L,:))
z = wpavg2(L,:);  % transparency map for nans
z(isnan(wpavg2(L,:))) = 0;
z(~isnan(wpavg2(L,:))) = 1;
alpha(z)
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, 'Color', [0.5 0.5 0.5], ...
         'XTick', [-170 -135 -90 -45 0 45 90 135 180], 'TickDir', 'in', ...
         'XTickLabel', {'-170', '', '', '', '0', '', '', '', '180'}, ...
         'YTick', [6 17 29 41], 'YTickLabel', {'8', '16', '32', '64'}, ...
         'clim', clim_)
%title('LIP_R', 'FontSize', tfs)
% plot annotation arrow 
[xb,yb] = axisXY2figXY(0,-1.5);
[xe,ye] = axisXY2figXY(0,0);
anh     = annotation('arrow', [xb xe], [yb ye]);
set(anh, 'Color', 'r', 'LineWidth', 0, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
          'HeadWidth', 1.5*ms, 'HeadLength', 1.5*ms)
[xb,yb] = axisXY2figXY(-180,-1.5);
[xe,ye] = axisXY2figXY(-180,0);
anh     = annotation('arrow', [xb xe], [yb ye]);
set(anh, 'Color', 'b', 'LineWidth', 0, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
          'HeadWidth', 1.5*ms, 'HeadLength', 1.5*ms)
xlabel('Direction tuning (degree)', 'FontSize', lfs)
ylabel('Threshold (% COH)', 'FontSize', lfs)


      

      
set(axes, 'Units', 'Normalized', 'Position', [0.2+2*(0.03+asize) 0.58 asize 1.3*asize*pw/ph], 'FontName', 'Helvetica')
imagesc(pdir, 1:sum(L), wnavg2(L,:))
z = wnavg2(L,:);  % transparency map for nans
z(isnan(wnavg2(L,:))) = 0;
z(~isnan(wnavg2(L,:))) = 1;
alpha(z)
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, 'Color', [0.5 0.5 0.5], ...
         'XTick', [-170 -135 -90 -45 0 45 90 135 180], 'TickDir', 'in', ...
         'XTickLabel', {'-170', '', '', '', '0', '', '', '', '180'}, ...
         'YTick', [6 17 29 41], 'YTickLabel', {}, ...
         'clim', clim_)
%title('LIP_R', 'FontSize', tfs)
% plot annotation arrow 
[xb,yb] = axisXY2figXY(0,-1.5);
[xe,ye] = axisXY2figXY(0,0);
anh     = annotation('arrow', [xb xe], [yb ye]);
set(anh, 'Color', 'r', 'LineWidth', 0, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
          'HeadWidth', 1.5*ms, 'HeadLength', 1.5*ms)
[xb,yb] = axisXY2figXY(-180,-1.5);
[xe,ye] = axisXY2figXY(-180,0);
anh     = annotation('arrow', [xb xe], [yb ye]);
set(anh, 'Color', 'b', 'LineWidth', 0, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
          'HeadWidth', 1.5*ms, 'HeadLength', 1.5*ms)
%xlabel('Direction tuning (degree)', 'FontSize', lfs)
%ylabel('Threshold (% COH)', 'FontSize', lfs)

      
set(axes, 'Units', 'Normalized', 'Position', [0.2+3*(0.03+asize) 0.58 asize 1.3*asize*pw/ph], 'FontName', 'Helvetica')
imagesc(pdir, 1:sum(L), wavg2(L,:))
z = wavg2(L,:);  % transparency map for nans
z(isnan(wavg2(L,:))) = 0;
z(~isnan(wavg2(L,:))) = 1;
alpha(z)
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, 'Color', [0.5 0.5 0.5], ...
         'XTick', [-170 -135 -90 -45 0 45 90 135 180], 'TickDir', 'in', ...
         'XTickLabel', {'-170', '', '', '', '0', '', '', '', '180'}, ...
         'YTick', [6 17 29 41], 'YTickLabel', {}, 'clim', clim_)
% xlabel('Direction tuning (degree)', 'FontSize', lfs)
% ylabel('Threshold (% COH)', 'FontSize', lfs)
%title('LIP_{R}-LIP_{L}', 'FontSize', tfs)
% plot annotation arrow 
[xb,yb] = axisXY2figXY(0,-1.5);
[xe,ye] = axisXY2figXY(0,0);
anh     = annotation('arrow', [xb xe], [yb ye]);
set(anh, 'Color', 'r', 'LineWidth', 0, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
          'HeadWidth', 1.5*ms, 'HeadLength', 1.5*ms)
[xb,yb] = axisXY2figXY(-180,-1.5);
[xe,ye] = axisXY2figXY(-180,0);
anh     = annotation('arrow', [xb xe], [yb ye]);
set(anh, 'Color', 'b', 'LineWidth', 0, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
          'HeadWidth', 1.5*ms, 'HeadLength', 1.5*ms)

      
% plot color scale
set(axes, 'Units', 'Normalized', 'Position', [0.17+4*(0.03+asize)+0.03 0.58+0.5*1.3*asize*pw/ph+0.02 0.2*asize 0.5*1.3*asize*pw/ph-0.02], 'FontName', 'Helvetica')
dat = clim_(2)*[1:-0.01:-1]';
imagesc(dat)
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, ...
         'XTick', [], 'TickDir', 'in', ...
         'XTickLabel', {}, ...
         'YTick', [1 101 201]-0.5, 'YTickLabel', {num2str(8) '0' num2str(-8) }, ...
         'clim', clim_)
text(0.5,-25, '\times10^{-3}', 'FontSize', afs)

      
 

% th = annotation('textbox', [0 0.7 0.05 0.05]);
% set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
%     'String', 'd', 'FontSize', ffs, 'FontWeight', 'bold', 'FontName', 'Helvetica') 

    
clim_ = [-0.004 0.004];
set(axes, 'Units', 'Normalized', 'Position', [0.07 0.58 asize 1.3*asize*pw/ph], 'FontName', 'Helvetica')
imagesc(pdir, 1:sum(L), wavg1(L,:))
z = wavg1(L,:);  % transparency map for nans
z(isnan(wavg1(L,:))) = 0;
z(~isnan(wavg1(L,:))) = 1;
alpha(z)
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, 'Color', [0.5 0.5 0.5], ...
         'XTickLabel', {'-170', '', '', '', '0', '', '', '', '180'}, ...
         'YTick', [6 17 29 41], 'YTickLabel', {'8', '16', '32', '64'}, ...
         'clim', clim_)
%title('Late', 'FontSize', tfs)
% plot annotation arrow 
[xb,yb] = axisXY2figXY(0,-1.5);
[xe,ye] = axisXY2figXY(0,0);
anh     = annotation('arrow', [xb xe], [yb ye]);
set(anh, 'Color', 'r', 'LineWidth', 0, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
          'HeadWidth', 1.5*ms, 'HeadLength', 1.5*ms)
[xb,yb] = axisXY2figXY(-180,-1.5);
[xe,ye] = axisXY2figXY(-180,0);
anh     = annotation('arrow', [xb xe], [yb ye]);
set(anh, 'Color', 'b', 'LineWidth', 0, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
          'HeadWidth', 1.5*ms, 'HeadLength', 1.5*ms)
xlabel('Direction tuning (degree)', 'FontSize', lfs)
ylabel('Threshold (% COH)', 'FontSize', lfs)

      
% plot color scale
set(axes, 'Units', 'Normalized', 'Position', [0.038+1*(0.03+asize)+0.03 0.58+0.5*1.3*asize*pw/ph+0.02 0.2*asize 0.5*1.3*asize*pw/ph-0.02], 'FontName', 'Helvetica')
dat = clim_(2)*[1:-0.01:-1]';
imagesc(dat)
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, ...
         'XTick', [], 'TickDir', 'in', ...
         'XTickLabel', {}, ...
         'YTick', [1 101 201]-0.5, 'YTickLabel', {num2str(4) '0' num2str(-4) }, ...
         'clim', clim_)
text(0.5,-25, '\times10^{-3}', 'FontSize', afs)

       
      

%% get LIP SNR
% get average SNR of pooled signal for both monkeys
% create figure
fh = figure; 
mm2inch = 25.4;
pw = 80/mm2inch;    % paper width, make it nature size of one column, 85mm x 85mm
ph = 80/mm2inch;   % paper height
wysifig(fh, pw, ph) % set it to US letter size

% define figure parameters
alw = 0.5;  % axis line width
afs = 8;    % axis font size
lfs = 8;    % label font size
ffs = 12;   % figure font size
tfs = 12;    % text font size
lw  = 1.5;  % line width
elw = 0.75;
ms  = 5;    % marker size

% define line color
lc1  = [0.2 0.2 0.2];
lc1e = [0.75 0.75 0.75];
lc2  = [1 0.2 0.2];
lc2e = [1 0.75 0.75];


asize = 0.7;


recompute = 0;

% model
Monk      = 'Cy';
NSEN      = 200;
NTUNE     = 36;
ELN       = [3.5]*1e-8;
DIRS      = [0 180];
SIMNUM    = [81];

COH  = 0.512;
TIME = 0.35;
[m1_c2, sd1_c2] = getML_RLFigS6PooledSNR(Monk, NSEN, NTUNE, ELN, SIMNUM, COH, TIME, recompute);

COH = 0.128;
[m2_c2, sd2_c2] = getML_RLFigS6PooledSNR(Monk, NSEN, NTUNE, ELN, SIMNUM, COH, TIME, recompute);

recompute = 0;

% model
Monk      = 'Cy';
NSEN      = 200;
NTUNE     = 36;
ELN       = [7]*1e-8;
DIRS      = [0 180];
SIMNUM    = [81];

COH  = 0.512;
TIME = 0.35;
[m1_c, sd1_c] = getML_RLFig5PooledSNR(Monk, NSEN, NTUNE, ELN, SIMNUM, COH, TIME, recompute);


COH = 0.128;
[m2_c, sd2_c] = getML_RLFig5PooledSNR(Monk, NSEN, NTUNE, ELN, SIMNUM, COH, TIME, recompute);

% data
fname     = 'CyTRain_LIP.txt';
crtf      = [0 1];
nvf       = [0 1];
dirf      = [0];
bb        = 200;
be        = 500;
bw        = 300;
bs        = bw;
bins      = [[bb:bs:be-bw]' [bb+bw:bs:be]'];
[LSNR_c, LTR_c] = getML_RLFig5LIPSNR(fname, crtf, nvf, dirf, bins, 0);
lsnr1_c   = squeeze(LSNR_c(:,6,:));
lsnr2_c   = squeeze(LSNR_c(:,4,:));

%% plot SNR of pooled signal for coh 1
Nc = 115985;
Nz = 77447;


set(axes, 'Units', 'Normalized', 'Position', [0.15 0.15 asize asize*pw/ph], 'FontName', 'Helvetica')
hold on
tr   = [1:1000:Nc Nc];
snr  = m1_c./(sd1_c+sqrt(mean(m1_c))+5);
h1 = plot(tr, snr, 'o', 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', 'none', 'MarkerSize', ms)
%plot(LTR_c, lsnr1_c, 'o', 'MarkerFaceColor', lc2, 'MarkerEdgeColor', 'none', 'MarkerSize', ms)

tr   = [1:1000:Nc];
snr  = m1_c2./(sd1_c2+sqrt(mean(m1_c2))+5);
h2 = plot(tr, snr, 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerSize', ms)
plot(LTR_c, lsnr1_c, 'o', 'MarkerFaceColor', lc2, 'MarkerEdgeColor', 'none', 'MarkerSize', ms)
hold off
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, ...
         'xlim', [0 Nc], 'ylim', [-0.2 3], 'YTick', [0 1.5 3],  ...
         'XTick', [0 0.5 1]*10e4, 'XTickLabel', {'0' '50' '100'}) 
%xlabel('Trial (\times 1000)', 'FontSize', lfs)
%ylabel('SNR ratio', 'FontSize', lfs)
ylabel('SNR', 'FontSize', lfs)
%title('Monkey C', 'FontSize', lfs)
%text(0.03*Nc,3, '51.2% coh', 'FontSize', lfs, 'HorizontalAlignment' , 'left', 'VerticalAlignment', 'top')

legend([h1 h2], 'One pool', 'Two pools', 'Location', 'NorthEast', 'FontSize', afs)
legend('boxoff')
