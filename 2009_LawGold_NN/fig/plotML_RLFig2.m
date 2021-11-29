% plot reinforcement learning figure 2 - behavior

% Matlab script for plotting Fig 2  (Comparison of behavioral performance of monkey and model)
% for submission to Nature Neuroscience
%
% created by jcl on 08/07/08

%%
% create figure
fh = figure; 
mm2inch = 25.4;
pw = 160/mm2inch;    % paper width, make it nature size of one column, 85mm x 85mm
ph = 85/mm2inch;   % paper height
wysifig(fh, pw, ph) % set it to US letter size

% define figure parameters
alw = 0.5;  % axis line width
afs = 7;    % axis font size
lfs = 7;    % label font size
ffs = 12;   % figure font size
tfs = 10;    % text font size
lw  = 1.5;  % line width
elw = 0.5;
ms  = 2;    % marker size

% define line color
lc1  = [0.2 0.2 0.2];
lc1e = [0.75 0.75 0.75];
lc2  = [1 0.2 0.2];
lc2e = [1 0.75 0.75];


asize = 0.38;



%% Figure 2a: Thresholds and lapses for monkeys as a function of training
% get fits
recompute = 0;
[fits_c, sems_c, tbins_c, lp_c, lpd_c, lN_c, lbins_c] = getML_RLFig2MonkAvg('Cy', recompute);
[fits_z, sems_z, tbins_z, lp_z, lpd_z, lN_z, lbins_z] = getML_RLFig2MonkAvg('ZZ', recompute);

% label
th = annotation('textbox', [0 0.95 0.05 0.05]);
set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'String', 'a', 'FontSize', ffs, 'FontWeight', 'bold', 'FontName', 'Helvetica') 


% plot Cy
set(axes, 'Units', 'Normalized', 'Position', [0.068 0.55 asize*ph/pw asize], 'FontName', 'Helvetica')
Monk = 'Cy';
load([Monk 'Combined.mat'])
N    = sum(data(:,3)>=0);
init = [0.08; 0.7; 10000];
bcon = [0.06 0.2; 0.50 1; 0 100000];
Lgd  = fits_c(1,:)>0.03;
err  = nanmean([fits_c(1,:)-shiftdim(sems_c(1,1,:),1); shiftdim(sems_c(1,2,:),1)-fits_c(1,:)]);
err(err==0) = nanmean(err);
[btmk_c, btmkd_c, gof, ymt_c] = exp_fitWd(tbins_c(Lgd), fits_c(1,Lgd)', err(Lgd)', init, bcon, {100,68});
hold on
%plot([tbins_c(Lgd) tbins_c(Lgd)]',100*[fits_c(1,Lgd)-sems_c(1,Lgd); fits_c(1,Lgd)+sems_c(1,Lgd)], '-', 'Color', lc1e, 'LineWidth', elw)
plot([tbins_c(Lgd) tbins_c(Lgd)]',100*shiftdim(sems_c(1,:,:),1), '-', 'Color', lc1e, 'LineWidth', elw)
plot(tbins_c(Lgd),100*fits_c(1,Lgd),'o', 'MarkerFaceColor', lc1, 'MarkerEdgeColor', 'none', 'MarkerSize', ms)
plot(tbins_c(Lgd), 100*ymt_c, 'k', 'LineWidth', 2)
hold off
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, ...
         'xlim', [0 N], 'ylim', [5 100], 'YScale', 'log', ...
         'XTick', [0 0.5 1]*10e4, 'XTickLabel', {'0' '50' '100'}, ...  
         'YTick', [5 10 20 40 80], 'YTickLabel', {'5' '10', '20', '40', '80'}) 
%xlabel('Trial (\times 1000)', 'FontSize', lfs)
ylabel('Threshold (% coh)', 'FontSize', lfs)
title('Monkey', 'FontSize', ffs-1)
text(0.98*N, 68, 'C', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'FontSize', tfs)


set(axes, 'Units', 'Normalized', 'Position', [0.068 0.55 asize*ph/pw asize], 'FontName', 'Helvetica')
Monk = 'Cy';
% load([Monk 'Combined.mat'])
% crt  = data(:,3);
% Lgd  = crt>=0;
% crt  = crt(Lgd);
% ddir = data(Lgd,4);
% coh  = data(Lgd,5);
% vt   = data(Lgd,6);
% N    = length(crt);
% L = coh/100>0.9 & vt/1000>1;
% [blmk_c, bld_c, ym_c] = exp_fitd_bino(find(L), 1-crt(L), [], [0 0; 0.5 0.5; 0 N]);
hold on
plot([lbins_c lbins_c]',1-[lp_c-lpd_c lp_c+lpd_c]', '-', 'Color', lc2e, 'LineWidth', elw)
plot(lbins_c,1-lp_c, 'v', 'MarkerFaceColor', lc2, 'MarkerEdgeColor', 'none', 'MarkerSize', ms)
init = [0; 0.5; 2000];
bcon = [0 0; 0.5 0.5; 100 5000];
Lgd  = ~isnan(lp_c);
lpd_c(lpd_c==0) = mynanmedian(lpd_c(lpd_c~=0));
[blmk_c, blmkd_c, gof, yml_c] = exp_fitWd(lbins_c(Lgd), 1-lp_c(Lgd), lpd_c(Lgd), init, bcon, {100,68});
plot(lbins_c, blmk_c(1)+(blmk_c(2)-blmk_c(1)).*exp(-1/blmk_c(3).*lbins_c), 'r', 'LineWidth', 2)
hold off
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, ...
         'Color', 'none', 'YColor', 'r', 'YAxisLocation', 'right', ...
         'xlim', [1 N], 'ylim', [0 0.5], 'XTick', [0 0.5 1]*10e4, 'YTick', [0 0.25 0.5], ...
         'XTickLabel', {}) 
%ylabel('Lapse', 'Color', 'r')





% plot them together
Monk = 'ZZ';
load([Monk 'Combined.mat'])
N    = sum(data(:,3)>=0);
set(axes, 'Units', 'Normalized', 'Position', [0.068 0.1 asize*ph/pw asize], 'FontName', 'Helvetica')
init = [0.08; 1; 80000];
bcon = [0.1 0.25; 0.8 1; 0 200000];
Lgd = fits_z(1,:)>0.03;
err  = nanmean([fits_z(1,:)-shiftdim(sems_z(1,1,:),1); shiftdim(sems_z(1,2,:),1)-fits_z(1,:)]);
err(err==0) = nanmean(err);
err(1:5)    = nanmean(err);
[btmk_z, btmkd_z, gof, ymt_z] = exp_fitWd(tbins_z(Lgd), fits_z(1,Lgd)', err(Lgd)', init, bcon, {100,68});
hold on
%plot([tbins_z(Lgd) tbins_z(Lgd)]',100*[fits_z(1,Lgd)-sems_z(1,Lgd); fits_z(1,Lgd)+sems_z(1,Lgd)], '-', 'Color', lc1e, 'LineWidth', elw)
plot([tbins_z(Lgd) tbins_z(Lgd)]',100*shiftdim(sems_z(1,:,:),1), '-', 'Color', lc1e, 'LineWidth', elw)
plot(tbins_z(Lgd),100*fits_z(1,Lgd),'o', 'MarkerFaceColor', lc1, 'MarkerEdgeColor', 'none', 'MarkerSize', ms)
plot(tbins_z(Lgd), 100*ymt_z, 'k', 'LineWidth', 2)
hold off
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, ...
         'xlim', [0 N], 'ylim', [5 100], 'YScale', 'log', ...
         'XTick', [0 0.35 0.7]*10e4, 'XTickLabel', {'0' '35' '70'}, ...  
         'YTick', [5 10 20 40 80], 'YTickLabel', {'5' '10', '20', '40', '80'}) 
xlabel('Trial (\times 1000)', 'FontSize', lfs)
ylabel('Threshold (% coh)', 'FontSize', lfs)
text(0.98*N, 68, 'Z', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'FontSize', tfs)


set(axes, 'Units', 'Normalized', 'Position', [0.068 0.1 asize*ph/pw asize], 'FontName', 'Helvetica')
Monk = 'ZZ';
% load([Monk 'Combined.mat'])
% crt  = data(:,3);
% Lgd  = crt>=0;
% crt  = crt(Lgd);
% ddir = data(Lgd,4);
% coh  = data(Lgd,5);
% vt   = data(Lgd,6);
% N    = length(crt);
% L = coh/100>0.9 & vt/1000>1;
% [blmk_z, bld_z, ym_z] = exp_fitd_bino(find(L), 1-crt(L), [], [0 0; 0.5 0.5; 0 N]);
hold on
plot([lbins_z lbins_z]',1-[lp_z-lpd_z lp_z+lpd_z]', '-', 'Color', lc2e, 'LineWidth', elw)
plot(lbins_z,1-lp_z, 'v', 'MarkerFaceColor', lc2, 'MarkerEdgeColor', 'none', 'MarkerSize', ms)
init = [0; 0.5; 1000];
bcon = [0 0; 0.5 0.5; 100 10000];
Lgd  = ~isnan(lp_z);
lpd_z(lpd_z==0) = nanmean(lpd_z(lpd_z~=0));
[blmk_z, blmkd_z, gof, yml_z] = exp_fitWd(lbins_z(Lgd), 1-lp_z(Lgd), lpd_z(Lgd), init, bcon, {100,68});
plot(lbins_z, blmk_z(1)+(blmk_z(2)-blmk_z(1)).*exp(-1/blmk_z(3).*lbins_z), 'r', 'LineWidth', 2)
hold off
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, ...
         'Color', 'none', 'YColor', 'r', 'YAxisLocation', 'right', ...
         'xlim', [1 N], 'ylim', [0 0.5], 'XTick', [0 0.35 0.7]*10e4, 'YTick', [0 0.25 0.5], ...
         'XTickLabel', {}) 
%ylabel('Lapse', 'Color', 'r')









%% Figure 2b: Thresholds and lapses for model with best-matched learning rate as a function of training
recompute = 0;

Monk      = 'Cy';
NSEN      = 200;
NTUNE     = 36;
ELN       = [7]*1e-8; %for correlated response
DIRS      = [-90 90];
SIMNUM    = [81:90];

[fits_c, sems_c, tbins_c, lp_c, lpd_c, lN_c, lbins_c, blmd_c, blmdd_c] = getML_RLFig2ModelAvg(Monk, NSEN, NTUNE, ELN, DIRS, SIMNUM, recompute);


% label
th = annotation('textbox', [0.325 0.95 0.05 0.05]);
set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'String', 'b', 'FontSize', ffs, 'FontWeight', 'bold', 'FontName', 'Helvetica') 


% plot Cy
set(axes, 'Units', 'Normalized', 'Position', [0.38 0.55 asize*ph/pw asize], 'FontName', 'Helvetica')
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
%xlabel('Trial (\times 1000)', 'FontSize', lfs)
%ylabel('Threshold (% coh)', 'FontSize', lfs)
title('Model', 'FontSize', ffs-1)
text(0.98*N, 68, 'C', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'FontSize', tfs)


set(axes, 'Units', 'Normalized', 'Position', [0.38 0.55 asize*ph/pw asize], 'FontName', 'Helvetica')
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
ylabel('Lapse', 'Color', 'r')

% plot early mid late arrows
tr      = [4000 40000 N];
[xb,yb] = axisXY2figXY(tr(1),0.45);
[xe,ye] = axisXY2figXY(tr(1),0.40);
anh     = annotation('arrow', [xb xe], [yb ye]);
set(anh, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
          'HeadWidth', 1.5*ms, 'HeadLength', 1.5*ms)
[xb,yb] = axisXY2figXY(tr(2),0.35);
[xe,ye] = axisXY2figXY(tr(2),0.30);
anh     = annotation('arrow', [xb xe], [yb ye]);
set(anh, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
          'HeadWidth', 1.5*ms, 'HeadLength', 1.5*ms)
[xb,yb] = axisXY2figXY(tr(3),0.25);
[xe,ye] = axisXY2figXY(tr(3),0.20);
anh     = annotation('arrow', [xb xe], [yb ye]);
set(anh, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
          'HeadWidth', 1.5*ms, 'HeadLength', 1.5*ms)





Monk      = 'ZZ';
NSEN      = 200;
NTUNE     = 36;
ELN       = [1]*1e-8; %for correlated response
DIRS      = [-90 90];
SIMNUM    = [81:90];

[fits_z, sems_z, tbins_z, lp_z, lpd_z, lN_z, lbins_z, blmd_z, blmdd_z] = getML_RLFig2ModelAvg(Monk, NSEN, NTUNE, ELN, DIRS, SIMNUM, recompute);


% plot them together
Monk = 'ZZ';
load([Monk 'Combined.mat'])
N    = sum(data(:,3)>=0);
set(axes, 'Units', 'Normalized', 'Position', [0.38 0.1 asize*ph/pw asize], 'FontName', 'Helvetica')
init = [0.2; 0.9; 15000];
bcon = [0.18 0.21; 0.75 1; 0 100000];
Lgd = fits_z(1,:)>0.03;
[bt_z, btd_z, gof, ymt_z] = exp_fitWd(tbins_z(Lgd), fits_z(1,Lgd)', sems_z(1,Lgd)', init, bcon, {100, 68});
hold on
plot([tbins_z(Lgd) tbins_z(Lgd)]',100*[fits_z(1,Lgd)-sems_z(1,Lgd); fits_z(1,Lgd)+sems_z(1,Lgd)], '-', 'Color', lc1e, 'LineWidth', elw)
plot(tbins_z(Lgd),100*fits_z(1,Lgd),'o', 'MarkerFaceColor', lc1, 'MarkerEdgeColor', 'none', 'MarkerSize', ms)
plot(tbins_z(Lgd), 100*ymt_z, 'k', 'LineWidth', lw)
hold off
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, ...
         'xlim', [0 N], 'ylim', [5 100], 'YScale', 'log', ...
         'XTick', [0 0.35 0.7]*10e4, 'XTickLabel', {'0' '35' '70'}, ...  
         'YTick', [5 10 20 40 80], 'YTickLabel', {'5' '10', '20', '40', '80'}) 
xlabel('Trial (\times 1000)', 'FontSize', lfs)
%ylabel('Threshold (% coh)', 'FontSize', lfs)
text(0.98*N, 68, 'Z', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'FontSize', tfs)


set(axes, 'Units', 'Normalized', 'Position', [0.38 0.1 asize*ph/pw asize], 'FontName', 'Helvetica')
Monk = 'ZZ';
hold on
plot([lbins_z(lN_z>15) lbins_z(lN_z>15)]',1-[lp_z(lN_z>15)-lpd_z(lN_z>15) lp_z(lN_z>15)+lpd_z(lN_z>15)]', '-', 'Color', lc2e, 'LineWidth', elw)
plot(lbins_z(lN_z>15),1-lp_z(lN_z>15), 'v', 'MarkerFaceColor', lc2, 'MarkerEdgeColor', 'none', 'MarkerSize', ms)
plot(lbins_z, blmd_z(1)+(blmd_z(2)-blmd_z(1)).*exp(-1/blmd_z(3).*lbins_z), 'r', 'LineWidth', lw)
hold off
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, ...
         'Color', 'none', 'YColor', 'r', 'YAxisLocation', 'right', ...
         'xlim', [1 N], 'ylim', [0 0.5], 'XTick', [0 0.35 0.7]*10e4, 'YTick', [0 0.25 0.5], ...
         'XTickLabel', {}) 
ylabel('Lapse', 'Color', 'r')




%% Figure 2c: Time constants for thresholds and lapses for different learning rate

% get data
Monk      = 'Cy';
NSEN      = 200;
NTUNE     = 36;
ELN       = [2.5 3 5 7 9]*1e-8; %for correlated response
DIRS      = [-90 90];
SIMNUM    = [81:85];


[hdir, ldir, cdir, tdir] = dirnames;
fname = ['/getML_RL_' Monk '_' int2str(NSEN) '_' int2str(NTUNE) '_' num2str(ELN(1)) '_v' num2str(SIMNUM(1)) '.mat'];
load([tdir fname])
fprintf([fname '\n']);
N      = length(FIRA);

% get thresholds
bb     = 0;
be     = N;
bw     = 1000;
bs     = bw;
bins   = [[bb:bs:be-bw]' [bb+bw:bs:be]'];
[fits_c, fitsd_c] = getML_RLFig2Th(Monk, NSEN, NTUNE, ELN, DIRS, SIMNUM, bins, 0);

% get lapse
bb     = 0;
be     = N;
bw     = 250;
bs     = bw;
bins  = [[bb:bs:be-bw]' [bb+bw:bs:be]'];
[lp_c, lpd_c, lN_c] = getML_RLFig2Lapse(Monk, NSEN, NTUNE, ELN, DIRS, SIMNUM, bins, 0);

% get time constants
[bt_c, btd_c, bl_c, bld_c] = getML_RLFig2Tau(Monk, NSEN, NTUNE, ELN, DIRS, SIMNUM, bins, 0);


Monk      = 'ZZ';
NSEN      = 200;
NTUNE     = 36;
ELN       = [0.5 0.75 1 1.5 2]*1e-8; %for correlated response
DIRS      = [-90 90];
SIMNUM    = [81:85];

[hdir, ldir, cdir, tdir] = dirnames;
fname = ['/getML_RL_' Monk '_' int2str(NSEN) '_' int2str(NTUNE) '_' num2str(ELN(1)) '_v' num2str(SIMNUM(1)) '.mat'];
load([tdir fname])
fprintf([fname '\n']);
N      = length(FIRA);

% get thresholds
bb     = 0;
be     = N;
bw     = 1000;
bs     = bw;
bins   = [[bb:bs:be-bw]' [bb+bw:bs:be]'];
[fits_z, fitsd_z] = getML_RLFig2Th(Monk, NSEN, NTUNE, ELN, DIRS, SIMNUM, bins, 0);

% get lapse
bb     = 0;
be     = N;
bw     = 250;
bs     = bw;
bins  = [[bb:bs:be-bw]' [bb+bw:bs:be]'];
[lp_z, lpd_z, lN_z] = getML_RLFig2Lapse(Monk, NSEN, NTUNE, ELN, DIRS, SIMNUM, bins, 0);

% get time constants
[bt_z, btd_z, bl_z, bld_z] = getML_RLFig2Tau(Monk, NSEN, NTUNE, ELN, DIRS, SIMNUM, bins, 0);




% plot
th = annotation('textbox', [0.68 0.95 0.05 0.05]);
set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'String', 'c', 'FontSize', ffs, 'FontWeight', 'bold', 'FontName', 'Helvetica') 


set(axes, 'Units', 'Normalized', 'Position', [0.76 0.55 asize*ph/pw asize], 'FontName', 'Helvetica')
Monk = 'Cy';
ml   = bl_c(3,:);
sel  = nanmean([bl_c(3,:)'-shiftdim(bld_c(3,1,:)) shiftdim(bld_c(3,2,:))-bl_c(3,:)'],2)';
mt   = bt_c(3,:);
sett = nanmean([bt_c(3,:)'-shiftdim(btd_c(3,1,:)) shiftdim(btd_c(3,2,:))-bt_c(3,:)'],2)';

hold on
plot([ml; ml], [mt-sett; mt+sett], '-', 'LineWidth', 0.5, 'Color', lc1e)
plot([ml-sel; ml+sel], [mt; mt], '-', 'LineWidth', 0.5, 'Color', lc1e)
h1 = plot(ml,mt, 'ok', 'MarkerSize', 2*ms, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none')
b       = regressXYW(ml',mt',ones(size(sel')),ones(size(sett')));
[ycb p] = regresscb(ml',mt', b, [-500:100:7000], [blmk_c(3) btmk_c(3)]);
plot([-500:100:7000], b(1)+b(2)*[-500:100:7000], '-k', 'LineWidth', lw)
plot([-500:100:7000]', ycb(:,1), ':', 'LineWidth', 0.5, 'Color', 'k')
plot([-500:100:7000]', ycb(:,2), ':', 'LineWidth', 0.5, 'Color', 'k')
h2 = plot(blmk_c(3), btmk_c(3), 'r*', 'MarkerSize', 3*ms)
hold off
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, ...
         'xlim', [-500 7000], 'ylim', [0 60000], 'XTick', [0 3.5 7]*1e3, 'YTick', [0 3 6]*1e4, ...
         'XTickLabel', {'0' '3.5' '7'}, 'YTickLabel', {'0' '30' '60'}) 
%xlabel('\tau, lapse (\times 1000 trials)', 'FontSize', lfs)
ylabel('\tau_{th} (\times 1000 trials)', 'FontSize', lfs)
text(-500+0.05*7500, 0.99*60000, 'C', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', tfs)
%text(-500+0.98*7500, 100, ['p=' sprintf('%.3f',p)], 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'FontSize', tfs-1)
hold on
plot(4000,5000, 'r*', 'MarkerSize', 3*ms)
text(6800, 5000, 'Monkey', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'Middle', 'FontSize', afs)
plot(4000,10000, 'ok', 'MarkerSize', 2*ms, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none')
text(6800,10000, 'Model', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'Middle', 'FontSize', afs)
hold off

set(axes, 'Units', 'Normalized', 'Position', [0.76 0.1 asize*ph/pw asize], 'FontName', 'Helvetica')
Monk = 'ZZ';
ml   = bl_z(3,:);
sel  = nanmean([bl_z(3,:)'-shiftdim(bld_z(3,1,:)) shiftdim(bld_z(3,2,:))-bl_z(3,:)'],2)';
mt   = bt_z(3,:);
sett = nanmean([bt_z(3,:)'-shiftdim(btd_z(3,1,:)) shiftdim(btd_z(3,2,:))-bt_z(3,:)'],2)';

hold on
plot([ml; ml], [mt-sett; mt+sett], '-', 'LineWidth', 0.5, 'Color', lc1e)
plot([ml-sel; ml+sel], [mt; mt], '-', 'LineWidth', 0.5, 'Color', lc1e)
plot(ml,mt, 'ok', 'MarkerSize', 2*ms, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none')
b = regressXYW(ml',mt',ones(size(sel')),ones(size(sett')));
[ycb p] = regresscb(ml',mt', b, [0:100:15000], [blmk_c(3) btmk_c(3)]);
plot([0:100:15000], b(1)+b(2)*[0:100:15000], '-k', 'LineWidth', lw)
plot([0:100:15000]', ycb(:,1), ':', 'LineWidth', 0.5, 'Color', 'k')
plot([0:100:15000]', ycb(:,2), ':', 'LineWidth', 0.5, 'Color', 'k')
plot(blmk_z(3), btmk_z(3), 'r*', 'MarkerSize', 3*ms)
hold off
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, ...
         'xlim', [0 15000], 'ylim', [0 50000], 'XTick', [0 0.75 1.5]*1e4, 'YTick', [0 2.5 5]*1e4, ...
         'XTickLabel', {'0' '7.5' '15'}, 'YTickLabel', {'0' '25' '50'}) 
xlabel('\tau_{la} (\times 1000 trials)', 'FontSize', lfs)
ylabel('\tau_{th} (\times 1000 trials)', 'FontSize', lfs)
text(0.05*15000, 0.99*50000, 'Z', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', tfs)
%text(0.98*15000, 100, ['p=' sprintf('%.3f',p)], 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'FontSize', tfs-1)







