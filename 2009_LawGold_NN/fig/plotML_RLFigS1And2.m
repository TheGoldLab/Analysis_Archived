% plot figures for supplementary documents
%%
% create figure
mm2inch = 25.4;
pw = 60/mm2inch;    % paper width, make it nature size of one column, 85mm x 85mm
ph = 60/mm2inch;   % paper height

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


asize = 0.6;



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
fh = figure; 
wysifig(fh, pw, ph) % set it to US letter size


% 
% plot time constants
% 

% get data
SIMNAME   = 'TwoLIPPools';
Monk      = 'Cy';
NSEN      = 200;
NTUNE     = 36;
ELN       = [1 2 3.5 5]*1e-8; %for correlated response
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
th = annotation('textbox', [0.45 0.9 0.05 0.05]);
set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'String', 'b', 'FontSize', ffs, 'FontWeight', 'bold', 'FontName', 'Helvetica') 

set(axes, 'Units', 'Normalized', 'Position', [0.58 0.2 asize*ph/pw asize], 'FontName', 'Helvetica')
Monk = 'Cy';
ml   = bl_c(3,:);
sel  = nanmean([bl_c(3,:)'-shiftdim(bld_c(3,1,:)) shiftdim(bld_c(3,2,:))-bl_c(3,:)'],2)';
mt   = bt_c(3,:);
sett = nanmean([bt_c(3,:)'-shiftdim(btd_c(3,1,:)) shiftdim(btd_c(3,2,:))-bt_c(3,:)'],2)';

hold on
plot([mlF2; mlF2], [mtF2-settF2; mtF2+settF2], '-', 'LineWidth', 0.5, 'Color', lc3e)
plot([mlF2-selF2; mlF2+selF2], [mtF2; mtF2], '-', 'LineWidth', 0.5, 'Color', lc3e)
plot(mlF2,mtF2, 'ok', 'MarkerSize', 2*ms, 'MarkerFaceColor', lc3e, 'MarkerEdgeColor', 'none')
plot([-500:100:7000], bF2(1)+bF2(2)*[-500:100:7000], '-', 'Color', lc3e', 'LineWidth', lw)

plot([ml; ml], [mt-sett; mt+sett], '-', 'LineWidth', 0.5, 'Color', lc1e)
plot([ml-sel; ml+sel], [mt; mt], '-', 'LineWidth', 0.5, 'Color', lc1e)
plot(ml,mt, 'ok', 'MarkerSize', 2*ms, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none')
%b = regressXYW(ml',mt',ones(size(sel')),ones(size(sett')));
b = regressXYW(ml',mt',ones(size(sel')),ones(size(sett')));
plot([-500:100:7000], b(1)+b(2)*[-500:100:7000], '-k', 'LineWidth', lw)
plot(blmk_c, btmk_c, 'r*', 'MarkerSize', 3*ms)
hold off
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, ...
         'xlim', [-500 7000], 'ylim', [0 60000], 'XTick', [0 3.5 7]*1e3, 'YTick', [0 3 6]*1e4, ...
         'XTickLabel', {'0' '3.5' '7'}, 'YTickLabel', {'0' '30' '60'}) 
%xlabel('\tau, lapse (\times 1000 trials)', 'FontSize', lfs)
xlabel('\tau_{la} (\times 1000 trials)', 'FontSize', lfs)
ylabel('\tau_{th} (\times 1000 trials)', 'FontSize', lfs)






% 
% plot time constants
% 
% get data
SIMNAME   = 'TwoLIPPools';
Monk      = 'Cy';
NSEN      = 200;
NTUNE     = 36;
ELN       = [1 2 3.5 5]*1e-8; %for correlated response
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

% label
th = annotation('textbox', [0 0.9 0.05 0.05]);
set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'String', 'a', 'FontSize', ffs, 'FontWeight', 'bold', 'FontName', 'Helvetica') 


% plot Cy
set(axes, 'Units', 'Normalized', 'Position', [0.1 0.2 asize*ph/pw asize], 'FontName', 'Helvetica')
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
%xlabel('Trial (\times 1000)', 'FontSize', lfs)
%ylabel('Threshold (% coh)', 'FontSize', lfs)
title('Model', 'FontSize', ffs-1)
text(0.98*N, 68, 'C', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'FontSize', tfs)


set(axes, 'Units', 'Normalized', 'Position', [0.1 0.2 asize*ph/pw asize], 'FontName', 'Helvetica')
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


%% Additive noise
% plot
fh = figure; 
wysifig(fh, pw, ph) % set it to US letter size


% 
% plot time constants
% 

% get data
SIMNAME   = 'AddNoise';
Monk      = 'Cy';
NSEN      = 200;
NTUNE     = 36;
ELN       = [1.25 1.5 2 4 5]*1e-8; %for correlated response (NOTE: add 2 more sim, 2 and 4)
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
% th = annotation('textbox', [0.45 0.9 0.05 0.05]);
% set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
%     'String', 'b', 'FontSize', ffs, 'FontWeight', 'bold', 'FontName', 'Helvetica') 

set(axes, 'Units', 'Normalized', 'Position', [0.2 0.2 asize*ph/pw asize], 'FontName', 'Helvetica')
Monk = 'Cy';
ml   = bl_c(3,:);
sel  = nanmean([bl_c(3,:)'-shiftdim(bld_c(3,1,:)) shiftdim(bld_c(3,2,:))-bl_c(3,:)'],2)';
mt   = bt_c(3,:);
sett = nanmean([bt_c(3,:)'-shiftdim(btd_c(3,1,:)) shiftdim(btd_c(3,2,:))-bt_c(3,:)'],2)';

hold on
plot([mlF2; mlF2], [mtF2-settF2; mtF2+settF2], '-', 'LineWidth', 0.5, 'Color', lc3e)
plot([mlF2-selF2; mlF2+selF2], [mtF2; mtF2], '-', 'LineWidth', 0.5, 'Color', lc3e)
plot(mlF2,mtF2, 'ok', 'MarkerSize', 2*ms, 'MarkerFaceColor', lc3e, 'MarkerEdgeColor', 'none')
plot([-500:100:7000], bF2(1)+bF2(2)*[-500:100:7000], '-', 'Color', lc3e', 'LineWidth', lw)

plot([ml; ml], [mt-sett; mt+sett], '-', 'LineWidth', 0.5, 'Color', lc1e)
plot([ml-sel; ml+sel], [mt; mt], '-', 'LineWidth', 0.5, 'Color', lc1e)
plot(ml,mt, 'ok', 'MarkerSize', 2*ms, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none')
b = regressXYW(ml',mt',ones(size(sel')),ones(size(sett')));
[ycb p] = regresscb(ml',mt', b, [-500:100:7000], [blmk_c btmk_c]);
plot([-500:100:7000], b(1)+b(2)*[-500:100:7000], '-k', 'LineWidth', lw)
plot([-500:100:7000]', ycb(:,1), ':', 'LineWidth', 0.5, 'Color', 'k')
plot([-500:100:7000]', ycb(:,2), ':', 'LineWidth', 0.5, 'Color', 'k')
plot(blmk_c, btmk_c, 'r*', 'MarkerSize', 3*ms)
hold off
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, ...
         'xlim', [-500 7000], 'ylim', [0 60000], 'XTick', [0 3.5 7]*1e3, 'YTick', [0 3 6]*1e4, ...
         'XTickLabel', {'0' '3.5' '7'}, 'YTickLabel', {'0' '30' '60'}) 
%xlabel('\tau, lapse (\times 1000 trials)', 'FontSize', lfs)
xlabel('\tau_{la} (\times 1000 trials)', 'FontSize', lfs)
ylabel('\tau_{th} (\times 1000 trials)', 'FontSize', lfs)



%% Multiplicative noise
% plot
fh = figure; 
wysifig(fh, pw, ph) % set it to US letter size


% 
% plot time constants
% 

% get data
SIMNAME   = 'MulNoise';
Monk      = 'Cy';
NSEN      = 200;
NTUNE     = 36;
ELN       = [0.5 0.7 1 3]*1e-8; %for correlated response 
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
[bt_c, btd_c, bl_c, bld_c] = getML_RLSuppDocTau(Monk, NSEN, NTUNE, ELN, DIRS, SIMNUM, lbins_c, SIMNAME,0);


% plot
% label
% th = annotation('textbox', [0.45 0.9 0.05 0.05]);
% set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
%     'String', 'b', 'FontSize', ffs, 'FontWeight', 'bold', 'FontName', 'Helvetica') 

set(axes, 'Units', 'Normalized', 'Position', [0.2 0.2 asize*ph/pw asize], 'FontName', 'Helvetica')
Monk = 'Cy';
ml   = bl_c(3,:);
sel  = nanmean([bl_c(3,:)'-shiftdim(bld_c(3,1,:)) shiftdim(bld_c(3,2,:))-bl_c(3,:)'],2)';
mt   = bt_c(3,:);
sett = nanmean([bt_c(3,:)'-shiftdim(btd_c(3,1,:)) shiftdim(btd_c(3,2,:))-bt_c(3,:)'],2)';

hold on
plot([mlF2; mlF2], [mtF2-settF2; mtF2+settF2], '-', 'LineWidth', 0.5, 'Color', lc3e)
plot([mlF2-selF2; mlF2+selF2], [mtF2; mtF2], '-', 'LineWidth', 0.5, 'Color', lc3e)
plot(mlF2,mtF2, 'ok', 'MarkerSize', 2*ms, 'MarkerFaceColor', lc3e, 'MarkerEdgeColor', 'none')
plot([-500:100:7000], bF2(1)+bF2(2)*[-500:100:7000], '-', 'Color', lc3e', 'LineWidth', lw)

plot([ml; ml], [mt-sett; mt+sett], '-', 'LineWidth', 0.5, 'Color', lc1e)
plot([ml-sel; ml+sel], [mt; mt], '-', 'LineWidth', 0.5, 'Color', lc1e)
plot(ml,mt, 'ok', 'MarkerSize', 2*ms, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none')
b = regressXYW(ml',mt',ones(size(sel')),ones(size(sett')));
[ycb p] = regresscb(ml',mt', b, [-500:100:7000], [blmk_c btmk_c]);
plot([-500:100:7000], b(1)+b(2)*[-500:100:7000], '-k', 'LineWidth', lw)
plot([-500:100:7000]', ycb(:,1), ':', 'LineWidth', 0.5, 'Color', 'k')
plot([-500:100:7000]', ycb(:,2), ':', 'LineWidth', 0.5, 'Color', 'k')
plot(blmk_c, btmk_c, 'r*', 'MarkerSize', 3*ms)
hold off
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, ...
         'xlim', [-500 7000], 'ylim', [0 60000], 'XTick', [0 3.5 7]*1e3, 'YTick', [0 3 6]*1e4, ...
         'XTickLabel', {'0' '3.5' '7'}, 'YTickLabel', {'0' '30' '60'}) 
%xlabel('\tau, lapse (\times 1000 trials)', 'FontSize', lfs)
xlabel('\tau_{la} (\times 1000 trials)', 'FontSize', lfs)
ylabel('\tau_{th} (\times 1000 trials)', 'FontSize', lfs)




%% No Normalization
% plot
fh = figure; 
wysifig(fh, pw, ph) % set it to US letter size


% 
% plot time constants
% 

% get data
SIMNAME   = 'NoNorm';
Monk      = 'Cy';
NSEN      = 200;
NTUNE     = 36;
ELN       = [1.5 2.5 3 5 7]*1e-8; %for correlated response 
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
% th = annotation('textbox', [0.45 0.9 0.05 0.05]);
% set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
%     'String', 'b', 'FontSize', ffs, 'FontWeight', 'bold', 'FontName', 'Helvetica') 

set(axes, 'Units', 'Normalized', 'Position', [0.2 0.2 asize*ph/pw asize], 'FontName', 'Helvetica')
Monk = 'Cy';
ml   = bl_c(3,:);
sel  = nanmean([bl_c(3,:)'-shiftdim(bld_c(3,1,:)) shiftdim(bld_c(3,2,:))-bl_c(3,:)'],2)';
mt   = bt_c(3,:);
sett = nanmean([bt_c(3,:)'-shiftdim(btd_c(3,1,:)) shiftdim(btd_c(3,2,:))-bt_c(3,:)'],2)';

hold on
plot([mlF2; mlF2], [mtF2-settF2; mtF2+settF2], '-', 'LineWidth', 0.5, 'Color', lc3e)
plot([mlF2-selF2; mlF2+selF2], [mtF2; mtF2], '-', 'LineWidth', 0.5, 'Color', lc3e)
plot(mlF2,mtF2, 'ok', 'MarkerSize', 2*ms, 'MarkerFaceColor', lc3e, 'MarkerEdgeColor', 'none')
plot([-500:100:16000], bF2(1)+bF2(2)*[-500:100:16000], '-', 'Color', lc3e', 'LineWidth', lw)

plot([ml; ml], [mt-sett; mt+sett], '-', 'LineWidth', 0.5, 'Color', lc1e)
plot([ml-sel; ml+sel], [mt; mt], '-', 'LineWidth', 0.5, 'Color', lc1e)
plot(ml,mt, 'ok', 'MarkerSize', 2*ms, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none')
b = regressXYW(ml',mt',ones(size(sel')),ones(size(sett')));
[ycb p] = regresscb(ml',mt', b, [-500:100:7000], [blmk_c btmk_c]);
plot([-500:100:7000], b(1)+b(2)*[-500:100:7000], '-k', 'LineWidth', lw)
plot([-500:100:7000]', ycb(:,1), ':', 'LineWidth', 0.5, 'Color', 'k')
plot([-500:100:7000]', ycb(:,2), ':', 'LineWidth', 0.5, 'Color', 'k')
plot(blmk_c, btmk_c, 'r*', 'MarkerSize', 3*ms)
hold off
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, ...
         'xlim', [-500 7000], 'ylim', [0 60000], 'XTick', [0 3.5 7]*1e3, 'YTick', [0 3 6]*1e4, ...
         'XTickLabel', {'0' '3.5' '7'}, 'YTickLabel', {'0' '30' '60'}) 
%xlabel('\tau, lapse (\times 1000 trials)', 'FontSize', lfs)
xlabel('\tau_{la} (\times 1000 trials)', 'FontSize', lfs)
ylabel('\tau_{th} (\times 1000 trials)', 'FontSize', lfs)








%% Subtractive Normalization
% plot
fh = figure; 
wysifig(fh, pw, ph) % set it to US letter size


% 
% plot time constants
% 

% get data
SIMNAME   = 'SubNorm';
Monk      = 'Cy';
NSEN      = 200;
NTUNE     = 36;
ELN       = [1.5 2.5 3 5 7]*1e-8; %for correlated response 
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
[bt_c, btd_c, bl_c, bld_c] = getML_RLSuppDocTau(Monk, NSEN, NTUNE, ELN, DIRS, SIMNUM, lbins_c, SIMNAME,0);


% plot
% label
% th = annotation('textbox', [0.45 0.9 0.05 0.05]);
% set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
%     'String', 'b', 'FontSize', ffs, 'FontWeight', 'bold', 'FontName', 'Helvetica') 

set(axes, 'Units', 'Normalized', 'Position', [0.2 0.2 asize*ph/pw asize], 'FontName', 'Helvetica')
Monk = 'Cy';
ml   = bl_c(3,:);
sel  = nanmean([bl_c(3,:)'-shiftdim(bld_c(3,1,:)) shiftdim(bld_c(3,2,:))-bl_c(3,:)'],2)';
mt   = bt_c(3,:);
sett = nanmean([bt_c(3,:)'-shiftdim(btd_c(3,1,:)) shiftdim(btd_c(3,2,:))-bt_c(3,:)'],2)';

hold on
plot([mlF2; mlF2], [mtF2-settF2; mtF2+settF2], '-', 'LineWidth', 0.5, 'Color', lc3e)
plot([mlF2-selF2; mlF2+selF2], [mtF2; mtF2], '-', 'LineWidth', 0.5, 'Color', lc3e)
plot(mlF2,mtF2, 'ok', 'MarkerSize', 2*ms, 'MarkerFaceColor', lc3e, 'MarkerEdgeColor', 'none')
plot([-500:100:11000], bF2(1)+bF2(2)*[-500:100:11000], '-', 'Color', lc3e', 'LineWidth', lw)

plot([ml; ml], [mt-sett; mt+sett], '-', 'LineWidth', 0.5, 'Color', lc1e)
plot([ml-sel; ml+sel], [mt; mt], '-', 'LineWidth', 0.5, 'Color', lc1e)
plot(ml,mt, 'ok', 'MarkerSize', 2*ms, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none')
b = regressXYW(ml',mt',ones(size(sel')),ones(size(sett')));
[ycb p] = regresscb(ml',mt', b, [-500:100:11000], [blmk_c btmk_c]);
plot([-500:100:11000], b(1)+b(2)*[-500:100:11000], '-k', 'LineWidth', lw)
plot([-500:100:11000]', ycb(:,1), ':', 'LineWidth', 0.5, 'Color', 'k')
plot([-500:100:11000]', ycb(:,2), ':', 'LineWidth', 0.5, 'Color', 'k')
plot(blmk_c, btmk_c, 'r*', 'MarkerSize', 3*ms)
hold off
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, ...
         'xlim', [-500 11000], 'ylim', [0 70000], 'XTick', [0 5 10]*1e3, 'YTick', [0 3.5 7]*1e4, ...
         'XTickLabel', {'0' '5' '10'}, 'YTickLabel', {'0' '35' '70'}) 
%xlabel('\tau, lapse (\times 1000 trials)', 'FontSize', lfs)
xlabel('\tau_{la} (\times 1000 trials)', 'FontSize', lfs)
ylabel('\tau_{th} (\times 1000 trials)', 'FontSize', lfs)







%% Fitting sigmodal function
% plot
fh = figure; 
wysifig(fh, pw, ph) % set it to US letter size


% 
% plot time constants
% 

% get data
SIMNAME   = 'FitSig';
Monk      = 'Cy';
NSEN      = 200;
NTUNE     = 36;
ELN       = [1 1.8 3 5]*1e-8; %for correlated response 
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
[bt_c, btd_c, bl_c, bld_c] = getML_RLSuppDocTau(Monk, NSEN, NTUNE, ELN, DIRS, SIMNUM, lbins_c, SIMNAME,0);


% plot
% label
% th = annotation('textbox', [0.45 0.9 0.05 0.05]);
% set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
%     'String', 'b', 'FontSize', ffs, 'FontWeight', 'bold', 'FontName', 'Helvetica') 

set(axes, 'Units', 'Normalized', 'Position', [0.2 0.2 asize*ph/pw asize], 'FontName', 'Helvetica')
Monk = 'Cy';
ml   = bl_c(3,:);
sel  = nanmean([bl_c(3,:)'-shiftdim(bld_c(3,1,:)) shiftdim(bld_c(3,2,:))-bl_c(3,:)'],2)';
mt   = bt_c(3,:);
sett = nanmean([bt_c(3,:)'-shiftdim(btd_c(3,1,:)) shiftdim(btd_c(3,2,:))-bt_c(3,:)'],2)';

hold on
plot([mlF2; mlF2], [mtF2-settF2; mtF2+settF2], '-', 'LineWidth', 0.5, 'Color', lc3e)
plot([mlF2-selF2; mlF2+selF2], [mtF2; mtF2], '-', 'LineWidth', 0.5, 'Color', lc3e)
plot(mlF2,mtF2, 'ok', 'MarkerSize', 2*ms, 'MarkerFaceColor', lc3e, 'MarkerEdgeColor', 'none')
plot([-500:100:7000], bF2(1)+bF2(2)*[-500:100:7000], '-', 'Color', lc3e', 'LineWidth', lw)

plot([ml; ml], [mt-sett; mt+sett], '-', 'LineWidth', 0.5, 'Color', lc1e)
plot([ml-sel; ml+sel], [mt; mt], '-', 'LineWidth', 0.5, 'Color', lc1e)
plot(ml,mt, 'ok', 'MarkerSize', 2*ms, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none')
b = regressXYW(ml',mt',ones(size(sel')),ones(size(sett')));
[ycb p] = regresscb(ml',mt', b, [-500:100:7000], [blmk_c btmk_c]);
plot([-500:100:7000], b(1)+b(2)*[-500:100:7000], '-k', 'LineWidth', lw)
plot([-500:100:7000]', ycb(:,1), ':', 'LineWidth', 0.5, 'Color', 'k')
plot([-500:100:7000]', ycb(:,2), ':', 'LineWidth', 0.5, 'Color', 'k')
plot(blmk_c, btmk_c, 'r*', 'MarkerSize', 3*ms)
hold off
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, ...
         'xlim', [-500 7000], 'ylim', [0 70000], 'XTick', [0 3.5 7]*1e3, 'YTick', [0 3 6]*1e4, ...
         'XTickLabel', {'0' '3.5' '7'}, 'YTickLabel', {'0' '30' '60'}) 
%xlabel('\tau, lapse (\times 1000 trials)', 'FontSize', lfs)
xlabel('\tau_{la} (\times 1000 trials)', 'FontSize', lfs)
ylabel('\tau_{th} (\times 1000 trials)', 'FontSize', lfs)











%% Oja learning rule
% plot
fh = figure; 
wysifig(fh, pw, ph) % set it to US letter size


% 
% plot time constants
% 

% get data
SIMNAME   = 'Oja';
Monk      = 'Cy';
NSEN      = 200;
NTUNE     = 36;
ELN       = [0.3 0.4 0.5 0.7 1]*1e-8; %for correlated response 
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
[bt_c, btd_c, bl_c, bld_c] = getML_RLSuppDocTau(Monk, NSEN, NTUNE, ELN, DIRS, SIMNUM, lbins_c, SIMNAME,0);


% plot
% label
% th = annotation('textbox', [0.45 0.9 0.05 0.05]);
% set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
%     'String', 'b', 'FontSize', ffs, 'FontWeight', 'bold', 'FontName', 'Helvetica') 

set(axes, 'Units', 'Normalized', 'Position', [0.2 0.2 asize*ph/pw asize], 'FontName', 'Helvetica')
Monk = 'Cy';
ml   = bl_c(3,:);
sel  = nanmean([bl_c(3,:)'-shiftdim(bld_c(3,1,:)) shiftdim(bld_c(3,2,:))-bl_c(3,:)'],2)';
mt   = bt_c(3,:);
sett = nanmean([bt_c(3,:)'-shiftdim(btd_c(3,1,:)) shiftdim(btd_c(3,2,:))-bt_c(3,:)'],2)';

hold on
plot([mlF2; mlF2], [mtF2-settF2; mtF2+settF2], '-', 'LineWidth', 0.5, 'Color', lc3e)
plot([mlF2-selF2; mlF2+selF2], [mtF2; mtF2], '-', 'LineWidth', 0.5, 'Color', lc3e)
plot(mlF2,mtF2, 'ok', 'MarkerSize', 2*ms, 'MarkerFaceColor', lc3e, 'MarkerEdgeColor', 'none')
plot([-500:100:9000], bF2(1)+bF2(2)*[-500:100:9000], '-', 'Color', lc3e', 'LineWidth', lw)

plot([ml; ml], [mt-sett; mt+sett], '-', 'LineWidth', 0.5, 'Color', lc1e)
plot([ml-sel; ml+sel], [mt; mt], '-', 'LineWidth', 0.5, 'Color', lc1e)
plot(ml,mt, 'ok', 'MarkerSize', 2*ms, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none')
b = regressXYW(ml',mt',ones(size(sel')),ones(size(sett')));
[ycb p] = regresscb(ml',mt', b, [-500:100:9000], [blmk_c btmk_c]);
plot([-500:100:9000], b(1)+b(2)*[-500:100:9000], '-k', 'LineWidth', lw)
plot([-500:100:9000]', ycb(:,1), ':', 'LineWidth', 0.5, 'Color', 'k')
plot([-500:100:9000]', ycb(:,2), ':', 'LineWidth', 0.5, 'Color', 'k')
plot(blmk_c, btmk_c, 'r*', 'MarkerSize', 3*ms)
hold off
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, ...
         'xlim', [-500 9000], 'ylim', [0 70000], 'XTick', [0 4.5 9]*1e3, 'YTick', [0 3.5 7]*1e4, ...
         'XTickLabel', {'0' '4.5' '9'}, 'YTickLabel', {'0' '35' '70'}) 
%xlabel('\tau, lapse (\times 1000 trials)', 'FontSize', lfs)
xlabel('\tau_{la} (\times 1000 trials)', 'FontSize', lfs)
ylabel('\tau_{th} (\times 1000 trials)', 'FontSize', lfs)











%% Nonlinear pooling
% plot
fh = figure; 
wysifig(fh, pw, ph) % set it to US letter size


% 
% plot time constants
% 

% get data
SIMNAME   = 'NonLinPool';
Monk      = 'Cy';
NSEN      = 200;
NTUNE     = 36;
ELN       = [0.5 0.55 0.6 0.8; 0.2 0.5 0.6 0.7; 0.3 0.5 0.7 2]*1e-8; %for correlated response 
NL        = [1.19; 1.41; 2];
DIRS      = [-90 90];
SIMNUM    = [81];

N         = 115985;


btall      = cell(1,length(NL));
btdall     = cell(1,length(NL));
blall      = cell(1,length(NL));
bldall     = cell(1,length(NL));


for i = 1:length(NL)
    % get thresholds
    bb     = 0;
    be     = N;
    bw     = 1000;
    bs     = bw;
    tbins_c   = [[bb:bs:be-bw]' [bb+bw:bs:be]'];
    [fits_c, fitsd_c] = getML_RLSuppDocTh(Monk, NSEN, NTUNE, ELN(i,:), DIRS, SIMNUM, tbins_c, SIMNAME, 0, NL(i));

    % get lapse
    bb     = 0;
    be     = N;
    bw     = 250;
    bs     = bw;
    lbins_c  = [[bb:bs:be-bw]' [bb+bw:bs:be]'];
    [lp_c, lpd_c, lN_c] = getML_RLSuppDocLapse(Monk, NSEN, NTUNE, ELN(i,:), DIRS, SIMNUM, lbins_c, SIMNAME, 0, NL(i));

    % get time constants
    [btall{i}, btdall{i}, blall{i}, bldall{i}] = getML_RLSuppDocTau(Monk, NSEN, NTUNE, ELN(i,:), DIRS, SIMNUM, lbins_c, SIMNAME, 0, NL(i));
end


% plot
% label
% th = annotation('textbox', [0.45 0.9 0.05 0.05]);
% set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
%     'String', 'b', 'FontSize', ffs, 'FontWeight', 'bold', 'FontName',
%     'Helvetica') 

lc = [0.8 0.8 0.8; 0.55 0.55 0.55; 0.3 0.3 0.3; 0.05 0.05 0.05];

set(axes, 'Units', 'Normalized', 'Position', [0.2 0.2 asize*ph/pw asize], 'FontName', 'Helvetica')
hold on
plot([mlF2; mlF2], [mtF2-settF2; mtF2+settF2], '-', 'LineWidth', 0.5, 'Color', lc(1,:))
plot([mlF2-selF2; mlF2+selF2], [mtF2; mtF2], '-', 'LineWidth', 0.5, 'Color', lc(1,:))
plot(mlF2,mtF2, 'ok', 'MarkerSize', 2*ms, 'MarkerFaceColor', lc(1,:), 'MarkerEdgeColor', 'none')
plot([-500:100:25000], bF2(1)+bF2(2)*[-500:100:25000], '-', 'Color', lc(1,:)', 'LineWidth', lw)
hold off

hold on
for i = 1:length(NL)
    Monk = 'Cy';
    bt_c  = btall{i};
    btd_c = btdall{i};
    bl_c  = blall{i};
    bld_c = bldall{i};
    ml   = bl_c(3,:);
    sel  = nanmean([bl_c(3,:)'-shiftdim(bld_c(3,1,:)) shiftdim(bld_c(3,2,:))-bl_c(3,:)'],2)';
    mt   = bt_c(3,:);
    sett = nanmean([bt_c(3,:)'-shiftdim(btd_c(3,1,:)) shiftdim(btd_c(3,2,:))-bt_c(3,:)'],2)';
   
    plot([ml; ml], [mt-sett; mt+sett], '-', 'LineWidth', 0.5, 'Color', lc(i+1,:))
    plot([ml-sel; ml+sel], [mt; mt], '-', 'LineWidth', 0.5, 'Color', lc(i+1,:))
    plot(ml,mt, 'ok', 'MarkerSize', 2*ms, 'MarkerFaceColor', lc(i+1,:), 'MarkerEdgeColor', 'none')
    b = regressXYW(ml',mt',ones(size(sel')),ones(size(sett')));
    [ycb p] = regresscb(ml',mt', b, [-500:100:25000], [blmk_c btmk_c]);
    plot([-500:100:25000], b(1)+b(2)*[-500:100:25000], '-', 'LineWidth', lw, 'Color', lc(i+1,:))
    % plot([-500:100:25000]', ycb(:,1), ':', 'LineWidth', 0.5, 'Color', lc(i+1,:))
    % plot([-500:100:25000]', ycb(:,2), ':', 'LineWidth', 0.5, 'Color', lc(i+1,:))
    plot(blmk_c, btmk_c, 'r*', 'MarkerSize', 3*ms)
end
hold off

set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, ...
         'xlim', [-500 25000], 'ylim', [0 80000], 'XTick', [0 12.5 25]*1e3, 'YTick', [0 4 8]*1e4, ...
         'XTickLabel', {'0' '12.5' '25'}, 'YTickLabel', {'0' '40' '80'}) 
%xlabel('\tau, lapse (\times 1000 trials)', 'FontSize', lfs)
xlabel('\tau_{la} (\times 1000 trials)', 'FontSize', lfs)
ylabel('\tau_{th} (\times 1000 trials)', 'FontSize', lfs)













%% Generalized RL learning models (11)
% plot
fh = figure; 
wysifig(fh, pw, ph) % set it to US letter size


% 
% plot time constants
% 

% get data
SIMNAME   = 'GenRL11';
Monk      = 'Cy';
NSEN      = 200;
NTUNE     = 36;
ELN       = [2 3 5 7 9]*1e-8; %for correlated response 
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
[bt_c, btd_c, bl_c, bld_c] = getML_RLSuppDocTau(Monk, NSEN, NTUNE, ELN, DIRS, SIMNUM, lbins_c, SIMNAME,0);


% plot
% label
% th = annotation('textbox', [0.45 0.9 0.05 0.05]);
% set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
%     'String', 'b', 'FontSize', ffs, 'FontWeight', 'bold', 'FontName', 'Helvetica') 

set(axes, 'Units', 'Normalized', 'Position', [0.2 0.2 asize*ph/pw asize], 'FontName', 'Helvetica')
Monk = 'Cy';
ml   = bl_c(3,:);
sel  = nanmean([bl_c(3,:)'-shiftdim(bld_c(3,1,:)) shiftdim(bld_c(3,2,:))-bl_c(3,:)'],2)';
mt   = bt_c(3,:);
sett = nanmean([bt_c(3,:)'-shiftdim(btd_c(3,1,:)) shiftdim(btd_c(3,2,:))-bt_c(3,:)'],2)';

hold on
plot([mlF2; mlF2], [mtF2-settF2; mtF2+settF2], '-', 'LineWidth', 0.5, 'Color', lc3e)
plot([mlF2-selF2; mlF2+selF2], [mtF2; mtF2], '-', 'LineWidth', 0.5, 'Color', lc3e)
plot(mlF2,mtF2, 'ok', 'MarkerSize', 2*ms, 'MarkerFaceColor', lc3e, 'MarkerEdgeColor', 'none')
plot([-500:100:7000], bF2(1)+bF2(2)*[-500:100:7000], '-', 'Color', lc3e', 'LineWidth', lw)

plot([ml; ml], [mt-sett; mt+sett], '-', 'LineWidth', 0.5, 'Color', lc1e)
plot([ml-sel; ml+sel], [mt; mt], '-', 'LineWidth', 0.5, 'Color', lc1e)
plot(ml,mt, 'ok', 'MarkerSize', 2*ms, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none')
b = regressXYW(ml',mt',ones(size(sel')),ones(size(sett')));
[ycb p] = regresscb(ml',mt', b, [-500:100:7000], [blmk_c btmk_c]);
plot([-500:100:7000], b(1)+b(2)*[-500:100:7000], '-k', 'LineWidth', lw)
plot([-500:100:7000]', ycb(:,1), ':', 'LineWidth', 0.5, 'Color', 'k')
plot([-500:100:7000]', ycb(:,2), ':', 'LineWidth', 0.5, 'Color', 'k')
plot(blmk_c, btmk_c, 'r*', 'MarkerSize', 3*ms)
hold off
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, ...
         'xlim', [-500 9000], 'ylim', [0 70000], 'XTick', [0 4.5 9]*1e3, 'YTick', [0 3.5 7]*1e4, ...
         'XTickLabel', {'0' '3.5' '7'}, 'YTickLabel', {'0' '30' '60'}) 
%xlabel('\tau, lapse (\times 1000 trials)', 'FontSize', lfs)
xlabel('\tau_{la} (\times 1000 trials)', 'FontSize', lfs)
ylabel('\tau_{th} (\times 1000 trials)', 'FontSize', lfs)





%% Generalized RL learning models (10)
% plot
fh = figure; 
wysifig(fh, pw, ph) % set it to US letter size


% 
% plot time constants
% 

% get data
SIMNAME   = 'GenRL10';
Monk      = 'Cy';
NSEN      = 200;
NTUNE     = 36;
ELN       = [2.75 3 5 7]*1e-8; %for correlated response 
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
[bt_c, btd_c, bl_c, bld_c] = getML_RLSuppDocTau(Monk, NSEN, NTUNE, ELN, DIRS, SIMNUM, lbins_c, SIMNAME,0);


% plot
% label
% th = annotation('textbox', [0.45 0.9 0.05 0.05]);
% set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
%     'String', 'b', 'FontSize', ffs, 'FontWeight', 'bold', 'FontName', 'Helvetica') 

set(axes, 'Units', 'Normalized', 'Position', [0.2 0.2 asize*ph/pw asize], 'FontName', 'Helvetica')
Monk = 'Cy';
ml   = bl_c(3,:);
sel  = nanmean([bl_c(3,:)'-shiftdim(bld_c(3,1,:)) shiftdim(bld_c(3,2,:))-bl_c(3,:)'],2)';
mt   = bt_c(3,:);
sett = nanmean([bt_c(3,:)'-shiftdim(btd_c(3,1,:)) shiftdim(btd_c(3,2,:))-bt_c(3,:)'],2)';

hold on
plot([mlF2; mlF2], [mtF2-settF2; mtF2+settF2], '-', 'LineWidth', 0.5, 'Color', lc3e)
plot([mlF2-selF2; mlF2+selF2], [mtF2; mtF2], '-', 'LineWidth', 0.5, 'Color', lc3e)
plot(mlF2,mtF2, 'ok', 'MarkerSize', 2*ms, 'MarkerFaceColor', lc3e, 'MarkerEdgeColor', 'none')
plot([-500:100:7000], bF2(1)+bF2(2)*[-500:100:7000], '-', 'Color', lc3e', 'LineWidth', lw)

plot([ml; ml], [mt-sett; mt+sett], '-', 'LineWidth', 0.5, 'Color', lc1e)
plot([ml-sel; ml+sel], [mt; mt], '-', 'LineWidth', 0.5, 'Color', lc1e)
plot(ml,mt, 'ok', 'MarkerSize', 2*ms, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none')
b = regressXYW(ml',mt',ones(size(sel')),ones(size(sett')));
[ycb p] = regresscb(ml',mt', b, [-500:100:7000], [blmk_c btmk_c]);
plot([-500:100:7000], b(1)+b(2)*[-500:100:7000], '-k', 'LineWidth', lw)
plot([-500:100:7000]', ycb(:,1), ':', 'LineWidth', 0.5, 'Color', 'k')
plot([-500:100:7000]', ycb(:,2), ':', 'LineWidth', 0.5, 'Color', 'k')
plot(blmk_c, btmk_c, 'r*', 'MarkerSize', 3*ms)
hold off
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, ...
         'xlim', [-500 7000], 'ylim', [0 60000], 'XTick', [0 3.5 7]*1e3, 'YTick', [0 3 6]*1e4, ...
         'XTickLabel', {'0' '3.5' '7'}, 'YTickLabel', {'0' '30' '60'}) 
%xlabel('\tau, lapse (\times 1000 trials)', 'FontSize', lfs)
xlabel('\tau_{la} (\times 1000 trials)', 'FontSize', lfs)
ylabel('\tau_{th} (\times 1000 trials)', 'FontSize', lfs)





%% Generalized RL learning models (01)
% plot
fh = figure; 
wysifig(fh, pw, ph) % set it to US letter size


% 
% plot time constants
% 

% get data
SIMNAME   = 'GenRL01';
Monk      = 'Cy';
NSEN      = 200;
NTUNE     = 36;
ELN       = [0.3 0.4 0.5 0.6]*1e-8; %for correlated response 
DIRS      = [-90 90];
SIMNUM    = [81];

N      = 115985;


% get thresholds
bb     = 0;
be     = N;
bw     = 1000;
bs     = bw;
tbins_c   = [[bb:bs:be-bw]' [bb+bw:bs:be]'];
[fits_c, fitsd_c] = getML_RLSuppDocTh(Monk, NSEN, NTUNE, ELN, DIRS, SIMNUM, tbins_c, SIMNAME,0);

% get lapse
bb     = 0;
be     = N;
bw     = 250;
bs     = bw;
lbins_c  = [[bb:bs:be-bw]' [bb+bw:bs:be]'];
[lp_c, lpd_c, lN_c] = getML_RLSuppDocLapse(Monk, NSEN, NTUNE, ELN, DIRS, SIMNUM, lbins_c, SIMNAME, 0);

% get time constants
[bt_c, btd_c, bl_c, bld_c] = getML_RLSuppDocTau(Monk, NSEN, NTUNE, ELN, DIRS, SIMNUM, lbins_c, SIMNAME,0);


% plot
% label
% th = annotation('textbox', [0.45 0.9 0.05 0.05]);
% set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
%     'String', 'b', 'FontSize', ffs, 'FontWeight', 'bold', 'FontName', 'Helvetica') 

set(axes, 'Units', 'Normalized', 'Position', [0.2 0.2 asize*ph/pw asize], 'FontName', 'Helvetica')
Monk = 'Cy';
ml   = bl_c(3,:);
sel  = nanmean([bl_c(3,:)'-shiftdim(bld_c(3,1,:)) shiftdim(bld_c(3,2,:))-bl_c(3,:)'],2)';
mt   = bt_c(3,:);
sett = nanmean([bt_c(3,:)'-shiftdim(btd_c(3,1,:)) shiftdim(btd_c(3,2,:))-bt_c(3,:)'],2)';

hold on
plot([mlF2; mlF2], [mtF2-settF2; mtF2+settF2], '-', 'LineWidth', 0.5, 'Color', lc3e)
plot([mlF2-selF2; mlF2+selF2], [mtF2; mtF2], '-', 'LineWidth', 0.5, 'Color', lc3e)
plot(mlF2,mtF2, 'ok', 'MarkerSize', 2*ms, 'MarkerFaceColor', lc3e, 'MarkerEdgeColor', 'none')
plot([-500:100:12000], bF2(1)+bF2(2)*[-500:100:12000], '-', 'Color', lc3e', 'LineWidth', lw)

plot([ml; ml], [mt-sett; mt+sett], '-', 'LineWidth', 0.5, 'Color', lc1e)
plot([ml-sel; ml+sel], [mt; mt], '-', 'LineWidth', 0.5, 'Color', lc1e)
plot(ml,mt, 'ok', 'MarkerSize', 2*ms, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none')
b = regressXYW(ml',mt',ones(size(sel')),ones(size(sett')));
[ycb p] = regresscb(ml',mt', b, [-500:100:12000], [blmk_c btmk_c]);
plot([-500:100:12000], b(1)+b(2)*[-500:100:12000], '-k', 'LineWidth', lw)
plot([-500:100:12000]', ycb(:,1), ':', 'LineWidth', 0.5, 'Color', 'k')
plot([-500:100:12000]', ycb(:,2), ':', 'LineWidth', 0.5, 'Color', 'k')
plot(blmk_c, btmk_c, 'r*', 'MarkerSize', 3*ms)
hold off
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, ...
         'xlim', [-500 12000], 'ylim', [0 80000], 'XTick', [0 6 12]*1e3, 'YTick', [0 4 8]*1e4, ...
         'XTickLabel', {'0' '6' '12'}, 'YTickLabel', {'0' '40' '80'}) 
%xlabel('\tau, lapse (\times 1000 trials)', 'FontSize', lfs)
xlabel('\tau_{la} (\times 1000 trials)', 'FontSize', lfs)
ylabel('\tau_{th} (\times 1000 trials)', 'FontSize', lfs)







%% Generalized RL learning models (00)
% plot
fh = figure; 
wysifig(fh, pw, ph) % set it to US letter size


% 
% plot time constants
% 

% get data
SIMNAME   = 'GenRL00';
Monk      = 'Cy';
NSEN      = 200;
NTUNE     = 36;
ELN       = [0.5 0.55 0.6 0.7]*1e-8; %for correlated response 
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
[bt_c, btd_c, bl_c, bld_c] = getML_RLSuppDocTau(Monk, NSEN, NTUNE, ELN, DIRS, SIMNUM, lbins_c, SIMNAME,0);


% plot
% label
% th = annotation('textbox', [0.45 0.9 0.05 0.05]);
% set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
%     'String', 'b', 'FontSize', ffs, 'FontWeight', 'bold', 'FontName', 'Helvetica') 

set(axes, 'Units', 'Normalized', 'Position', [0.2 0.2 asize*ph/pw asize], 'FontName', 'Helvetica')
Monk = 'Cy';
ml   = bl_c(3,:);
sel  = nanmean([bl_c(3,:)'-shiftdim(bld_c(3,1,:)) shiftdim(bld_c(3,2,:))-bl_c(3,:)'],2)';
mt   = bt_c(3,:);
sett = nanmean([bt_c(3,:)'-shiftdim(btd_c(3,1,:)) shiftdim(btd_c(3,2,:))-bt_c(3,:)'],2)';

hold on
plot([mlF2; mlF2], [mtF2-settF2; mtF2+settF2], '-', 'LineWidth', 0.5, 'Color', lc3e)
plot([mlF2-selF2; mlF2+selF2], [mtF2; mtF2], '-', 'LineWidth', 0.5, 'Color', lc3e)
plot(mlF2,mtF2, 'ok', 'MarkerSize', 2*ms, 'MarkerFaceColor', lc3e, 'MarkerEdgeColor', 'none')
plot([-500:100:7000], bF2(1)+bF2(2)*[-500:100:7000], '-', 'Color', lc3e', 'LineWidth', lw)

plot([ml; ml], [mt-sett; mt+sett], '-', 'LineWidth', 0.5, 'Color', lc1e)
plot([ml-sel; ml+sel], [mt; mt], '-', 'LineWidth', 0.5, 'Color', lc1e)
plot(ml,mt, 'ok', 'MarkerSize', 2*ms, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none')
b = regressXYW(ml',mt',ones(size(sel')),ones(size(sett')));
[ycb p] = regresscb(ml',mt', b, [-500:100:7000], [blmk_c btmk_c]);
plot([-500:100:7000], b(1)+b(2)*[-500:100:7000], '-k', 'LineWidth', lw)
plot([-500:100:7000]', ycb(:,1), ':', 'LineWidth', 0.5, 'Color', 'k')
plot([-500:100:7000]', ycb(:,2), ':', 'LineWidth', 0.5, 'Color', 'k')
plot(blmk_c, btmk_c, 'r*', 'MarkerSize', 3*ms)
hold off
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, ...
         'xlim', [-500 7000], 'ylim', [0 100000], 'XTick', [0 3.5 7]*1e3, 'YTick', [0 5 10]*1e4, ...
         'XTickLabel', {'0' '3.5' '7'}, 'YTickLabel', {'0' '50' '100'}) 
%xlabel('\tau, lapse (\times 1000 trials)', 'FontSize', lfs)
xlabel('\tau_{la} (\times 1000 trials)', 'FontSize', lfs)
ylabel('\tau_{th} (\times 1000 trials)', 'FontSize', lfs)





%% Bound
% plot
fh = figure; 
wysifig(fh, pw, ph) % set it to US letter size


% 
% plot time constants
% 

% get data
SIMNAME   = 'Bound';
Monk      = 'Cy';
NSEN      = 200;
NTUNE     = 36;
ELN       = [20 40 60 80]; %for correlated response 
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
[bt_c, btd_c, bl_c, bld_c] = getML_RLSuppDocTau(Monk, NSEN, NTUNE, ELN, DIRS, SIMNUM, lbins_c, SIMNAME,0);


% plot example
fh = figure; 
ind = 2;
wysifig(fh, pw, ph) % set it to US letter size
set(axes, 'Units', 'Normalized', 'Position', [0.2 0.2 asize*ph/pw asize], 'FontName', 'Helvetica')
plot(nanmean(tbins_c,2),100*fits_c(1,:,ind),'o', 'MarkerFaceColor', lc1, 'MarkerEdgeColor', 'none', 'MarkerSize', ms)
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, ...
         'xlim', [0 N], 'ylim', [5 100], 'YScale', 'log', ...
         'XTick', [0 0.5 1]*10e4, 'XTickLabel', {'0' '50' '100'}, ...  
         'YTick', [5 10 20 40 80], 'YTickLabel', {'5' '10', '20', '40', '80'}) 
xlabel('Trials')
ylabel('Threshold')

set(axes, 'Units', 'Normalized', 'Position', [0.2 0.2 asize*ph/pw asize], 'FontName', 'Helvetica')
hold on
plot(nanmean(lbins_c,2),1-lp_c(:,ind), 'v', 'MarkerFaceColor', lc2, 'MarkerEdgeColor', 'none', 'MarkerSize', ms)
hold off
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, ...
         'Color', 'none', 'YColor', 'r', 'YAxisLocation', 'right', ...
         'xlim', [1 N], 'ylim', [0 0.5], 'XTick', [0 0.5 1]*10e4, 'YTick', [0 0.25 0.5], ...
         'XTickLabel', {}) 
ylabel('Lapse', 'Color', 'r')



% plot
% label
% th = annotation('textbox', [0.45 0.9 0.05 0.05]);
% set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
%     'String', 'b', 'FontSize', ffs, 'FontWeight', 'bold', 'FontName', 'Helvetica') 

fh = figure; 
wysifig(fh, pw, ph) % set it to US letter size

set(axes, 'Units', 'Normalized', 'Position', [0.2 0.2 asize*ph/pw asize], 'FontName', 'Helvetica')
Monk = 'Cy';
ml   = bl_c(3,:);
sel  = nanmean([bl_c(3,:)'-shiftdim(bld_c(3,1,:)) shiftdim(bld_c(3,2,:))-bl_c(3,:)'],2)';
mt   = bt_c(3,:);
sett = nanmean([bt_c(3,:)'-shiftdim(btd_c(3,1,:)) shiftdim(btd_c(3,2,:))-bt_c(3,:)'],2)';

hold on
plot([mlF2; mlF2], [mtF2-settF2; mtF2+settF2], '-', 'LineWidth', 0.5, 'Color', lc3e)
plot([mlF2-selF2; mlF2+selF2], [mtF2; mtF2], '-', 'LineWidth', 0.5, 'Color', lc3e)
plot(mlF2,mtF2, 'ok', 'MarkerSize', 2*ms, 'MarkerFaceColor', lc3e, 'MarkerEdgeColor', 'none')
plot([-500:100:52000], bF2(1)+bF2(2)*[-500:100:52000], '-', 'Color', lc3e', 'LineWidth', lw)

plot([ml; ml], [mt-sett; mt+sett], '-', 'LineWidth', 0.5, 'Color', lc1e)
plot([ml-sel; ml+sel], [mt; mt], '-', 'LineWidth', 0.5, 'Color', lc1e)
plot(ml,mt, 'ok', 'MarkerSize', 2*ms, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none')
b = regressXYW(ml',mt',ones(size(sel')),ones(size(sett')));
[ycb p] = regresscb(ml',mt', b, [-500:100:52000], [blmk_c btmk_c]);
plot([-500:100:52000], b(1)+b(2)*[-500:100:52000], '-k', 'LineWidth', lw)
plot([-500:100:52000]', ycb(:,1), ':', 'LineWidth', 0.5, 'Color', 'k')
plot([-500:100:52000]', ycb(:,2), ':', 'LineWidth', 0.5, 'Color', 'k')
plot(blmk_c, btmk_c, 'r*', 'MarkerSize', 3*ms)
hold off
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, ...
         'xlim', [-500 52000], 'ylim', [0 40000], 'XTick', [0 26 52]*1e3, 'YTick', [0 16 32]*1e4, ...
         'XTickLabel', {'0' '26' '52'}, 'YTickLabel', {'0' '16' '32'}) 
%xlabel('\tau, lapse (\times 1000 trials)', 'FontSize', lfs)
xlabel('\tau_{la} (\times 1000 trials)', 'FontSize', lfs)
ylabel('\tau_{th} (\times 1000 trials)', 'FontSize', lfs)

