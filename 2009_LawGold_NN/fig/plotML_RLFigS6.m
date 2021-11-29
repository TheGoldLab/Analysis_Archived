% plot reinforcement learning figure 6 - PCA
%
% Matlab script for plotting Fig 6
% for submission to Nature Neuroscience
%
% created by jcl on 08/07/08

%%
% create figure
fh = figure; 
mm2inch = 25.4;
pw = 85/mm2inch;    % paper width, make it nature size of one column, 85mm x 85mm
ph = 150/mm2inch;   % paper height
wysifig(fh, pw, ph) % set it to US letter size

% define figure parameters
alw = 0.5;  % axis line width
afs = 5;    % axis font size
lfs = 6;    % label font size
ffs = 12;   % figure font size
tfs = 12;    % text font size
lw  = 1;  % line width
elw = 0.75;
ms  = 2;    % marker size

% define line color
lc11  = hsv2rgb([0,0.8,1]);
lc12  = hsv2rgb([0,0.4,1]);
lc21  = hsv2rgb([0.667,0.8,1]);
lc22  = hsv2rgb([0.667,0.4,1]);

asize = 0.68;




%% get data
% simulate performance of the model as a function of motion directions for
% coarse discriminatin task
Monk      = 'Cy';
NSEN      = 200;
NTUNE     = 36;
ELN       = [7]*1e-8; %for correlated response
DIRS      = [-90 90];
cohS      = [6.4 99.9]/100;
TIME      = 1;
N         = 100;
SIMNUM    = [51];
recompute = 0;

[coeff_c, score_c, latent_c] = getML_RLFig6SimResponses(Monk, NSEN, NTUNE, ELN, SIMNUM, DIRS, cohS, TIME, N, recompute);

PC1_c = reshape(score_c(:,1),N,length(cohS),length(DIRS));
PC2_c = reshape(score_c(:,2),N,length(cohS),length(DIRS));


% get weights
[hdir, ldir, cdir, tdir] = dirnames;
fname = ['/getML_RL_' Monk '_' int2str(NSEN) '_' int2str(NTUNE) '_' num2str(ELN) '_v' num2str(SIMNUM(1)) '.mat'];
load([tdir fname])
fprintf([fname '\n']);
N      = length(FIRA);

TRIALS    = [0 1000 N];
trialnum  = floor(TRIALS/1000)+1;

w_c = [];
for i = 1:length(trialnum)
    w_  = W(:,:,trialnum(i));
    w_c = [w_c, w_(:)];
end
wpca_c = w_c'*coeff_c;


%% get data
% simulate performance of the model as a function of motion directions for
% fine discriminatin task
Monk      = 'Cy';
NSEN      = 200;
NTUNE     = 36;
ELN       = [15]*1e-8; %for correlated response
DIRS      = [-10 10];
cohS      = [25.6 99.9]/100;
TIME      = 1;
N         = 100;
SIMNUM    = [72];
recompute = 0;

[coeff_f, score_f, latent_f] = getML_RLFig6SimResponses(Monk, NSEN, NTUNE, ELN, SIMNUM, DIRS, cohS, TIME, N, recompute);

PC1_f = reshape(score_f(:,1),N,length(cohS),length(DIRS));
PC2_f = reshape(score_f(:,2),N,length(cohS),length(DIRS));


% get weights
[hdir, ldir, cdir, tdir] = dirnames;
fname = ['/getML_RL_' Monk '_' int2str(NSEN) '_' int2str(NTUNE) '_' num2str(ELN) '_v' num2str(SIMNUM(1)) '.mat'];
load([tdir fname])
fprintf([fname '\n']);
N      = length(FIRA);

TRIALS    = [0 1000 N];
trialnum  = floor(TRIALS/1000)+1;

w_f = [];
for i = 1:length(trialnum)
    w_  = W(:,:,trialnum(i));
    w_f = [w_f, w_(:)];
end
wpca_f = w_f'*coeff_f;

%% plot
th = annotation('textbox', [0 0.95 0.05 0.05]);
set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'String', 'a', 'FontSize', ffs, 'FontWeight', 'bold', 'FontName', 'Helvetica') 



% course discrimination
set(axes, 'Units', 'Normalized', 'Position', [0.18 0.56 asize asize*pw/ph], 'FontName', 'Helvetica')
hold on
plot(PC1_c(:,2,1), PC2_c(:,2,1), 'o', 'MarkerFaceColor', lc21, 'MarkerEdgeColor', lc21, 'MarkerSize', ms)
plot(PC1_c(:,2,2), PC2_c(:,2,2), 'o', 'MarkerFaceColor', lc11, 'MarkerEdgeColor', lc11, 'MarkerSize', ms)
plot(PC1_c(:,1,1), PC2_c(:,1,1), 'o', 'MarkerFaceColor', lc22, 'MarkerEdgeColor', lc22, 'MarkerSize', ms)
plot(PC1_c(:,1,2), PC2_c(:,1,2), 'o', 'MarkerFaceColor', lc12, 'MarkerEdgeColor', lc12, 'MarkerSize', ms)
for i = 1:size(wpca_c,1)
    x = linspace(-1500,1500,100);
    h=plot(x,-wpca_c(i,1)/wpca_c(i,2)*x, 'k');
    if     i==1
        set(h, 'linestyle', ':', 'LineWidth', lw)
    elseif i==2
        set(h, 'linestyle', '--', 'LineWidth', lw)
    else   i==3
        set(h, 'linestyle', '-', 'LineWidth', lw)
    end
    xlim([-1200 1200])
    ylim([-800 800])
end
hold off
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, ...
         'XLim', [-1200 1200], 'XTick', [-1200 -600 0 600 1200], ...
         'YLim', [-800 800], 'YTick', [-800 -400 0 400 800])
xlabel('PC1', 'FontSize', lfs)
ylabel('PC2', 'FontSize', lfs)
title('Coarse direction discrimination', 'FontSize', tfs)        





% fine discrimination
th = annotation('textbox', [0 0.45 0.05 0.05]);
set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'String', 'b', 'FontSize', ffs, 'FontWeight', 'bold', 'FontName', 'Helvetica') 

set(axes, 'Units', 'Normalized', 'Position', [0.18 0.05 asize asize*pw/ph], 'FontName', 'Helvetica')
hold on
plot(PC2_f(:,2,1), PC1_f(:,2,1), 'o', 'MarkerFaceColor', lc21, 'MarkerEdgeColor', lc21, 'MarkerSize', ms)
plot(PC2_f(:,2,2), PC1_f(:,2,2), 'o', 'MarkerFaceColor', lc11, 'MarkerEdgeColor', lc11, 'MarkerSize', ms)
plot(PC2_f(:,1,1), PC1_f(:,1,1), 'o', 'MarkerFaceColor', lc22, 'MarkerEdgeColor', lc22, 'MarkerSize', ms)
plot(PC2_f(:,1,2), PC1_f(:,1,2), 'o', 'MarkerFaceColor', lc12, 'MarkerEdgeColor', lc12, 'MarkerSize', ms)
for i = 1:size(wpca_f,1)
    x = linspace(-1500,1500,100);
    h=plot(x,-wpca_f(i,2)/wpca_f(i,1)*x, 'k');
    if     i==1
        set(h, 'linestyle', ':', 'LineWidth', lw)
    elseif i==2
        set(h, 'linestyle', '--', 'LineWidth', lw)
    else   i==3
        set(h, 'linestyle', '-', 'LineWidth', lw)
    end
end
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, ...
         'XLim', [-450 450], 'XTick', [-400 -200 0 200 400], ...
         'YLim', [-800 800], 'YTick', [-800 -400 0 400 800])
xlabel('PC1', 'FontSize', lfs)
ylabel('PC2', 'FontSize', lfs)
title('Fine direction discrimination', 'FontSize', tfs)        
        














