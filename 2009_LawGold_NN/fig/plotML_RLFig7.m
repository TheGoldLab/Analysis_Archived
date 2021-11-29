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
afs = 7;    % axis font size
lfs = 7;    % label font size
ffs = 12;   % figure font size
tfs = 8;    % text font size
lw  = 1; % line width
elw = 0.5;
ms  = 3.8;    % marker size

% define line color
lc1  = [0 0 0];
lc1e = [0.5 0.5 0.5];
lc2  = [1 0.15 0.15];
lc2e = [1 0.6 0.6];
lc3  = [0 1 1];
lc3e = [0.5 1 1];


asize = 0.33;


%% load weights and simulate performance for coarse task
Monk      = 'Cy';
NSEN      = 200;
NTUNE     = 36;
eln_      = [9]*1e-8; %for correlated response
DIRS      = [-90 90];
SIMNUM    = [51];
RECOMPUTE = 0;
[FIRA, W, MTind] = getML_RLSimPerceptualLearning11(Monk, NSEN, NTUNE, eln_, DIRS, DIRS, SIMNUM, RECOMPUTE);

recompute = 0;
DIRS      = -80:5:80;
DIRSC      = [DIRS; DIRS+180];
COHS      = [0 3.2 6.4 12.8 25.6 51.2 99.9]/100;
w         = circshift(W(:,:,end),[0,-10]);  % use weight from simulation for previous submission (trained fixed angle but [-90 90] degrees, so shift weights
N         = 500;
TAG       = 'coarse';
[pc]       = getML_RLPerformanceDirCoh(Monk, NSEN, NTUNE, DIRSC, COHS, w, MTind, N, TAG, recompute);



%% load weights and simulate performance for fine task
Monk      = 'Cy';
NSEN      = 200;
NTUNE     = 36;
eln_      = [15]*1e-8; %for correlated response
DIRS      = [-10 10];
SIMNUM    = [71];
RECOMPUTE = 0;
[FIRA, W, MTind] = getML_RLSimPerceptualLearningFine(Monk, NSEN, NTUNE, eln_, DIRS, DIRS, SIMNUM, RECOMPUTE);

recompute = 0;
DIRS      = -80:1:80;
DIRSF      = [DIRS-10; DIRS+10];
COHS      = [0 3.2 6.4 12.8 25.6 51.2 99.9]/100;
w         = W(:,:,end);
N         = 500;
TAG       = 'fine';
[pf]       = getML_RLPerformanceDirCohFine(Monk, NSEN, NTUNE, DIRSF, COHS, w, MTind, N, TAG, recompute);



%%
% get thresholds
thc  = nans(length(DIRSC), 1);
thcd = nans(length(DIRSC), 2); 
for i = 1:length(DIRSC)
    fits = ctPsych_fit(@quickTsFixLapse, [COHS(COHS~=0)' ones(size(COHS(COHS~=0)'))],  pc(i,(COHS~=0))', [], {}, [], [1,0,0,0], 0);
    thc(i)  = fits(1);
    %thcd(i,:) = sems(1,:);
end

% get thresholds
thf = [];
for i = 1:length(DIRSF)
    fits = ctPsych_fit(@quickTsFixLapse, [COHS(COHS~=0)' ones(size(COHS(COHS~=0)'))],  pf(i,(COHS~=0))', [], {}, [], [1,0,0,0], 0);
    thf(i) = fits(1);
end

% initial values for thresholds and lapse from exponential fits in Figure 2
thc_bl = 0.45;
lpc_bl = 0.5;

thf_bl = 1;
lpf_bl = 0.5;


%% plot
% coarse
th = annotation('textbox', [0.07-0.07 0.58+asize+0.02 0.05 0.05]);
set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'String', 'a', 'FontSize', ffs, 'FontWeight', 'bold', 'FontName', 'Helvetica') 

set(axes, 'Units', 'Normalized', 'Position', [0.11 0.525 asize asize*pw/ph], 'FontName', 'Helvetica')
hold on
plot(-90:5:90, lpc_bl*ones(size(-90:5:90)), ':', 'Color', lc2, 'LineWidth', lw);
plot(DIRSC(1,:), pc(:,end), 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerSize', ms);
hold off
set(gca, 'LineWidth', alw, 'FontSize', afs, ...
    'XLim', [-90 90], 'XTick', [-90 0 90], ...
    'YLim', [0.4 1.02], 'YTick', [0.5 0.75 1], 'YTickLabel', {'50' '75' '100'})
ylabel('Percentage correct', 'FontSize', lfs)
title(sprintf('Associative'), 'FontSize', 8)


th = annotation('textbox', [0.07-0.07 0.42 0.05 0.05]);
set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'String', 'c', 'FontSize', ffs, 'FontWeight', 'bold', 'FontName', 'Helvetica') 

set(axes, 'Units', 'Normalized', 'Position', [0.11 0.1 asize asize*pw/ph], 'FontName', 'Helvetica')
hold on
plot(-90:5:90, thc_bl*ones(size(-90:5:90)), ':', 'Color', lc2, 'LineWidth', lw);
plot(DIRSC(1,:), thc, 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerSize', ms);
hold off
set(gca, 'LineWidth', alw, 'FontSize', afs, ...
    'XLim', [-90 90], 'XTick', [-90 0 90], ...
    'yscale', 'log', 'YLim', [0.05 1], 'YTick', [0.05 0.1 0.2 0.4 0.8], 'YTickLabel', {'5', '10', '20', '40', '80'})
xlabel('Direction axis', 'FontSize', lfs)
ylabel('Threshold (% COH)', 'FontSize', lfs)
title(sprintf('Perceputal'), 'FontSize', 8)

th = annotation('textbox', [0.238 0.538+asize+0.02 0.1 0.1]);
set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'String', 'Coarse discrimination', 'FontSize', ffs-2, 'FontName', 'Helvetica') 



% fine
th = annotation('textbox', [0.53 0.58+asize+0.02 0.05 0.05]);
set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'String', 'b', 'FontSize', ffs, 'FontWeight', 'bold', 'FontName', 'Helvetica') 

set(axes, 'Units', 'Normalized', 'Position', [0.65 0.525 asize asize*pw/ph], 'FontName', 'Helvetica')
hold on
plot(-90:5:90, lpf_bl*ones(size(-90:5:90)), ':', 'Color', lc2, 'LineWidth', lw);
plot(mean(DIRSF), pf(:,end), 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerSize', ms);
hold off
set(gca, 'LineWidth', alw, 'FontSize', afs, ...
    'XLim', [-90 90], 'XTick', [-90 0 90], ...
    'YLim', [0.4 1.02], 'YTick', [0.5 0.75 1], 'YTickLabel', {'50' '75' '100'})
title('Fine discrimination task', 'FontSize', 8)
ylabel('Percentage correct', 'FontSize', lfs)
title(sprintf('Associative'), 'FontSize', 8)


th = annotation('textbox', [0.53 0.42 0.05 0.05]);
set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'String', 'd', 'FontSize', ffs, 'FontWeight', 'bold', 'FontName', 'Helvetica') 

set(axes, 'Units', 'Normalized', 'Position', [0.65 0.1 asize asize*pw/ph], 'FontName', 'Helvetica')
hold on
plot(-90:5:90, thf_bl*ones(size(-90:5:90)), ':', 'Color', lc2, 'LineWidth', lw);
plot(mean(DIRSF), thf, 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerSize', ms);
hold off
set(gca, 'LineWidth', alw, 'FontSize', afs, ...
    'XLim', [-90 90], 'XTick', [-90 0 90], ...
    'yscale', 'log', 'YLim', [0.2 1.05], 'YTick', [0.2 0.4 0.8], 'YTickLabel', {'20', '40', '80'})
xlabel('Direction axis', 'FontSize', lfs)
ylabel('Threshold (% COH)', 'FontSize', lfs)
title(sprintf('Perceputal'), 'FontSize', 8)


th = annotation('textbox', [0.77 0.538+asize+0.02 0.1 0.1]);
set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'String', 'Fine discrimination', 'FontSize', ffs-2, 'FontName', 'Helvetica') 


% 
% 
% subplot(2,2,2)
% plot(DIRSC(1,:), thc, '.');
% set(gca, 'ylim', [0.05 1.4], 'yscale', 'log')
% 
% subplot(2,2,3)
% plot(mean(DIRSF), pf(:,end), '.');
% set(gca, 'ylim', [0.5 1])
% 
% subplot(2,2,4)
% plot(mean(DIRSF), thf, '.');
% set(gca, 'ylim', [0.2 1], 'yscale', 'log')



%% 
%% plot data
%%%%%%%%%%%%%%%%%%%%%%%%%% plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create figure
fh = figure; 
mm2inch = 25.4;
pw = 50/mm2inch;    % paper width, make it nature size of one column, 85mm x 85mm
ph = 50/mm2inch;    % paper height
wysifig(fh, pw, ph) % set it to US letter size

% define figure parameters
alw = 0.5;  % axis line width
afs = 7;    % axis font size
lfs = 7;    % label font size
ffs = 12;   % figure font size
tfs = 8;    % text font size
lw  = 1; % line width
elw = 0.5;
ms  = 3.8;    % marker size

% define line color
lc1  = [0 0 0];
lc1e = [0.5 0.5 0.5];
lc2  = [1 0.15 0.15];
lc2e = [1 0.6 0.6];
lc3  = [0 1 1];
lc3e = [0.5 1 1];


asize = 0.75;
%%
recompute = 0;
Monk      = 'ZZ';
NSEN      = 200;
NTUNE     = 36;
eln_      = [7]*1e-8; %for correlated response
DIRS      = [-90 90];
SIMNUM    = [81:85];
[fits, sems, lp, lpd, lN] = getML_RLFig8ModelAvgPerSes(Monk, NSEN, NTUNE, eln_, DIRS, SIMNUM, recompute);



plotflag = 0;
[rzmo, dmo, stat, plotdat] = getML_RLFig8psySpecificityDuringTraining(Monk, plotflag);



% %% 
% set(axes, 'Units', 'Normalized', 'Position', [0.22 0.22 asize asize*pw/ph], 'FontName', 'Helvetica')
% 
% hold on
% plot(abs(rz), d, 'o', 'MarkerFaceColor', [0.3 0.3 0.3], 'MarkerEdgeColor', 'none', 'MarkerSize', ms)
% b = regress(d,[ones(size(rz)) abs(rz)]);
% plot([0:0.1:ceil(max(abs(rz)))], b(1)+b(2).*[0:0.1:ceil(max(abs(rz)))], '-k', 'LineWidth', lw);
% hold off
% set(gca, 'LineWidth', alw, 'FontSize', afs, ...
%     'XLim', [0 2], 'XTick', [0 1 2], ...
%     'YLim', [-0.1 0.2], 'YTick', [-0.1 0 0.1 0.2], 'YTickLabel', {'-10' '0' '10' '20'})
% ylabel('Difference in threshold (% COH)', 'FontSize', lfs)
% xlabel(sprintf('Absolute z score\nof motion direction'), 'FontSize', lfs)
% 
% 
% 
% 
% 
% 
% %%
% recompute = 0;
% Monk      = 'ZZ';
% NSEN      = 200;
% NTUNE     = 36;
% eln_      = [1]*1e-8; %for correlated response
% DIRS      = [-90 90];
% SIMNUM    = [81:85];
% [fits, sems, lp, lpd, lN] = getML_RLFig8ModelAvgPerSes(Monk, NSEN, NTUNE, eln_, DIRS, SIMNUM, recompute);
% 
% 
% monk = 'ZZ';
% plotflag = 0;
% [rzmo, dmo, stat, plotdat] = getML_RLFig8psySpecificityDuringTraining(monk, plotflag);
% 
% 
% [rMonk, dMonk, sp_cy]       = getML_psySpecificity('ZZ', 0);
% 
% 




%%
% plot behavioral specificity to target angles

% get data differently for different groups of monkeys
monk = 'ZZ';
% get task info
fn    = [monk 'TRain_psy.txt'];
a     = getML_txt(fn);
ddir  = a.data{strcmp(a.name,'ddir')};
tdir  = a.data{strcmp(a.name,'trg_dir')};
ses   = a.data{strcmp(a.name,'session')};
bs    = 25;     % include data from this session and after, this is b/c early training data are noisy

% transform data from its "direction axis" and to -90 - 90 degree
L       = ddir>90 & ddir<270;
ddir(L) = mod(ddir(L)+180,360);
L       = ddir>=180;
ddir(L) = -1*(360-ddir(L));

L       = tdir>90 & tdir<270;
tdir(L) = mod(tdir(L)+180,360);
L       = tdir>=180;
tdir(L) = -1*(360-tdir(L));

testdir = ddir;

% get psychometric performance
[fits, sems, th] = getML_psyPerformanceCT(fn,0);


fits = real(fits);
sems = real(sems);
th   = 100*real(th);
% get rid of bad fit
Lbd = th(:)<3 | ses <= bs;
fits(:,Lbd) = [];
sems(:,Lbd) = [];
ses(Lbd)    = [];
testdir(Lbd)= [];
th(Lbd)     = [];



%%% plot

 

 
%% plot against zscore for all session
% figure propreties
% lw    = 4;                      % line width
% alw   = 2;                      % axis line width
rmw   = 4;                      % running mean width

% compute values
dth = (th(:)'-nanrungeomean(th(:)',rmw));      % percentage change
dth = dth(:);

rz = nans(size(testdir));           % new way, plot against running zscore
for i = 1:length(testdir)
    tmp   = zscore(testdir([1:i]));
    rz(i) = tmp(i);
end

[R,p] = corr(dth(~isnan(dth)), abs(rz(~isnan(dth))));
[R,p] 
% % plot
% figure
% set(gcf, 'Units', 'inches', 'Position', [0 0 6 4])
% hold on
% plot(abs(rz), dth, 'o', 'MarkerFaceColor', [0.3 0.3 0.3], 'MarkerEdgeColor', 'none', 'MarkerSize', 5)
% b = regress(dth,[ones(size(testdir)) abs(rz)]);
% plot([0:0.1:ceil(max(abs(rz)))], b(1)+b(2).*[0:0.1:ceil(max(abs(rz)))], '-k', 'LineWidth', 4);
% hold off
% set(gca, 'Box', 'on', 'LineWidth', alw, 'FontSize', fsa)  % general
% XLabel('Absolute z score', 'FontSize', fsl)
% YLabel('% change in threshold', 'FontSize', fsl)
% xx = get(gca, 'XLim');
% yy = get(gca, 'YLim');
% text('Position', [0.1 0.98*yy(2)], 'String', monk, 'VerticalAlignment', 'Top', ...
%       'FontSize', fsm, 'FontWeight', 'bold') 
% text('Position', [0.98*xx(2) 0.95*yy(1)], 'String', sprintf('R=%.2f, p=%.4f', R, p), ...
%       'HorizontalAlignment', 'right', 'VerticalAlignment', 'Bottom', 'FontSize', fsa) 
% if strcmp(monk, 'Av')
%     ylim([-100 150])    % force ylim to be smaller because of a few outliers
% end




% %%%% plot relationship between R and variance of target angles
% fsa   = 18;                     % font size for axis
% fsl   = 20;                     % font size for axis labels
% fsm   = 24;                     % font size for monkey name
% ms    = 5;                      % marker size    
% 
% M  = {'Cy', 'ZZ', 'At', 'Av'};
% R  = [0.3795 0.2886 0.4337 0.1107]';
% sd = [34.0792 23.0434 40.4697 54.7698]';
% r  = [0.8259 0.9220 0.8052 0.6381]';  % concentration computed by vector average
% 
% % plot
% figure
% set(gcf, 'Units', 'inches', 'Position', [0 0 7 5])
% hold on
% plot(r, R, 'o', 'MarkerSize', 12, 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', 'none')
% [b bi xxx xxx stat] = regress(R,[ones(size(r)) r]);
% plot([0.5:0.1:1], b(1)+b(2).*[0.5:0.1:1], '-k', 'LineWidth', 3);
% hold off
% set(gca, 'Box', 'on', 'LineWidth', 2, 'FontSize', fsa, 'XLim', [0.5 1], 'YLim', [0 0.5])  % general
% XLabel('Target angles vector strength', 'FontSize', fsl)
% YLabel('Correlation')
% % plot monkey names next to markers
% for i = 1:4
%     text('Position', [r(i)+0.008 R(i)], 'String', M{i}, 'VerticalAlignment', 'Bottom', ...
%           'FontSize', fsm, 'FontWeight', 'bold') 
% end
% [aa bb] = corr(r,R);
% text('Position', [0.98 0.02], 'String', sprintf('R=%.2f, p=%.4f', aa, bb), ...
%       'HorizontalAlignment', 'right', 'VerticalAlignment', 'Bottom', 'FontSize', fsa) 
%   
%   
  
  
  
  

%%
Monk = 'ZZ';
if strcmp(Monk, 'Cy')
    set(axes, 'Units', 'Normalized', 'Position', [0.22 0.22 asize asize*pw/ph], 'FontName', 'Helvetica')

    hold on
    plot(abs(rz), dth, 'o', 'MarkerFaceColor', [1 0.3 0.3], 'MarkerEdgeColor', 'none', 'MarkerSize', ms)
    b = regress(dth,[ones(size(rz)) abs(rz)]);
    plot([0:0.1:ceil(max(abs(rz)))], b(1)+b(2).*[0:0.1:ceil(max(abs(rz)))], '-r', 'LineWidth', lw);
    
    plot(abs(rzmo), 100*dmo, 'o', 'MarkerFaceColor', [0.3 0.3 0.3], 'MarkerEdgeColor', 'none', 'MarkerSize', ms)
    b = regress(100*dmo,[ones(size(rzmo)) abs(rzmo)]);
    plot([0:0.1:ceil(max(abs(rzmo)))], b(1)+b(2).*[0:0.1:ceil(max(abs(rzmo)))], '-k', 'LineWidth', lw);
    hold off
    set(gca, 'LineWidth', alw, 'FontSize', afs, ...
        'XLim', [0 2], 'XTick', [0 1 2], ...
        'YLim', [-10 20], 'YTick', [-10 0 10 20], 'YTickLabel', {'-10' '0' '10' '20'})
    ylabel('Difference in threshold (% COH)', 'FontSize', lfs)
    xlabel(sprintf('Absolute z score\nof motion direction'), 'FontSize', lfs)
else
    set(axes, 'Units', 'Normalized', 'Position', [0.22 0.22 asize asize*pw/ph], 'FontName', 'Helvetica')

    hold on
    plot(abs(rz), dth, 'o', 'MarkerFaceColor', [1 0.3 0.3], 'MarkerEdgeColor', 'none', 'MarkerSize', ms)
    b = regress(dth,[ones(size(rz)) abs(rz)]);
    plot([0:0.1:ceil(max(abs(rz)))], b(1)+b(2).*[0:0.1:ceil(max(abs(rz)))], '-r', 'LineWidth', lw);
    
    plot(abs(rzmo), 100*dmo, 'o', 'MarkerFaceColor', [0.3 0.3 0.3], 'MarkerEdgeColor', 'none', 'MarkerSize', ms)
    b = regress(100*dmo,[ones(size(rzmo)) abs(rzmo)]);
    plot([0:0.1:ceil(max(abs(rzmo)))], b(1)+b(2).*[0:0.1:ceil(max(abs(rzmo)))], '-k', 'LineWidth', lw);
    hold off
    set(gca, 'LineWidth', alw, 'FontSize', afs, ...
        'XLim', [0 3], 'XTick', [0 1.5 3], ...
        'YLim', [-20 40], 'YTick', [-20 0 20 40])
    ylabel('Difference in threshold (% COH)', 'FontSize', lfs)
    xlabel(sprintf('Absolute z score\nof motion direction'), 'FontSize', lfs)
end

