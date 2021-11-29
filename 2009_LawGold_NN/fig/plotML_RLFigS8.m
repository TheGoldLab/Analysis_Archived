%% Figure S11
% time courses of learning for learning diff versus easy first

%% simulate performance
Monk      = 'Cy';
NSEN      = 200;
NTUNE     = 36;
eln_      = [7]*1e-8; %for correlated response
DIRS      = [-90 90];
SIMNUM    = [81];
DFLAG     = 0;
RECOMPUTE = 0;
[FIRA, W, MTind] = getML_RLSimPerceptualLearningEasyDiff(Monk, NSEN, NTUNE, eln_, DIRS, DIRS, SIMNUM, DFLAG, RECOMPUTE);


ind       = [1:15 16:3:60 61:4:117];
recompute = 0;
DIRSC     = [0; 180];
COHS      = [3.2 6.4 12.8 25.6 51.2 99.9]/100;
N         = 100;
pe        = nans(length(ind),length(COHS));
for i = 1:length(ind)
    ind(i)
    TAG       = ['easyfirst' '_' int2str(ind(i))];
    [pe(i,:)] = getML_RLPerformanceDirCoh(Monk, NSEN, NTUNE, DIRSC, COHS, W(:,:,ind(i)), MTind, N, TAG, recompute);
end


%% simulate performance
Monk      = 'Cy';
NSEN      = 200;
NTUNE     = 36;
eln_      = [7]*1e-8; %for correlated response
DIRS      = [-90 90];
SIMNUM    = [81];
DFLAG     = 1;
RECOMPUTE = 0;
[FIRA, W, MTind] = getML_RLSimPerceptualLearningEasyDiff(Monk, NSEN, NTUNE, eln_, DIRS, DIRS, SIMNUM, DFLAG, RECOMPUTE);


ind       = [1:15 16:3:60 61:4:117];
recompute = 0;
DIRSC     = [0; 180];
COHS      = [3.2 6.4 12.8 25.6 51.2 99.9]/100;
N         = 100;
pd        = nans(length(ind),length(COHS));
for i = 1:length(ind)
    ind(i)
    TAG       = ['difffirst' '_' int2str(ind(i))];
    [pd(i,:)] = getML_RLPerformanceDirCoh(Monk, NSEN, NTUNE, DIRSC, COHS, W(:,:,ind(i)), MTind, N, TAG, recompute);
end


%% get COH
Monk      = 'Cy';
NSEN      = 200;
NTUNE     = 36;
eln_      = [7]*1e-8; %for correlated response
DIRS      = [-90 90];
SIMNUM    = [81:85];
DFLAG     = 0;
RECOMPUTE = 0;

acohe = []; % all coh
for i = 1:length(SIMNUM)
    [FIRA, W, MTind] = getML_RLSimPerceptualLearningEasyDiff(Monk, NSEN, NTUNE, eln_, DIRS, DIRS, SIMNUM(i), DFLAG, RECOMPUTE);
    acohe = [acohe FIRA(:,2)];
end


Monk      = 'Cy';
NSEN      = 200;
NTUNE     = 36;
eln_      = [7]*1e-8; %for correlated response
DIRS      = [-90 90];
SIMNUM    = [81:85];
DFLAG     = 1;
RECOMPUTE = 0;

acohn = []; % all coh
for i = 1:length(SIMNUM)
    [FIRA, W, MTind] = getML_RLSimPerceptualLearningEasyDiff(Monk, NSEN, NTUNE, eln_, DIRS, DIRS, SIMNUM(i), DFLAG, RECOMPUTE);
    acohn = [acohn FIRA(:,2)];
end


%% compute thresholds
% get thresholds
the  = nans(length(ind), 1);
thed = nans(length(ind), 2); 
for i = 1:length(ind)
    fits = ctPsych_fit(@quickTsFixLapse, [COHS(COHS~=0)' ones(size(COHS(COHS~=0)'))],  pe(i,(COHS~=0))', [], {}, [], [1,0,1,0], 0);
    the(i)  = fits(1);
    %thcd(i,:) = sems(1,:);
end

% get thresholds
thd  = nans(length(ind), 1);
thdd = nans(length(ind), 2); 
for i = 1:length(ind)
    fits = ctPsych_fit(@quickTsFixLapse, [COHS(COHS~=0)' ones(size(COHS(COHS~=0)'))],  pd(i,(COHS~=0))', [], {}, [], [1,0,1,0], 0);
    thd(i)  = fits(1);
    %thcd(i,:) = sems(1,:);
end



%%
fh = figure; 
mm2inch = 25.4;
pw = 150/mm2inch;    % paper width, make it nature size of one column, 85mm x 85mm
ph = 60/mm2inch;   % paper height
wysifig(fh, pw, ph) % set it to US letter size

% define figure parameters
alw = 0.5;  % axis line width
afs = 7;    % axis font size
lfs = 7;    % label font size
ffs = 12;   % figure font size
tfs = 10;    % text font size
lw  = 1.5;  % line width
elw = 0.5;
ms  = 5;    % marker size

% define line color
lc1  = [0.5 0.5 1];
lc1e = [0.75 0.75 0.75];
lc2  = [1 0.5 0.5];
lc2e = [1 0.75 0.75];


asize = 0.6;

%%
set(axes, 'Units', 'Normalized', 'Position', [0.1 0.2 asize*ph/pw asize], 'FontName', 'Helvetica')
hold on
N = 115985;
plot(1:115985, 100*mean(acohe,2), 'o', 'MarkerFaceColor', lc1, 'MarkerEdgeColor', 'none', 'MarkerSize', 0.1)
plot(1:115985, 100*mean(acohn,2), 'o', 'MarkerFaceColor', lc2, 'MarkerEdgeColor', 'none', 'MarkerSize', 0.1)
hold off
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, ...
         'xlim', [0 N], 'ylim', [0 100], 'YScale', 'Linear', ...
         'XTick', [0 0.5 1]*10e4, 'XTickLabel', {'0' '50' '100'}, ...
         'YTick', [0 25 50 75 100], 'YTickLabel', {'0' '25' '50' '75' '100'}) 
xlabel('Trial (\times 1000)', 'FontSize', lfs)
ylabel('Mean stimulus strength (% COH)', 'FontSize', lfs)



set(axes, 'Units', 'Normalized', 'Position', [0.41 0.2 asize*ph/pw asize], 'FontName', 'Helvetica')
hold on
plot(ind*1000,1-pe(:,end), 'o', 'MarkerFaceColor', lc1, 'MarkerEdgeColor', 'none', 'MarkerSize', ms)
plot(ind*1000,1-pd(:,end), 'o', 'MarkerFaceColor', lc2, 'MarkerEdgeColor', 'none', 'MarkerSize', ms)
init = [0; 0.5; 500];
bcon = [0 0; 0.5 0.5; 100 5000];
Lgd  = ~isnan(1-pe(:,end));
[ble, bled, gof, ymt_c] = exp_fitWd(ind(Lgd)'*1000, 1-pe(Lgd,end), 0.05*ones(size(1-pe(Lgd,end))), init, bcon, {100,68});
plot(ind*1000, ble(1)+(ble(2)-ble(1)).*exp(-1/ble(3).*ind*1000), 'b', 'LineWidth', 2)
Lgd  = ~isnan(1-pd(:,end));
[bld, bldd, gof, ymt_c] = exp_fitWd(ind(Lgd)'*1000, 1-pd(Lgd,end), 0.05*ones(size(1-pd(Lgd,end))), init, bcon, {100,68});
plot(ind*1000, bld(1)+(bld(2)-bld(1)).*exp(-1/bld(3).*ind*1000), 'r', 'LineWidth', 2)
hold off
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, ...
         'xlim', [0 N], 'ylim', [0 0.5], 'YScale', 'linear', ...
         'XTick', [0 0.5 1]*10e4, 'XTickLabel', {'0' '50' '100'}, ...
         'YTick', [0 0.25 0.5], 'YTickLabel', {'0' '0.25' '0.5'}) 
xlabel('Trial (\times 1000)', 'FontSize', lfs)
ylabel('Lapse rate', 'FontSize', lfs)





set(axes, 'Units', 'Normalized', 'Position', [0.72 0.2 asize*ph/pw asize], 'FontName', 'Helvetica')
hold on
plot(ind*1000,100*the, 'o', 'MarkerFaceColor', lc1, 'MarkerEdgeColor', 'none', 'MarkerSize', ms)
plot(ind*1000,100*thd, 'o', 'MarkerFaceColor', lc2, 'MarkerEdgeColor', 'none', 'MarkerSize', ms)
init = [0.08; 0.7; 10000];
bcon = [0.06 0.2; 0.70 1; 0 100000];
Lgd  = the>0.03;
[bte, bte, gof, ymt_c] = exp_fitWd(ind(Lgd)'*1000, the(Lgd), 0.05*ones(size(the(Lgd))), init, bcon, {100,68});
plot(ind*1000, 100*(bte(1)+(bte(2)-bte(1)).*exp(-1/bte(3).*ind*1000)), 'b', 'LineWidth', 2)
init = [0.08; 0.7; 10000];
bcon = [0.06 0.2; 0.70 1; 0 100000];
Lgd  = thd>0.03;
[btd, btd, gof, ymt_c] = exp_fitWd(ind(Lgd)'*1000, thd(Lgd), 0.05*ones(size(thd(Lgd))), init, bcon, {100,68});
plot(ind*1000, 100*(btd(1)+(btd(2)-btd(1)).*exp(-1/btd(3).*ind*1000)), 'r', 'LineWidth', 2)


hold off
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, ...
         'xlim', [0 N], 'ylim', [5 100], 'YScale', 'log', ...
         'XTick', [0 0.5 1]*10e4, 'XTickLabel', {'0' '50' '100'}, ...
         'YTick', [5 10 20 40 80], 'YTickLabel', {'5' '10', '20', '40', '80'}) 
xlabel('Trial (\times 1000)', 'FontSize', lfs)
ylabel('Threshold (% coh)', 'FontSize', lfs)

     
     
     





%% plot pooling weights
recompute = 0;

Nc = 115985;
Nz = 77447;

TRIALS = [0 4000 40000 Nc];
bb     = 2.32;
be     = 6.53;
bw     = (6.53-2.32)/50;
THBINS = [2.^[bb:bw:be-bw]' 2.^[bb+bw:bw:be]'];
   
Monk      = 'Cy';
NSEN      = 200;
NTUNE     = 36;
ELN       = [7]*1e-8;
DIRS      = [180 0];
SIMNUM    = [81:85];
DFLAG     = 0;
[wavge_c, wsee_c] = getML_RLFigS11AvgWSortedByTh(Monk, NSEN, NTUNE, ELN, DIRS, SIMNUM, TRIALS, THBINS, DFLAG, recompute);


recompute = 0;

Nc = 115985;
Nz = 77447;

TRIALS = [0 4000 40000 Nc];
bb     = 2.32;
be     = 6.53;
bw     = (6.53-2.32)/50;
THBINS = [2.^[bb:bw:be-bw]' 2.^[bb+bw:bw:be]'];
   
Monk      = 'Cy';
NSEN      = 200;
NTUNE     = 36;
ELN       = [7]*1e-8;
DIRS      = [180 0];
SIMNUM    = [81:85];
DFLAG     = 1;
[wavgd_c, wsed_c] = getML_RLFigS11AvgWSortedByTh(Monk, NSEN, NTUNE, ELN, DIRS, SIMNUM, TRIALS, THBINS, DFLAG, recompute);


% % shift weights
wavg_c  = circshift(wavg_c,[0,-1]);
wse_c   = circshift(wse_c,[0,-1]);
woavg_c = circshift(woavg_c,[0,-1]);
DTUNE   = 360/NTUNE;
pdir    = linspace(-170,180,NTUNE);

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
%% plot weights (binned by thresholds)
% label
% th = annotation('textbox', [0 0.95 0.05 0.05]);
% set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
%     'String', 'a', 'FontSize', ffs, 'FontWeight', 'bold', 'FontName', 'Helvetica') 



clim_ = [-0.004 0.004];
thm   = nangeomean(THBINS,2);
L     = thm>6 & thm<80;


% plot weights
set(axes, 'Units', 'Normalized', 'Position', [0.07 0.6 asize 1.3*asize*pw/ph], 'FontName', 'Helvetica')
imagesc(pdir, 1:sum(L), wavge_c(L,:,1))
z = wavge_c(L,:,1);  % transparency map for nans
z(isnan(wavge_c(L,:,1))) = 0;
z(~isnan(wavge_c(L,:,1))) = 1;
alpha(z)
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, 'Color', [0.5 0.5 0.5], ...
         'XTick', [-170 -135 -90 -45 0 45 90 135 180], 'TickDir', 'in', ...
         'XTickLabel', {'-170', '', '', '', '0', '', '', '', '180'}, ...
         'YTick', [6 17 29 41], 'YTickLabel', {'8', '16', '32', '64'}, 'clim', clim_)
%xlabel('Direction tuning (degree)', 'FontSize', lfs)
ylabel('Threshold (% COH)', 'FontSize', lfs)
title('Initial', 'FontSize', tfs)
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



set(axes, 'Units', 'Normalized', 'Position', [0.07+1*(0.03+asize) 0.6 asize 1.3*asize*pw/ph], 'FontName', 'Helvetica')
imagesc(pdir, 1:sum(L), wavge_c(L,:,2))
z = wavge_c(L,:,2);  % transparency map for nans
z(isnan(wavge_c(L,:,2))) = 0;
z(~isnan(wavge_c(L,:,2))) = 1;
alpha(z)
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, 'Color', [0.5 0.5 0.5], ...
         'XTick', [-170 -135 -90 -45 0 45 90 135 180], 'TickDir', 'in', ...
         'XTickLabel', {'-170', '', '', '', '0', '', '', '', '180'}, ...
         'YTick', [6 17 29 41], 'YTickLabel', {}, ...
         'clim', clim_)
title('Early', 'FontSize', tfs)
hla = xlabel('Direction tuning (degree)', 'FontSize', lfs)
lpos = get(hla, 'position');
lpos(1) = 220;
set(hla, 'position', lpos)

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

      
set(axes, 'Units', 'Normalized', 'Position', [0.07+2*(0.03+asize) 0.6 asize 1.3*asize*pw/ph], 'FontName', 'Helvetica')
imagesc(pdir, 1:sum(L), wavge_c(L,:,3))
z = wavge_c(L,:,3);  % transparency map for nans
z(isnan(wavge_c(L,:,3))) = 0;
z(~isnan(wavge_c(L,:,3))) = 1;
alpha(z)
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, 'Color', [0.5 0.5 0.5], ...
         'XTick', [-170 -135 -90 -45 0 45 90 135 180], 'TickDir', 'in', ...
         'XTickLabel', {'-170', '', '', '', '0', '', '', '', '180'}, ...
         'YTick', [6 17 29 41], 'YTickLabel', {}, ...
         'clim', clim_)
title('Mid', 'FontSize', tfs)
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

      

set(axes, 'Units', 'Normalized', 'Position', [0.07+3*(0.03+asize) 0.6 asize 1.3*asize*pw/ph], 'FontName', 'Helvetica')
imagesc(pdir, 1:sum(L), wavge_c(L,:,4))
z = wavge_c(L,:,4);  % transparency map for nans
z(isnan(wavge_c(L,:,4))) = 0;
z(~isnan(wavge_c(L,:,4))) = 1;
alpha(z)
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, 'Color', [0.5 0.5 0.5], ...
         'XTick', [-170 -135 -90 -45 0 45 90 135 180], 'TickDir', 'in', ...
         'XTickLabel', {'-170', '', '', '', '0', '', '', '', '180'}, ...
         'YTick', [6 17 29 41], 'YTickLabel', {}, ...
         'clim', clim_)
title('Late', 'FontSize', tfs)
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
set(axes, 'Units', 'Normalized', 'Position', [0.07+4*(0.03+asize)+0.03 0.6+0.5*1.3*asize*pw/ph+0.02 0.2*asize 0.5*1.3*asize*pw/ph-0.02], 'FontName', 'Helvetica')
dat = clim_(2)*[1:-0.01:-1]';
imagesc(dat)
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, ...
         'XTick', [], 'TickDir', 'in', ...
         'XTickLabel', {}, ...
         'YTick', [1 101 201]-0.5, 'YTickLabel', {num2str(4) '0' num2str(-4) }, ...
         'clim', clim_)
text(0.5,-25, '\times10^{-3}', 'FontSize', afs)





%% plot weights (binned by thresholds)
% label
% th = annotation('textbox', [0 0.95 0.05 0.05]);
% set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
%     'String', 'a', 'FontSize', ffs, 'FontWeight', 'bold', 'FontName', 'Helvetica') 



clim_ = [-0.004 0.004];
thm   = nangeomean(THBINS,2);
L     = thm>6 & thm<80;


% plot weights
set(axes, 'Units', 'Normalized', 'Position', [0.07 0.1 asize 1.3*asize*pw/ph], 'FontName', 'Helvetica')
imagesc(pdir, 1:sum(L), wavgd_c(L,:,1))
z = wavgd_c(L,:,1);  % transparency map for nans
z(isnan(wavgd_c(L,:,1))) = 0;
z(~isnan(wavgd_c(L,:,1))) = 1;
alpha(z)
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, 'Color', [0.5 0.5 0.5], ...
         'XTick', [-170 -135 -90 -45 0 45 90 135 180], 'TickDir', 'in', ...
         'XTickLabel', {'-170', '', '', '', '0', '', '', '', '180'}, ...
         'YTick', [6 17 29 41], 'YTickLabel', {'8', '16', '32', '64'}, 'clim', clim_)
%xlabel('Direction tuning (degree)', 'FontSize', lfs)
ylabel('Threshold (% COH)', 'FontSize', lfs)
title('Initial', 'FontSize', tfs)
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



set(axes, 'Units', 'Normalized', 'Position', [0.07+1*(0.03+asize) 0.1 asize 1.3*asize*pw/ph], 'FontName', 'Helvetica')
imagesc(pdir, 1:sum(L), wavgd_c(L,:,2))
z = wavgd_c(L,:,2);  % transparency map for nans
z(isnan(wavgd_c(L,:,2))) = 0;
z(~isnan(wavgd_c(L,:,2))) = 1;
alpha(z)
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, 'Color', [0.5 0.5 0.5], ...
         'XTick', [-170 -135 -90 -45 0 45 90 135 180], 'TickDir', 'in', ...
         'XTickLabel', {'-170', '', '', '', '0', '', '', '', '180'}, ...
         'YTick', [6 17 29 41], 'YTickLabel', {}, ...
         'clim', clim_)
title('Early', 'FontSize', tfs)
hla = xlabel('Direction tuning (degree)', 'FontSize', lfs)
lpos = get(hla, 'position');
lpos(1) = 220;
set(hla, 'position', lpos)

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

      
set(axes, 'Units', 'Normalized', 'Position', [0.07+2*(0.03+asize) 0.1 asize 1.3*asize*pw/ph], 'FontName', 'Helvetica')
imagesc(pdir, 1:sum(L), wavgd_c(L,:,3))
z = wavgd_c(L,:,3);  % transparency map for nans
z(isnan(wavgd_c(L,:,3))) = 0;
z(~isnan(wavgd_c(L,:,3))) = 1;
alpha(z)
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, 'Color', [0.5 0.5 0.5], ...
         'XTick', [-170 -135 -90 -45 0 45 90 135 180], 'TickDir', 'in', ...
         'XTickLabel', {'-170', '', '', '', '0', '', '', '', '180'}, ...
         'YTick', [6 17 29 41], 'YTickLabel', {}, ...
         'clim', clim_)
title('Mid', 'FontSize', tfs)
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

      

set(axes, 'Units', 'Normalized', 'Position', [0.07+3*(0.03+asize) 0.1 asize 1.3*asize*pw/ph], 'FontName', 'Helvetica')
imagesc(pdir, 1:sum(L), wavgd_c(L,:,4))
z = wavgd_c(L,:,4);  % transparency map for nans
z(isnan(wavgd_c(L,:,4))) = 0;
z(~isnan(wavgd_c(L,:,4))) = 1;
alpha(z)
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, 'Color', [0.5 0.5 0.5], ...
         'XTick', [-170 -135 -90 -45 0 45 90 135 180], 'TickDir', 'in', ...
         'XTickLabel', {'-170', '', '', '', '0', '', '', '', '180'}, ...
         'YTick', [6 17 29 41], 'YTickLabel', {}, ...
         'clim', clim_)
title('Late', 'FontSize', tfs)
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
set(axes, 'Units', 'Normalized', 'Position', [0.07+4*(0.03+asize)+0.03 0.6+0.5*1.3*asize*pw/ph+0.02 0.2*asize 0.5*1.3*asize*pw/ph-0.02], 'FontName', 'Helvetica')
dat = clim_(2)*[1:-0.01:-1]';
imagesc(dat)
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, ...
         'XTick', [], 'TickDir', 'in', ...
         'XTickLabel', {}, ...
         'YTick', [1 101 201]-0.5, 'YTickLabel', {num2str(4) '0' num2str(-4) }, ...
         'clim', clim_)
text(0.5,-25, '\times10^{-3}', 'FontSize', afs)
