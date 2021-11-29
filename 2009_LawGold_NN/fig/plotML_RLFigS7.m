%% optimal weight as a function of tuning width

%% simulation parameters
recompute = 0;

Nc = 115985;
Nz = 77447;

TRIALS = [0 8000 40000 Nc];
bb     = 2.32;
be     = 6.53;
bw     = (6.53-2.32)/50;
THBINS = [2.^[bb:bw:be-bw]' 2.^[bb+bw:bw:be]'];

Monk      = 'Cy';
NSEN      = 200;
NTUNE     = 36;
ELN       = [15]*1e-8;
DIRS      = [-10,10];
SIMNUM    = [71:80];
WIDTH     = [10 20 40 60 80];

% optimal weights
NORMSCALE = sqrt(0.02);

[woavg_c, wod] = getML_RLFigS10OptimalWSortedByTh(Monk, NSEN, NTUNE, ELN, DIRS, SIMNUM, THBINS, NORMSCALE, WIDTH, recompute);


%% estimate peak of tuning width from data
% estimate from theshold bin 6 to 17, correspond to 6.9 to 13.1% COH
% thresholds
tbini = 1:length(THBINS);%[6 8:17]; % bin 7 doesn't get any data
pdir = -180:10:170;
pwp  = nans(size(WIDTH));
pwpd = nans(size(WIDTH));
pwn  = nans(size(WIDTH));
pwnd = nans(size(WIDTH));

for i = 1:length(WIDTH)
    wmax = nans(size(tbini));
    wmin = nans(size(tbini));
    for j = 1:length(tbini)
        w = woavg_c(tbini(j),:,i);
        if sum(isnan(w))==0
            % population of cell with maximum slope at 0 degree
            [kkk iii] = max(w);
            wmax(j) = pdir(iii);
            % population of cell with minimum slope at 0 degree
            [kkk iii] = min(w);
            wmin(j) = pdir(iii);
        end
    end
    pwp(i)  = nanmean(wmax);
    pwpd(i) = nanstd(wmax);
    pwn(i)  = nanmean(wmin);
    pwnd(i) = nanstd(wmin);
end


%% Estimate tuning function with the steepest slope at 0 degree

tn   = 10:1:80;
ptnp = nans(1,size(tn));
ptnn = nans(1,size(tn));
for i = 1:length(tn)
    PH = -180:1:180;
    x  = -180:1:180;
    stn = nans(1,length(PH));
    for j = 1:length(stn)
        y  = exp(-(x-PH(j)).^2/tn(i).^2/2);
        xd = x(1:end-1)+0.5;
        yd = diff(y);
        stn(j) = mean(yd(180:181));
    end
    % population of cell with maximum slope at 0 degree
    [kkk iii] = max(stn);
    ptnp(i) = PH(iii);
    % population of cell with minimum slope at 0 degree
    [kkk iii] = min(stn);
    ptnn(i) = PH(iii);
end

%% plot
% create figure
fh = figure; 
mm2inch = 25.4;
pw = 150/mm2inch;    % paper width, make it nature size of one column, 85mm x 85mm
ph = 100/mm2inch;   % paper height
wysifig(fh, pw, ph) % set it to US letter size

% define figure parameters
alw = 0.5;  % axis line width
afs = 5;    % axis font size
lfs = 6;    % label font size
ffs = 12;   % figure font size
tfs = 10;    % text font size
lw  = 1.5;  % line width
elw = 0.75;
ms  = 5;    % marker size

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

wavg_c = woavg_c;

clim_ = [-0.004 0.004];
thm   = nangeomean(THBINS,2);
L     = thm>6 & thm<80;


% plot weights
set(axes, 'Units', 'Normalized', 'Position', [0.07 0.58 asize 1.3*asize*pw/ph], 'FontName', 'Helvetica')
imagesc(pdir, 1:sum(L), wavg_c(L,:,1))
z = wavg_c(L,:,1);  % transparency map for nans
z(isnan(wavg_c(L,:,1))) = 0;
z(~isnan(wavg_c(L,:,1))) = 1;
alpha(z)
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, 'Color', [0.5 0.5 0.5], ...
         'XTick', [-170 -135 -90 -45 0 45 90 135 180], 'TickDir', 'in', ...
         'XTickLabel', {'-170', '', '', '', '0', '', '', '', '180'}, ...
         'YTick', [6 17 29 41], 'YTickLabel', {'8', '16', '32', '64'}, 'clim', clim_)
xlabel('Direction tuning (degree)', 'FontSize', lfs)
ylabel('Threshold (% COH)', 'FontSize', lfs)
%title('Initial', 'FontSize', tfs)
% plot annotation arrow 
[xb,yb] = axisXY2figXY(10,1.5);
[xe,ye] = axisXY2figXY(10,0.2);
anh     = annotation('arrow', [xb xe], [yb ye]);
set(anh, 'Color', 'r', 'LineWidth', 0, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
          'HeadWidth', 1.5*ms, 'HeadLength', 1.5*ms)
[xb,yb] = axisXY2figXY(-10,1.5);
[xe,ye] = axisXY2figXY(-10,0.2);
anh     = annotation('arrow', [xb xe], [yb ye]);
set(anh, 'Color', 'b', 'LineWidth', 0, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
          'HeadWidth', 1.5*ms, 'HeadLength', 1.5*ms)


 

set(axes, 'Units', 'Normalized', 'Position', [0.07+1*(0.03+asize) 0.58 asize 1.3*asize*pw/ph], 'FontName', 'Helvetica')
imagesc(pdir, 1:sum(L), wavg_c(L,:,2))
z = wavg_c(L,:,2);  % transparency map for nans
z(isnan(wavg_c(L,:,2))) = 0;
z(~isnan(wavg_c(L,:,2))) = 1;
alpha(z)
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, 'Color', [0.5 0.5 0.5], ...
         'XTick', [-170 -135 -90 -45 0 45 90 135 180], 'TickDir', 'in', ...
         'XTickLabel', {'-170', '', '', '', '0', '', '', '', '180'}, ...
         'YTick', [6 17 29 41], 'YTickLabel', {}, 'clim', clim_)
%xlabel('Direction tuning (degree)', 'FontSize', lfs)
%title('Early', 'FontSize', tfs)
% plot annotation arrow 
[xb,yb] = axisXY2figXY(10,1.5);
[xe,ye] = axisXY2figXY(10,0.2);
anh     = annotation('arrow', [xb xe], [yb ye]);
set(anh, 'Color', 'r', 'LineWidth', 0, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
          'HeadWidth', 1.5*ms, 'HeadLength', 1.5*ms)
[xb,yb] = axisXY2figXY(-10,1.5);
[xe,ye] = axisXY2figXY(-10,0.2);
anh     = annotation('arrow', [xb xe], [yb ye]);
set(anh, 'Color', 'b', 'LineWidth', 0, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
          'HeadWidth', 1.5*ms, 'HeadLength', 1.5*ms)

 
 

set(axes, 'Units', 'Normalized', 'Position', [0.07+2*(0.03+asize) 0.58 asize 1.3*asize*pw/ph], 'FontName', 'Helvetica')
imagesc(pdir, 1:sum(L), wavg_c(L,:,3))
z = wavg_c(L,:,3);  % transparency map for nans
z(isnan(wavg_c(L,:,3))) = 0;
z(~isnan(wavg_c(L,:,3))) = 1;
alpha(z)
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, 'Color', [0.5 0.5 0.5], ...
         'XTick', [-170 -135 -90 -45 0 45 90 135 180], 'TickDir', 'in', ...
         'XTickLabel', {'-170', '', '', '', '0', '', '', '', '180'}, ...
         'YTick', [6 17 29 41], 'YTickLabel', {}, 'clim', clim_)
%xlabel('Direction tuning (degree)', 'FontSize', lfs)
%title('Mid', 'FontSize', tfs)
% plot annotation arrow 
[xb,yb] = axisXY2figXY(10,1.5);
[xe,ye] = axisXY2figXY(10,0.2);
anh     = annotation('arrow', [xb xe], [yb ye]);
set(anh, 'Color', 'r', 'LineWidth', 0, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
          'HeadWidth', 1.5*ms, 'HeadLength', 1.5*ms)
[xb,yb] = axisXY2figXY(-10,1.5);
[xe,ye] = axisXY2figXY(-10,0.2);
anh     = annotation('arrow', [xb xe], [yb ye]);
set(anh, 'Color', 'b', 'LineWidth', 0, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
          'HeadWidth', 1.5*ms, 'HeadLength', 1.5*ms)

 

set(axes, 'Units', 'Normalized', 'Position', [0.07+3*(0.03+asize) 0.58 asize 1.3*asize*pw/ph], 'FontName', 'Helvetica')
imagesc(pdir, 1:sum(L), wavg_c(L,:,4))
z = wavg_c(L,:,4);  % transparency map for nans
z(isnan(wavg_c(L,:,4))) = 0;
z(~isnan(wavg_c(L,:,4))) = 1;
alpha(z)
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, 'Color', [0.5 0.5 0.5], ...
         'XTick', [-170 -135 -90 -45 0 45 90 135 180], 'TickDir', 'in', ...
         'XTickLabel', {'-170', '', '', '', '0', '', '', '', '180'}, ...
         'YTick', [6 17 29 41], 'YTickLabel', {}, 'clim', clim_)
%%xlabel('Direction tuning (degree)', 'FontSize', lfs)
%title('\sigma_{\theta}=80\circ', 'FontSize', tfs)
% plot annotation arrow 
[xb,yb] = axisXY2figXY(10,1.5);
[xe,ye] = axisXY2figXY(10,0.2);
anh     = annotation('arrow', [xb xe], [yb ye]);
set(anh, 'Color', 'r', 'LineWidth', 0, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
          'HeadWidth', 1.5*ms, 'HeadLength', 1.5*ms)
[xb,yb] = axisXY2figXY(-10,1.5);
[xe,ye] = axisXY2figXY(-10,0.2);
anh     = annotation('arrow', [xb xe], [yb ye]);
set(anh, 'Color', 'b', 'LineWidth', 0, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
          'HeadWidth', 1.5*ms, 'HeadLength', 1.5*ms)

      
%% plot color scale
set(axes, 'Units', 'Normalized', 'Position', [0.07+4*(0.03+asize)+0.03 0.58+0.5*1.3*asize*pw/ph+0.02 0.2*asize 0.5*1.3*asize*pw/ph-0.02], 'FontName', 'Helvetica')
dat = clim_(2)*[1:-0.01:-1]';
imagesc(dat)
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, ...
         'XTick', [], 'TickDir', 'in', ...
         'XTickLabel', {}, ...
         'YTick', [1 101 201]-0.5, 'YTickLabel', {num2str(4) '0' num2str(-4) }, ...
         'clim', clim_)
text(0.5,-25, '\times10^{-3}', 'FontSize', afs)
 


%% plot
% create figure
fh = figure; 
mm2inch = 25.4;
pw = 150/mm2inch;    % paper width, make it nature size of one column, 85mm x 85mm
ph = 100/mm2inch;   % paper height
wysifig(fh, pw, ph) % set it to US letter size

% define figure parameters
alw = 0.5;  % axis line width
afs = 5;    % axis font size
lfs = 6;    % label font size
ffs = 12;   % figure font size
tfs = 10;    % text font size
lw  = 1.5;  % line width
elw = 0.75;
ms  = 5;    % marker size

% define line color
lc1  = [1 0.2 0.2];
lc1e = [1 0.75 0.75];
lc2  = [0.2 0.2 1];
lc2e = [0.75 0.75 1];
lc3  = [0.2 0.2 0.2];
lc3e = [0.75 0.75 0.75];


asize = 0.18;

%% Location of maximum weight
set(axes, 'Units', 'Normalized', 'Position', [0.07+0.23*asize*pw/ph 0.1 asize asize*pw/ph], 'FontName', 'Helvetica')
hold on
plot(WIDTH, pwp, 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'none', 'MarkerSize', ms)
plot(WIDTH, pwn, 'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'none', 'MarkerSize', ms)
plot(0:90, zeros(size(0:90)), 'k:', 'LineWidth', 0.5)
hold off
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, ...
         'XLIM', [0 90], 'YLIM', [-180,180], ...
         'XTick', [0 20 40 60 80], 'TickDir', 'in', ...
         'YTick', [-180 -90 0 90 180])
xlabel('\sigma_{\theta} (degree)', 'FontSize', lfs)
ylabel('Location of extremum weight (degree)', 'FontSize', lfs)




%% tuning curves
set(axes, 'Units', 'Normalized', 'Position', [0.07+0.3*asize*pw/ph+(0.07+asize) 0.1 asize asize*pw/ph], 'FontName', 'Helvetica')
hold on
plot(zeros(size(0:0.1:1)), 0:0.1:1, 'k:', 'LineWidth', 0.5)
plot(-180:1:180, exp(-([-180:1:180]-80).^2/40.^2/2), 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5)
plot([-10 10], exp(-([-10 10]-80).^2/40.^2/2), 'o', 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', 'none', 'MarkerSize', ms+1)
plot(-180:1:180, exp(-([-180:1:180]-40).^2/40.^2/2), 'Color', 'k', 'LineWidth', 1.5)
plot([-10 10], exp(-([-10 10]-40).^2/40.^2/2), 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerSize', ms+1)
hold off
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, ...
         'XLIM', [-180 180], 'YLIM', [0 1], ...
         'YTick', [0 0.5 1], 'TickDir', 'in', ...
         'XTick', [-180 -90 0 90 180])
xlabel('Motion direction (degree)', 'FontSize', lfs)
ylabel('Normalized response', 'FontSize', lfs)


%% Location of maximum slopes
set(axes, 'Units', 'Normalized', 'Position', [0.07+0.42*asize*pw/ph+2*(0.07+asize) 0.1 asize asize*pw/ph], 'FontName', 'Helvetica')
hold on
plot(tn, ptnp, 'r-', 'LineWidth', 1.5)
plot(tn, ptnn, 'b-', 'LineWidth', 1.5)
plot(0:90, zeros(size(0:90)), 'k:', 'LineWidth', 0.5)
hold off
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, ...
         'XLIM', [0 90], 'YLIM', [-180,180], ...
         'XTick', [0 20 40 60 80], 'TickDir', 'in', ...
         'YTick', [-180 -90 0 90 180])
xlabel('\sigma_{\theta} (degree)', 'FontSize', lfs)
ylabel('Direction tuning of neurons with the steepest slope at 0\circ (degree)', 'FontSize', lfs)


