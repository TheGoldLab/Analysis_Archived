%% plotML_fig5.m
% Matlab script for plotting Fig 5  (Relationship of LIP DV with behavioral, motor and 
% motivational parameters)
% Final version for publication to Nature Neuroscience
%
%


%% create figure
% create figure
fh = figure;
mm2inch = 25.4;
pw = 180/mm2inch;   % paper width
ph = 60/mm2inch;   % paper height
wysifig(fh, pw, ph) % set it to US letter size

% define figure parameters
alw = 0.5;  % axis line width
afs = 5;    % axis font size
lfs = 5;    % label font size
ffs = 12;   % figure font size
tfs = 6;    % text font size
lw  = 1.5;  % line width
elw = 0.5;
ms  = 3;    % marker size

% define line color
lc1  = [0 0 0];
lc1e = [0.75 0.75 0.75];
lc2  = [1 0.3 0.3];
lc2e = [1 0.6 0.6];


%% plot correlation of LIP threshold with different parameters
% LIP
a                    = getML_txt('CyTRain_LIP.txt');
sl_cy                = a.data{strcmp(a.name, 'session')};   
a                    = getML_txt('ZZTRain_LIP.txt');
sl_zz                = a.data{strcmp(a.name, 'session')};   
A_cy                 = getML_LIPModelroc('Cy', 0);
A_zz                 = getML_LIPModelroc('ZZ', 0);

[pt_cy nt_cy ptci_cy]   = getML_performanceAsso('CyTRain_psy.txt', 0);
[pt_zz nt_zz ptci_zz]   = getML_performanceAsso('ZZTRain_psy.txt', 0);
[Be_cy, Bed_cy]   = getML_psyPerformanceCTFixLapse('CyTRain_psy.txt', 0, 1-pt_cy);
[Be_zz, Bed_zz]   = getML_psyPerformanceCTFixLapse('ZZTRain_psy.txt', 0, 1-pt_zz);
%[Be_cy, Bed_cy] = getML_psyPerformanceGSDV('CyTRain_psy.txt', 0);
%[Be_zz, Bed_zz] = getML_psyPerformanceGSDV('ZZTRain_psy.txt', 0);
[MT_cy, MTd_cy]   = getML_neurometricROCT('CyTRain_MT.txt', [], 0, [0 1], 0, 0);
[MT_zz, MTd_zz]   = getML_neurometricROCT('CyTRain_MT.txt', [], 0, [0 1], 0, 0);
%[LIP_cy, LIPd_cy, p] = getML_LIPModel('Cy', 0, 1000*A_cy(3,:)');
%[LIP_zz, LIPd_zz, p] = getML_LIPModel('ZZ', 0, 1000*A_zz(3,:)');
[L99_cy, Ld99_cy] = getML_LIPModelroc99('Cy', 0);
[L99_zz, Ld99_zz] = getML_LIPModelroc99('ZZ', 0);
[La99_cy, Lad99_cy] = getML_LIPModelroc99asso('Cy', 0);
[La99_zz, Lad99_zz] = getML_LIPModelroc99asso('ZZ', 0);
[LIP_cy, LIPd_cy] = getML_LIPCTdependence('Cy', 0);
[LIP_zz, LIPd_zz] = getML_LIPCTdependence('ZZ', 0);
Ll_cy = getML_LapseSelectionArray('Cy', 'LIP');
Ll_zz = getML_LapseSelectionArray('ZZ', 'LIP');


% th = annotation('textbox', [0 0.55 0.05 0.05]);
% set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
%     'String', 'd', 'FontSize', ffs) 


asize = 0.08;

th = annotation('textbox', [0.18 0.95 0.05 0.05]);
set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'String', 'Behaviour', 'FontSize', tfs, 'FontName', 'Helvetica') 
annotation('line', [0.075 0.335],   [0.93 0.93])
annotation('line', [0.075 0.075],   [0.92 0.93])
annotation('line', [0.335 0.335],   [0.92 0.93])
annotation('line', [0.205 0.205],   [0.93 0.935])

th = annotation('textbox', [0.425 0.95 0.2 0.05]);
set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'String', 'Motor', 'FontSize', tfs, 'FontName', 'Helvetica') 
annotation('line', [0.375 0.675],   [0.93 0.93])
annotation('line', [0.375 0.375],   [0.92 0.93])
annotation('line', [0.675 0.675],   [0.92 0.93])
annotation('line', [0.525 0.525],   [0.93 0.935])

th = annotation('textbox', [0.82 0.95 0.05 0.05]);
set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'String', 'Motivation', 'FontSize', tfs, 'FontName', 'Helvetica') 
annotation('line', [0.695 0.99],    [0.93 0.93])
annotation('line', [0.695 0.695],   [0.92 0.93])
annotation('line', [0.99 0.99],     [0.92 0.93])
annotation('line', [0.8425 0.8425], [0.93 0.935])



set(axes, 'Units', 'Normalized', 'Position', [0.03 0.92 asize asize*pw/ph], 'FontName', 'Helvetica')
axis off
text(0,0,'Monkey C','HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
     'FontSize', tfs, 'Rotation', 0)

set(axes, 'Units', 'Normalized', 'Position', [0.03 0.51 asize asize*pw/ph], 'FontName', 'Helvetica')
axis off
text(0,0,'Monkey Z','HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
     'FontSize', tfs, 'Rotation', 0)



%% beta 1
set(axes, 'Units', 'Normalized', 'Position', [0.08 0.6 asize asize*pw/ph], 'FontName', 'Helvetica')
LIP = LIP_cy(4,:)';
Be  = 100*Be_cy(1,sl_cy(:))';
Bed = 100*shiftdim(Bed_cy(1,:,sl_cy(:)),1);
Lp  = getML_performanceAsso('CyTRain_psy.txt', 0);
Lp  = 1-Lp(sl_cy(:))';
% get residual of LIP given lapse
[b bi LIP_Lp LIPd_Lp] = regress(LIP,[ones(size(Lp)) Lp], 0.68);
[b bi Be_Lp Bed_Lp]   = regress(Be,[ones(size(Lp)) Lp], 0.68);
Lgd     = Lp<=0.05;
LIP     = LIP(Lgd);
LIP_Lp  = LIP_Lp(Lgd);
LIPd_Lp = LIPd_Lp(Lgd,:);
Bed     = Bed(:,Lgd);
Be      = Be(Lgd);
Lp      = Lp(Lgd);
%plot
hold on
plot([Be'; Be'], LIPd_Lp', 'LineWidth', elw, 'Color', lc1e)
plot(Bed, [LIP_Lp'; LIP_Lp'], 'LineWidth', elw, 'Color', lc1e)
plot(Be, LIP_Lp, 'ko', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerSize', ms-1)
hold off
Lbd = isnan(LIP_Lp) | isnan(Be);
[rlbc, plbc] = corr(LIP_Lp(~Lbd),log(Be(~Lbd)));
title(['\itr\rm = ' num2str(rlbc, '%.2f') '^\clubsuit' '*' ], 'FontSize', tfs)
% set(gca, 'YLim', [-200 400], 'YTick', [-200 0 200 400], ...
%            'XLim', [0 28], 'XTick', [0 14 28], 'XTickLabel', [0 14 28], ...
%            'FontSize', afs)
set(gca, 'YLim', [-6 13], 'YTick', [-6 0 6 12], ...
           'XLim', [8 60], 'XTick', [10 20 40], 'XTickLabel', [10 20 40], ...
           'XScale', 'log', 'FontSize', afs)
ylabel(['k_3 (LIP)' sprintf('\nafter accounting for\nerrors at 99%% coh')], 'FontSize', lfs)
% normal corr between beta 1 and threshold
Lbd = isnan(LIP) | isnan(Be);
[r, p] = corr(LIP(~Lbd),log(Be(~Lbd)));



set(axes, 'Units', 'Normalized', 'Position', [0.08 0.2 asize asize*pw/ph], 'FontName', 'Helvetica')
LIP = LIP_zz(4,:)';
Be  = 100*Be_zz(1,sl_zz(:))';
Bed = 100*shiftdim(Bed_zz(1,:,sl_zz(:)),1);
Lp  = getML_performanceAsso('ZZTRain_psy.txt', 0);
Lp  = 1-Lp(sl_zz(:))';
% get residual of LIP given lapse
[b bi LIP_Lp LIPd_Lp] = regress(LIP,[ones(size(Lp)) Lp], 0.68);
[b bi Be_Lp Bed_Lp]   = regress(Be,[ones(size(Lp)) Lp], 0.68);
Lgd    = Lp<=0.15;
LIP    = LIP(Lgd);
LIP_Lp = LIP_Lp(Lgd);
LIPd_Lp = LIPd_Lp(Lgd,:);
Bed    = Bed(:,Lgd);
Be     = Be(Lgd);
Lp     = Lp(Lgd);
%plot
hold on
plot([Be'; Be'], LIPd_Lp', 'LineWidth', elw, 'Color', lc1e)
plot(Bed, [LIP_Lp'; LIP_Lp'], 'LineWidth', elw, 'Color', lc1e)
plot(Be, LIP_Lp, 'ko', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerSize', ms-1)
hold off
Lbd = isnan(LIP_Lp) | isnan(Be);
[rlbz, plbz] = corr(LIP_Lp(~Lbd),log(Be(~Lbd)));
title(['\itr\rm = ' num2str(rlbz, '%.2f') '^\clubsuit' '*' ], 'FontSize', tfs)
set(gca, 'YLim', [-3 6], 'YTick', [-3 0 3 6], ...
           'XLim', [12 65], 'XTick', [12 24 48], 'XTickLabel', [12 24 48], ...
           'XScale', 'log', 'FontSize', afs)
xlabel('Threshold (% coh)', 'FontSize', lfs)
ylabel(['k_3 (LIP)' sprintf('\nafter accounting for\nerrors at 99%% coh')], 'FontSize', lfs)
% normal corr between beta 1 and threshold
Lbd = isnan(LIP) | isnan(Be);
[r, p] = corr(LIP(~Lbd),log(Be(~Lbd)));



% lapse 
set(axes, 'Units', 'Normalized', 'Position', [0.08+1*(asize+0.09) 0.6 asize asize*pw/ph], 'FontName', 'Helvetica')
LIP = LIP_cy(4,:)';
Be  = 100*Be_cy(1,sl_cy(:))';
[Lp n Lpd] = getML_performanceAsso('CyTRain_psy.txt', 0);
Lp  = 1-Lp(sl_cy(:))';
Lpd = 1-Lpd(:,sl_cy(:))';
[b bi LIP_Be LIPd_Be] = regress(LIP,[ones(size(Be)) log(Be)], 0.68);
[b bi Lp_Be Lpd_Be] = regress(Lp,[ones(size(Be)) log(Be)], 0.68);
hold on
plot([Lp'; Lp'], LIPd_Be', 'LineWidth', elw, 'Color', lc1e)
plot(Lpd', [LIP_Be'; LIP_Be'], 'LineWidth', elw, 'Color', lc1e)
plot(Lp, LIP_Be, 'ko', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerSize', ms-1)
hold off
Lbd = isnan(LIP_Be) | isnan(Lp_Be);
[rllc, pllc] = corr(LIP_Be(~Lbd),Lp(~Lbd));
title(['\itr\rm = ' num2str(rllc, '%.2f') '  *'], 'FontSize', tfs)
ylabel(['k_3 (LIP)' sprintf('\nafter accounting for\n') 'behavioural threshold'], 'FontSize', lfs)
set(gca, 'YLim', [-6 13], 'YTick', [-6 0 6 12], ...
           'XLim', [0 0.4], 'XTick', [0 0.2 0.4], 'FontSize', afs)
% normal corr between beta 1 and threshold
Lbd = isnan(LIP) | isnan(Lp);
[r, p] = corr(LIP(~Lbd),Lp(~Lbd));


       
set(axes, 'Units', 'Normalized', 'Position', [0.08+1*(asize+0.09) 0.2 asize asize*pw/ph], 'FontName', 'Helvetica')
LIP = LIP_zz(4,:)';
Be  = 100*Be_zz(1,sl_zz(:))';
[Lp n Lpd] = getML_performanceAsso('ZZTRain_psy.txt', 0);
Lp  = 1-Lp(sl_zz(:))';
Lpd = 1-Lpd(:,sl_zz(:))';
[b bi LIP_Be LIPd_Be] = regress(LIP,[ones(size(Be)) log(Be)], 0.68);
[b bi Lp_Be Lpd_Be] = regress(Lp,[ones(size(Be)) log(Be)], 0.68);
hold on
plot([Lp'; Lp'], LIPd_Be', 'LineWidth', elw, 'Color', lc1e)
plot(Lpd', [LIP_Be'; LIP_Be'], 'LineWidth', elw, 'Color', lc1e)
plot(Lp, LIP_Be, 'ko', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerSize', ms-1)
hold off
Lbd = isnan(LIP_Be) | isnan(Lp_Be);
[rllz, pllz] = corr(LIP_Be(~Lbd),Lp(~Lbd));
title(['\itr\rm = ' num2str(rllz, '%.2f') '  *'], 'FontSize', tfs)
set(gca, 'YLim', [-3 6], 'YTick', [-3 0 3 6], ...
           'XLim', [0 0.45], 'XTick', [0 0.2 0.4], 'FontSize', afs)
ylabel(['k_3 (LIP)' sprintf('\nafter accounting for\n') 'behavioural threshold'], 'FontSize', lfs)
xlabel([sprintf('Error rate\nat 99%% coh (') '\lambda' ')'], 'FontSize', lfs)
% normal corr between beta 1 and threshold
Lbd = isnan(LIP) | isnan(Lp);
[r, p] = corr(LIP(~Lbd),Lp(~Lbd));





%%         motor           %%%
% Saccade latency
set(axes, 'Units', 'Normalized', 'Position', [0.1+2*(asize+0.02)+0.084 0.6 asize asize*pw/ph], 'FontName', 'Helvetica')
[s sd] = getML_saccade('CyTRain_psy.txt', 0);
L  = Ll_cy;
x  = LIP_cy(4,L);
xs = LIPd_cy(4,L);
xs(xs==0) = nanmean(xs);
xd = [x-LIPd_cy(4,L); LIPd_cy(4,L)+x];
y  = s(sl_cy(L),7)'/1000;
ys = sd(sl_cy(L),7)'/1000;
yd = sd(sl_cy(L),7)'/1000;
[b, bint, p, r, rp] = getLinearDependence(x,y);
hold on
plot([y-yd; y+yd], [x; x], 'LineWidth', elw, 'Color', lc1e)
plot([y; y], xd, 'LineWidth', elw, 'Color', lc1e)
plot(y, x, 'ko', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerSize', ms-1)
hold off
ylabel('k_3 (LIP)', 'FontSize', lfs)
set(gca, 'YLim', [-4 18], 'YTick', [0 9 18], ...
           'XLim', [0.140 0.220], 'XTick', [0.140 0.180 0.220], 'FontSize', afs)
title(['\itr\rm = ' num2str(r,'%.2f') '^\clubsuit'], 'FontSize', tfs)


set(axes, 'Units', 'Normalized', 'Position', [0.1+2*(asize+0.02)+0.084 0.2 asize asize*pw/ph], 'FontName', 'Helvetica')
[s sd] = getML_saccade('ZZTRain_psy.txt', 0);
L  = Ll_zz;
x  = LIP_zz(4,L);
xs = LIPd_zz(4,L);
xs(xs==0) = nanmean(xs);
xd = [x-LIPd_zz(4,L); LIPd_zz(4,L)+x];
y  = s(sl_zz(L),7)'/1000;
ys = sd(sl_zz(L),7)'/1000;
yd = sd(sl_zz(L),7)'/1000;
[b, bint, p, r, rp] = getLinearDependence(x,y);
hold on
plot([y-yd; y+yd], [x; x], 'LineWidth', elw, 'Color', lc1e)
plot([y; y], xd, 'LineWidth', elw, 'Color', lc1e)
plot(y, x, 'ko', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerSize', ms-1)
hold off
ylabel('k_3 (LIP)', 'FontSize', lfs)
set(gca, 'YLim', [-2 10], 'YTick', [0 5 10], ...
           'XLim', [0.170 0.220], 'XTick', [0.170 0.195 0.220], 'FontSize', afs)
title(['\itr\rm = ' num2str(r,'%.2f') '  *'], 'FontSize', tfs)
xlabel(sprintf('Saccade latency\n(s)'), 'FontSize', lfs)



% avg vel
set(axes, 'Units', 'Normalized', 'Position', [0.1+3*(asize+0.02)+0.084 0.6 asize asize*pw/ph], 'FontName', 'Helvetica')
[s sd] = getML_saccade('CyTRain_psy.txt', 0);
L  = Ll_cy;
x  = LIP_cy(4,L);
xs = LIPd_cy(4,L);
xs(xs==0) = nanmean(xs);
xd = [x-LIPd_cy(4,L); LIPd_cy(4,L)+x];
y  = 1000*s(sl_cy(L),4)';
ys = 1000*sd(sl_cy(L),4)';
yd = 1000*sd(sl_cy(L),4)';
[b, bint, p, r, rp] = getLinearDependence(x,y);
hold on
plot([y-yd; y+yd], [x; x], 'LineWidth', elw, 'Color', lc1e)
plot([y; y], xd, 'LineWidth', elw, 'Color', lc1e)
plot(y, x, 'ko', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerSize', ms-1)
hold off
set(gca, 'YLim', [-4 18], 'YTick', [0 9 18], 'YTickLabel', [], ...
          'XLim', [120 240], 'XTick', [120 180 240], 'FontSize', afs)
title(['\itr\rm = ' num2str(r,'%.2f') '  *'], 'FontSize', tfs)


set(axes, 'Units', 'Normalized', 'Position', [0.1+3*(asize+0.02)+0.084 0.2 asize asize*pw/ph], 'FontName', 'Helvetica')
[s sd] = getML_saccade('ZZTRain_psy.txt', 0);
L  = Ll_zz;
x  = LIP_zz(4,L);
xs = LIPd_zz(4,L);
xs(xs==0) = nanmean(xs);
xd = [x-LIPd_zz(4,L); LIPd_zz(4,L)+x];
y  = 1000*s(sl_zz(L),4)';
ys = 1000*sd(sl_zz(L),4)';
yd = 1000*sd(sl_zz(L),4)';
[b, bint, p, r, rp] = getLinearDependence(x,y);
hold on
plot([y-yd; y+yd], [x; x], 'LineWidth', elw, 'Color', lc1e)
plot([y; y], xd, 'LineWidth', elw, 'Color', lc1e)
plot(y, x, 'ko', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerSize', ms-1)
hold off
set(gca, 'YLim', [-2 10], 'YTick', [0 5 10], 'YTickLabel', [], ...
          'XLim', [150 310], 'XTick', [150 225 300], 'FontSize', afs)
title(['\itr\rm = ' num2str(r,'%.2f')], 'FontSize', tfs)
xlabel(sprintf('Saccade\nvelocity (deg s^{-1})'), 'FontSize', lfs)



% accuracy
set(axes, 'Units', 'Normalized', 'Position', [0.1+4*(asize+0.02)+0.084 0.6 asize asize*pw/ph], 'FontName', 'Helvetica')
[s sd] = getML_saccade('CyTRain_psy.txt', 0);
L  = Ll_cy;
x  = LIP_cy(4,L);
xs = LIPd_cy(4,L);
xs(xs==0) = nanmean(xs);
xd = [x-LIPd_cy(4,L); LIPd_cy(4,L)+x];
y  = s(sl_cy(L),6)';
ys = sd(sl_cy(L),6)';
yd = sd(sl_cy(L),6)';
[b, bint, p, r, rp] = getLinearDependence(x,y);
hold on
plot([y-yd; y+yd], [x; x], 'LineWidth', elw, 'Color', lc1e)
plot([y; y], xd, 'LineWidth', elw, 'Color', lc1e)
plot(y, x, 'ko', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerSize', ms-1)
hold off
set(gca, 'YLim', [-4 18], 'YTick', [0 9 18], 'YTickLabel', [], ...
          'XLim', [-2 2], 'XTick', [-2 0 2], 'FontSize', afs)
title(['\itr\rm = ' num2str(r,'%.2f')], 'FontSize', tfs)


set(axes, 'Units', 'Normalized', 'Position', [0.1+4*(asize+0.02)+0.084 0.2 asize asize*pw/ph], 'FontName', 'Helvetica')
[s sd] = getML_saccade('ZZTRain_psy.txt', 0);
L  = true(size(Ll_zz));
x  = LIP_zz(4,L);
xs = LIPd_zz(4,L);
xs(xs==0) = nanmean(xs);
xd = [x-LIPd_zz(4,L); LIPd_zz(4,L)+x];
y  = s(sl_zz(L),6)';
ys = sd(sl_zz(L),6)';
yd = sd(sl_zz(L),6)';
[b, bint, p, r, rp] = getLinearDependence(x,y);
hold on
plot([y-yd; y+yd], [x; x], 'LineWidth', elw, 'Color', lc1e)
plot([y; y], xd, 'LineWidth', elw, 'Color', lc1e)
plot(y, x, 'ko', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerSize', ms-1)
hold off
set(gca, 'YLim', [-2 10], 'YTick', [0 5 10], 'YTickLabel', [], ...
          'XLim', [-2 2], 'XTick', [-2 0 2], 'FontSize', afs)
title(['\itr\rm = ' num2str(r,'%.2f')], 'FontSize', tfs)
xlabel(sprintf('Saccade\naccuracy (deg)'), 'FontSize', lfs)



%% Motivation
% time to attain fixation
set(axes, 'Units', 'Normalized', 'Position', [0.1+5*(asize+0.02)+0.1 0.6 asize asize*pw/ph], 'FontName', 'Helvetica')
[s sd] = getML_saccade('CyTRain_psy.txt', 0);
L  = Ll_cy;
x  = LIP_cy(4,L);
xs = LIPd_cy(4,L);
xs(xs==0) = nanmean(xs);
xd = [x-LIPd_cy(4,L); LIPd_cy(4,L)+x];
y  = s(sl_cy(L),9)'/1000;
ys = sd(sl_cy(L),9)'/1000;
yd = sd(sl_cy(L),9)'/1000;
[b, bint, p, r, rp] = getLinearDependence(x,y);
hold on
plot([y-yd; y+yd], [x; x], 'LineWidth', elw, 'Color', lc1e)
plot([y; y], xd, 'LineWidth', elw, 'Color', lc1e)
plot(y, x, 'ko', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerSize', ms-1)
hold off
set(gca, 'YLim', [-4 18], 'YTick', [0 9 18], 'YTickLabel', [], ...
          'XLim', [0.5 1.5], 'XTick', [0.5 1 1.5], 'XTickLabel', {'0.5', '1.0', '1.5'}, 'FontSize', afs)
title(['\itr\rm = ' num2str(r,'%.2f') '  *'], 'FontSize', tfs)


set(axes, 'Units', 'Normalized', 'Position', [0.1+5*(asize+0.02)+0.1 0.2 asize asize*pw/ph], 'FontName', 'Helvetica')
[s sd] = getML_saccade('ZZTRain_psy.txt', 0);
L  = Ll_zz;
x  = LIP_zz(4,L);
xs = LIPd_zz(4,L);
xs(xs==0) = nanmean(xs);
xd = [x-LIPd_zz(4,L); LIPd_zz(4,L)+x];
y  = s(sl_zz(L),9)'/1000;
ys = sd(sl_zz(L),9)'/1000;
yd = sd(sl_zz(L),9)'/1000;
[b, bint, p, r, rp] = getLinearDependence(x,y);
hold on
plot([y-yd; y+yd], [x; x], 'LineWidth', elw, 'Color', lc1e)
plot([y; y], xd, 'LineWidth', elw, 'Color', lc1e)
plot(y, x, 'ko', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerSize', ms-1)
hold off
set(gca, 'YLim', [-2 10], 'YTick', [0 5 10], 'YTickLabel', [], ...
          'XLim', [0.5 1.5], 'XTick', [0.5 1 1.5], 'XTickLabel', {'0.5', '1.0', '1.5'}, 'FontSize', afs)
title(['\itr\rm = ' num2str(r,'%.2f') '  *'], 'FontSize', tfs)
xlabel(sprintf('Time to attain\nfixation (s)'), 'FontSize', lfs)



% percent bf or nc
set(axes, 'Units', 'Normalized', 'Position', [0.1+6*(asize+0.02)+0.1 0.6 asize asize*pw/ph], 'FontName', 'Helvetica')
[ts, ps, tse, pse] = getML_BehResponseStat('Cy',0);
L  = Ll_cy;
x  = LIP_cy(4,L);
xs = LIPd_cy(4,L);
xs(xs==0) = nanmean(xs);
xd = [x-LIPd_cy(4,L); LIPd_cy(4,L)+x];
y  = ts(sl_cy(L),4)';
yd = tse(sl_cy(L),4)';
[b, bint, p, r, rp] = getLinearDependence(x,y);
hold on
plot([y-yd; y+yd], [x; x], 'LineWidth', elw, 'Color', lc1e)
plot([y; y], xd, 'LineWidth', elw, 'Color', lc1e)
plot(y, x, 'ko', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerSize', ms-1)
hold off
title(['\itr\rm = ' num2str(r, '%.2f') '  *'], 'FontSize', tfs)
set(gca, 'YLim', [-4 18], 'YTick', [0 9 18], 'YTickLabel', [], ...
           'XLim', [0 0.4], 'XTick', [0 0.2 0.4], 'FontSize', afs)
       
set(axes, 'Units', 'Normalized', 'Position', [0.1+6*(asize+0.02)+0.1 0.2 asize asize*pw/ph], 'FontName', 'Helvetica')
[ts, ps, tse, pse] = getML_BehResponseStat('ZZ',0);
L  = true(size(Ll_zz));
x  = LIP_zz(4,L);
xs = LIPd_zz(4,L);
xs(xs==0) = nanmean(xs);
xd = [x-LIPd_zz(4,L); LIPd_zz(4,L)+x];
y  = ts(sl_zz(L),4)';
yd = tse(sl_zz(L),4)';
[b, bint, p, r, rp] = getLinearDependence(x,y);
hold on
plot([y-yd; y+yd], [x; x], 'LineWidth', elw, 'Color', lc1e)
plot([y; y], xd, 'LineWidth', elw, 'Color', lc1e)
plot(y, x, 'ko', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerSize', ms-1)
hold off
title(['\itr\rm = ' num2str(r, '%.2f') '^\clubsuit' '*'], 'FontSize', tfs)
set(gca, 'YLim', [-2 10], 'YTick', [0 5 10], 'YTickLabel', [], ...
           'XLim', [0 0.4], 'XTick', [0 0.2 0.4], 'FontSize', afs)
xlabel(sprintf('Percentage\nbroken fixation'), 'FontSize', lfs)



% # of reward
set(axes, 'Units', 'Normalized', 'Position', [0.1+7*(asize+0.02)+0.1 0.6 asize asize*pw/ph], 'FontName', 'Helvetica')
[re red] = getML_numReward('CyTRain_psy.txt', 0);
L  = Ll_cy;
x  = LIP_cy(4,L);
xs = LIPd_cy(4,L);
xs(xs==0) = nanmean(xs);
xd = [x-LIPd_cy(4,L); LIPd_cy(4,L)+x];
y  = re(sl_cy(L))';
ys = red(sl_cy(L)')';
yd = red(sl_cy(L)')';
[b, bint, p, r, rp] = getLinearDependence(x,y);
hold on
plot([y-yd; y+yd], [x; x], 'LineWidth', elw, 'Color', lc1e)
plot([y; y], xd, 'LineWidth', elw, 'Color', lc1e)
plot(y, x, 'ko', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerSize', ms-1)
hold off
set(gca, 'YLim', [-4 18], 'YTick', [0 9 18], 'YTickLabel', [], ...
          'XLim', [1 3], 'XTick', [1 2 3], 'FontSize', afs)
title(['\itr\rm = ' num2str(r,'%.2f') '  *'], 'FontSize', tfs)


set(axes, 'Units', 'Normalized', 'Position', [0.1+7*(asize+0.02)+0.1 0.2 asize asize*pw/ph], 'FontName', 'Helvetica')
[re red] = getML_numReward('ZZTRain_psy.txt', 0);
L  = Ll_zz;
x  = LIP_zz(4,L);
xs = LIPd_zz(4,L);
xs(xs==0) = nanmean(xs);
xd = [x-LIPd_zz(4,L); LIPd_zz(4,L)+x];
y  = re(sl_zz(L))';
ys = red(sl_zz(L)')';
yd = red(sl_zz(L)')';
[b, bint, p, r, rp] = getLinearDependence(x,y);
hold on
plot([y-yd; y+yd], [x; x], 'LineWidth', elw, 'Color', lc1e)
plot([y; y], xd, 'LineWidth', elw, 'Color', lc1e)
plot(y, x, 'ko', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerSize', ms-1)
hold off
set(gca, 'YLim', [-2 10], 'YTick', [0 5 10], 'YTickLabel', [], ...
          'XLim', [1 3], 'XTick', [1 2 3], 'FontSize', afs)
title(['\itr\rm = ' num2str(r,'%.2f') '  *'], 'FontSize', tfs)
xlabel(sprintf('Number of\nreward per trial'), 'FontSize', lfs)




