%% plotML_fig6.m
% Matlab script for plotting Fig 6  (Decision variable extracted from behavioral, MT and LIP data)
% Final version for publication to Nature Neuroscience
%
%   a) Model schematic (not plotted in Matlab)
%   b) DVs for behavioral, MT and LIP data
%   c) correlations of behavioral DV with MT and LIP DVs
%
%

%% get decision variables for behavior, MT and LIP responses
% using the Gold Shadlen diffusion model

% behavior
a               = getML_txt('CyTRain_psy.txt');
s_cy            = a.data{strcmp(a.name, 'session')};   
a               = getML_txt('ZZTRain_psy.txt');
s_zz            = a.data{strcmp(a.name, 'session')};   
[DV_cy, DVd_cy] = getML_psyPerformanceGSDV('CyTRain_psy.txt', 0);
[DV_zz, DVd_zz] = getML_psyPerformanceGSDV('ZZTRain_psy.txt', 0);
% remove sessions with high behavioral lapse
Lp_cy = true(size(getML_LapseSelectionArray('Cy', 'psy')));
Lp_zz = true(size(getML_LapseSelectionArray('ZZ', 'psy')));
[pt_cy nt_cy ptci_cy]     = getML_performanceAsso('CyTRain_psy.txt', 0);
[pt_zz nt_zz ptci_zz]     = getML_performanceAsso('ZZTRain_psy.txt', 0);



% MT
a               = getML_txt('CyTRain_MT.txt');
sm_cy           = a.data{strcmp(a.name, 'session')};   
a               = getML_txt('ZZTRain_MT.txt');
sm_zz           = a.data{strcmp(a.name, 'session')};   
[MT_cy, fano]   = getML_SigNoise('CyTRain_MT.txt', 0, [0 1], 0, 0);
[MT_zz, fano]   = getML_SigNoise('ZZTRain_MT.txt', 0, [0 1], 0, 0);
% remove sessions with high behavioral lapse
Lm_cy = true(size(getML_LapseSelectionArray('Cy', 'MT')));
Lm_zz = true(size(getML_LapseSelectionArray('ZZ', 'MT')));



% LIP
a                    = getML_txt('CyTRain_LIP.txt');
sl_cy                = a.data{strcmp(a.name, 'session')};   
a                    = getML_txt('ZZTRain_LIP.txt');
sl_zz                = a.data{strcmp(a.name, 'session')};   
A_cy                 = getML_LIPModelroc('Cy', 0);
A_zz                 = getML_LIPModelroc('ZZ', 0);
[LIP_cy, LIPd_cy, p] = getML_LIPModel('Cy', 0, 1000*A_cy(3,:)');
[LIP_zz, LIPd_zz, p] = getML_LIPModel('ZZ', 0, 1000*A_zz(3,:)');
% remove sessions with high behavioral lapse
Ll_cy = true(size(getML_LapseSelectionArray('Cy', 'LIP')));
Ll_zz = true(size(getML_LapseSelectionArray('ZZ', 'LIP')));


%%
%%%%%%%%%%%%%%%%%%%%%%%%%% plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create figure
fh = figure; 
mm2inch = 25.4;
pw = 85/mm2inch;    % paper width, make it nature size of one column, 85mm x 85mm
ph = 110/mm2inch;    % paper height
wysifig(fh, pw, ph) % set it to US letter size

% define figure parameters
alw = 0.5;  % axis line width
afs = 5;    % axis font size
lfs = 5;    % label font size
ffs = 12;   % figure font size
tfs = 6;    % text font size
lw  = 1.5; % line width
elw = 0.5;
ms  = 2;    % marker size

% define line color
lc1  = [0 0 0];
lc1e = [0.75 0.75 0.75];
lc2  = [1 0.3 0.3];
lc2e = [1 0.75 0.75];
lc3  = [0 1 1];
lc3e = [0.5 1 1];

asize = 0.33;
th = annotation('textbox', [0    0.95 0.05 0.05]);
set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'String', 'a', 'FontSize', ffs, 'FontWeight', 'bold', 'FontName', 'Helvetica') 

th = annotation('textbox', [0 1-pw/2/ph-0.06 0.05 0.05]);
set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'String', 'b', 'FontSize', ffs, 'FontWeight', 'bold', 'FontName', 'Helvetica') 

set(axes, 'Units', 'Normalized', 'Position', [0.105 1-pw/2/ph-asize*pw/ph-0.06 asize asize*pw/ph], 'FontName', 'Helvetica')
hold on
% plot MT
line(repmat(sm_cy',2,1), [(MT_cy(1,:)-MT_cy(2,:)); (MT_cy(1,:)+MT_cy(2,:))],...
                                'Color', lc3e, 'LineWidth', elw)
plot(sm_cy(~Lm_cy), MT_cy(1,~Lm_cy), 'o', 'MarkerSize', ms, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', lc3)
plot(sm_cy(Lm_cy), MT_cy(1,Lm_cy), 'o', 'MarkerSize', ms, 'MarkerFaceColor', lc3, 'MarkerEdgeColor', 'none')
MT_cy(2,MT_cy(2,:)==0) = nanmean(MT_cy(2,MT_cy(2,:)~=0));
L = ~isnan(MT_cy(1,:))' & Lm_cy; % plot regression lines
[b be h pmc] = nestedFW(MT_cy(1,L)', ones(size(MT_cy(2,L)')), [ones(sum(L),1) sm_cy(L)]);
%[b be r ri stats] = regress(MT_cy(1,L)', [ones(sum(L),1) sm_cy(L)]);
%pmc = stats(3);
line([1:nanmax(s_cy)], b(1)+b(2)*[1:nanmax(s_cy)], 'Color', lc3, 'LineWidth', lw)

% plot LIP
line(repmat(sl_cy',2,1), shiftdim(1000*LIPd_cy(3,:,:),1), 'Color', lc2e, 'LineWidth', elw)
plot(sl_cy(~Ll_cy), 1000*LIP_cy(3,~Ll_cy), 'o', 'MarkerSize', ms, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', lc2)
plot(sl_cy(Ll_cy), 1000*LIP_cy(3,Ll_cy), 'o', 'MarkerSize', ms, 'MarkerFaceColor', lc2, 'MarkerEdgeColor', 'none')
L = ~isnan(1000*LIP_cy(3,:)) & Ll_cy'; % plot regression lines
LIPs_cy = mean([LIP_cy(3,:)-shiftdim(LIPd_cy(3,1,:))'; shiftdim(LIPd_cy(3,2,:))'-LIP_cy(3,:)]);
LIPs_cy(LIPs_cy==0) = nanmean(LIPs_cy);
[b be h plc] = nestedFW(1000*LIP_cy(3,L)', 1000*ones(size(LIPs_cy(L)')), [ones(sum(L),1) sl_cy(L)]);
[b be r ri stats] = regress(1000*LIP_cy(3,L)', [ones(sum(L),1) sl_cy(L)]);
plc = stats(3);
line([1:nanmax(s_cy)], b(1)+b(2)*[1:nanmax(s_cy)], 'Color', lc2, 'LineWidth', lw)

set(gca, 'XLim', [1, max(s_cy)], 'YAxisLocation', 'right', ...
          'ylim', [-100 500], 'YTick', [0 250 500], 'FontSize', afs) 
hold off



% plot behavior
set(axes, 'Units', 'Normalized', 'Position', [0.105 1-pw/2/ph-asize*pw/ph-0.06 asize asize*pw/ph], 'FontName', 'Helvetica')
hold on
line(repmat(s_cy',2,1), shiftdim(DVd_cy(1,:,:),1), 'Color', lc1e, 'LineWidth', elw)
plot(s_cy(~Lp_cy), DV_cy(1,~Lp_cy), 'o', 'MarkerSize', ms, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', lc1)
plot(s_cy(Lp_cy), DV_cy(1,Lp_cy), 'o', 'MarkerSize', ms, 'MarkerFaceColor', lc1, 'MarkerEdgeColor', 'none')
L = ~isnan(DV_cy(1,:))' & Lp_cy; % plot regression lines
DVs_cy = mean([DV_cy(1,:)-shiftdim(DVd_cy(1,1,:))'; shiftdim(DVd_cy(1,2,:))'-DV_cy(1,:)]);
DVs_cy(DVs_cy==0) = nanmean(DVs_cy);
[b be h pbc] = nestedFW(DV_cy(1,L)', ones(size(DVs_cy(L)')), [ones(sum(L),1) s_cy(L)]);
line([1:nanmax(s_cy)], b(1)+b(2)*[1:nanmax(s_cy)], 'Color', lc1, 'LineWidth', lw)
text(5, 36, 'Monkey C', 'FontSize', tfs, ...
     'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')
hold off
set(gca, 'Color', 'none', 'XLim', [1, max(s_cy)], 'XTick', [1 50 100 150], ...
          'YAxisLocation', 'left', 'YLim', [0 36], 'YTick', [0 18 36], 'FontSize', afs) 
xlabel('Session (d)', 'FontSize', lfs)
ylabel('\ita\rm_{be} (spikes per s^2 per coh)', 'FontSize', lfs)

% text(-40, -2, '\ita\rm (Behavior, spikes/s^2/coh)', ...
%       'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
%       'FontSize', lfs, 'Rotation', 90)
% text(200, -2, '\ita\rm (\color[rgb]{0.3 0.3 1}MT\color[rgb]{0 0 0} & \color[rgb]{1 0.3 0.3}LIP\color[rgb]{0 0 0}, spikes/s^2/coh)', ...
%       'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
%       'FontSize', lfs, 'Rotation', 90)



% ZZ
set(axes, 'Units', 'Normalized', 'Position', [0.55 1-pw/2/ph-asize*pw/ph-0.06 asize asize*pw/ph], 'FontName', 'Helvetica')
hold on
% plot MT
line(repmat(sm_zz',2,1), [(MT_zz(1,:)-MT_zz(2,:)); (MT_zz(1,:)+MT_zz(2,:))],...
                                'Color', lc3e, 'LineWidth', elw)
plot(sm_zz(~Lm_zz), MT_zz(1,~Lm_zz), 'o', 'MarkerSize', ms, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', lc3)
plot(sm_zz(Lm_zz), MT_zz(1,Lm_zz), 'o', 'MarkerSize', ms, 'MarkerFaceColor', lc3, 'MarkerEdgeColor', 'none')
MT_zz(2,MT_zz(2,:)==0) = nanmean(MT_zz(2,MT_zz(2,:)~=0));
L = ~isnan(MT_zz(1,:))' & Lm_zz; % plot regression lines
[b be h pmz] = nestedFW(MT_zz(1,L)', ones(size(MT_zz(2,L)')), [ones(sum(L),1) sm_zz(L)]);
%[b be r ri stats] = regress(MT_zz(1,L)', [ones(sum(L),1) sm_zz(L)]);
%pmz = stats(3);
line([1:nanmax(s_zz)], b(1)+b(2)*[1:nanmax(s_zz)], 'Color', lc3, 'LineWidth', lw)

% plot LIP
line(repmat(sl_zz',2,1), shiftdim(1000*LIPd_zz(3,:,:),1), 'Color', lc2e, 'LineWidth', elw)
plot(sl_zz(~Ll_zz), 1000*LIP_zz(3,~Ll_zz), 'o', 'MarkerSize', ms, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', lc2)
plot(sl_zz(Ll_zz), 1000*LIP_zz(3,Ll_zz), 'o', 'MarkerSize', ms, 'MarkerFaceColor', lc2, 'MarkerEdgeColor', 'none')
L = ~isnan(1000*LIP_zz(3,:)) & Ll_zz'; % plot regression lines
LIPs_zz = mean([LIP_zz(3,:)-shiftdim(LIPd_zz(3,1,:))'; shiftdim(LIPd_zz(3,2,:))'-LIP_zz(3,:)]);
LIPs_zz(LIPs_zz==0) = nanmean(LIPs_zz);
[b be h plz] = nestedFW(1000*LIP_zz(3,L)', 1000*ones(size(LIPs_zz(L)')), [ones(sum(L),1) sl_zz(L)]);
%[b be r ri stats] = regress(1000*LIP_zz(3,L)', [ones(sum(L),1) sl_zz(L)]);
%plz = stats(3);
line([1:nanmax(s_zz)], b(1)+b(2)*[1:nanmax(s_zz)], 'Color', lc2, 'LineWidth', lw)

set(gca, 'XLim', [1, max(s_zz)], 'XTick', [], 'YAxisLocation', 'right', ...
            'YLim', [-50 350], 'YTick', [0 175 350], 'FontSize', afs) 
hold off
ylabel('\color[rgb]{0 1 1}\ita\rm_{MT}\color[rgb]{0 0 0} & \color[rgb]{1 0.3 0.3}\ita\rm_{LIP}\color[rgb]{0 0 0} (spikes per s^2 per coh)', 'FontSize', lfs)


% plot behavior
set(axes, 'Units', 'Normalized', 'Position', [0.55 1-pw/2/ph-asize*pw/ph-0.06 asize asize*pw/ph], 'FontName', 'Helvetica')
hold on
line(repmat(s_zz',2,1), shiftdim(DVd_zz(1,:,:),1), 'Color', lc1e, 'LineWidth', elw)
plot(s_zz(~Lp_zz), DV_zz(1,~Lp_zz), 'o', 'MarkerSize', ms, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', [0 0 0])
plot(s_zz(Lp_zz), DV_zz(1,Lp_zz), 'o', 'MarkerSize', ms, 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', 'none')
L = ~isnan(DV_zz(1,:))' & Lp_zz; % plot regression lines
DVs_zz = mean([DV_zz(1,:)-shiftdim(DVd_zz(1,1,:))'; shiftdim(DVd_zz(1,2,:))'-DV_zz(1,:)]);
DVs_zz(DVs_zz==0) = nanmean(DVs_zz);
[b be h pbz] = nestedFW(DV_zz(1,L)', ones(size(DVs_zz(L)')), [ones(sum(L),1) s_zz(L)]);
line([1:nanmax(s_zz)], b(1)+b(2)*[1:nanmax(s_zz)], 'Color', [0 0 0], 'LineWidth', lw)
text(5, 20, 'Monkey Z', 'FontSize', tfs, ...
     'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')
hold off
set(gca, 'Color', 'none', 'XLim', [1, max(s_zz)], 'XTick', [1 40 80 120], ...
          'YAxisLocation', 'left', 'YLim', [0 20], 'YTick', [0 10 20], 'FontSize', afs) 
xlabel('Session (d)', 'FontSize', lfs)



% get correlation between MT, LIP with behavior
% cyrus
% align MT and LIP DV with DV from behavior
L  = ~isnan(MT_cy(1,:)) & Lm_cy';
MT = MT_cy(1,L);
Be = DV_cy(1,sm_cy(L));
[RBMc, PBMc] = corr(MT', Be'); 


L   = ~isnan(LIP_cy(3,:)) & LIP_cy(3,:)<1 & Ll_cy';
L(1) = 0;
LIP = LIP_cy(3,L);
Be  = DV_cy(1,sl_cy(L));
[RBLc, PBLc] = corr(LIP', Be'); 



% zsa zsa
% align MT and LIP DV with DV from behavior
L  = ~isnan(MT_zz(1,:)) & Lm_zz';
MT = MT_zz(1,L);
Be = DV_zz(1,sm_zz(L));
[RBMz, PBMz] = corr(MT', Be'); 


L   = ~isnan(LIP_zz(3,:)) & Ll_zz';
LIP = LIP_zz(3,L);
Be  = DV_zz(1,sl_zz(L));
[RBLz, PBLz] = corr(LIP', Be'); 



%% plot correlation of MT, LIP DV with behavior
th = annotation('textbox', [0 0.2 0.05 0.05]);
set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'String', 'c', 'FontSize', ffs, 'FontWeight', 'bold', 'FontName', 'Helvetica') 


% cy
set(axes, 'Units', 'Normalized', 'Position', [0.105 0.08 0.4*asize 0.4*asize*pw/ph], 'FontName', 'Helvetica')
LIP  = 1000*LIP_cy(3,:)';
LIPd = 1000*shiftdim(LIPd_cy(3,:,:));
Be   = DV_cy(1,sl_cy(:))';
Bed  = shiftdim(DVd_cy(1,:,sl_cy(:)),1);
%plot
hold on
plot(LIPd, [Be'; Be'], 'LineWidth', elw, 'Color', lc1e)
plot([LIP'; LIP'], Bed, 'LineWidth', elw, 'Color', lc1e)
plot(LIP, Be, 'ko', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerSize', ms-0.5)
hold off
Lbd = isnan(LIP) | isnan(Be);
[rbl, pbl] = corr(LIP(~Lbd),Be(~Lbd));
title(['\itr\rm = ' num2str(rbl, '%.2f') '^\clubsuit'], 'FontSize', tfs)
set(gca, 'YLim', [0 28], 'YTick', [0 14 28], ...
           'XLim', [-50 500], 'XTick', [0 250 500], 'XTickLabel', [0 250 500], ...
           'FontSize', afs)
xlabel('\ita\rm_{LIP}', 'FontSize', lfs)
ylabel('\ita\rm_{be}', 'FontSize', lfs)


set(axes, 'Units', 'Normalized', 'Position', [0.105+0.4*asize+0.2*asize 0.08 0.4*asize 0.4*asize*pw/ph], 'FontName', 'Helvetica')
MT   = MT_cy(1,:)';
MTd  = [(MT_cy(1,:)-MT_cy(2,:)); (MT_cy(1,:)+MT_cy(2,:))];
Be   = DV_cy(1,sm_cy(:))';
Bed  = shiftdim(DVd_cy(1,:,sm_cy(:)),1);
%plot
hold on
plot(MTd, [Be'; Be'], 'LineWidth', elw, 'Color', lc1e)
plot([MT'; MT'], Bed, 'LineWidth', elw, 'Color', lc1e)
plot(MT, Be, 'ko', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerSize', ms-0.5)
hold off
Lbd = isnan(MT) | isnan(Be);
[rbm, pbm] = corr(MT(~Lbd),Be(~Lbd));
title(['\itr\rm = ' num2str(rbm, '%.2f')], 'FontSize', tfs)
set(gca, 'YLim', [0 28], 'YTick', [0 14 28], 'YTickLabel', [],...
           'XLim', [-50 100], 'XTick', [-50 0 50 100], 'XTickLabel', [-50 0 50 100], ...
           'FontSize', afs)
xlabel('\ita\rm_{MT}', 'FontSize', lfs)




% zz
set(axes, 'Units', 'Normalized', 'Position', [0.55 0.08 0.4*asize 0.4*asize*pw/ph], 'FontName', 'Helvetica')
LIP  = 1000*LIP_zz(3,:)';
LIPd = 1000*shiftdim(LIPd_zz(3,:,:));
Be   = DV_zz(1,sl_zz(:))';
Bed  = shiftdim(DVd_zz(1,:,sl_zz(:)),1);
%plot
hold on
plot(LIPd, [Be'; Be'], 'LineWidth', elw, 'Color', lc1e)
plot([LIP'; LIP'], Bed, 'LineWidth', elw, 'Color', lc1e)
plot(LIP, Be, 'ko', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerSize', ms-0.5)
hold off
Lbd = isnan(LIP) | isnan(Be);
[rbl, pbl] = corr(LIP(~Lbd),Be(~Lbd));
title(['\itr\rm = ' num2str(rbl, '%.2f') '^\clubsuit'], 'FontSize', tfs)
set(gca, 'YLim', [0 20], 'YTick', [0 10 20], ...
           'XLim', [-50 350], 'XTick', [0 175 350], 'XTickLabel', [0 175 350], ...
           'FontSize', afs)
xlabel('\ita\rm_{LIP}', 'FontSize', lfs)
ylabel('\ita\rm_{be}', 'FontSize', lfs)


set(axes, 'Units', 'Normalized', 'Position', [0.55+0.4*asize+0.2*asize 0.08 0.4*asize 0.4*asize*pw/ph], 'FontName', 'Helvetica')
MT   = MT_zz(1,:)';
MTd  = [(MT_zz(1,:)-MT_zz(2,:)); (MT_zz(1,:)+MT_zz(2,:))];
Be   = DV_zz(1,sm_zz(:))';
Bed  = shiftdim(DVd_zz(1,:,sm_zz(:)),1);
%plot
hold on
plot(MTd, [Be'; Be'], 'LineWidth', elw, 'Color', lc1e)
plot([MT'; MT'], Bed, 'LineWidth', elw, 'Color', lc1e)
plot(MT, Be, 'ko', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerSize', ms-0.5)
hold off
Lbd = isnan(MT) | isnan(Be);
[rbm, pbm] = corr(MT(~Lbd),Be(~Lbd));
title(['\itr\rm = ' num2str(rbm, '%.2f')], 'FontSize', tfs)
set(gca, 'YLim', [0 20], 'YTick', [0 10 20], 'YTickLabel', [],...
           'XLim', [-50 100], 'XTick', [-50 0 50 100], 'XTickLabel', [-50 0 50 100], ...
           'FontSize', afs)
xlabel('\ita\rm_{MT}', 'FontSize', lfs)








