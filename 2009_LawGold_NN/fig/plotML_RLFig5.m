% plot reinforcement learning figure 5 - pooled SNR
%
% Matlab script for plotting Fig 5
% for submission to Nature Neuroscience
%
% created by jcl on 08/07/08

%%
% create figure
fh = figure; 
mm2inch = 25.4;
pw = 120/mm2inch;    % paper width, make it nature size of one column, 85mm x 85mm
ph = 120/mm2inch;   % paper height
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
lc1  = [0.2 0.2 0.2];
lc1e = [0.75 0.75 0.75];
lc2  = [1 0.2 0.2];
lc2e = [1 0.75 0.75];


asize = 0.23;

%%
% get data
Monk      = 'Cy';
NSEN      = 200;
NTUNE     = 36;
ELN       = [7]*1e-8;
DIRS      = [180 0];
SIMNUM    = [81];
COHS      = [0 3.2 6.4 12.8 25.6 51.2 99.9]/100;
TIMES     = [0:0.1:1];
recompute = 0;

[mpool] = getML_RLFig5PooledCT(Monk, NSEN, NTUNE, ELN, SIMNUM, COHS, TIMES, DIRS, recompute);
mpool = nanmean(mpool,5);


%% plot
th = annotation('textbox', [0.1-0.1 0.96 0.05 0.05]);
set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'String', 'a', 'FontSize', ffs, 'FontWeight', 'bold', 'FontName', 'Helvetica') 

for i = 1:3
    % align to dots onset
    ah(i) = axes;
    set(ah(i), 'Units', 'Normalized', 'Position', [0.15+(i-1)*(asize+0.03) 0.71 0.9*asize*ph/pw asize], 'FontName', 'Helvetica')
    mp = mpool(:,:,1,i);
    mn = mpool(:,:,2,i);
    mp(1,:) = nan;
    mn(1,:) = nan;
    
    plotML_rateByDirCoh(TIMES, [mn; mp], gca, lw)
    set(ah(i), 'LineWidth', alw, 'FontSize', afs, 'Box', 'off', ...
        'XLim', [0 1], 'XTickLabel', [0 0.5 1])
    xlabel([])
    ylabel([])
    
    if i==1
        set(gca, 'YLim', [-100 100], 'YTick', [-100 0 100], ...
            'XLim', [0 1], 'XTick', [0 0.5 1])
        ylabel('Response (spikes per s)', 'FontSize', lfs)
        title('Early', 'FontSize', tfs)
    elseif i==2       
        set(gca, 'YLim', [-100 100], 'YTick', [-100 0 100], 'YTickLabel', [], ...
            'XLim', [0 1], 'XTick', [0 0.5 1])
        title('Mid', 'FontSize', tfs)
        xlabel('Time (s)', 'FontSize', lfs)
       
    else       
        set(gca, 'YLim', [-100 100], 'YTick', [-100 0 100], 'YTickLabel', [], ...
            'XLim', [0 1], 'XTick', [0 0.5 1])
        title('Late', 'FontSize', tfs)  
    end
end



%% get data
% get average SNR of pooled signal for both monkeys
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
[mo1_c, sdo1_c] = getML_RLFig5PooledSNROptimal(Monk, NSEN, NTUNE, ELN, SIMNUM, COH, TIME, recompute);

COH = 0.128;
[m2_c, sd2_c] = getML_RLFig5PooledSNR(Monk, NSEN, NTUNE, ELN, SIMNUM, COH, TIME, recompute);
[mo2_c, sdo2_c] = getML_RLFig5PooledSNROptimal(Monk, NSEN, NTUNE, ELN, SIMNUM, COH, TIME, recompute);


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





% model
Monk      = 'ZZ';
NSEN      = 200;
NTUNE     = 36;
ELN       = [0.75]*1e-8;
DIRS      = [0 180];
SIMNUM    = [81];

COH = 0.999;
TIME = 0.6;
[m1_z, sd1_z] = getML_RLFig5PooledSNR(Monk, NSEN, NTUNE, ELN, SIMNUM, COH, TIME, recompute);
[mo1_z, sdo1_z] = getML_RLFig5PooledSNROptimal(Monk, NSEN, NTUNE, ELN, SIMNUM, COH, TIME, recompute);

COH = 0.256;
[m2_z, sd2_z] = getML_RLFig5PooledSNR(Monk, NSEN, NTUNE, ELN, SIMNUM, COH, TIME, recompute);
[mo2_z, sdo2_z] = getML_RLFig5PooledSNROptimal(Monk, NSEN, NTUNE, ELN, SIMNUM, COH, TIME, recompute);


% data
fname     = 'ZZTRain_LIP.txt';
crtf      = [0 1];
nvf       = [0 1];
dirf      = [0];
bb        = 600;
be        = 900;
bw        = 300;
bs        = bw;
bins      = [[bb:bs:be-bw]' [bb+bw:bs:be]'];
[LSNR_z, LTR_z] = getML_RLFig5LIPSNR(fname, crtf, nvf, dirf, bins, 0);
lsnr1_z   = squeeze(LSNR_z(:,7,:));
lsnr2_z   = squeeze(LSNR_z(:,5,:));


        
    

%% plot SNR of pooled signal for coh 1
th = annotation('textbox', [0.1-0.1 0.58 0.05 0.05]);
set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'String', 'b', 'FontSize', ffs, 'FontWeight', 'bold', 'FontName', 'Helvetica') 


Nc = 115985;
Nz = 77447;


set(axes, 'Units', 'Normalized', 'Position', [0.22 0.08 asize asize*pw/ph], 'FontName', 'Helvetica')
tr   = [1:1000:Nc Nc];
snr  = m2_c./(sd2_c+sqrt(mean(m2_c))+5);
snro = mo2_c./(sdo2_c+sqrt(mo2_c)+5);
hold on
plot(tr, snr, 'o', 'MarkerFaceColor', lc1, 'MarkerEdgeColor', 'none', 'MarkerSize', ms)
plot(LTR_c, lsnr2_c, 'o', 'MarkerFaceColor', lc2, 'MarkerEdgeColor', 'none', 'MarkerSize', ms)
plot(tr, snro*ones(size(tr)), ':', 'Color', 'k', 'LineWidth', 1)
hold off
% set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, ...
%          'xlim', [0 Nc], 'ylim', [-0.5 1.2], 'YTick', [-0.5 0 0.5 1], ...
%          'XTick', [0 0.5 1]*10e4, 'XTickLabel', {'0' '50' '100'}) 
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, ...
         'xlim', [0 Nc], 'ylim', [-0.2 3.2], 'YTick', [0 1.5 3], 'YTickLabel', {'0', '1.5', '3'},...
         'XTick', [0 0.5 1]*10e4, 'XTickLabel', {'0' '50' '100'}) 
%xlabel('Trial (\times 1000)', 'FontSize', lfs)
ylabel('SNR', 'FontSize', lfs)
xlabel('Trial (\times 1000)', 'FontSize', lfs)
text(0.03*Nc,3.1, '12.8% coh', 'FontSize', lfs, 'HorizontalAlignment' , 'left', 'VerticalAlignment', 'top')


set(axes, 'Units', 'Normalized', 'Position', [0.22 0.08+0.05+asize asize asize*pw/ph], 'FontName', 'Helvetica')
tr   = [1:1000:Nc Nc];
snr  = m1_c./(sd1_c+sqrt(mean(m1_c))+5);
snro = mo1_c./(sdo1_c+sqrt(mo1_c)+5);
hold on
plot(tr, snr, 'o', 'MarkerFaceColor', lc1, 'MarkerEdgeColor', 'none', 'MarkerSize', ms)
plot(LTR_c, lsnr1_c, 'o', 'MarkerFaceColor', lc2, 'MarkerEdgeColor', 'none', 'MarkerSize', ms)
plot(tr, snro*ones(size(tr)), ':', 'Color', 'k', 'LineWidth', 1)
hold off
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, ...
         'xlim', [0 Nc], 'ylim', [-0.2 3.2], 'YTick', [0 1.5 3], 'YTickLabel', {'0', '1.5', '3'},...
         'XTick', [0 0.5 1]*10e4, 'XTickLabel', {'0' '50' '100'}) 
%xlabel('Trial (\times 1000)', 'FontSize', lfs)
ylabel('SNR', 'FontSize', lfs)
title('Monkey C', 'FontSize', lfs)
text(0.03*Nc,3.1, '51.2% coh', 'FontSize', lfs, 'HorizontalAlignment' , 'left', 'VerticalAlignment', 'top')


set(axes, 'Units', 'Normalized', 'Position', [0.22+asize+0.08 0.08 asize asize*pw/ph], 'FontName', 'Helvetica')
tr   = [1:1000:Nz Nz];
snr  = m2_z./(sd2_z+sqrt(2*mean(m2_z))+5);
snro = mo2_z./(sdo2_z+sqrt(2*mo2_z)+5);
hold on
plot(tr, snr, 'o', 'MarkerFaceColor', lc1, 'MarkerEdgeColor', 'none', 'MarkerSize', ms)
plot(LTR_z, lsnr2_z, 'o', 'MarkerFaceColor', lc2, 'MarkerEdgeColor', 'none', 'MarkerSize', ms)
plot(tr, snro*ones(size(tr)), ':', 'Color', 'k', 'LineWidth', 1)
hold off
% set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, ...
%          'xlim', [0 Nz], 'ylim', [-0.3 1.22], 'YTick', [0 0.5 1], ...
%          'XTick', [0 0.35 0.7]*10e4, 'XTickLabel', {'0' '35' '70'})
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, ...
         'xlim', [0 Nz], 'ylim', [-0.2 3.2], 'YTick', [0 1.5 3], 'YTickLabel', {'0', '1.5', '3'},...
         'XTick', [0 0.35 0.7]*10e4, 'XTickLabel', {'0' '35' '70'})       
xlabel('Trial (\times 1000)', 'FontSize', lfs)
text(0.03*Nz,3.1, '25.6% coh', 'FontSize', lfs, 'HorizontalAlignment' , 'left', 'VerticalAlignment', 'top')


set(axes, 'Units', 'Normalized', 'Position', [0.22+asize+0.08 0.08+0.05+asize asize asize*pw/ph], 'FontName', 'Helvetica')
tr   = [1:1000:Nz Nz];
snr  = m1_z./(sd1_z+sqrt(2*mean(m1_z))+5);
snro = mo1_z./(sdo1_z+sqrt(2*mo1_z)+5);
hold on
plot(tr, snr, 'o', 'MarkerFaceColor', lc1, 'MarkerEdgeColor', 'none', 'MarkerSize', ms)
plot(LTR_z, lsnr1_z, 'o', 'MarkerFaceColor', lc2, 'MarkerEdgeColor', 'none', 'MarkerSize', ms)
%plot(tr, snro*ones(size(tr)), ':', 'Color', 'k', 'LineWidth', 1)
plot(tr, (3.2-0.05)*ones(size(tr)), ':', 'Color', 'k', 'LineWidth', 1)
hold off
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, ...
         'xlim', [0 Nz], 'ylim', [-0.2 3.2], 'YTick', [0 1.5 4], 'YTickLabel', {'0', '1.5', '3'},...
         'XTick', [0 0.35 0.7]*10e4, 'XTickLabel', {'0' '35' '70'})       
%ylabel('SNR ratio', 'FontSize', lfs)
title('Monkey Z', 'FontSize', lfs)
text(0.03*Nz,3.1, '99.9% coh', 'FontSize', lfs, 'HorizontalAlignment' , 'left', 'VerticalAlignment', 'top')


%% get r-square values for fits
% find model SNR corresponding to the session in which neural data was
% collected
tr_c    = [1:1000:Nc Nc];
snr1_c  = m1_c./(sd1_c+sqrt(2*mean(m1_c))+5);
snr2_c  = m2_c./(sd2_c+sqrt(2*mean(m2_c))+5);

snrm1_c = nans(size(LTR_c));
snrm2_c = nans(size(LTR_c));
for i = 1:length(LTR_c)
    if ~isnan(lsnr1_c(i))
        [xxx, iii] = nanmin(abs(tr_c-LTR_c(i)));
        snrm1_c(i) = snr1_c(iii);
    end
    
    if ~isnan(lsnr2_c(i))
        [xxx, iii] = nanmin(abs(tr_c-LTR_c(i)));
        snrm2_c(i) = snr2_c(iii);
    end          
end


tr_z    = [1:1000:Nz Nz];
snr1_z  = m1_z./(sd1_z+sqrt(2*mean(m1_z))+5);
snr2_z  = m2_z./(sd2_z+sqrt(2*mean(m2_z))+5);

snrm1_z = nans(size(LTR_z));
snrm2_z = nans(size(LTR_z));
for i = 1:length(LTR_z)
    if ~isnan(lsnr1_z(i))
        [xxx, iii] = nanmin(abs(tr_z-LTR_z(i)));
        snrm1_z(i) = snr1_z(iii);
    end
    
    if ~isnan(lsnr2_z(i))
        [xxx, iii] = nanmin(abs(tr_z-LTR_z(i)));
        snrm2_z(i) = snr2_z(iii);
    end      
end




figure
subplot(2,2,1)
hold on
plot(snrm1_c, lsnr1_c, '.k')
Lgd = ~isnan(snrm1_c);
[rrr ppp] = corr(snrm1_c(Lgd)', lsnr1_c(Lgd))
[h ppp] = ttest(snrm1_c(Lgd)', lsnr1_c(Lgd))
lsline
title(sprintf('r=%.2f, p=%.4f', rrr, ppp))
hold off

subplot(2,2,2)
plot(snrm2_c, lsnr2_c, '.k')
Lgd = ~isnan(snrm2_c);
[rrr ppp] = corr(snrm2_c(Lgd)', lsnr2_c(Lgd))
[h ppp] = ttest(snrm2_c(Lgd)', lsnr2_c(Lgd))

title(sprintf('r=%.2f, p=%.4f', rrr, ppp))
lsline  
    
subplot(2,2,3)
hold on
plot(snrm1_z, lsnr1_z, '.k')
Lgd = ~isnan(snrm1_z) & ~isinf(lsnr1_z');
[rrr ppp] = corr(snrm1_z(Lgd)', lsnr1_z(Lgd))
[h ppp] = ttest(snrm1_z(Lgd)', lsnr1_z(Lgd))

lsline
title(sprintf('r=%.2f, p=%.4f', rrr, ppp))
hold off

subplot(2,2,4)
plot(snrm2_z, lsnr2_z, '.k')
lsline  
Lgd = ~isnan(snrm2_z);
[rrr ppp] = corr(snrm2_z(Lgd)', lsnr2_z(Lgd))
[h ppp] = ttest(snrm2_z(Lgd)', lsnr2_z(Lgd))

title(sprintf('r=%.2f, p=%.4f', rrr, ppp))













