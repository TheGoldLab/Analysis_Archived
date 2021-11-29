%% plotML_figS1.m
% Matlab script for plotting Fig S1  (Distributions of MT and LIP receptive/response field)
% Final version for publication to Nature Neuroscience
%
%


monk     = 'Cy';
a        = getML_txt([monk 'TRain_psy.txt']);
a_cy     = a.data{strcmp(a.name,'a')};
th_cy    = a.data{strcmp(a.name,'th')};
d_cy     = a.data{strcmp(a.name,'d')};
trga_cy  = a.data{strcmp(a.name,'trg_a')};
trgth_cy = a.data{strcmp(a.name,'trg_dir')};
dd_cy    = a.data{strcmp(a.name,'ddir')};
s_cy     = a.data{strcmp(a.name,'s')};


a        = getML_txt([monk 'TRain_MT.txt']);
use      = a.data{strcmp(a.name,'usable')};
a_cy     = a.data{strcmp(a.name,'a')};  a_cy(use~=1)  = [];
th_cy    = a.data{strcmp(a.name,'th')}; th_cy(use~=1) = [];
d_cy     = a.data{strcmp(a.name,'d')};  d_cy(use~=1)  = [];
dd_cy    = a.data{strcmp(a.name,'ddir')};  dd_cy(use~=1)  = [];

a        = getML_txt([monk 'TRain_LIP.txt']);
use      = a.data{strcmp(a.name,'usable')};
trga_cy  = a.data{strcmp(a.name,'trg_a')};  trga_cy(use~=1)  = [];
trgth_cy = a.data{strcmp(a.name,'trg_dir')}; trgth_cy(use~=1) = [];

a        = getML_txt([monk 'PRe_MT.txt']);
use      = a.data{strcmp(a.name,'usable')};
a_cyp    = a.data{strcmp(a.name,'a')};  a_cyp(use~=1)  = [];
th_cyp   = a.data{strcmp(a.name,'th')}; th_cyp(use~=1) = [];
d_cyp    = a.data{strcmp(a.name,'d')};  d_cyp(use~=1)  = [];
dd_cyp   = a.data{strcmp(a.name,'ddir')};  dd_cyp(use~=1) = [];
s_cyp    = a.data{strcmp(a.name,'s')};


% transform dot direction to its direction axis, a value between -90 to 90
L        = dd_cy>90 & dd_cy<270;
dd_cy(L) = mod(dd_cy(L)+180,360);
L        = dd_cy>=180;
dd_cy(L) = -1*(360-dd_cy(L));

L        = trgth_cy>90 & trgth_cy<270;
trgth_cy(L) = mod(trgth_cy(L)+180,360);
L        = trgth_cy>=180;
trgth_cy(L) = -1*(360-trgth_cy(L));

L         = dd_cyp>90 & dd_cyp<270;
dd_cyp(L) = mod(dd_cyp(L)+180,360);
L         = dd_cyp>=180;
dd_cyp(L) = -1*(360-dd_cyp(L));


monk     = 'ZZ';
a        = getML_txt([monk 'TRain_psy.txt']);
a_zz     = a.data{strcmp(a.name,'a')};
th_zz    = a.data{strcmp(a.name,'th')};
d_zz     = a.data{strcmp(a.name,'d')};
trga_zz  = a.data{strcmp(a.name,'trg_a')};
trgth_zz = a.data{strcmp(a.name,'trg_dir')};
dd_zz    = a.data{strcmp(a.name,'ddir')};
s_zz     = a.data{strcmp(a.name,'s')};

a        = getML_txt([monk 'TRain_MT.txt']);
use      = a.data{strcmp(a.name,'usable')};
a_zz     = a.data{strcmp(a.name,'a')};  a_zz(use~=1)  = [];
th_zz    = a.data{strcmp(a.name,'th')}; th_zz(use~=1) = [];
d_zz     = a.data{strcmp(a.name,'d')};  d_zz(use~=1)  = [];
dd_zz    = a.data{strcmp(a.name,'ddir')};  dd_zz(use~=1) = [];

a        = getML_txt([monk 'TRain_LIP.txt']);
use      = a.data{strcmp(a.name,'usable')};
trga_zz  = a.data{strcmp(a.name,'trg_a')};  trga_zz(use~=1)  = [];
trgth_zz = a.data{strcmp(a.name,'trg_dir')}; trgth_zz(use~=1) = [];

a        = getML_txt([monk 'PRe_MT.txt']);
use      = a.data{strcmp(a.name,'usable')};
a_zzp    = a.data{strcmp(a.name,'a')};  a_zzp(use~=1)  = [];
th_zzp   = a.data{strcmp(a.name,'th')}; th_zzp(use~=1) = [];
d_zzp    = a.data{strcmp(a.name,'d')};  d_zzp(use~=1)  = [];
dd_zzp   = a.data{strcmp(a.name,'ddir')};  dd_zzp(use~=1)  = [];
s_zzp    = a.data{strcmp(a.name,'s')};

% transform dot direction to its direction axis, a value between -90 to 90
L        = dd_zz>90 & dd_zz<270;
dd_zz(L) = mod(dd_zz(L)+180,360);
L        = dd_zz>=180;
dd_zz(L) = -1*(360-dd_zz(L));

L        = trgth_zz>90 & trgth_zz<270;
trgth_zz(L) = mod(trgth_zz(L)+180,360);
L        = trgth_zz>=180;
trgth_zz(L) = -1*(360-trgth_zz(L));

L         = dd_zzp>90 & dd_zzp<270;
dd_zzp(L) = mod(dd_zzp(L)+180,360);
L         = dd_zzp>=180;
dd_zzp(L) = -1*(360-dd_zzp(L));



%%% plot MT RF, LIP response field distribution for both monkeys %%%
fh = figure;
pw = 8;  % paper width
ph = 10;   % paper height
wysifig(fh, pw, ph) % set it to US letter size

% define figure parameters
alw = 1;    % axis line width
afs = 14;   % axis font size
lfs = 18;   % label font size
ffs = 30;   % figure font size
tfs = 14;   % text font size
lw  = 2;    % line width
ms  = 5;    % marker size

asize = 0.24;

% th = annotation('textbox', [0.1-0.1 0.58+asize+0.018 0.05 0.05]);
% set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
%     'String', 'a', 'FontSize', ffs, 'FontWeight', 'bold') 

% cyrus
% plot RF distribution
ah = axes;
set(ah, 'Units', 'Normalized', 'Position', [0.18 0.4 asize asize*pw/ph]);
plotML_RFdist(a_cyp, th_cyp, d_cyp, 0, ah);
set(ah, 'FontSize', afs)
xl = get(gca, 'XLim');
yl = get(gca, 'YLim');
xlabel('Degree visual angle', 'FontSize', afs)
ylabel('Degree visual angle', 'FontSize', afs)
% plot direction distribution
apos  = get(gca, 'Position');
apos2 = [apos(1)+12*apos(3)/20 apos(2)+12.5*apos(4)/20 7*apos(3)/20 6*apos(4)/20];
set(axes, 'Units', 'Normalized', 'Position', apos2);
dd       = dd_cyp;
bins     = [-90:180/5:90];   % bin edges
n        = histc(dd,bins);
bh       = bar(bins, n, 0.99, 'histc');
set(bh,'LineWidth', 1, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', [1 1 1])
set(gca, 'XLim', [min(bins) max(bins)], 'XTick', [], 'YLim', [0 20], 'YTick', [0 10 20], ...
            'Layer', 'top', 'box', 'off')
yl = get(gca, 'YLim');
for i=1:length(bins)-1
    angle = mean(bins(i:i+1)); 
    text(angle, -0.15*yl(2), '\rightarrow', 'FontSize', 11, 'FontWeight', 'bold', 'Rotation', angle, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle') 
end
text(-200, 10, sprintf('# of\nsession'), 'FontSize', 9, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle')


ah = axes;
set(ah, 'Units', 'Normalized', 'Position', [0.18+asize+0.05+0.01 0.4 asize asize*pw/ph]);
plotML_RFdist(a_cy, th_cy, d_cy, 0, ah);
set(ah, 'FontSize', afs)
set(ah, 'YTickLabel', [])
% plot direction distribution
apos  = get(gca, 'Position');
apos2 = [apos(1)+12*apos(3)/20 apos(2)+12.5*apos(4)/20 7*apos(3)/20 6*apos(4)/20];
set(axes, 'Units', 'Normalized', 'Position', apos2);
dd       = dd_cy;
bins     = [-90:180/5:90];   % bin edges
n        = histc(dd,bins);
bh       = bar(bins, n, 0.99, 'histc');
set(bh,'LineWidth', 1, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', [1 1 1])
set(gca, 'XLim', [min(bins) max(bins)], 'XTick', [], 'YLim', [0 40], 'YTick', [0 20 40], ...
            'Layer', 'top', 'box', 'off')
yl = get(gca, 'YLim');
for i=1:length(bins)-1
    angle = mean(bins(i:i+1)); 
    text(angle, -0.15*yl(2), '\rightarrow', 'FontSize', 11, 'FontWeight', 'bold', 'Rotation', angle, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle') 
end


ah = axes;
set(ah, 'Units', 'Normalized', 'Position', [0.18+2*asize+0.08+0.01 0.4 asize asize*pw/ph]);
plotML_RFdist(trga_cy, trgth_cy, 20*ones(size(trga_cy)), 1, ah);
set(ah, 'FontSize', afs)
set(ah, 'YTickLabel', [])




% zsa zsa
ah = axes;
set(ah, 'Units', 'Normalized', 'Position', [0.18 0.38+asize asize asize*pw/ph]);
plotML_RFdist(a_zzp, th_zzp, d_zzp, 0, ah);
set(ah, 'FontSize', afs)
set(ah, 'XTickLabel', [])
% plot direction distribution
apos  = get(gca, 'Position');
apos2 = [apos(1)+12*apos(3)/20 apos(2)+12.5*apos(4)/20 7*apos(3)/20 6*apos(4)/20];
set(axes, 'Units', 'Normalized', 'Position', apos2);
dd       = dd_zzp;
bins     = [-90:180/5:90];   % bin edges
n        = histc(dd,bins);
bh       = bar(bins, n, 0.99, 'histc');
set(bh,'LineWidth', 1, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', [1 1 1])
set(gca, 'XLim', [min(bins) max(bins)], 'XTick', [], 'YLim', [0 20], 'YTick', [0 10 20], ...
            'Layer', 'top', 'box', 'off')
yl = get(gca, 'YLim');
for i=1:length(bins)-1
    angle = mean(bins(i:i+1)); 
    text(angle, -0.15*yl(2), '\rightarrow', 'FontSize', 11, 'FontWeight', 'bold', 'Rotation', angle, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle') 
end


ah = axes;
set(ah, 'Units', 'Normalized', 'Position', [0.18+asize+0.05+0.01 0.38+asize asize asize*pw/ph]);
plotML_RFdist(a_zz, th_zz, d_zz, 0, ah);
set(ah, 'FontSize', afs)
set(ah, 'XTickLabel', [], 'YTickLabel', [])
% plot direction distribution
apos  = get(gca, 'Position');
apos2 = [apos(1)+12*apos(3)/20 apos(2)+12.5*apos(4)/20 7*apos(3)/20 6*apos(4)/20];
set(axes, 'Units', 'Normalized', 'Position', apos2);
dd       = dd_zz;
bins     = [-90:180/5:90];   % bin edges
n        = histc(dd,bins);
bh       = bar(bins, n, 0.99, 'histc');
set(bh,'LineWidth', 1, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', [1 1 1])
set(gca, 'XLim', [min(bins) max(bins)], 'XTick', [], 'YLim', [0 30], 'YTick', [0 15 30], ...
            'Layer', 'top', 'box', 'off')
yl = get(gca, 'YLim');
for i=1:length(bins)-1
    angle = mean(bins(i:i+1)); 
    text(angle, -0.15*yl(2), '\rightarrow', 'FontSize', 11, 'FontWeight', 'bold', 'Rotation', angle, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle') 
end


ah = axes;
set(ah, 'Units', 'Normalized', 'Position', [0.18+2*asize+0.08+0.01 0.38+asize asize asize*pw/ph]);
plotML_RFdist(trga_zz, trgth_zz, 20*ones(size(trga_zz)), 1, ah);
set(ah, 'FontSize', afs)
set(ah, 'XTickLabel', [], 'YTickLabel', [])



%%% report statistics %%%
% cyrus
cybe4 = [nanmean(a_cyp) meanangle(th_cyp) nanmean(d_cyp) meanangle(dd_cyp)  nanmean(s_cyp); ...
         nanstd(a_cyp)  stdangle(th_cyp)  nanstd(d_cyp)  stdangle(dd_cyp)   nanstd(s_cyp)];

cydur = [nanmean(a_cy) meanangle(th_cy) nanmean(d_cy) meanangle(dd_cy) nanmean(s_cy) nanmean(trga_cy) meanangle(trgth_cy); ...
         nanstd(a_cy)  stdangle(th_cy)  nanstd(d_cy)  stdangle(dd_cy)  nanstd(s_cy)  nanstd(trga_cy)  stdangle(trgth_cy)];

     
zzbe4 = [nanmean(a_zzp) meanangle(th_zzp) nanmean(d_zzp) meanangle(dd_zzp)  nanmean(s_zzp); ...
         nanstd(a_zzp)  stdangle(th_zzp)  nanstd(d_zzp)  stdangle(dd_zzp)   nanstd(s_zzp)];

zzdur = [nanmean(a_zz) meanangle(th_zz) nanmean(d_zz) meanangle(dd_zz) nanmean(s_zz) nanmean(trga_zz) meanangle(trgth_zz); ...
         nanstd(a_zz)  stdangle(th_zz)  nanstd(d_zz)  stdangle(dd_zz)  nanstd(s_zz)  nanstd(trga_zz)  stdangle(trgth_zz)];
     
     
% subplot(2,2,1)
% %% plot direction distribution
% % figure propreties
% dd    = dd_cyp;
% lw    = 3;                      % line width
% binc  = pi/180*[0:30:360-30];   % bin center
% 
% rose(dd*pi/180, binc);
% set(findall(gcf,'Type', 'text'),'FontSize', 18)             % set tick labels font size
% set(findall(gcf,'Type', 'line'),'LineWidth', 1)             % set grid line line widths
% set(findobj('Type','Line'),'LineWidth', lw, 'Color', 'k')   % set rose plot line width
% 
% 
% subplot(2,2,2)
% %% plot direction distribution
% % figure propreties
% dd    = dd_cy;
% lw    = 3;                      % line width
% binc  = pi/180*[0:30:360-30];   % bin center
% 
% rose(dd*pi/180, binc);
% set(findall(gcf,'Type', 'text'),'FontSize', 18)             % set tick labels font size
% set(findall(gcf,'Type', 'line'),'LineWidth', 1)             % set grid line line widths
% set(findobj('Type','Line'),'LineWidth', lw, 'Color', 'k')   % set rose plot line width
% 
% subplot(2,2,3)
% %% plot direction distribution
% % figure propreties
% dd    = dd_zzp;
% lw    = 3;                      % line width
% binc  = pi/180*[0:30:360-30];   % bin center
% 
% rose(dd*pi/180, binc);
% set(findall(gcf,'Type', 'text'),'FontSize', 18)             % set tick labels font size
% set(findall(gcf,'Type', 'line'),'LineWidth', 1)             % set grid line line widths
% set(findobj('Type','Line'),'LineWidth', lw, 'Color', 'k')   % set rose plot line width
% 
% 
% subplot(2,2,4)
% %% plot direction distribution
% % figure propreties
% dd    = dd_zz;
% lw    = 3;                      % line width
% binc  = pi/180*[0:30:360-30];   % bin center
% 
% rose(dd*pi/180, binc);
% set(findall(gcf,'Type', 'text'),'FontSize', 18)             % set tick labels font size
% set(findall(gcf,'Type', 'line'),'LineWidth', 1)             % set grid line line widths
% set(findobj('Type','Line'),'LineWidth', lw, 'Color', 'k')   % set rose plot line width
% 
