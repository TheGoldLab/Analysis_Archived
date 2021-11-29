%% plotML_fig3.m
% Matlab script for plotting Fig 3  (MT data)
% Final version for publication to Nature Neuroscience
%
%   a) average responses at difference training epochs
%   b) fitted parameters to simple linear model
%   c) choice probability vs thresholds sensitivity of MT neuron
%
%
% created by jcl on 1/16/07

%%
% load data
[m_cy, ms_cy, p_cy, mp_cy, msp_cy, pp_cy] = getML_MTCTdependence('Cy', 0);
[m_zz, ms_zz, p_zz, mp_zz, msp_zz, pp_zz] = getML_MTCTdependence('ZZ', 0);
a             = getML_txt('CyTRain_MT.txt');
mses_cy       = a.data{strcmp(a.name,'session')};
fn_cy         = a.data{strcmp(a.name,'dat_fn')};
a             = getML_txt('ZZTRain_MT.txt');
mses_zz       = a.data{strcmp(a.name,'session')};
bses_cy       = [1:160]';
bses_zz       = [1:130]';



fh = figure; 
mm2inch = 25.4;
pw = 100/mm2inch;   % paper width, make it nature size of one column, 85mm x 85mm
ph = 120/mm2inch;   % paper height
wysifig(fh, pw, ph) % set it to US letter size

% figure parameters
alw = 0.5;  % axis line width
afs = 5;    % axis font size
lfs = 5;    % label font size
ffs = 12;   % figure font size
tfs = 5;    % text font size
lw  = 1;    % line width
elw = 0.5;
ms  = 2;    % marker size
lcc = [0 0 0];
lce = [0.75 0.75 0.75];

asize = 0.1558;
%% 1a: plot example fits
th = annotation('textbox', [0.1-0.1 0.96 0.05 0.05]);
set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'String', 'a', 'FontSize', ffs, 'FontWeight', 'bold', 'FontName', 'Helvetica') 

EgInd = [43 6 62 105]; % example fits (pre, post-early, post-late)

           
% get data
% pre-training
Monk = 'Cy';
a      = getML_txt([Monk 'PRe_MT.txt']);
fn     = a.data{strcmp(a.name,'dat_fn')};
usable = a.data{strcmp(a.name,'usable')};
    
crtf = [1];
nvf  = [0 1];
dirf = 0;

bb     = -200;
be     = 1500;
bs     = 25;
bw     = 100;
bins   = [[bb:bs:be-bw]' [bb+bw:bs:be]'];
mbins  = mean(bins,2);

[rMp_cy, rSDp_cy, rNp_cy, rMnp_cy] = getML_mRate([Monk 'PRe_MT.txt'], bins, crtf, nvf, dirf, 0);


% post-training
Monk = 'Cy';
a      = getML_txt([Monk 'TRain_MT.txt']);
fn     = a.data{strcmp(a.name,'dat_fn')};
ses    = a.data{strcmp(a.name,'session')};
usable = a.data{strcmp(a.name,'usable')};
    
crtf = [1];
nvf  = [0 1];
dirf = 0;

bb     = -200;
be     = 1500;
bs     = 25;
bw     = 100;
bins   = [[bb:bs:be-bw]' [bb+bw:bs:be]'];
mbins  = mean(bins,2);

[rM_cy, rSD_cy, rN_cy, rMn_cy] = getML_mRate([Monk 'TRain_MT.txt'], bins, crtf, nvf, dirf, 0);


% get average
rMp_avg = nanmean(rMp_cy,3);
a       = getML_txt([Monk 'PRe_MT.txt']);
rNp     = sum(a.data{strcmp(a.name,'usable')}==1);


N    = 3;
Ms   = max(ses);
Bses = [1:floor(Ms/N):Ms]';
Bses = [Bses(1:N) [Bses(2:N)-1; Ms]];
Bses = [1 50; 51 100; 101 160];


rM_avg     = nans(14,length(mbins),N);
rMn_avg    = nans(14,length(mbins),N);

a          = getML_txt([Monk 'TRain_MT.txt']);
use        = a.data{strcmp(a.name,'usable')};

for i = 1:N
    Lses = ses>=Bses(i,1) & ses<Bses(i,2);
    rN(i)          = sum(usable==1&Lses);
    rM_avg(:,:,i)  = nanmean(rM_cy(:,:,Lses),3);
    rMn_avg(:,:,i) = nanmean(rMn_cy(:,:,Lses),3);
end



% plot
ah(1) = axes;
set(ah(1), 'Units', 'Normalized', 'Position', [0.1 0.82 asize*ph/pw asize], 'FontName', 'Helvetica')
hold on 

lc  = [];
for j = 1:7
    lc = [lc; [0.9 0.9 0.9]-(j-1)/7];
end

bb     = -200;
be     = 1500;
bs     = 25;
bw     = 100;
bins   = [[bb:bs:be-bw]' [bb+bw:bs:be]'];
mbins  = mean(bins,2);
             
hold on
plotML_rateByDirCoh(mbins/1000, rMp_avg, gca, lw)
set(ah(1), 'LineWidth', alw, 'FontSize', afs, 'Box', 'off', ...
           'XLim', [0 1], 'XTick', [0 0.5 1], 'XTickLabel', {'0', '0.5', '1.0'} )
ylabel('Response (spikes per s)', 'FontSize', lfs)
xlabel('Time (s)', 'FontSize', lfs)
text(0.99, 60, [sprintf('Pre-training\n') '\itn\rm' sprintf(' = %.0f',rNp)], ...
        'FontSize', tfs-1, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top')
set(gca, 'YTick', [0 30 60])


% post-training
for i = 1:3
    ah(i+1) = axes;
    set(ah(i+1), 'Units', 'Normalized', 'Position', [0.1+i*(asize+0.02)*ph/pw 0.82 asize*ph/pw asize], 'FontName', 'Helvetica')
    plotML_rateByDirCoh(mbins/1000, rM_avg(:,:,i), gca, lw)
    set(ah(i+1), 'LineWidth', alw, 'FontSize', afs, 'Box', 'off', ...
        'XLim', [0 1], 'XTick', [0 0.5 1], 'XTickLabel', {'0', '0.5', '1.0'})
    xlabel([])
    ylabel([])
  
        text(0.99, 60, [sprintf('Sessions %.0f-%.0f\n',Bses(i,1), Bses(i,2)) '\itn\rm' sprintf(' = %.0f', rN(i))], ...
            'FontSize', tfs-1, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top')
        set(gca, 'YTick', [0 30 60], 'YTickLabel', [])
        ylim([0 60])
end
   



%% MT
th = annotation('textbox', [0 0.72 0.05 0.05]);
set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'String', 'b', 'FontSize', ffs, 'FontWeight', 'bold', 'FontName', 'Helvetica') 

% cyrus
% coh
set(axes, 'Units', 'Normalized', 'Position', [0.15 0.62 0.5*asize*ph/pw 0.75*asize], 'FontName', 'Helvetica')
hold on
% cy
L     = ~isnan(mp_cy(1,:));
N     = sum(L);
mpses = 1:N;
plot([1:N],repmat(0,1,N),':k')
plot([mpses; mpses],[mp_cy(2,L)-msp_cy(2,L); mp_cy(2,L)+msp_cy(2,L)],'-', 'LineWidth', elw, 'Color', lce)
plot(mpses,mp_cy(2,L),'ok', 'MarkerSize', ms, 'MarkerFaceColor', lcc, 'MarkerEdgeColor', 'none')
[b bi h p] = nestedFW(mp_cy(2,L)', ones(size(msp_cy(2,L)')), [ones(N,1) mpses']);
if p<0.05
    plot(1:N, b(1)+b(2)*[1:N], '-k', 'LineWidth', 1)
else
%    plot(1:N, b(1)+b(2)*[1:N], '--b', 'LineWidth', 1)
end
hold off
set(gca, 'xlim', [1 50], 'XTick', [1 25 50], 'XTickLabel', [], ...
         'ylim', [-3 3], 'YTick', [-3 0 3], 'YTickLabel', [-3 0 3], ...
         'FontSize', afs, 'LineWidth', alw)
text(5, 2.9, 'Pre-training', 'FontSize', tfs-1, ...
     'HorizontalAlignment', 'left', 'VerticalAlignment', 'Top')
text(5, -2.9,  'Monkey C', 'FontSize', tfs-1, ...
     'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom')
ylabel('% coh^{-1}', 'FontSize', lfs)
text(-45, 0, 'Coh (k_1)', 'FontSize', tfs+2, ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
     'Rotation', 90)

 
 
set(axes, 'Units', 'Normalized', 'Position', [0.15+1*(0.5*asize+0.025)*ph/pw 0.62 asize*ph/pw 0.75*asize], 'FontName', 'Helvetica')
hold on
% cy
plot([1:max(bses_cy)],repmat(0,1,max(bses_cy)),':k')
plot([mses_cy'; mses_cy'],[m_cy(2,:)-ms_cy(2,:); m_cy(2,:)+ms_cy(2,:)],'-', 'LineWidth', elw, 'Color', lce)
plot(mses_cy,m_cy(2,:),'ok', 'MarkerSize', ms, 'MarkerFaceColor', lcc, 'MarkerEdgeColor', 'none')
Lgd = ~isnan(m_cy(2,:)');
[b bi h p] = nestedFW(m_cy(2,Lgd)', ones(size(ms_cy(2,Lgd)')), [ones(sum(Lgd),1) mses_cy(Lgd)]);
if p<0.05
    plot(1:max(bses_cy), b(1)+b(2)*[1:max(bses_cy)], '-k', 'LineWidth', 1)
else
%    plot(1:max(bses_cy), b(1)+b(2)*[1:max(bses_cy)], '--b', 'LineWidth', 1)
end
hold off
text(5, 2.9,  'Training', 'FontSize', tfs-1, ...
     'HorizontalAlignment', 'left', 'VerticalAlignment', 'Top')
set(gca, 'xlim', [1 max(bses_cy)], 'XTick', [1 50 100 150], 'XTickLabel', [], ...
         'ylim', [-3 3], 'YTick', [-3 0 3], 'YTickLabel', [], 'YColor', [1 1 1], ...
         'FontSize', afs, 'LineWidth', alw)
 

% time 
set(axes, 'Units', 'Normalized', 'Position', [0.15 0.62-1*(0.75*asize+0.025) 0.5*asize*ph/pw 0.75*asize])
hold on
% cy
L     = ~isnan(mp_cy(1,:));
N     = sum(L);
mpses = 1:N;
plot([1:N],repmat(0,1,N),':k')
plot([mpses; mpses],[mp_cy(3,L)-msp_cy(3,L); mp_cy(3,L)+msp_cy(3,L)],'-', 'LineWidth', elw, 'Color', lce)
plot(mpses,mp_cy(3,L),'ok', 'MarkerSize', ms, 'MarkerFaceColor', lcc, 'MarkerEdgeColor', 'none')
[b bi h p] = nestedFW(mp_cy(3,L)', ones(size(msp_cy(3,L)')), [ones(N,1) mpses']);
if p<0.05
    plot(1:N, b(1)+b(2)*[1:N], '-k', 'LineWidth', 1)
else
%    plot(1:N, b(1)+b(2)*[1:N], '--b', 'LineWidth', 1)
end
hold off
set(gca, 'xlim', [1 50], 'XTick', [1 25 50], 'XTickLabel', [], ...
         'ylim', [-3 3], 'YTick', [-3 0 3], 'YTickLabel', [-3 0 3], ...
         'FontSize', afs, 'LineWidth', alw)
ylabel('s^{-1}', 'FontSize', lfs)
text(-45, 0, 'Time (k_2)', 'FontSize', tfs+2, ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
     'Rotation', 90)

set(axes, 'Units', 'Normalized', 'Position', [0.15+1*(0.5*asize+0.025)*ph/pw 0.62-1*(0.75*asize+0.025) asize*ph/pw 0.75*asize], 'FontName', 'Helvetica')
hold on
% cy
plot([1:max(bses_cy)],repmat(0,1,max(bses_cy)),':k')
plot([mses_cy'; mses_cy'],[m_cy(3,:)-ms_cy(3,:); m_cy(3,:)+ms_cy(3,:)],'-', 'LineWidth', elw, 'Color', lce)
plot(mses_cy,m_cy(3,:),'ok', 'MarkerSize', ms, 'MarkerFaceColor', lcc, 'MarkerEdgeColor', 'none')
Lgd = ~isnan(m_cy(3,:)');
[b bi h p] = nestedFW(m_cy(3,Lgd)', ones(size(ms_cy(3,Lgd)')), [ones(sum(Lgd),1) mses_cy(Lgd)]);
if p<0.05
    plot(1:max(bses_cy), b(1)+b(2)*[1:max(bses_cy)], '-k', 'LineWidth', 1)
else
%    plot(1:max(bses_cy), b(1)+b(2)*[1:max(bses_cy)], '--b', 'LineWidth',
%    1)
end
hold off
set(gca, 'xlim', [1 max(bses_cy)], 'XTick', [1 50 100 150], 'XTickLabel', [], ...
         'ylim', [-3 3], 'YTick', [-3 0 3], 'YTickLabel', [], 'YColor', [1 1 1], ...
         'FontSize', afs, 'LineWidth', alw)

 

% coh*time
set(axes, 'Units', 'Normalized', 'Position', [0.15 0.62-2*(0.75*asize+0.025) 0.5*asize*ph/pw 0.75*asize])
hold on
% cy
L     = ~isnan(mp_cy(1,:));
N     = sum(L);
mpses = 1:N;
plot([1:N],repmat(0,1,N),':k')
plot([mpses; mpses],[mp_cy(4,L)-msp_cy(4,L); mp_cy(4,L)+msp_cy(4,L)],'-', 'LineWidth', elw, 'Color', lce)
plot(mpses,mp_cy(4,L),'ok', 'MarkerSize', ms, 'MarkerFaceColor', lcc, 'MarkerEdgeColor', 'none')
[b bi h p] = nestedFW(mp_cy(4,L)', ones(size(msp_cy(4,L)')), [ones(N,1) mpses']);
if p<0.05
    plot(1:N, b(1)+b(2)*[1:N], '-k', 'LineWidth', 1)
else
%    plot(1:N, b(1)+b(2)*[1:N], '--b', 'LineWidth', 1)
end
hold off
set(gca, 'xlim', [1 50], 'XTick', [1 25 50], 'XTickLabel', [1 25 50], ...
         'ylim', [-3 3], 'YTick', [-3 0 3], 'YTickLabel', [-3 0 3], ...
         'FontSize', afs, 'LineWidth', alw)
ylabel('(coh \times s)^{-1}', 'FontSize', lfs)
text(-45, 0, 'Coh \times time (k_3)', 'FontSize', tfs+2, ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
     'Rotation', 90)
xlabel('Session (d)', 'FontSize', lfs)

 
set(axes, 'Units', 'Normalized', 'Position', [0.15+1*(0.5*asize+0.025)*ph/pw 0.62-2*(0.75*asize+0.025) asize*ph/pw 0.75*asize], 'FontName', 'Helvetica')
hold on
% cy
plot([1:max(bses_cy)],repmat(0,1,max(bses_cy)),':k')
plot([mses_cy'; mses_cy'],[m_cy(4,:)-ms_cy(4,:); m_cy(4,:)+ms_cy(4,:)],'-', 'LineWidth', elw, 'Color', lce)
plot(mses_cy,m_cy(4,:),'ok', 'MarkerSize', ms, 'MarkerFaceColor', lcc, 'MarkerEdgeColor', 'none')
Lgd = ~isnan(m_cy(4,:)');
[b bi h p] = nestedFW(m_cy(4,Lgd)', ones(size(ms_cy(4,Lgd)')), [ones(sum(Lgd),1) mses_cy(Lgd)]);
if p<0.05
    plot(1:max(bses_cy), b(1)+b(2)*[1:max(bses_cy)], '-k', 'LineWidth', 1)
else
%    plot(1:max(bses_cy), b(1)+b(2)*[1:max(bses_cy)], '--b', 'LineWidth', 1)
end
hold off
set(gca, 'xlim', [1 max(bses_cy)], 'XTick', [1 50 100 150], 'XTickLabel', [1 50 100 150], ...
         'ylim', [-3 3], 'YTick', [-3 0 3], 'YTickLabel', [], 'YColor', [1 1 1], ...
         'FontSize', afs, 'LineWidth', alw)
xlabel('Session (d)', 'FontSize', lfs)


 
 
%% ZsaZsa
% coh
set(axes, 'Units', 'Normalized', 'Position', [0.58 0.62 0.5*asize*ph/pw 0.75*asize], 'FontName', 'Helvetica')
hold on
% zz
plot([1:max(bses_cy)],repmat(0,1,max(bses_cy)),':k')
L     = ~isnan(mp_zz(1,:));
N     = sum(L);
mpses = 1:N;
plot([mpses; mpses],[mp_zz(2,L)-msp_zz(2,L); mp_zz(2,L)+msp_zz(2,L)],'-', 'LineWidth', elw, 'Color', lce)
plot(mpses,mp_zz(2,L),'ok', 'MarkerSize', ms, 'MarkerFaceColor', lcc, 'MarkerEdgeColor', 'none')
[b bi h p] = nestedFW(mp_zz(2,L)', ones(size(msp_zz(2,L)')), [ones(N,1) mpses']);
if p<0.05
    plot(1:N, b(1)+b(2)*[1:N], '-k', 'LineWidth', 1)
else 
%    plot(1:N, b(1)+b(2)*[1:N], '--r', 'LineWidth', 1)
end
hold off
set(gca, 'xlim', [1 50], 'XTick', [1 25 50], 'XTickLabel', [], ...
         'ylim', [-3 3], 'YTick', [-3 0 3], 'YTickLabel', [-3 0 3], ...
         'FontSize', afs, 'LineWidth', alw)
text(5, -2.9, 'Monkey Z', 'FontSize', tfs-1, ...
     'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom')
text(5, 2.9, 'Pre-training', 'FontSize', tfs-1, ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'Top')



set(axes, 'Units', 'Normalized', 'Position', [0.58+1*(0.5*asize+0.025)*ph/pw 0.62 asize*ph/pw 0.75*asize], 'FontName', 'Helvetica')
hold on
% zz
plot([1:max(bses_cy)],repmat(0,1,max(bses_cy)),':k')
plot([mses_zz'; mses_zz'],[m_zz(2,:)-ms_zz(2,:); m_zz(2,:)+ms_zz(2,:)],'-', 'LineWidth', elw, 'Color', lce)
plot(mses_zz,m_zz(2,:),'ok', 'MarkerSize', ms, 'MarkerFaceColor', lcc, 'MarkerEdgeColor', 'none')
Lgd = ~isnan(m_zz(2,:)');
[b bi h p] = nestedFW(m_zz(2,Lgd)', ones(size(ms_zz(2,Lgd)')), [ones(sum(Lgd),1) mses_zz(Lgd)]);
if p<0.05
    plot(1:max(bses_zz), b(1)+b(2)*[1:max(bses_zz)], '-k', 'LineWidth', 1)
else
%    plot(1:max(bses_zz), b(1)+b(2)*[1:max(bses_zz)], '--r', 'LineWidth', 1)
end
hold off
set(gca, 'xlim', [1 max(bses_zz)], 'XTick', [1 40 80 120], 'XTickLabel', [], ...
         'ylim', [-3 3], 'YTick', [-3 0 3], 'YTickLabel', [], 'YColor', [1 1 1], ...
         'FontSize', afs, 'LineWidth', alw)
text(5, 2.9,  'Training', 'FontSize', tfs-1, ...
     'HorizontalAlignment', 'left', 'VerticalAlignment', 'Top')
 

% time 
set(axes, 'Units', 'Normalized', 'Position', [0.58 0.62-1*(0.75*asize+0.025) 0.5*asize*ph/pw 0.75*asize], 'FontName', 'Helvetica')
hold on
plot([1:max(bses_cy)],repmat(0,1,max(bses_cy)),':k')
% zz
L     = ~isnan(mp_zz(1,:));
N     = sum(L);
mpses = 1:N;
plot([mpses; mpses],[mp_zz(3,L)-msp_zz(3,L); mp_zz(3,L)+msp_zz(3,L)],'-', 'LineWidth', elw, 'Color', lce)
plot(mpses,mp_zz(3,L),'ok', 'MarkerSize', ms, 'MarkerFaceColor', lcc, 'MarkerEdgeColor', 'none')
[b bi h p] = nestedFW(mp_zz(3,L)', ones(size(msp_zz(3,L)')), [ones(N,1) mpses']);
if p<0.05
    plot(1:N, b(1)+b(2)*[1:N], '-.k', 'LineWidth', 1)
else
%    plot(1:N, b(1)+b(2)*[1:N], '--r', 'LineWidth', 1)
end
hold off
set(gca, 'xlim', [1 50], 'XTick', [1 25 50], 'XTickLabel', [], ...
         'ylim', [-3 3], 'YTick', [-3 0 3], 'YTickLabel', [-3 0 3], ...
         'FontSize', afs, 'LineWidth', alw)

 
set(axes, 'Units', 'Normalized', 'Position', [0.58+1*(0.5*asize+0.025)*ph/pw 0.62-1*(0.75*asize+0.025) asize*ph/pw 0.75*asize], 'FontName', 'Helvetica')
hold on
% zz
plot([1:max(bses_cy)],repmat(0,1,max(bses_cy)),':k')
plot([mses_zz'; mses_zz'],[m_zz(3,:)-ms_zz(3,:); m_zz(3,:)+ms_zz(3,:)],'-', 'LineWidth', elw, 'Color', lce)
plot(mses_zz,m_zz(3,:),'ok', 'MarkerSize', ms, 'MarkerFaceColor', lcc, 'MarkerEdgeColor', 'none')
Lgd = ~isnan(m_zz(3,:)');
[b bi h p] = nestedFW(m_zz(3,Lgd)', ones(size(ms_zz(3,Lgd)')), [ones(sum(Lgd),1) mses_zz(Lgd)]);
if p<0.05
    plot(1:max(bses_zz), b(1)+b(2)*[1:max(bses_zz)], '-k', 'LineWidth', 1)
else
%    plot(1:max(bses_zz), b(1)+b(2)*[1:max(bses_zz)], '--r', 'LineWidth', 1)
end
hold off
set(gca, 'xlim', [1 max(bses_zz)], 'XTick', [1 40 80 120], 'XTickLabel', [], ...
         'ylim', [-3 3], 'YTick', [-3 0 3], 'YTickLabel', [], 'YColor', [1 1 1], ...
         'FontSize', afs, 'LineWidth', alw)

 

% coh*time
set(axes, 'Units', 'Normalized', 'Position', [0.58 0.62-2*(0.75*asize+0.025) 0.5*asize*ph/pw 0.75*asize], 'FontName', 'Helvetica')
hold on
% zz
plot([1:max(bses_cy)],repmat(0,1,max(bses_cy)),':k')
L     = ~isnan(mp_zz(1,:));
N     = sum(L);
mpses = 1:N;
plot([mpses; mpses],[mp_zz(4,L)-msp_zz(4,L); mp_zz(4,L)+msp_zz(4,L)],'-', 'LineWidth', elw, 'Color', lce)
plot(mpses,mp_zz(4,L),'ok', 'MarkerSize', ms, 'MarkerFaceColor', lcc, 'MarkerEdgeColor', 'none')
[b bi h p] = nestedFW(mp_zz(4,L)', ones(size(msp_zz(4,L)')), [ones(N,1) mpses']);
if p<0.05
    plot(1:N, b(1)+b(2)*[1:N], '-k', 'LineWidth', 1)
else
%    plot(1:N, b(1)+b(2)*[1:N], '--r', 'LineWidth', 1)
end
hold off
set(gca, 'xlim', [1 50], 'XTick', [1 25 50], 'XTickLabel', [1 25 50], ...
         'ylim', [-3 3], 'YTick', [-3 0 3], 'YTickLabel', [-3 0 3], ...
         'FontSize', afs, 'LineWidth', alw)
xlabel('Session (d)', 'FontSize', lfs)

 
set(axes, 'Units', 'Normalized', 'Position', [0.58+1*(0.5*asize+0.025)*ph/pw 0.62-2*(0.75*asize+0.025) asize*ph/pw 0.75*asize], 'FontName', 'Helvetica')
hold on
plot([1:max(bses_cy)],repmat(0,1,max(bses_cy)),':k')
% zz
plot([mses_zz'; mses_zz'],[m_zz(4,:)-ms_zz(4,:); m_zz(4,:)+ms_zz(4,:)],'-', 'LineWidth', elw, 'Color', lce)
plot(mses_zz,m_zz(4,:),'ok', 'MarkerSize', ms, 'MarkerFaceColor', lcc, 'MarkerEdgeColor', 'none')
Lgd = ~isnan(m_zz(4,:)');
[b bi h p] = nestedFW(m_zz(4,Lgd)', ones(size(ms_zz(4,Lgd)')), [ones(sum(Lgd),1) mses_zz(Lgd)]);
if p<0.05
    plot(1:max(bses_zz), b(1)+b(2)*[1:max(bses_zz)], '-k', 'LineWidth', 1)
else
%    plot(1:max(bses_zz), b(1)+b(2)*[1:max(bses_zz)], '--r', 'LineWidth', 1)
end
hold off
set(gca, 'xlim', [1 max(bses_zz)], 'XTick', [1 40 80 120], 'XTickLabel', [1 40 80 120], ...
         'ylim', [-3 3], 'YTick', [-3 0 3], 'YTickLabel', [], 'YColor', [1 1 1], ...
         'FontSize', afs, 'LineWidth', alw)
xlabel('Session (d)', 'FontSize', lfs)



% %% plot neurometric thresholds
% th = annotation('textbox', [0.1-0.1 0.35 0.05 0.05]);
% set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
%     'String', 'c', 'FontSize', ffs, 'FontWeight', 'bold') 
% 
% % load data
% Monk = 'Cy';
% [fitsp_cy, semsp_cy] = getML_neurometricROCT([Monk 'PRe_MT.txt'], bins, 1, [0 1], 0, 0);
% [fits_cy, sems_cy]   = getML_neurometricROCT([Monk 'TRain_MT.txt'], bins, 1, [0 1], 0, 0);
% 
% Monk = 'ZZ';
% [fitsp_zz, semsp_zz] = getML_neurometricROCT([Monk 'PRe_MT.txt'], bins, 1, [0 1], 0, 0);
% [fits_zz, sems_zz]   = getML_neurometricROCT([Monk 'TRain_MT.txt'], bins, 1, [0 1], 0, 0);
% 
% fitsp_cy(1,:)   = 100*fitsp_cy(1,:);
% fitsp_zz(1,:)   = 100*fitsp_zz(1,:);
% fits_cy(1,:)    = 100*fits_cy(1,:);
% fits_zz(1,:)    = 100*fits_zz(1,:);
% 
% semsp_cy(1,:,:) = 100*semsp_cy(1,:,:);
% semsp_zz(1,:,:) = 100*semsp_zz(1,:,:);
% sems_cy(1,:,:)  = 100*sems_cy(1,:,:);
% sems_zz(1,:,:)  = 100*sems_zz(1,:,:);
% 
% 
% 
% % cy
% set(axes, 'Units', 'Normalized', 'Position', [0.15 0.27 0.5*asize*ph/pw 0.75*asize])
% hold on
% L     = ~isnan(fitsp_cy(1,:));
% N     = sum(L);
% mpses = 1:N;
% plot([mpses; mpses],shiftdim(semsp_cy(1,:,L),1),'-', 'LineWidth', elw, 'Color', lce)
% plot(mpses,fitsp_cy(1,L),'ok', 'MarkerSize', ms, 'MarkerFaceColor', lcc, 'MarkerEdgeColor', 'none')
% plot([1:N],repmat(nangeomean(fitsp_cy(1,:)),1,N),'k', 'LineWidth', lw)
% hold off
% set(gca, 'LineWidth', alw, 'FontSize', afs, 'Box', 'off', ...
%          'XLim', [1 50], 'YLim', [3 180], ...
%          'XTick', [1 25 50], ...
%          'YScale', 'log', 'YTick', [5 20 80])
% ylabel('Threshold (% coh)', 'FontSize', lfs)
% % text(-51, nangeomean([3 180]), 'Threshold', 'FontSize', tfs+3, ...
% %      'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
% %      'Rotation', 90)
% xlabel('Session (d)', 'FontSize', lfs)
% text(5, 180, 'Pre-training', 'FontSize', tfs-1, ...
%      'HorizontalAlignment', 'left', 'VerticalAlignment', 'Top')
% text(5, 3.2,  'Monkey C', 'FontSize', tfs-1, ...
%      'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom')
% 
%  
% set(axes, 'Units', 'Normalized', 'Position', [0.15+1*(0.5*asize+0.025)*ph/pw 0.27 asize*ph/pw 0.75*asize])
% hold on
% plot([mses_cy'; mses_cy'],shiftdim(sems_cy(1,:,:),1),'-', 'LineWidth', elw, 'Color', lce)
% plot(mses_cy,fits_cy(1,:),'ok', 'MarkerSize', ms, 'MarkerFaceColor', lcc, 'MarkerEdgeColor', 'none')
% plot([1:160],repmat(nangeomean(fits_cy(1,:)),1,160),'k', 'LineWidth', lw)
% hold off
% set(gca, 'LineWidth', alw, 'FontSize', afs, 'Box', 'off', ...
%          'XLim', [1 160], 'YLim', [3 180], ...
%          'XTick', [1 50 100 150], ...
%          'YScale', 'log', 'YTick', [5 20 80], 'YTickLabel', [])
% xlabel('Session (d)', 'FontSize', lfs)
% text(5, 180, 'Pre-training', 'FontSize', tfs-1, ...
%      'HorizontalAlignment', 'left', 'VerticalAlignment', 'Top')
% 
% 
% % zz
% set(axes, 'Units', 'Normalized', 'Position', [0.58 0.27 0.5*asize*ph/pw 0.75*asize])
% hold on
% L     = ~isnan(fitsp_zz(1,:));
% N     = sum(L);
% mpses = 1:N;
% plot([mpses; mpses],shiftdim(semsp_zz(1,:,L),1),'-', 'LineWidth', elw, 'Color', lce)
% plot(mpses,fitsp_zz(1,L),'ok', 'MarkerSize', ms, 'MarkerFaceColor', lcc, 'MarkerEdgeColor', 'none')
% plot([1:N],repmat(nangeomean(fitsp_zz(1,:)),1,N),'k', 'LineWidth', lw)
% hold off
% set(gca, 'LineWidth', alw, 'FontSize', afs, 'Box', 'off', ...
%          'XLim', [1 50], 'YLim', [3 180], ...
%          'XTick', [1 25 50], ...
%          'YScale', 'log', 'YTick', [5 20 80])
% xlabel('Session (d)', 'FontSize', lfs)
% text(5, 180, 'Pre-training', 'FontSize', tfs-1, ...
%      'HorizontalAlignment', 'left', 'VerticalAlignment', 'Top')
% text(5, 3.2,  'Monkey Z', 'FontSize', tfs-1, ...
%      'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom')
% 
%  
% set(axes, 'Units', 'Normalized', 'Position', [0.58+1*(0.5*asize+0.025)*ph/pw 0.27 asize*ph/pw 0.75*asize])
% hold on
% plot([mses_zz'; mses_zz'],shiftdim(sems_zz(1,:,:),1),'-', 'LineWidth', elw, 'Color', lce)
% plot(mses_zz,fits_zz(1,:),'ok', 'MarkerSize', ms, 'MarkerFaceColor', lcc, 'MarkerEdgeColor', 'none')
% plot([1:130],repmat(nangeomean(fits_zz(1,:)),1,130),'k', 'LineWidth', lw)
% hold off
% set(gca, 'LineWidth', alw, 'FontSize', afs, 'Box', 'off', ...
%          'XLim', [1 130], 'YLim', [3 180], ...
%          'XTick', [1 40 80 120], ...
%          'YScale', 'log', 'YTick', [5 20 80], 'YTickLabel', [])
% xlabel('Session (d)', 'FontSize', lfs)
% text(5, 180, 'Pre-training', 'FontSize', tfs-1, ...
%      'HorizontalAlignment', 'left', 'VerticalAlignment', 'Top')
% 
% 
% 
% 
% % % test whether threshold change systematically with recording/training sessions
% % L = ~isnan(fitsp_zz(1,:));
% % [b bi h pvzp] = regressW(log(100*fitsp_zz(1,L)'), log(100*shiftdim(mean(semsp_zz(1,:,L)),2)), [ones(sum(L),1) [1:sum(L)]']);
% % % test whether lowest 50% threshold change systematically with recording/training sessions
% % L = ~isnan(fits_zz(1,:))' & [fits_zz(1,:)<myprctile(fits_zz(1,:),25)]';
% % [b bi h pvzL50] = regressW(log(100*fits_zz(1,L&Lp_zz)'), log(100*shiftdim(mean(sems_zz(1,:,L&Lp_zz)),2)), [ones(sum(L&Lp_zz),1) ses_zz(L&Lp_zz)]);
% % 
% % 
% % 
% % % t-test for pre-post thres
% % prz   = 100*fitsp_zz(1,:); prz(isnan(prz)) = [];
% % trz   = 100*fits_zz(1,Lp_zz);  trz(isnan(trz)) = [];
% % przm  = nangeomean(prz);
% % trzm  = nangeomean(trz);
% % pvalz = ranksum(prz,trz);
% 

%%     plot choice probability vs threshold     %%%
% Fig 1d, choice probability vs threshold
% define line color
lc1  = [0 0 0];
lc1e = [0.75 0.75 0.75];
lc2  = [1 0.3 0.3];
lc2e = [1 0.6 0.6];

th = annotation('textbox', [0 0.22 0.05 0.05]);
set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'String', 'c', 'FontSize', ffs, 'FontWeight', 'bold', 'FontName', 'Helvetica') 

% plot data from both monkeys together
Monk = 'Cy';
[fits_cy, sems_cy]   = getML_neurometricROCT([Monk 'TRain_MT.txt'], bins, 1, [0 1], 0, 0);
[roc_cy, ns_cy, rocd_cy, p_cy] = getML_MTChoiceProbability('Cy', 0);
Lp_cy = true(size(getML_LapseSelectionArray('Cy', 'MT')));

Monk = 'ZZ';
[fits_zz, sems_zz]   = getML_neurometricROCT([Monk 'TRain_MT.txt'], bins, 1, [0 1], 0, 0);
[roc_zz, ns_zz, rocd_zz, p_zz] = getML_MTChoiceProbability('ZZ', 0);
Lp_zz = true(size(getML_LapseSelectionArray('ZZ', 'MT')));



Monk    = 'Cy';
a       = getML_txt([Monk, 'TRain_MT.txt']);
ses_cy  = a.data{strcmp(a.name,'session')};
N       = 3;
Ms      = 160;
Bses_cy = [1:floor(Ms/N):Ms]';
Bses_cy = [Bses_cy(1:N) [Bses_cy(2:N)-1; Ms]];
Bses_cy = [1 50; 51 100; 101 160];

Monk    = 'ZZ';
a       = getML_txt([Monk, 'TRain_MT.txt']);
ses_zz  = a.data{strcmp(a.name,'session')};
N       = 3;
Ms      = 130;
Bses_zz = [1:floor(Ms/N):Ms]';
Bses_zz = [Bses_zz(1:N) [Bses_zz(2:N)-1; Ms]];
%Bses_zz = [1 49; 50 93; 94 130];
Bses_zz = [1 43; 44 87; 88 130];

for i = 1:N
    set(axes, 'Units', 'Normalized', 'Position', [0.13+(i-1)*0.55975/2 0.065 asize*ph/pw asize], 'FontName', 'Helvetica')
    Lses_cy = ses_cy>=Bses_cy(i,1) & ses_cy<Bses_cy(i,2);
    Lses_zz = ses_zz>=Bses_zz(i,1) & ses_zz<Bses_zz(i,2);   
    hold on
    plot([0:100], 0.5*ones(size([0:100])), 'k:')
    
    th      = [100*fits_cy(1,Lses_cy)'; 100*fits_zz(1,Lses_zz)'];
    thd     = [shiftdim(100*sems_cy(1,:,Lses_cy),1)'; shiftdim(100*sems_zz(1,:,Lses_zz),1)'];
    roc     = [roc_cy(Lses_cy,1); roc_zz(Lses_zz,1)];
    rocd    = [rocd_cy(Lses_cy,1); rocd_zz(Lses_zz,1)];
    
    % plot errorbars
    plot([th th]', [roc-rocd roc+rocd]', 'Color', lc1e, 'LineWidth', elw)
    plot(thd', [roc roc]', 'Color', lc1e, 'LineWidth', elw)
   
    plot(100*fits_zz(1,Lses_zz&Lp_zz), roc_zz(Lses_zz&Lp_zz,1), 'kv', 'MarkerFaceColor', 'k', 'MarkerSize', ms)
    plot(100*fits_cy(1,Lses_cy&Lp_cy), roc_cy(Lses_cy&Lp_cy,1), 'ks', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerSize', ms)
    
    plot(100*fits_zz(1,Lses_zz&~Lp_zz), roc_zz(Lses_zz&~Lp_zz,1), 'kv', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'k', 'MarkerSize', ms)
    plot(100*fits_cy(1,Lses_cy&~Lp_cy), roc_cy(Lses_cy&~Lp_cy,1), 'ks', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'k', 'MarkerSize', ms)
    
    % get regression
    L       = ~isnan(th) & ~isnan(roc);
    ths     = mean([th-thd(:,1) thd(:,2)-th],2);
    %[b, bi, h, p, F, df]  = regressXYW(log(th(L)), roc(L));
    [b, bi]  = regressXYW(log(th(L)), roc(L));
    %[b bi h p] = nestedFW(roc(L), 0.1*ones(size(rocd(L))), [ones(sum(L),1) log(th(L))]);
    plot(2:100, b(1)+b(2)*log(2:100), 'k', 'LineWidth', lw)
    hold off
    set(gca, 'LineWidth', alw, 'FontSize', afs, 'Box', 'off')
    if i==1
        set(gca, 'xlim', [3 100], 'XScale', 'log', 'XTick', [3 5 10 20 40 80], ...
                   'ylim', [0.35 0.75], 'YTick', [0.4 0.5 0.6 0.7])
        xlabel('Threshold (% coh)')
        ylabel('Choice probability')
%         text(100, 0.75, sprintf('First 1/3\nsessions'), 'FontSize', tfs, ...
%                 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top')
        title(sprintf('1^{st} 1/3 sessions'), 'FontSize', tfs+1)
    elseif i==2
        set(gca, 'xlim', [5 100], 'XScale', 'log', 'XTick', [5 10 20 40 80], ...
                   'ylim', [0.35 0.75], 'YTick', [0.4 0.5 0.6 0.7], 'YTickLabel', [])
%         text(100, 0.75, sprintf('Second 1/3\nsessions'), 'FontSize', tfs, ...
%                 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top')
        title(sprintf('2^{nd} 1/3 sessions'), 'FontSize', tfs+1)
    else
        set(gca, 'xlim', [5 100], 'XScale', 'log', 'XTick', [5 10 20 40 80], ...
                   'ylim', [0.35 0.75], 'YTick', [0.4 0.5 0.6 0.7], 'YTickLabel', [])
%         text(100, 0.75, sprintf('Third 1/3\nsessions'), 'FontSize', tfs, ...
%                 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top')
        title(sprintf('3^{rd} 1/3 sessions'), 'FontSize', tfs+1)
    end    
    [rcpth(i) pcpth(i)] = corr(roc(L), log(th(L)), 'tail', 'lt');
end
