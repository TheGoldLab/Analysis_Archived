%% plotML_fig4.m
% Matlab script for plotting Fig 4  (LIP data)
% Final version for publication to Nature Neuroscience
%
%   a) average responses at difference training epochs, align to motion and
%   sac onset
%   b) fitted parameters to simple linear model
%   c) changes of LIP activity at each coherence
%
%


%% plot LIP fig
% create figure
% load data
[l_cy, ls_cy] = getML_LIPCTdependence('Cy', 0);
[l_zz, ls_zz] = getML_LIPCTdependence('ZZ', 0);
a             = getML_txt('CyTRain_LIP.txt');
lses_cy       = a.data{strcmp(a.name,'session')};
a             = getML_txt('ZZTRain_LIP.txt');
lses_zz       = a.data{strcmp(a.name,'session')};
bses_cy       = [1:160]';
bses_zz       = [1:130]';



fh = figure; 
mm2inch = 25.4;
pw = 110/mm2inch;   % paper width, make it nature size of one column, 85mm x 85mm
ph = 88/mm2inch;   % paper height
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
lcz = [0 0 0];
lce = [0.75 0.75 0.75];
asize = 0.1558*120/88;
           
%% get data
% get data for LIP neurons
Monk = 'Cy';
a      = getML_txt([Monk 'TRain_LIP.txt']);
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

[rM_cy, rSD_cy, rN_cy, rMn_cy] = getML_mRate([Monk 'TRain_LIP.txt'], bins, crtf, nvf, dirf, 0);

N    = 3;
Ms   = max(ses);
Bses = [1:floor(Ms/N):Ms]';
Bses = [Bses(1:N) [Bses(2:N)-1; Ms]];
Bses = [1 50; 51 100; 101 160];

  

r_avg     = nans(14,length(mbins),N);
rn_avg    = nans(14,length(mbins),N);

a          = getML_txt([Monk 'TRain_LIP.txt']);
use        = a.data{strcmp(a.name,'usable')};

for i = 1:N
    Lses = ses>=Bses(i,1) & ses<Bses(i,2);
    rN(i)         = sum(usable==1&Lses);
    r_avg(:,:,i)  = nanmean(rM_cy(:,:,Lses),3);
    rn_avg(:,:,i) = nanmean(rMn_cy(:,:,Lses),3);
end


fname = 'CyTRain_LIP.txt';
crtf = [1];
nvf  = [0 1];
dirf = 0;

bb     = -500;
be     = 300;
bs     = 25;
bw     = 100;
bins   = [[bb:bs:be-bw]' [bb+bw:bs:be]'];
mbinss  = mean(bins,2);

[rM_cy, rSD_cy, rN_cy, rMn_cy] = getML_mRateSac(fname, bins, crtf, nvf, dirf, 0);


Monk = 'Cy';
a      = getML_txt([Monk 'TRain_LIP.txt']);
fn     = a.data{strcmp(a.name,'dat_fn')};
ses    = a.data{strcmp(a.name,'session')};
usable = a.data{strcmp(a.name,'usable')};


N    = 3;
Ms   = max(ses);
Bses = [1:floor(Ms/N):Ms]';
Bses = [Bses(1:N) [Bses(2:N)-1; Ms]];
Bses = [1 50; 51 100; 101 160];
  

rs_avg     = nans(14,length(mbinss),N);
rns_avg    = nans(14,length(mbinss),N);

a          = getML_txt([Monk 'TRain_LIP.txt']);
use        = a.data{strcmp(a.name,'usable')};

for i = 1:N
    Lses = ses>=Bses(i,1) & ses<Bses(i,2);
    rN(i)         = sum(usable==1&Lses);
   rs_avg(:,:,i)  = nanmean(rM_cy(:,:,Lses),3);
   rns_avg(:,:,i) = nanmean(rMn_cy(:,:,Lses),3); 
end


%%
% test whether activity during the delay period changed significantly as a
% function of session
[pt_cy nt_cy ptci_cy]   = getML_performanceAsso('CyTRain_psy.txt', 0);
[pt_zz nt_zz ptci_zz]   = getML_performanceAsso('ZZTRain_psy.txt', 0);
[Be_cy, Bed_cy]   = getML_psyPerformanceCTFixLapse('CyTRain_psy.txt', 0, 1-pt_cy);
[Be_zz, Bed_zz]   = getML_psyPerformanceCTFixLapse('ZZTRain_psy.txt', 0, 1-pt_zz);

ind    = find(mbinss>=-500 & mbinss<=-100);
sesbn = 15;
Monk   = 'Cy';
[rM, rSD] = getML_mRateSac([Monk 'TRain_LIP.txt'], bins, crtf, nvf, dirf, 0);
r_cy   = [];
for i = ind'
    r_cy   = cat(3,r_cy,squeeze(rM(:,i,:)));
end
r_cy   = nanmean(r_cy,3);
a      = getML_txt([Monk 'TRain_LIP.txt']);
fn     = a.data{strcmp(a.name,'dat_fn')};
ses_cy = a.data{strcmp(a.name,'session')};
[b bi p_cy0 r rp] = getLinearDependence(ses_cy, nanmean(r_cy(1:7,:))); 
[b bi p_cy1 r rp] = getLinearDependence(ses_cy(ses_cy<=sesbn), nanmean(r_cy(1:7,(ses_cy<=sesbn)))); 
[b bi p_cy2 r rp] = getLinearDependence(ses_cy(ses_cy>sesbn), nanmean(r_cy(1:7,(ses_cy>sesbn)))); 
r_cy = nanmean(r_cy(1:7,:));
L    = ses_cy<15 & ~isnan(r_cy');    
[h p_cy3] = ttest2(r_cy(L),r_cy(~L));


Monk   = 'ZZ';
[rM, rSD] = getML_mRateSac([Monk 'TRain_LIP.txt'], bins, crtf, nvf, dirf, 0);
r_zz   = [];
for i = ind'
    r_zz   = cat(3,r_zz,squeeze(rM(:,i,:)));
end
r_zz   = nanmean(r_zz,3);
a      = getML_txt([Monk 'TRain_LIP.txt']);
fn     = a.data{strcmp(a.name,'dat_fn')};
ses_zz = a.data{strcmp(a.name,'session')};
[b bi p_zz0 r rp] = getLinearDependence(ses_zz, nanmean(r_zz(1:7,:))'); 
[b bi p_zz1 r rp] = getLinearDependence(ses_zz(ses_zz<=sesbn), nanmean(r_zz(1:7,(ses_zz<=sesbn)))); 
[b bi p_zz2 r rp] = getLinearDependence(ses_zz(ses_zz>sesbn), nanmean(r_zz(1:7,(ses_zz>sesbn)))); 
r_zz = nanmean(r_zz(1:7,:));
L    = ses_zz<15 & ~isnan(r_zz');    
[h p_zz3] = ttest2(r_zz(L),r_zz(~L));


[p_cy0 p_cy1 p_cy2 nanmean(r_cy(L)) nanmean(r_cy(~L)) p_cy3]
[p_zz0 p_zz1 p_zz2 nanmean(r_zz(L)) nanmean(r_zz(~L)) p_zz3]





%% 1a: plot example fits
th = annotation('textbox', [0.1-0.1 0.96 0.05 0.05]);
set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'String', 'a', 'FontSize', ffs, 'FontWeight', 'bold', 'FontName', 'Helvetica') 

for i = 1:3
    % align to dots onset
    ah(2*(i-1)+1) = axes;
    set(ah(2*(i-1)+1), 'Units', 'Normalized', 'Position', [0.08+(i-1)*(0.9*asize+0.025)*ph/pw+(i-1)*(0.5*asize+0.05)*ph/pw 0.76 0.9*asize*ph/pw asize], 'FontName', 'Helvetica')
    plotML_rateByDirCoh(mbins/1000, r_avg(:,:,i), gca, lw)
    set(ah(2*(i-1)+1), 'LineWidth', alw, 'FontSize', afs, 'Box', 'off', ...
        'XLim', [0 0.8], 'XTickLabel', [0 0.2 0.4 0.6])
    xlabel([])
    ylabel([])
    text(0.05, 42, [sprintf('Sessions %.0f-%.0f\n',Bses(i,1), Bses(i,2)) '\itn\rm' sprintf(' = %.0f', rN(i))], ...
        'FontSize', tfs-1, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')
    if 2*(i-1)+1>1
        set(gca, 'YLim', [10 42], 'YTick', [10 20 30 40], 'YTickLabel', [], ...
            'XLim', [0 0.8], 'XTick', [0 0.2 0.4 0.6])
    else
        set(gca, 'YLim', [10 42], 'YTick', [10 20 30 40], ...
            'XLim', [0 0.8], 'XTick', [0 0.2 0.4 0.6])
        ylabel('Response (spikes per s)', 'FontSize', lfs)
        xlabel('Time (s)', 'FontSize', lfs)
    end
    
    
    % align to sac onset
    ah(2*(i-1)+2) = axes;
    set(ah(2*(i-1)+2), 'Units', 'Normalized', 'Position', [0.08+i*(0.9*asize+0.025)*ph/pw+(i-1)*(0.5*asize+0.05)*ph/pw 0.76 0.45*asize*ph/pw asize], 'FontName', 'Helvetica')
    hold on
    plot([0 0], [10 42], 'k', 'LineWidth', alw)
    plot([-0.0075 0.0075], [10 10], 'k', 'LineWidth', alw)
    plot([-0.0075 0.0075], [20 20], 'k', 'LineWidth', alw)
    plot([-0.0075 0.0075], [30 30], 'k', 'LineWidth', alw)
    plot([-0.0075 0.0075], [40 40], 'k', 'LineWidth', alw)
    plotML_rateByDirCoh(mbinss/1000, rs_avg(:,:,i), gca, lw)
    hold off
    set(ah(2*(i-1)+2), 'LineWidth', alw, 'FontSize', afs, 'Box', 'off', ...
        'XLim', [-0.2 0.2], 'XTickLabel', [-0.2 0 0.2])
    xlabel([])
    ylabel([])
    set(gca, 'YLim', [10 42], 'YTick', [], 'YTickLabel', [], ...
        'XLim', [-0.2 0.2], 'XTick', [-0.2 0 0.2], 'YColor', [1 1 1])
    
    
end

%% LIP
% cyrus
th = annotation('textbox', [0 0.62 0.05 0.05]);
set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'String', 'b', 'FontSize', ffs, 'FontWeight', 'bold', 'FontName', 'Helvetica') 

set(axes, 'Units', 'Normalized', 'Position', [0.08+0.3*(asize+0.071)*ph/pw 0.48 asize*ph/pw 0.75*asize], 'FontName', 'Helvetica')
hold on
% cy
plot([1:max(bses_cy)],repmat(0,1,max(bses_cy)),':k')
plot([lses_cy'; lses_cy'],[l_cy(2,:)-ls_cy(2,:); l_cy(2,:)+ls_cy(2,:)],'-', 'LineWidth', elw, 'Color', lce)
plot(lses_cy,l_cy(2,:),'ok', 'MarkerSize', ms, 'MarkerFaceColor', lcc, 'MarkerEdgeColor', 'none')
Lgd = ~isnan(l_cy(2,:)');
[b bi h p] = nestedFW(l_cy(2,Lgd)', ls_cy(2,Lgd)', [ones(sum(Lgd),1) lses_cy(Lgd)]);
if p<0.05
    plot(1:max(bses_cy), b(1)+b(2)*[1:max(bses_cy)], '-k', 'LineWidth', lw)
else
%    plot(1:max(bses_cy), b(1)+b(2)*[1:max(bses_cy)], '--b', 'LineWidth', lw)
end
hold off
set(gca, 'xlim', [1 max(bses_cy)], 'XTick', [1 50 100 150], 'XTickLabel', [], ...
         'ylim', [-5 20], 'YTick', [0 10 20], 'YTickLabel', [0 10 20], ...
         'FontSize', afs, 'LineWidth', alw)
text(10, 20, 'Monkey C', 'FontSize', tfs-1, ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'Top')
ylabel('% coh^{-1}', 'FontSize', lfs)
text(-90, 7.5, 'Coh (k_1)', 'FontSize', lfs+2, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Rotation', 90)


set(axes, 'Units', 'Normalized', 'Position', [0.08+0.3*(asize+0.071)*ph/pw 0.48-1*(0.75*asize+0.02) asize*ph/pw 0.75*asize], 'FontName', 'Helvetica')
hold on
% cy
plot([1:max(bses_cy)],repmat(0,1,max(bses_cy)),':k')
plot([lses_cy'; lses_cy'],[l_cy(3,:)-ls_cy(3,:); l_cy(3,:)+ls_cy(3,:)],'-', 'LineWidth', elw, 'Color', lce)
plot(lses_cy,l_cy(3,:),'ok', 'MarkerSize', ms, 'MarkerFaceColor', lcc, 'MarkerEdgeColor', 'none')
Lgd = ~isnan(l_cy(3,:)');
[b bi h p] = nestedFW(l_cy(3,Lgd)', ones(size(ls_cy(3,Lgd)))', [ones(sum(Lgd),1) lses_cy(Lgd)]);
if p<0.05
    plot(1:max(bses_cy), b(1)+b(2)*[1:max(bses_cy)], '-k', 'LineWidth', lw)
else
%    plot(1:max(bses_cy), b(1)+b(2)*[1:max(bses_cy)], '--b', 'LineWidth', lw)
end
hold off
set(gca, 'xlim', [1 max(bses_cy)], 'XTick', [1 50 100 150], 'XTickLabel', [], ...
         'ylim', [-5 20], 'YTick', [0 10 20], 'YTickLabel', [0 10 20], ...
         'FontSize', afs, 'LineWidth', alw)
ylabel('s^{-1}', 'FontSize', lfs)
text(-90, 7.5, 'Time (k_2)', 'FontSize', lfs+2, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Rotation', 90)



set(axes, 'Units', 'Normalized', 'Position', [0.08+0.3*(asize+0.071)*ph/pw 0.48-2*(0.75*asize+0.02) asize*ph/pw 0.75*asize], 'FontName', 'Helvetica')
hold on
% cy
plot([1:max(bses_cy)],repmat(0,1,max(bses_cy)),':k')
plot([lses_cy'; lses_cy'],[l_cy(4,:)-ls_cy(4,:); l_cy(4,:)+ls_cy(4,:)],'-', 'LineWidth', elw, 'Color', lce)
plot(lses_cy,l_cy(4,:),'ok', 'MarkerSize', ms, 'MarkerFaceColor', lcc, 'MarkerEdgeColor', 'none')
Lgd = ~isnan(l_cy(4,:)');
[b bi h p] = nestedFW(l_cy(4,Lgd)', ones(size(ls_cy(4,Lgd)')), [ones(sum(Lgd),1) lses_cy(Lgd)]);
if p<0.05
    plot(1:max(bses_cy), b(1)+b(2)*[1:max(bses_cy)], '-k', 'LineWidth', lw)
else
%    plot(1:max(bses_cy), b(1)+b(2)*[1:max(bses_cy)], '--b', 'LineWidth', lw)
end
hold off
set(gca, 'xlim', [1 max(bses_cy)], 'XTick', [1 50 100 150], 'XTickLabel', [1 50 100 150], ...
         'ylim', [-5 20], 'YTick', [0 10 20], 'YTickLabel', [0 10 20], ...
         'FontSize', afs, 'LineWidth', alw)
text(-90, 7.5, 'Coh \times time (k_3)', 'FontSize', lfs+2, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Rotation', 90)
xlabel('Session (d)', 'FontSize', lfs)
ylabel('(% coh \times s)^{-1}', 'FontSize', lfs)


 
 
 
%% zsa zsa 
set(axes, 'Units', 'Normalized', 'Position', [0.08+1.4*(asize+0.071)*ph/pw 0.48-0*(0.75*asize+0.02) asize*ph/pw 0.75*asize], 'FontName', 'Helvetica')
hold on
plot([1:max(bses_cy)],repmat(0,1,max(bses_cy)),':k')
% zz
plot([lses_zz'; lses_zz'],[l_zz(2,:)-ls_zz(2,:); l_zz(2,:)+ls_zz(2,:)],'-', 'LineWidth', elw, 'Color', lce)
plot(lses_zz,l_zz(2,:),'ok', 'MarkerSize', ms, 'MarkerFaceColor', lcz, 'MarkerEdgeColor', 'none')
Lgd = ~isnan(l_zz(2,:)');
[b bi h p] = nestedFW(l_zz(2,Lgd)', ones(size(ls_zz(2,Lgd)))', [ones(sum(Lgd),1) lses_zz(Lgd)]);
if p<0.05
    plot(1:max(bses_zz), b(1)+b(2)*[1:max(bses_zz)], '-k', 'LineWidth', lw)
else
%    plot(1:max(bses_zz), b(1)+b(2)*[1:max(bses_zz)], '--r', 'LineWidth', lw)
end
hold off
set(gca, 'xlim', [1 max(bses_zz)], 'XTick', [1 40 80 120], 'XTickLabel', [], ...
         'ylim', [-5 10], 'YTick', [0 5 10], 'YTickLabel', [0 5 10], ...
         'FontSize', afs, 'LineWidth', alw)
text(10, 10, 'Monkey Z', 'FontSize', tfs-1, ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'Top')


set(axes, 'Units', 'Normalized', 'Position', [0.08+1.4*(asize+0.071)*ph/pw 0.48-1*(0.75*asize+0.02) asize*ph/pw 0.75*asize], 'FontName', 'Helvetica')
hold on
plot([1:max(bses_cy)],repmat(0,1,max(bses_cy)),':k')
% zz
plot([lses_zz'; lses_zz'],[l_zz(3,:)-ls_zz(3,:); l_zz(3,:)+ls_zz(3,:)],'-', 'LineWidth', elw, 'Color', lce)
plot(lses_zz,l_zz(3,:),'ok', 'MarkerSize', ms, 'MarkerFaceColor', lcz, 'MarkerEdgeColor', 'none')
Lgd = ~isnan(l_zz(3,:)');
[b bi h p] = nestedFW(l_zz(3,Lgd)', ones(size(ls_zz(3,Lgd)))', [ones(sum(Lgd),1) lses_zz(Lgd)]);
if p<0.05
    plot(1:max(bses_zz), b(1)+b(2)*[1:max(bses_zz)], '-k', 'LineWidth', lw)
else
%    plot(1:max(bses_zz), b(1)+b(2)*[1:max(bses_zz)], '--r', 'LineWidth', lw)
end
hold off
set(gca, 'xlim', [1 max(bses_zz)], 'XTick', [1 40 80 120], 'XTickLabel', [], ...
         'ylim', [-5 10], 'YTick', [0 5 10], 'YTickLabel', [0 5 10], ...
         'FontSize', afs, 'LineWidth', alw)


set(axes, 'Units', 'Normalized', 'Position', [0.08+1.4*(asize+0.071)*ph/pw 0.48-2*(0.75*asize+0.02) asize*ph/pw 0.75*asize], 'FontName', 'Helvetica')
hold on
plot([1:max(bses_cy)],repmat(0,1,max(bses_cy)),':k')
% zz
plot([lses_zz'; lses_zz'],[l_zz(4,:)-ls_zz(4,:); l_zz(4,:)+ls_zz(4,:)],'-', 'LineWidth', elw, 'Color', lce)
plot(lses_zz,l_zz(4,:),'ok', 'MarkerSize', ms, 'MarkerFaceColor', lcz, 'MarkerEdgeColor', 'none')
Lgd = ~isnan(l_zz(4,:)');
[b bi h p] = nestedFW(l_zz(4,Lgd)', ones(size(ls_zz(4,Lgd)')), [ones(sum(Lgd),1) lses_zz(Lgd)]);
if p<0.05
    plot(1:max(bses_zz), b(1)+b(2)*[1:max(bses_zz)], '-k', 'LineWidth', lw)
else
%    plot(1:max(bses_zz), b(1)+b(2)*[1:max(bses_zz)], '--r', 'LineWidth', lw)
end
hold off
set(gca, 'xlim', [1 max(bses_zz)], 'XTick', [1 40 80 120], 'XTickLabel', [1 40 80 120], ...
         'ylim', [-5 10], 'YTick', [0 5 10], 'YTickLabel', [0 5 10], ...
         'FontSize', afs, 'LineWidth', alw)
xlabel('Session (d)', 'FontSize', lfs)


 
 

% %%
% %%%% plot LIP slope as a function of session, and LIP mean vs variance and model
% %
% th = annotation('textbox', [0 0.22 0.05 0.05]);
% set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
%     'String', 'c', 'FontSize', ffs, 'FontWeight', 'bold') 
% 
% % cy slope vs ses
% set(axes, 'Units', 'Normalized', 'Position', [0.1+0.2*(asize+0.071)*ph/pw 0.225 asize*ph/pw 0.75*asize])
% [f fi s] = getML_LIPslope('Cy', 0);
% c = [1 3 6 13 26 51 99];
% hold on
% plot([0:100], zeros(1,101), 'k:', 'LineWidth', elw)
% plot([c; c],[f-fi; f+fi], 'k-', 'LineWidth', elw, 'Color', 'k')
% plot(c, f, 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerSize', ms+2)
% plot(c(s<0.05), f(s<0.05)+fi(s<0.05)+1e-3, '*', 'MarkerEdgeColor', 'k', 'MarkerSize', ms+1)
% hold off
% set(gca, 'xscale', 'log', 'xlim', [0.7 99], 'XTick', c, 'XTickLabel', cellstr(num2str([0 c(2:end)]','%.0f')), ...
%          'ylim', [-0.01 0.07], 'YTick', [0 0.035 0.07], 'FontSize', afs)
% text(1, 0.07, 'Monkey C', 'FontSize', tfs, ...
%     'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')
% ylabel('s^{-1} Session^{-1}', 'FontSize', lfs)
% xlabel(sprintf('%% coherence'), 'FontSize', lfs)
% 
% 
% % zz slope vs ses
% set(axes, 'Units', 'Normalized', 'Position', [0.1+1.4*(asize+0.071)*ph/pw 0.225 asize*ph/pw 0.75*asize])
% [f fi s] = getML_LIPslope('ZZ', 0);
% c = [1 3 6 13 26 51 99];
% hold on
% plot([0:100], zeros(1,101), 'k:', 'LineWidth', elw)
% plot([c; c],[f-fi; f+fi], 'k-', 'LineWidth', elw, 'Color', 'k')
% plot(c, f, 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerSize', ms+2)
% plot(c(s<0.05), f(s<0.05)+fi(s<0.05)+5e-4, '*', 'MarkerEdgeColor', 'k', 'MarkerSize', ms+1)
% hold off
% set(gca, 'xscale', 'log', 'xlim', [0.7 99], 'XTick', c, 'XTickLabel', cellstr(num2str([0 c(2:end)]','%.0f')), ...
%          'ylim', [-0.02 0.045], 'YTick', [0 0.02 0.04], 'FontSize', afs)
% text(1, 0.045, 'Monkey Z', 'FontSize', tfs, ...
%     'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')
% xlabel(sprintf('%% coherence'), 'FontSize', lfs)
% 



% %% plot LIP mean vs var
% th = annotation('textbox', [0.62 0.62 0.05 0.05]);
% set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
%     'String', 'c', 'FontSize', ffs, 'FontWeight', 'bold') 
% 
% 
% % cy LIP mean vs var
% set(axes, 'Units', 'Normalized', 'Position', [0.75 0.41 asize*ph/pw asize])
% [R m v]      = getML_LIPSlopeDist('Cy',0);
% [xsim, ysim] = getML_simLIPSlopeMvsV('Cy', 1000, 0);
% 
% lc  = [];
% for i = 1:7
%     lc = [lc; [0.9 0.9 0.9]-(i-1)/7];
% end
% 
% W = 20;
% hold on
% for i=1:7
%     plot(m(:,i), v(:,i), 'o', 'MarkerFaceColor', lc(i,:), ...
%                                 'MarkerEdgeColor', 'none', 'MarkerSize', ms-1)
%     plot(nanrunmean(xsim(:,i),W), nanrunmean(ysim(:,i),W), '-', 'Color', lc(i,:), 'LineWidth', lw-0.5)
%     plot(m(:,i+7), v(:,i+7), 'o', 'MarkerFaceColor', lc(i,:), ...
%                                 'MarkerEdgeColor', 'none', 'MarkerSize', ms-1)
%     plot(nanrunmean(xsim(:,i+7),W), nanrunmean(ysim(:,i+7),W), '-', 'Color', lc(i,:), 'LineWidth', lw-0.5)
% end
% x   = m(:);
% y   = v(:);
% Lgd = ~isnan(x) & ~isnan(y);
% b   = regress(y(Lgd),[ones(size(x(Lgd))) x(Lgd)]);
% plot([1:120],b(1)+b(2)*[1:120], 'LineWidth', lw/2, 'Color', 'r', 'LineStyle', '--')
% hold off
% set(gca, 'FontSize', afs, ...
%           'xlim', [-20 160], 'XTick', [0 80 160] , 'XTickLabel', [0 80 160], ...
%           'ylim', [0 12e3], 'YTick', [0 6e3 12e3] , 'YTickLabel', [0 6e3 12e3])
% text(-10, 12000, 'Monkey C', 'FontSize', tfs-1, ...
%     'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')
% 
% 
% 
%       
% % zz LIP mean vs var
% set(axes, 'Units', 'Normalized', 'Position', [0.75 0.12132 asize*ph/pw asize])
% [R m v]      = getML_LIPSlopeDist('ZZ',0);
% [xsim, ysim] = getML_simLIPSlopeMvsV('ZZ', 1000, 0);
% 
% lc  = [];
% for i = 1:7
%     lc = [lc; [0.9 0.9 0.9]-(i-1)/7];
% end
% 
% hold on
% for i=1:7
%     plot(m(:,i), v(:,i), 'o', 'MarkerFaceColor', lc(i,:), ...
%                                 'MarkerEdgeColor', 'none', 'MarkerSize', ms-1)
%     plot(nanrunmean(xsim(:,i),W), nanrunmean(ysim(:,i),W), '-', 'Color', lc(i,:), 'LineWidth', lw-0.5)
%     plot(m(:,i+7), v(:,i+7), 'o', 'MarkerFaceColor', lc(i,:), ...
%                                 'MarkerEdgeColor', 'none', 'MarkerSize', ms-1)
%     plot(nanrunmean(xsim(:,i+7),W), nanrunmean(ysim(:,i+7),W), '-', 'Color', lc(i,:), 'LineWidth', lw-0.5)
% end
% x   = m(:);
% y   = v(:);
% y   = y(x>=myprctile(x,3) & x<=myprctile(x,97));
% x   = x(x>=myprctile(x,3) & x<=myprctile(x,97));
% Lgd = ~isnan(x) & ~isnan(y);
% b   = regress(y(Lgd),[ones(size(x(Lgd))) x(Lgd)]);
% plot([1:120],b(1)+b(2)*[1:120], 'LineWidth', lw/2, 'Color', 'r', 'LineStyle', '--')
% hold off
% set(gca, 'FontSize', afs, ...
%           'xlim', [-20 160], 'XTick', [0 80 160] , 'XTickLabel', [0 80 160], ...
%           'ylim', [0 12e3], 'YTick', [0 6e3 12e3] , 'YTickLabel', [0 6e3 12e3])
% text(-10, 12000, 'Monkey Z', 'FontSize', tfs-1, ...
%     'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')
% xlabel(sprintf('slope, mean (spikes/s^2)'), 'FontSize', lfs)
% ylabel(sprintf('slope, variance\n(spikes^2/s^4)'), 'FontSize', lfs)
% 

%% plot LIP mean vs var
th = annotation('textbox', [0.62 0.62 0.05 0.05]);
set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'String', 'c', 'FontSize', ffs, 'FontWeight', 'bold', 'FontName', 'Helvetica') 


% cy slope vs ses
set(axes, 'Units', 'Normalized', 'Position', [0.75 0.41 asize*ph/pw asize], 'FontName', 'Helvetica')
[f fi s] = getML_LIPslope('Cy', 0);
c = [1 3 6 13 26 51 99];
hold on
plot([0:100], zeros(1,101), 'k:', 'LineWidth', elw)
plot([c; c],[f-fi; f+fi], '-', 'LineWidth', elw, 'Color', lce)
plot(c, f, 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerSize', ms+2)
plot(c(s<0.05), f(s<0.05)+fi(s<0.05)+1e-3, '*', 'MarkerEdgeColor', 'k', 'MarkerSize', ms+1)
hold off
set(gca, 'xscale', 'log', 'xlim', [0.7 99], 'XTick', c, 'XTickLabel', [], ...
         'ylim', [-1e-3 40e-3], 'YTick', [0 20e-3 40e-3], 'FontSize', afs)
text(0.8, 40e-3, 'Monkey C', 'FontSize', tfs, ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')


% zz slope vs ses
set(axes, 'Units', 'Normalized', 'Position', [0.75 0.12132 asize*ph/pw asize], 'FontName', 'Helvetica')
[f fi s] = getML_LIPslope('ZZ', 0);
c = [1 3 6 13 26 51 99];
hold on
plot([0:100], zeros(1,101), 'k:', 'LineWidth', elw)
plot([c; c],[f-fi; f+fi], '-', 'LineWidth', elw, 'Color', lce)
plot(c, f, 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerSize', ms+2)
plot(c(s<0.05), f(s<0.05)+fi(s<0.05)+5e-4, '*', 'MarkerEdgeColor', 'k', 'MarkerSize', ms+1)
hold off
set(gca, 'xscale', 'log', 'xlim', [0.7 99], 'XTick', c, 'XTickLabel', cellstr(num2str([0 c(2:end)]','%.0f')), ...
         'ylim', [-15e-3 45e-3], 'YTick', [0 20e-3 40e-3], 'FontSize', afs)
text(0.8, 45e-3, 'Monkey Z', 'FontSize', tfs, ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')
%ylabel('spikes/s/session', 'FontSize', lfs)
xlabel(sprintf('Percentage coherence'), 'FontSize', lfs)
text(0.08, 0.06, 'Spikes per s per session', 'FontSize', lfs+2, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Rotation', 90)
    


