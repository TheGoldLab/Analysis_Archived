%% plotML_fig2.m
% Matlab script for plotting Fig 2  (Behavioral data during training)
% Final version for publication to Nature Neuroscience
%
%   a) behavioral performance as a function of time, example sessions
%   b) behavioral thresholds as a function of time, example sessions
%   c) behavioral thresholds as a function of training sessions
%   d) changes of thresholds as a function of sessions
%
%
% created by jcl on 1/16/07

%%
% create figure
fh = figure; 
mm2inch = 25.4;
pw = 85/mm2inch;   % paper width, make it nature size of one column, 85mm x 85mm
ph = 168/mm2inch;   % paper height
wysifig(fh, pw, ph) % set it to US letter size

% define figure parameters
alw = 0.5;  % axis line width
afs = 5;    % axis font size
lfs = 6;    % label font size
ffs = 12;   % figure font size
tfs = 6;    % text font size
lw  = 1.5;  % line width
elw = 0.75;
ms  = 2;    % marker size

% define line color
lc1  = [0.1 0.1 0.1];
lc1e = [0.75 0.75 0.75];
lc2  = [1 0.1 0.1];
lc2e = [1 0.75 0.75];


asize = 0.35;


% load data
% get threshold fitted by power function
[p_cy n_cy pci_cy]        = getML_performanceAsso('CyTRain_asso.txt', 0);
[pt_cy nt_cy ptci_cy]     = getML_performanceAsso('CyTRain_psy.txt', 0);
N_cy                      = length(p_cy);
[p_zz n_zz pci_zz]        = getML_performanceAsso('ZZTRain_asso.txt', 0);
[pt_zz nt_zz ptci_zz]     = getML_performanceAsso('ZZTRain_psy.txt', 0);
N_zz                      = length(p_zz);
[fits_cy, sems_cy, th_cy] = getML_psyPerformanceCTFixLapse('CyTRain_psy.txt', 0, 1-pt_cy);
[fits_zz, sems_zz, th_zz] = getML_psyPerformanceCTFixLapse('ZZTRain_psy.txt', 0, 1-pt_zz);
% Lp_cy = getML_LapseSelectionArray('Cy');
% Lp_zz = getML_LapseSelectionArray('ZZ');
Lp_cy = logical(ones(size(getML_LapseSelectionArray('Cy'))));
Lp_zz = logical(ones(size(getML_LapseSelectionArray('ZZ'))));



%% Fig 1a: Example fits - performance vs time
th = annotation('textbox', [0 0.96 0.05 0.05]);
set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'String', 'a', 'FontSize', ffs, 'FontWeight', 'bold', 'FontName', 'Helvetica') 

EgInd  = [11 142];

fn     = 'CyTRain_psy.txt';
a      = getML_txt(fn);
fname  = a.data{strcmp(a.name,'dat_fn')};
global FIRA

% early
set(axes, 'Units', 'Normalized', 'Position', [0.115 0.79 asize asize*pw/ph], 'FontName', 'Helvetica')
i = EgInd(1);
fprintf('%d: %s\n', i, fname{i})
openFIRA(fname{i})

warning off
bb   = 0;
be   = 1500;
bw   = 300;
bs   = 150;
bins = [[bb:bs:be-bw]' [bb+bw:bs:be]'];
mbin = mean(bins,2)./1000;
coh  = [0 3.2 6.4 12.8 25.6 51.2 99.9]';

% get performance for each coh and time
p = nans(length(coh),size(bins,1));
N = nans(length(coh),size(bins,1));
for j = 1:size(bins,1)
    Lt = FIRA.ecodes.data(:,getFIRA_ecodeColumnByName('dot_dur'))>=bins(j,1) ...
        & FIRA.ecodes.data(:,getFIRA_ecodeColumnByName('dot_dur'))<bins(j,2);
    for k = 1:length(coh)
        Lcoh = FIRA.ecodes.data(:,getFIRA_ecodeColumnByName('dot_coh'))==coh(k);
        N(k,j) = nansum(FIRA.ecodes.data(Lcoh&Lt,getFIRA_ecodeColumnByName('correct'))>=0);
        p(k,j) = nansum(FIRA.ecodes.data(Lcoh&Lt,getFIRA_ecodeColumnByName('correct'))==1)/N(k,j);
    end
end

% prepare plot parameters
lc  = [];
for j = 1:7
    lc = [lc; [0.9 0.9 0.9]-(j-1)/7];
end
% lc = colormap(bone(8));
% lc(end,:) = [];
% lc = flipud(lc);

hold on
for j = 1:7
    plot(mbin, p(j,:)', 'o', 'MarkerSize', ms+2, 'MarkerFaceColor', lc(j,:), 'MarkerEdgeColor', 'none')
    ym = 0.5 + (1 - 0.5 - (1-pt_cy(i))).*(1-exp(-((coh(j)/100)./(fits_cy(1,i).*mbin.^fits_cy(2,i))).^fits_cy(3,i)));
    plot(mbin, ym, 'LineWidth', lw, 'Color', lc(j,:));
end
hold off
set(gca, 'LineWidth', alw, 'FontSize', afs, 'Box', 'off', ...
         'XLim', [0 1], 'YLim', [0.38 1.02], ...
         'XTick', [0 0.5 1.0], 'XTickLabel', {'0', '0.5', '1.0'}, ...
         'YTick', [0.4 0.6 0.8 1], 'YTickLabel', [40 60 80 100])
warning on
xlabel('Time (s)', 'FontSize', lfs)
ylabel('Percentage correct', 'FontSize', lfs)
title('Early', 'FontSize', tfs+2)
text(0.98, 0.38, ['Session ' int2str(i)], 'FontSize', tfs-1, ...
     'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom')
% plot annotation arrow 
[xb,yb] = axisXY2figXYslog(0.6,1);
[xe,ye] = axisXY2figXYslog(0,ym(end));
anh     = annotation('arrow', [xb xe], [yb ye]);
set(anh, 'Color', 'k', 'LineWidth', lw/2, 'LineStyle', '--', 'HeadStyle', 'cback2', ...
          'HeadWidth', 3*ms, 'HeadLength', 3*ms)

 
 
% late
set(axes, 'Units', 'Normalized', 'Position', [0.53 0.79 asize asize*pw/ph], 'FontName', 'Helvetica')
i = EgInd(2);
fprintf('%d: %s\n', i, fname{i})
openFIRA(fname{i})

warning off
bb   = 0;
be   = 1500;
bw   = 300;
bs   = 150;
bins = [[bb:bs:be-bw]' [bb+bw:bs:be]'];
mbin = mean(bins,2)./1000;
coh  = [0 3.2 6.4 12.8 25.6 51.2 99.9]';

% get performance for each coh and time
p = nans(length(coh),size(bins,1));
N = nans(length(coh),size(bins,1));
for j = 1:size(bins,1)
    Lt = FIRA.ecodes.data(:,getFIRA_ecodeColumnByName('dot_dur'))>=bins(j,1) ...
        & FIRA.ecodes.data(:,getFIRA_ecodeColumnByName('dot_dur'))<bins(j,2);
    for k = 1:length(coh)
        Lcoh = FIRA.ecodes.data(:,getFIRA_ecodeColumnByName('dot_coh'))==coh(k);
        N(k,j) = nansum(FIRA.ecodes.data(Lcoh&Lt,getFIRA_ecodeColumnByName('correct'))>=0);
        p(k,j) = nansum(FIRA.ecodes.data(Lcoh&Lt,getFIRA_ecodeColumnByName('correct'))==1)/N(k,j);
    end
end



hold on
for j = 1:7
    plot(mbin, p(j,:)', 'o', 'MarkerSize', ms+2, 'MarkerFaceColor', lc(j,:), 'MarkerEdgeColor', 'none')
    ym = 0.5 + (1 - 0.5 - (1-pt_cy(i))).*(1-exp(-((coh(j)/100)./(fits_cy(1,i).*mbin.^fits_cy(2,i))).^fits_cy(3,i)));
    plot(mbin, ym, 'LineWidth', lw, 'Color', lc(j,:));
end
hold off
set(gca, 'LineWidth', alw, 'FontSize', afs, 'Box', 'off', ...
         'XLim', [0 1], 'YLim', [0.38 1.02], ...
         'XTick', [0 0.5 1], 'XTickLabel', {'0', '0.5', '1.0'}, ...
         'YTick', [0.4 0.6 0.8 1], 'YTickLabel', [])
warning on
title('Late', 'FontSize', tfs+2)
text(0.98, 0.38, ['Session ' int2str(i)], 'FontSize', tfs-1, ...
     'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom')

% plot annotation arrow 
[xb,yb] = axisXY2figXYslog(0.15,1);
[xe,ye] = axisXY2figXYslog(0,ym(end));
anh     = annotation('arrow', [xb xe], [yb ye]);
set(anh, 'Color', 'k', 'LineWidth', lw/2, 'LineStyle', '--', 'HeadStyle', 'cback2', ...
          'HeadWidth', 3*ms, 'HeadLength', 3*ms)



%% Fig 1b: Example fits - threshold vs time
th = annotation('textbox', [0.1-0.1 0.69 0.05 0.05]);
set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'String', 'b', 'FontSize', ffs, 'FontWeight', 'bold', 'FontName', 'Helvetica') 


% prepare data
% get threshold for each bin
[fits_cy, sems_cy, th_cy] = getML_psyPerformanceCT('CyTRain_psy.txt', 0);
a   = getML_txt('CyTRain_psy.txt');
fn  = a.data{strcmp(a.name,'dat_fn')};
ses = a.data{strcmp(a.name,'session')}; % choose session 142 or 138


% early
ind = EgInd(1);

bb   = 0;
be   = 1500;
bw   = 300;
bs   = 150;
bins = [[bb:bs:be-bw]' [bb+bw:bs:be]'];
mbin = mean(bins,2);
coh  = [0 3.2 6.4 12.8 25.6 51.2 99.9]';

% get performance for each coh and time
global FIRA
openFIRA(fn{ind})
fprintf('%s\n', fn{ind})
p = nans(length(coh),size(bins,1));
N = nans(length(coh),size(bins,1));
for i = 1:size(bins,1)
    Lt = FIRA.ecodes.data(:,getFIRA_ecodeColumnByName('dot_dur'))>=bins(i,1) ...
        & FIRA.ecodes.data(:,getFIRA_ecodeColumnByName('dot_dur'))<bins(i,2);         
    for j = 1:length(coh)
        Lcoh = FIRA.ecodes.data(:,getFIRA_ecodeColumnByName('dot_coh'))==coh(j);
        N(j,i) = sum(FIRA.ecodes.data(Lcoh&Lt,getFIRA_ecodeColumnByName('correct'))>=0);
        p(j,i) = sum(FIRA.ecodes.data(Lcoh&Lt,getFIRA_ecodeColumnByName('correct'))==1)/N(j,i);
    end
end


recompute = 0;
if recompute
    fits  = nans(2,size(bins,1));
    sems  = nans(2,2,size(bins,1));

    for i = 1:size(bins,1)
        i
        Lt  = FIRA.ecodes.data(:,getFIRA_ecodeColumnByName('dot_dur'))>=bins(i,1) ...
            & FIRA.ecodes.data(:,getFIRA_ecodeColumnByName('dot_dur'))<bins(i,2);
        Lcr = FIRA.ecodes.data(:,getFIRA_ecodeColumnByName('correct'))>=0;
        
        if sum(Lt&Lcr)>5
            [fits(:,i), sems(:,:,i)] = ctPsych_fit(@quick2, FIRA.ecodes.data(Lt&Lcr,getFIRA_ecodeColumnByName('dot_coh'))/100,...
                FIRA.ecodes.data(Lt&Lcr,getFIRA_ecodeColumnByName('correct')), [], {100,68,81.6}, []);
        end
    end
    
    [home_dir, lab_dir, current_dir, tmat_dir] = dirnames;
    save([tmat_dir '/plotML_fig1b-early.mat'], 'fits', 'sems')

else
    [home_dir, lab_dir, current_dir, tmat_dir] = dirnames;
    load([tmat_dir '/plotML_fig1b-early.mat'])
end



% plot thresholds against time
set(axes, 'Units', 'Normalized', 'Position', [0.115 0.54 asize asize*pw/ph], 'FontName', 'Helvetica')
lc  = [];
for j = 1:7
    lc = [lc; [0.9 0.9 0.9]-(j-1)/7];
end

% plot raw data and fitted data
et = 1;
hold on
alpha = fits_cy(1,ind).*[0.01:0.01:et].^fits_cy(2,ind);
%plot([mbin mbin]'/1000, [fits(1,:)-sems(1,:); fits(1,:)+sems(1,:)], 'k')
plot([mbin mbin]'/1000, 100*shiftdim(sems(1,:,:),1), 'Color', lc1e)
plot(mbin/1000, 100*fits(1,:), 'o', 'MarkerSize', ms+2, 'MarkerFaceColor', lc1, 'MarkerEdgeColor', 'none');
plot([0.01:0.01:et], 100*alpha, 'LineWidth', lw, 'Color', 'k');
plot(ones(size([0:0.1:100*alpha(end)])), [0:0.1:100*alpha(end)], 'LineWidth', lw/2, ...
        'Color', 'k', 'LineStyle', '--');
hold off
set(gca, 'LineWidth', alw, 'FontSize', afs, 'Box', 'off', ...
         'XLim', [0 et], 'XTick', [0 0.5 1], 'XTickLabel', {'0', '0.5', '1.0'}, ...
         'YLim', [7.5 100], 'YTick', [10 20 40 80], 'YScale', 'log')   
% plot annotation arrow
[xb,yb] = axisXY2figXYslog(et,100*alpha(end));
[xe,ye] = axisXY2figXYslog(0,100*alpha(end));
anh     = annotation('arrow', [xb xe], [yb ye]);
set(anh, 'Color', 'k', 'LineWidth', lw/2, 'LineStyle', '--', 'HeadStyle', 'cback2', ...
          'HeadWidth', 3*ms, 'HeadLength', 3*ms)
xlabel('Time (s)', 'FontSize', lfs)
ylabel('Threshold, \alpha (% coh)', 'FontSize', lfs)
text(0.95, 7.55, ['Session ' int2str(ind)], 'FontSize', tfs-1,...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom')

 

% late
ind = EgInd(2);

bb   = 0;
be   = 1500;
bw   = 300;
bs   = 150;
bins = [[bb:bs:be-bw]' [bb+bw:bs:be]'];
mbin = mean(bins,2);
coh  = [0 3.2 6.4 12.8 25.6 51.2 99.9]';

% get performance for each coh and time
global FIRA
openFIRA(fn{ind})
fprintf('%s\n', fn{ind})
p = nans(length(coh),size(bins,1));
N = nans(length(coh),size(bins,1));
for i = 1:size(bins,1)
    Lt = FIRA.ecodes.data(:,getFIRA_ecodeColumnByName('dot_dur'))>=bins(i,1) ...
        & FIRA.ecodes.data(:,getFIRA_ecodeColumnByName('dot_dur'))<bins(i,2);         
    for j = 1:length(coh)
        Lcoh = FIRA.ecodes.data(:,getFIRA_ecodeColumnByName('dot_coh'))==coh(j);
        N(j,i) = sum(FIRA.ecodes.data(Lcoh&Lt,getFIRA_ecodeColumnByName('correct'))>=0);
        p(j,i) = sum(FIRA.ecodes.data(Lcoh&Lt,getFIRA_ecodeColumnByName('correct'))==1)/N(j,i);
    end
end


recompute = 0;
if recompute
    fits  = nans(2,size(bins,1));
    sems  = nans(2,2,size(bins,1));

    for i = 1:size(bins,1)
        i
        Lt  = FIRA.ecodes.data(:,getFIRA_ecodeColumnByName('dot_dur'))>=bins(i,1) ...
            & FIRA.ecodes.data(:,getFIRA_ecodeColumnByName('dot_dur'))<bins(i,2);
        Lcr = FIRA.ecodes.data(:,getFIRA_ecodeColumnByName('correct'))>=0;
        
        if sum(Lt&Lcr)>5
            [fits(:,i), sems(:,:,i)] = ctPsych_fit(@quick2, FIRA.ecodes.data(Lt&Lcr,getFIRA_ecodeColumnByName('dot_coh'))/100,...
                FIRA.ecodes.data(Lt&Lcr,getFIRA_ecodeColumnByName('correct')), [], {100,68,81.6}, []);
        end
    end
    [home_dir, lab_dir, current_dir, tmat_dir] = dirnames;
    save([tmat_dir '/plotML_fig1b-late.mat'], 'fits', 'sems')

else
    [home_dir, lab_dir, current_dir, tmat_dir] = dirnames;
    load([tmat_dir '/plotML_fig1b-late.mat'])
end


% plot thresholds against time
set(axes, 'Units', 'Normalized', 'Position', [0.53 0.54 asize asize*pw/ph], 'FontName', 'Helvetica')
lc  = [];
for j = 1:7
    lc = [lc; [0.9 0.9 0.9]-(j-1)/7];
end

% plot raw data and fitted data
et = 1;
hold on
alpha = fits_cy(1,ind).*[0.01:0.01:et].^fits_cy(2,ind);
%plot([mbin mbin]'/1000, [fits(1,:)-sems(1,:); fits(1,:)+sems(1,:)], 'k')
plot([mbin mbin]'/1000, 100*shiftdim(sems(1,:,:),1), 'Color', lc1e)
plot(mbin/1000, 100*fits(1,:), 'o', 'MarkerSize', ms+2, 'MarkerFaceColor', lc1, 'MarkerEdgeColor', 'none');
plot([0.01:0.01:et], 100*alpha, 'LineWidth', lw, 'Color', 'k');
plot(ones(size([0:0.1:100*alpha(end)])), [0:0.1:100*alpha(end)], 'LineWidth', lw/2, ...
        'Color', 'k', 'LineStyle', '--');
hold off
set(gca, 'LineWidth', alw, 'FontSize', afs, 'Box', 'off', ...
         'XLim', [0 et], 'XTick', [0 0.5 1], 'XTickLabel', {'0', '0.5', '1.0'}, ...
         'YLim', [7.5 100], 'YTick', [10 20 40 80], 'YTickLabel', [], 'YScale', 'log')   
% plot annotation arrow
[xb,yb] = axisXY2figXYslog(et,100*alpha(end));
[xe,ye] = axisXY2figXYslog(0,100*alpha(end));
anh     = annotation('arrow', [xb xe], [yb ye]);
set(anh, 'Color', 'k', 'LineWidth', lw/2, 'LineStyle', '--', 'HeadStyle', 'cback2', ...
          'HeadWidth', 3*ms, 'HeadLength', 3*ms)
%xlabel('Time (s)', 'FontSize', lfs)
%ylabel('Threshold (% coh)', 'FontSize', lfs)
text(0.95, 7.55, ['Session ' int2str(ind)], 'FontSize', tfs-1,...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom')

%% Fig 1c and 1d: plot threshold and lapse vs session
th = annotation('textbox', [0.1-0.1 0.439 0.05 0.05]);
set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'String', 'c', 'FontSize', ffs, 'FontWeight', 'bold', 'FontName', 'Helvetica') 


% get threshold fitted by power function
[p_cy n_cy pci_cy]        = getML_performanceAsso('CyTRain_asso.txt', 0);
[pt_cy nt_cy ptci_cy]     = getML_performanceAsso('CyTRain_psy.txt', 0);
N_cy                      = length(p_cy);
[p_zz n_zz pci_zz]        = getML_performanceAsso('ZZTRain_asso.txt', 0);
[pt_zz nt_zz ptci_zz]     = getML_performanceAsso('ZZTRain_psy.txt', 0);
N_zz                      = length(p_zz);
[fits_cy, sems_cy, th_cy] = getML_psyPerformanceCTFixLapse('CyTRain_psy.txt', 0, 1-pt_cy);
[fits_zz, sems_zz, th_zz] = getML_psyPerformanceCTFixLapse('ZZTRain_psy.txt', 0, 1-pt_zz);
% Lp_cy = getML_LapseSelectionArray('Cy');
% Lp_zz = getML_LapseSelectionArray('ZZ');
Lp_cy = logical(ones(size(getML_LapseSelectionArray('Cy'))));
Lp_zz = logical(ones(size(getML_LapseSelectionArray('ZZ'))));


% plot threshold and error bars
set(axes, 'Units', 'Normalized', 'Position', [0.115 0.27 asize asize*pw/ph], 'FontName', 'Helvetica')
hold on
line(repmat([1:size(fits_cy,2)],2,1), shiftdim(sems_cy(1,:,:),1)*100,...
                                'Color', lc1e, 'LineWidth', elw)
plot(find(Lp_cy), fits_cy(1,Lp_cy)*100, 'o', 'MarkerSize', ms, ...
                'MarkerFaceColor', lc1, 'MarkerEdgeColor', 'none')
plot(find(~Lp_cy), fits_cy(1,~Lp_cy)*100, 'o', 'MarkerSize', ms, ...
                'MarkerFaceColor', 'w', 'MarkerEdgeColor', lc1)
% fit and plot threshold with exponential function
ys = shiftdim(sems_cy(1,2,:))-shiftdim(sems_cy(1,1,:));
[b, bsem, gof, ym] = exp_fitW([1:size(fits_cy,2)]', 100*real(fits_cy(1,:))', ones(size(100*ys)), [10; 80; 40], [0 20; 0 100; 0 160]);
[b(1)+b(2)                 b(1)              b(3); ...
 b(1)+b(2)-bsem(1)-bsem(2) b(1)-bsem(1) b(3)-bsem(3); ...
 b(1)+b(2)+bsem(1)+bsem(2) b(1)+bsem(1) b(3)+bsem(3); ...
 ym(end)                   bsem(1)+bsem(2)*exp(-1/bsem(3)*ym(end)) nan]
plot([1:size(fits_cy,2)], ym, 'Color', 'k', 'LineWidth', lw)
hold off

set(gca, 'LineWidth', alw, 'FontSize', afs, 'Box', 'off', ...
         'XLim', [-N_cy, size(fits_cy,2)], 'YLim', [5 100], ...
         'XTick', [1 50 100 150], 'YScale', 'log', 'YTick', [5 10 20 40 80])
%xlabel('Session (d)', 'FontSize', lfs)
ylabel('Threshold, \alpha (% coh)', 'FontSize', lfs)
title('Monkey C', 'FontSize', tfs+2)
% text(0.98*size(fits_cy,2), 98, 'Monkey C', 'FontSize', tfs, ...
%      'HorizontalAlignment', 'right', 'VerticalAlignment', 'top')


set(axes, 'Units', 'Normalized', 'Position', [0.53 0.27 asize asize*pw/ph])
hold on
line(repmat([1:size(fits_zz,2)],2,1), shiftdim(sems_zz(1,:,:),1)*100,...
                'Color', lc1e, 'LineWidth', elw)
plot(find(Lp_zz), fits_zz(1,Lp_zz)*100, 'o', 'MarkerSize', ms, ...
                'MarkerFaceColor', lc1, 'MarkerEdgeColor', 'none')
plot(find(~Lp_zz), fits_zz(1,~Lp_zz)*100, 'o', 'MarkerSize', ms, ...
                'MarkerFaceColor', 'w', 'MarkerEdgeColor', lc1)

% fit and plot threshold with exponential function
ys = shiftdim(sems_cy(1,2,:))-shiftdim(sems_cy(1,1,:));
L  = [1:86 88:130];
[b, bsem, gof, ym] = exp_fitW([1:size(fits_zz(1,L),2)], fits_zz(1,L)*100, ones(size(ys(L))), [10; 80; 100], [0 60; 0 100; 0 100]);
[b(1)+b(2)                 b(1)              b(3); ...
 b(1)+b(2)-bsem(1)-bsem(2) b(1)-bsem(1) b(3)-bsem(3); ...
 b(1)+b(2)+bsem(1)+bsem(2) b(1)+bsem(1) b(3)+bsem(3); ...
 ym(end)                   bsem(1)+bsem(2)*exp(-1/bsem(3)*ym(end)) nan]
plot(L, ym, 'Color', 'k', 'LineWidth', lw)
hold off
set(gca, 'LineWidth', alw, 'FontSize', afs, 'Box', 'off', ...
         'XLim', [-N_zz, size(fits_zz,2)], 'YLim', [5 100], ...
         'XTick', [-N_zz 1 40 80 120], 'YScale', 'log', 'YTick', [5 10 20 40 80], 'YTickLabel', {})
%xlabel('Session (d)', 'FontSize', lfs)
title('Monkey Z', 'FontSize', tfs+2)
% text(0.98*size(fits_zz,2), 98, 'Monkey Z', 'FontSize', tfs, ...
%      'HorizontalAlignment', 'right', 'VerticalAlignment', 'top')

 
% % plot a patch to indicate the pre-training period (99.9% only)
% % cy
% set(axes, 'Units', 'Normalized', 'Position', [0.115 0.1 asize asize])
% p1 = patch([-N_cy 0 0 -N_cy], [0 0 1 1], 'k');  
% %set(p1, 'FaceColor', 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.25)
% set(gca, 'XLim', [-N_cy, size(fits_cy,2)], 'YLim', [0 1], 'Color', 'none')
% axis('off')
% 
% % zz
% set(axes, 'Units', 'Normalized', 'Position', [0.55 0.1 asize asize])
% p1 = patch([-N_zz 0 0 -N_zz], [0 0 1 1], 'k');  
% set(p1, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.25)
% set(gca, 'XLim', [-N_zz, size(fits_zz,2)], 'YLim', [0 1], 'Color', 'none')
% axis('off')

 
% plot error rate at 99.9%
% cyrus
set(axes, 'Units', 'Normalized', 'Position', [0.115 0.27 asize asize*pw/ph], 'FontName', 'Helvetica')
hold on
lx   = [-N_cy:1:-1 1:size(pt_cy,2)];
ly   = 1-[p_cy pt_cy];
lyci = 1-[pci_cy ptci_cy];
lys  = mean(abs([lyci(1,:)-ly; lyci(2,:)-ly]));
lys(lys==0) = 0.001;   % make it a small number
L    = ~isnan(lx) & ~isnan(ly) & ~isnan(lys);
[blc, bsemlc, gof, lym] = exp_fit2W(lx(L), ly(L), lys(L), [0.5; 5], [0.4 0.6;0.01 100]);
line([lx; lx], lyci,'Color', lc2e, 'LineWidth', elw)
plot(lx, ly, '^', 'MarkerSize', ms, 'MarkerFaceColor', lc2, 'MarkerEdgeColor', 'none')
plot(lx(L), lym, 'LineWidth', lw, 'Color', 'r')
%plot([-N_cy:1:160], 0.15*ones(size([-N_cy:1:160])), 'LineWidth', 0.5, 'LineStyle', ':', 'Color', 'r')
hold off
set(gca, 'LineWidth', alw, 'FontSize', afs, 'Box', 'off', ...
         'XLim', [-N_cy, size(fits_cy,2)], 'XTick', [1 50 100 150], ...
         'YLim', [0 1], 'YTick', [0 0.25 0.5 0.75 1], 'YTickLabel', [], 'YColor', 'r', 'YAxisLocation', 'right', ...
         'Color', 'none')

% fit a line to lapse rate after 1.61 time constant of the session to see if
% it changed significantly
tau = round(1.61*blc(2));
L   = ~isnan(lx) & ~isnan(ly) & ~isnan(lys) & logical([zeros(1,tau), ones(1,length(ly)-tau)]);
[b, bint, h, p] = nestedFW(ly(L)', lys(L)', [ones(size(ly(L)')) lx(L)']);



% zsa zsa
set(axes, 'Units', 'Normalized', 'Position', [0.53 0.27 asize asize*pw/ph])
hold on
lx   = [-N_zz:1:-1 1:size(pt_zz,2)];
ly   = 1-[p_zz pt_zz];
lyci = 1-[pci_zz ptci_zz];
lys  = mean(abs([lyci(1,:)-ly; lyci(2,:)-ly]));
lys(lys==0) = 0.001;   % make it a small number
L    = ~isnan(lx) & ~isnan(ly) & ~isnan(lys);
[blz, bsemlz, gof, lym] = exp_fit2W(lx(L), ly(L), lys(L), [0.5; 5], [0.4 0.6;0.01 100]);
line([lx; lx], lyci,'Color', lc2e, 'LineWidth', elw)
plot(lx, ly, '^', 'MarkerSize', ms, 'MarkerFaceColor', lc2, 'MarkerEdgeColor', 'none')
plot(lx(L), lym, 'LineWidth', lw, 'Color', 'r')
%plot([-N_zz:1:130], 0.15*ones(size([-N_zz:1:130])), 'LineWidth', 0.5, 'LineStyle', ':', 'Color', 'r')
hold off
set(gca, 'LineWidth', alw, 'FontSize', afs, 'Box', 'off', ...
         'XLim', [-N_zz, size(fits_zz,2)], 'XTick', [-N_zz 1 40 80 120], ...
         'YLim', [0 1], 'YTick', [0 0.25 0.5 0.75 1], 'YTickLabel', {'0', '0.25', '0.5', '0.75', '1.00'}, ...
         'YColor', 'r', 'YAxisLocation', 'right', ...
         'Color', 'none')
ylabel(['Error rate at 99% coh, \lambda'] , 'FontSize', lfs)





% fit a line to lapse rate after 1.61 time constant of the session to see if
% it changed significantly
tau = round(1.61*blz(2));
L   = ~isnan(lx) & ~isnan(ly) & ~isnan(lys) & logical([zeros(1,tau), ones(1,length(ly)-tau)]);
[b, bint, h, p] = nestedFW(ly(L)', lys(L)', [ones(size(ly(L)')) lx(L)']);


%% see if time exponents change significantly with training
Lgd = fits_cy(2,:)>-3;
ses = 1:160;
y  = fits_cy(2,Lgd)';
ys = shiftdim(sems_cy(1,2,Lgd))-shiftdim(sems_cy(1,1,Lgd));
[b, bint, h, p] = nestedFW(y, ys, [ones(size(y)) [ses(Lgd)]']);

Lgd = fits_zz(2,:)>-3;
ses = 1:160;
y  = fits_zz(2,Lgd)';
ys = shiftdim(sems_zz(1,2,Lgd))-shiftdim(sems_zz(1,1,Lgd));
[b, bint, h, p] = nestedFW(y, ys, [ones(size(y)) [ses(Lgd)]']);



%%
% %%% old code, plot lapse from the fit %%%
% % plot lapse
% set(axes, 'Units', 'Normalized', 'Position', [0.115 0.1 asize asize])
% hold on
% line(repmat([1:size(fits_cy,2)],2,1), shiftdim(sems_cy(4,:,:),1)*100,...
%     'Color', lc2e, 'LineWidth', elw)
% plot([1:size(fits_cy,2)], 100*fits_cy(4,:), 'v', 'MarkerSize', ms, ...
%                                 'MarkerFaceColor', lc2, 'MarkerEdgeColor', 'none')
% plot([1:size(fits_cy,2)], mynanrunmedian(100*fits_cy(4,:),4), 'LineWidth', lw, 'Color', lc2)
% hold off
% set(gca, 'LineWidth', alw, 'FontSize', afs, 'Box', 'off', ...
%          'XLim', [1, size(fits_cy,2)], 'XTick', [1 50 100 150], ...
%          'YLim', [0 50], 'YTick', [0 10 20 30 40 50], 'YTickLabel', [], 'YColor', lc2, 'YAxisLocation', 'right', ...
%          'Color', 'none')
% %ylabel('Lapse rate', 'FontSize', lfs)
% 
% 
% set(axes, 'Units', 'Normalized', 'Position', [0.55 0.1 asize asize])
% hold on
% line(repmat([1:size(fits_zz,2)],2,1), shiftdim(sems_zz(4,:,:),1)*100,...
%     'Color', lc2e, 'LineWidth', elw)
% plot([1:size(fits_zz,2)], 100*fits_zz(4,:), 'v', 'MarkerSize', ms, ...
%                                 'MarkerFaceColor', lc2, 'MarkerEdgeColor', 'none')
% plot([1:size(fits_zz,2)], mynanrunmedian(100*fits_zz(4,:),4), 'LineWidth', lw, 'Color', lc2)
% hold off
% set(gca, 'LineWidth', alw, 'FontSize', afs, 'Box', 'off', ...
%          'XLim', [1, size(fits_zz,2)], 'XTick', [1 40 80 120], ...
%          'YLim', [0 50], 'YColor', lc2, 'YTick', [0 10 20 30 40 50], 'YAxisLocation', 'right', ...
%          'Color', 'none')
% ylabel('Lapse rate', 'FontSize', lfs, 'Color', lc2)



%% Fig 1d: learning rate
th = annotation('textbox', [0.1-0.1 0.2 0.05 0.05]);
set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'String', 'd', 'FontSize', ffs, 'FontWeight', 'bold', 'FontName', 'Helvetica') 

% cy
[rp, rpci, rpp, rpm, rth, rthci, rthp, rthm] = getML_learningRate('Cy', 40, 0);
for i = 1:length(rth)
    if mod(i,10)
        rp(i)    = nan;
        rpci(i)  = nan;
        rth(i)   = nan;
        rthci(i) = nan;
    end
end

set(axes, 'Units', 'Normalized', 'Position', [0.115 0.05 asize asize*pw/ph], 'FontName', 'Helvetica')
% plot threshold and error bars
hold on
plot(-N_cy:length(rth), zeros(1,length(rth)+N_cy+1), ':k')

line(repmat([1:length(rth)],2,1), [rth(1:end)'-rthci(1:end)'; rth(1:end)'+rthci(1:end)'], ...
                                'Color', lc1e, 'LineWidth', elw)
plot(1:length(rth), rth(1:end), 'o', 'MarkerSize', ms, ...
                'MarkerFaceColor', lc1, 'MarkerEdgeColor', 'none')
hold off
set(gca, 'LineWidth', alw, 'FontSize', afs, 'Box', 'off', ...
         'XLim', [-N_cy, length(rth)], 'XTick', [0 50 100 150], ...
         'YLim', [-0.04 0.02], 'YTick', [-0.04 -0.02 0], ...
         'YTickLabel', [-0.04 -0.02 0])
xlabel('Session (d)', 'FontSize', lfs)
ylabel('\Delta\alpha (log(coh) session^{-1})', 'FontSize', lfs)


set(axes, 'Units', 'Normalized', 'Position', [0.115 0.05 asize asize*pw/ph])
hold on
lx   = [-N_cy:1:-1 1:length(rth)];
ly   = rp;
lym  = rpm;
lyci = [rp'-rpci'; rp'+rpci'];
line([lx(1:end); lx(1:end)], lyci(:,1:end),'Color', lc2e, 'LineWidth', elw)
plot(lx(1:end), ly(1:end), '^', 'MarkerSize', ms, 'MarkerFaceColor', lc2, 'MarkerEdgeColor', 'none')
%plot(lx, lym, 'LineWidth', lw, 'Color', 'r')
hold off
set(gca, 'LineWidth', alw, 'FontSize', afs, 'Box', 'off', ...
         'XLim', [-N_cy, length(rth)], 'XTick', [0 50 100 150], ...
         'YLim', [-0.04 0.02], 'YTick', [-0.04 -0.02 0 0.02], ...
         'YTickLabel', [], ...
         'Color', 'none', 'YColor', 'r', 'YAxisLocation', 'right')
% text(0, 1e-2, 'Monkey C', 'FontSize', tfs, ...
%      'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')


     
% zz     
[rp, rpci, rpp, rpm, rth, rthci, rthp, rthm] = getML_learningRate('ZZ', 40, 0);
for i = 1:length(rth)
    if mod(i,10)
        rp(i)    = nan;
        rpci(i)  = nan;
        rth(i)   = nan;
        rthci(i) = nan;
    end
end

% plot threshold and error bars
set(axes, 'Units', 'Normalized', 'Position', [0.53 0.05 asize asize*pw/ph], 'FontName', 'Helvetica')
hold on
plot(-N_zz:length(rth), zeros(1,length(rth)+N_zz+1), ':k')
line(repmat([1:length(rth)],2,1), [rth(1:end)'-rthci(1:end)'; rth(1:end)'+rthci(1:end)'], ...
                                'Color', lc1e, 'LineWidth', elw)
plot(1:length(rth), rth(1:end), 'o', 'MarkerSize', ms, ...
                'MarkerFaceColor', lc1, 'MarkerEdgeColor', 'none')
hold off
set(gca, 'LineWidth', alw, 'FontSize', afs, 'Box', 'off', ...
         'XLim', [-N_zz, length(rth)], 'XTick', [-25 0 40 80 120], ...
         'YLim', [-0.04 0.02], 'YTick', [-0.04 -0.02 0], 'YTickLabel', [])
     
     
set(axes, 'Units', 'Normalized', 'Position', [0.53 0.05 asize asize*pw/ph])
hold on
lx   = [-N_zz:1:-1 1:length(rth)];
ly   = rp;
lym  = rpm;
lyci = [rp'-rpci'; rp'+rpci'];
line([lx(1:end); lx(1:end)], lyci(:,1:end),'Color', lc2e, 'LineWidth', elw)
plot(lx(1:end), ly(1:end), '^', 'MarkerSize', ms, 'MarkerFaceColor', lc2, 'MarkerEdgeColor', 'none')
%plot(lx, lym, 'LineWidth', lw, 'Color', 'r')
hold off
set(gca, 'LineWidth', alw, 'FontSize', afs, 'Box', 'off', ...
         'XLim', [-N_zz, length(rth)], 'XTick', [-25 0 40 80 120], ...
         'YLim', [-0.04 0.02], 'YTick', [-0.04 -0.02 0], 'YTickLabel', [-0.04 -0.02 0 0.02], ...
         'Color', 'none', 'YColor', 'r', 'YAxisLocation', 'right')
xlabel('Session (d)', 'FontSize', lfs)
ylabel('\Delta\lambda (session^{-1})', 'FontSize', lfs)






