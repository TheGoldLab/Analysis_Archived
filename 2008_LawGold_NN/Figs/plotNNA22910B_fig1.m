% plotML_fig1.m
% 

% created by jcl on 1/16/07

% Panel 1: task (not plotted in this script)
% Panel 2: plot example of behavioral performance in one session, and fits
% Panel 3: plot Cy's threshold against session
% Panel 4: plot ZZ's threshold against session

% create figure
fh = figure; 
mm2inch = 25.4;
pw = 85/mm2inch;   % paper width, make it nature size of one column, 85mm x 85mm
ph = 126/mm2inch;   % paper height
wysifig(fh, pw, ph) % set it to US letter size

% define figure parameters
alw = 0.5;  % axis line width
afs = 6;    % axis font size
lfs = 6;    % label font size
ffs = 12;   % figure font size
tfs = 6;    % text font size
lw  = 1.5; % line width
elw = 0.5;
ms  = 2;    % marker size

% define line color
lc1  = [0.3 0.3 0.3];
lc1e = [0.75 0.75 0.75];
lc2  = [1 0.3 0.3];
lc2e = [1 0.75 0.75];


asize = 0.35;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fig 1a: task
th = annotation('textbox', [0.1-0.1 0.58+asize+0.018 0.05 0.05]);
set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'String', 'a', 'FontSize', ffs, 'FontWeight', 'bold') 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fig 1c and 1d: plot threshold and lapse vs session
th = annotation('textbox', [0.1-0.1 0.1+asize+0.013+0.16 0.05 0.05]);
set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'String', 'c', 'FontSize', ffs, 'FontWeight', 'bold') 


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
set(axes, 'Units', 'Normalized', 'Position', [0.105 0.1+0.28 asize asize*pw/ph])
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
xlabel('Session (d)', 'FontSize', lfs)
ylabel('Threshold (% coh)', 'FontSize', lfs)
text(0.98*size(fits_cy,2), 98, 'Monkey C', 'FontSize', tfs, ...
     'HorizontalAlignment', 'right', 'VerticalAlignment', 'top')


set(axes, 'Units', 'Normalized', 'Position', [0.52 0.1+0.28 asize asize*pw/ph])
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
xlabel('Session (d)', 'FontSize', lfs)
text(0.98*size(fits_zz,2), 98, 'Monkey Z', 'FontSize', tfs, ...
     'HorizontalAlignment', 'right', 'VerticalAlignment', 'top')

 
% % plot a patch to indicate the pre-training period (99.9% only)
% % cy
% set(axes, 'Units', 'Normalized', 'Position', [0.105 0.1 asize asize])
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
set(axes, 'Units', 'Normalized', 'Position', [0.105 0.1+0.28 asize asize*pw/ph])
hold on
lx   = [-N_cy:1:-1 1:size(pt_cy,2)];
ly   = 1-[p_cy pt_cy];
lyci = 1-[pci_cy ptci_cy];
lys  = mean(abs([lyci(1,:)-ly; lyci(2,:)-ly]));
lys(lys==0) = 0.001;   % make it a small number
L    = ~isnan(lx) & ~isnan(ly) & ~isnan(lys);
[blc, bsemlc, gof, lym] = exp_fit2W(lx(L), ly(L), lys(L), [0.5; 5], [0.4 0.6;0.01 100]);
line([lx; lx], lyci,'Color', lc2e, 'LineWidth', elw)
plot(lx, ly, 'o', 'MarkerSize', ms, 'MarkerFaceColor', lc2, 'MarkerEdgeColor', 'none')
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
set(axes, 'Units', 'Normalized', 'Position', [0.52 0.1+0.28 asize asize*pw/ph])
hold on
lx   = [-N_zz:1:-1 1:size(pt_zz,2)];
ly   = 1-[p_zz pt_zz];
lyci = 1-[pci_zz ptci_zz];
lys  = mean(abs([lyci(1,:)-ly; lyci(2,:)-ly]));
lys(lys==0) = 0.001;   % make it a small number
L    = ~isnan(lx) & ~isnan(ly) & ~isnan(lys);
[blz, bsemlz, gof, lym] = exp_fit2W(lx(L), ly(L), lys(L), [0.5; 5], [0.4 0.6;0.01 100]);
line([lx; lx], lyci,'Color', lc2e, 'LineWidth', elw)
plot(lx, ly, 'o', 'MarkerSize', ms, 'MarkerFaceColor', lc2, 'MarkerEdgeColor', 'none')
plot(lx(L), lym, 'LineWidth', lw, 'Color', 'r')
%plot([-N_zz:1:130], 0.15*ones(size([-N_zz:1:130])), 'LineWidth', 0.5, 'LineStyle', ':', 'Color', 'r')
hold off
set(gca, 'LineWidth', alw, 'FontSize', afs, 'Box', 'off', ...
         'XLim', [-N_zz, size(fits_zz,2)], 'XTick', [-N_zz 1 40 80 120], ...
         'YLim', [0 1], 'YTick', [0 0.25 0.5 0.75 1], 'YColor', 'r', 'YAxisLocation', 'right', ...
         'Color', 'none')
ylabel('Error rate at 99% coh', 'FontSize', lfs)





% fit a line to lapse rate after 1.61 time constant of the session to see if
% it changed significantly
tau = round(1.61*blz(2));
L   = ~isnan(lx) & ~isnan(ly) & ~isnan(lys) & logical([zeros(1,tau), ones(1,length(ly)-tau)]);
[b, bint, h, p] = nestedFW(ly(L)', lys(L)', [ones(size(ly(L)')) lx(L)']);


% %%% old code, plot lapse from the fit %%%
% % plot lapse
% set(axes, 'Units', 'Normalized', 'Position', [0.105 0.1 asize asize])
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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Fig 1b: plot behavioral performance of an example session, and the fit
th = annotation('textbox', [0.54-0.1 0.58+asize+0.018 0.05 0.05]);
set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'String', 'b', 'FontSize', ffs, 'FontWeight', 'bold') 

% get threshold for each bin
[fits_cy, sems_cy, th_cy] = getML_psyPerformanceCT('CyTRain_psy.txt', 0);
a   = getML_txt('CyTRain_psy.txt');
fn  = a.data{strcmp(a.name,'dat_fn')};
ses = a.data{strcmp(a.name,'session')}; % choose session 142 or 138

ind   = find(ses==142);

bb   = 0;
be   = 1000;
bw   = 150;
bs   = 75;
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

% fits  = nans(3,size(bins,1));
% sems  = nans(3,size(bins,1));
% 
% for i = 1:size(bins,1)
%     [fits(:,i) sems(:,i)] = quick_fit([coh p(:,i) N(:,i)]);
% end

recompute = 0;
if recompute
    fits  = nans(2,size(bins,1));
    sems  = nans(2,2,size(bins,1));

    for i = 1:size(bins,1)
        i
        Lt  = FIRA.ecodes.data(:,getFIRA_ecodeColumnByName('dot_dur'))>=bins(i,1) ...
            & FIRA.ecodes.data(:,getFIRA_ecodeColumnByName('dot_dur'))<bins(i,2);
        Lcr = FIRA.ecodes.data(:,getFIRA_ecodeColumnByName('correct'))>=0;

        [fits(:,i), sems(:,:,i)] = ctPsych_fit(@quick2, FIRA.ecodes.data(Lt&Lcr,getFIRA_ecodeColumnByName('dot_coh'))/100,...
            FIRA.ecodes.data(Lt&Lcr,getFIRA_ecodeColumnByName('correct')), [], {100,68,81.6}, []);
    end
    save /work/mirror_jeff/code/matlab/tmp-mat/plotML_fig1b.mat fits sems

else

    load /work/mirror_jeff/code/matlab/tmp-mat/plotML_fig1b.mat
end



% plot thresholds against time
set(axes, 'Units', 'Normalized', 'Position', [0.52 0.71 asize asize*pw/ph])
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
plot(mbin/1000, 100*fits(1,:), 'o', 'MarkerSize', ms, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');
plot([0.01:0.01:et], 100*alpha, 'LineWidth', lw, 'Color', 'k');
plot(ones(size([0:0.1:100*alpha(end)])), [0:0.1:100*alpha(end)], 'LineWidth', lw/2, ...
        'Color', 'k', 'LineStyle', '--');
hold off
set(gca, 'LineWidth', alw, 'FontSize', afs, 'Box', 'off', ...
         'XLim', [0 et], 'YLim', [7.5 40], 'YTick', [10 20 40], 'YScale', 'log')   
% plot annotation arrow
[xb,yb] = axisXY2figXYslog(et,100*alpha(end));
[xe,ye] = axisXY2figXYslog(0,100*alpha(end));
anh     = annotation('arrow', [xb xe], [yb ye]);
set(anh, 'Color', 'k', 'LineWidth', lw/2, 'LineStyle', '--', 'HeadStyle', 'cback2', ...
          'HeadWidth', 3*ms, 'HeadLength', 3*ms)
xlabel('Time (s)', 'FontSize', lfs)
ylabel('Threshold (% coh)', 'FontSize', lfs)
text(1.18, 7.55, sprintf('%s', [fn{ind}(1:10) '\_' fn{ind}(12:end-4)]), 'FontSize', tfs-1,...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom')



% figure
% subplot(1,2,1)
% L = fits_cy(4,:)<=0.01;
% plot(fits_cy(4,~L),'ko')
% hold on
% plot(fits_cy(4,L),'ko', 'MarkerFaceColor', 'k')
% hold off
% 
% subplot(1,2,2)
% L = fits_zz(4,:)<=0.01;
% plot(fits_zz(4,~L),'ko')
% hold on
% plot(fits_zz(4,L),'ko', 'MarkerFaceColor', 'k')
% hold off
% 
% 
% 





% old code

% [fits_cy, sems_cy, th_cy] = getML_psyPerformanceCT('CyTRain_psy.txt', 0);
% a   = getML_txt('CyTRain_psy.txt');
% fn  = a.data{strcmp(a.name,'dat_fn')};
% ses = a.data{strcmp(a.name,'session')}; % choose session 142
% ind   = find(ses==142); %142
% 
% bb   = 0;
% be   = 1000;
% bw   = 250;
% bs   = 50;
% bins = [[bb:bs:be-bw]' [bb+bw:bs:be]'];
% mbin = mean(bins,2);
% coh  = [0 3.2 6.4 12.8 25.6 51.2 99.9]';
%  
% % get performance for each coh and time
% global FIRA
% openFIRA(fn{ind})
% fprintf('%s\n', fn{ind})
% p = nans(length(coh),size(bins,1));
% N = nans(length(coh),size(bins,1));
% for i = 1:size(bins,1)
%     Lt = FIRA.ecodes.data(:,getFIRA_ecodeColumnByName('dot_dur'))>=bins(i,1) ...
%         & FIRA.ecodes.data(:,getFIRA_ecodeColumnByName('dot_dur'))<bins(i,2);         
%     for j = 1:length(coh)
%         Lcoh = FIRA.ecodes.data(:,getFIRA_ecodeColumnByName('dot_coh'))==coh(j);
%         N(j,i) = sum(FIRA.ecodes.data(Lcoh&Lt,getFIRA_ecodeColumnByName('correct'))>=0);
%         p(j,i) = sum(FIRA.ecodes.data(Lcoh&Lt,getFIRA_ecodeColumnByName('correct'))==1)/N(j,i);
%     end
% end
% 
% 
% % plot threshold with time
% set(axes, 'Units', 'Normalized', 'Position', [0.55 0.55 asize asize])
% lc  = [];
% for j = 1:7
%     lc = [lc; [0.9 0.9 0.9]-(j-1)/7];
% end
% 
% % plot raw data and fitted data
% hold on
% for i = 1:7
%     alpha = fits_cy(1,ind).*(mbin/1000).^fits_cy(2,ind);
%     val = 0.5 + (0.5 - fits_cy(4,ind)).*(1-exp(-((coh(i)/100)./alpha).^fits_cy(3,ind)));
%     plot(mbin/1000, p(i,:), 'o', 'MarkerSize', ms, 'MarkerFaceColor', lc(i,:), 'MarkerEdgeColor', 'none');
%     plot(mbin/1000, val, 'LineWidth', lw, 'Color', lc(i,:));
% end
% hold off
% 
% set(gca, 'LineWidth', alw, 'FontSize', afs, 'Box', 'off', ...
%          'XLim', [0 0.8], 'YLim', [0.35 1])
% xlabel('Time (second)', 'FontSize', lfs)
% ylabel('Proportion correct', 'FontSize', lfs)
% text(0.78, 0.4, sprintf('%s', fn{ind}(1:end-4)), 'FontSize', tfs,...
%     'HorizontalAlignment', 'right', 'VerticalAlignment', 'top')
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
th = annotation('textbox', [0.1-0.1 0.04+asize*pw/ph 0.05 0.05]);
set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'String', 'd', 'FontSize', ffs, 'FontWeight', 'bold') 

% get data
[pp, ppp, bc, bci, pvc] = getML_performancePerCoh('Cy', 0);
[pp, ppp, bz, bzi, pvz] = getML_performancePerCoh('ZZ', 0);

% cy
set(axes, 'Units', 'Normalized', 'Position', [0.105 0.07 asize asize*pw/ph*0.9])
c = [1 3 6 13 26 51 99];
hold on
plot([0:100], zeros(1,101), 'k:', 'LineWidth', elw)
plot([c; c],[bc-bci; bc+bci], 'k-', 'LineWidth', elw, 'Color', 'k')
plot(c, bc, 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerSize', ms+2)
plot(c(pvc<0.05), bc(pvc<0.05)+bci(pvc<0.05)+5e-4, '*', 'MarkerEdgeColor', 'k', 'MarkerSize', ms+1)
hold off
set(gca, 'LineWidth', alw, 'FontSize', afs, 'Box', 'off', ...
         'ylim', [-1e-3 4e-3], 'YTick', [0 2e-3 4e-3], ...
         'xscale', 'log', 'xlim', [0.7 99], 'XTick', c, 'XTickLabel', [0 c(2:end)])   
xlabel('% coh', 'FontSize', lfs)
ylabel('% correct/session', 'FontSize', lfs)
text(1, 4e-3, 'Monkey C', 'FontSize', tfs, ...
     'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')


% zz
set(axes, 'Units', 'Normalized', 'Position', [0.52 0.07 asize asize*pw/ph*0.9])
c = [1 3 6 13 26 51 99];
hold on
plot([0:100], zeros(1,101), 'k:', 'LineWidth', elw)
plot([c; c],[bz-bzi; bz+bzi], 'k-', 'LineWidth', elw, 'Color', 'k')
plot(c, bz, 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerSize', ms+2)
plot(c(pvz<0.05), bz(pvz<0.05)+bzi(pvz<0.05)+5e-4, '*', 'MarkerEdgeColor', 'k', 'MarkerSize', ms+1)
hold off
set(gca, 'LineWidth', alw, 'FontSize', afs, 'Box', 'off', ...
         'ylim', [-1e-3 4e-3], 'YTick', [0 2e-3 4e-3], 'YTickLabel', [],...
         'xscale', 'log', 'xlim', [0.7 99], 'XTick', c, 'XTickLabel', [0 c(2:end)])   
xlabel('% coh', 'FontSize', lfs)
text(1, 4e-3, 'Monkey Z', 'FontSize', tfs, ...
     'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')

 
 



