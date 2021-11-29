% plot supplemental figure 3
%  plot average responses (spike rate) for MT, and ROC
%  areas during passive fixation
%%%          plot           %%%
fh = figure;
pw = 7;    % paper width
ph = 5;   % paper height
wysifig(fh, pw, ph) % set it to US letter size

% define figure parameters
alw = 1;    % axis line width
afs = 10;   % axis font size
lfs = 10;   % label font size
tfs = 8;   % text font size
ffs = 20;
lw  = 2;    % line width
ms  = 5;    % marker size
elw = 0.5;

asize = 0.18;
%% a: plot example fits
th = annotation('textbox', [0.1-0.1 0.96 0.05 0.05]);
set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'String', 'a', 'FontSize', ffs, 'FontWeight', 'bold') 

EgInd = [43 6 62 105]; % example fits (pre, post-early, post-late)

ah(1) = axes;
set(ah(1), 'Units', 'Normalized', 'Position', [0.09 0.69 asize asize*pw/ph])
hold on 
% get data
Monk = 'Cy';
[xx1, xx2, xx3, fits] = getML_MTCTdependence('Cy', 0);
i      = EgInd(1); % example fit
a      = getML_txt([Monk 'PRe_MT.txt']);
fn     = a.data{strcmp(a.name,'dat_fn')};
usable = a.data{strcmp(a.name,'usable')};

bb     = -200;
be     = 1500;
bs     = 100;
bw     = 200;
bins   = [[bb:bs:be-bw]' [bb+bw:bs:be]'];
mbins  = bins(:,2)/1000;

L      = mbins>=0.1; % assume neural delay is 150ms

[rM, rSD, rN, rMn] = getML_mRate([Monk 'PRe_MT.txt'], bins, [0 1], [], 0, 0);

m    = rMn(1:7,:,:)-rMn(8:14,:,:);            % select between normalized and unnormalized spikes
sd   = 1./(rN(1:7,:,:)+rN(8:14,:,:));         % weight by # of trials because I don't have sd for normalized spike rate
coh  = [0 0.032 0.064 0.128 0.256 0.512 0.999]';
Lcoh = logical([1 1 1 1 1 1 1]);
          

% plot
lc  = [];
for j = 1:7
    lc = [lc; [0.9 0.9 0.9]-(j-1)/7];
end

             
hold on
for j = 1:7
    plot(mbins, m(j,:,i), 'o', 'MarkerSize', ms+1, ...
        'MarkerFaceColor', lc(j,:), 'MarkerEdgeColor', 'none')
    mbP = mbins(L);
    t   = [mbP(1):((mbP(end)-mbP(1))/100):mbP(end)]';
    c   = repmat(coh(j),length(t),1);
    A   = [ones(length(t),1), c, t, c.*t];

    plot(t, A*fits(:,i), '-', 'LineWidth', lw, 'Color', lc(j,:))
end
hold off
set(gca, 'xlim', [0 1], 'XTick', [0 0.5 1], 'ylim', [-0.2 1.5], 'YTick', [0 0.5 1], ...
          'FontSize', afs, 'LineWidth', alw)
xlabel('Time (s)', 'FontSize', lfs)
ylabel([sprintf('Normalized rate\n') '(pref-null)'], 'FontSize', lfs)
title('Pre-training', 'FontSize', tfs+1)

   

i=1; 
ah(i+1) = axes;
set(ah(i+1), 'Units', 'Normalized', 'Position', [0.09+i*(asize+0.035) 0.69 asize asize*pw/ph])
hold on
% get data
Monk = 'Cy';
[fits] = getML_MTCTdependence('Cy', 0);
i      = EgInd(2); % example fit
a      = getML_txt([Monk 'TRain_MT.txt']);
fn     = a.data{strcmp(a.name,'dat_fn')};
ses    = a.data{strcmp(a.name,'session')};
usable = a.data{strcmp(a.name,'usable')};
a             = getML_txt('CyTRain_MT.txt');
mses_cy       = a.data{strcmp(a.name,'session')};


bb     = -200;
be     = 1500;
bs     = 100;
bw     = 200;
bins   = [[bb:bs:be-bw]' [bb+bw:bs:be]'];
mbins  = bins(:,2)/1000;

[rM, rSD, rN, rMn] = getML_mRate([Monk 'TRain_MT.txt'], bins, [0 1], [], 0, 0);

m    = rMn(1:7,:,:)-rMn(8:14,:,:);            % select between normalized and unnormalized spikes
sd   = 1./(rN(1:7,:,:)+rN(8:14,:,:));         % weight by # of trials because I don't have sd for normalized spike rate
coh  = [0 0.032 0.064 0.128 0.256 0.512 0.999]';
Lcoh = logical([1 1 1 1 1 1 1]);
               

% plot
lc  = [];
for j = 1:7
    lc = [lc; [0.9 0.9 0.9]-(j-1)/7];
end

hold on
for j = 1:7
    plot(mbins, m(j,:,i), 'o', 'MarkerSize', ms+1, ...
        'MarkerFaceColor', lc(j,:), 'MarkerEdgeColor', 'none')
    mbP = mbins(L);
    t   = [mbP(1):((mbP(end)-mbP(1))/100):mbP(end)]';
    c   = repmat(coh(j),length(t),1);
    A   = [ones(length(t),1), c, t, c.*t];

    plot(t, A*fits(:,i), '-', 'LineWidth', lw, 'Color', lc(j,:))
end
hold off
set(gca, 'xlim', [0 1], 'XTick', [0 0.5 1], 'ylim', [-0.2 1.5], 'YTick', [0 0.5 1], 'YTickLabel', [], ...
          'FontSize', afs, 'LineWidth', alw)
text(1, -0.198,  ['Session ' int2str(mses_cy(i))], 'FontSize', tfs-1, ...
     'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom')
title('Early', 'FontSize', tfs+1)



   


i=2; 
ah(i+1) = axes;
set(ah(i+1), 'Units', 'Normalized', 'Position', [0.09+i*(asize+0.035) 0.69 asize asize*pw/ph])
hold on
% get data
Monk = 'Cy';
[fits] = getML_MTCTdependence('Cy', 0);
i      = EgInd(3); % example fit
a      = getML_txt([Monk 'TRain_MT.txt']);
fn     = a.data{strcmp(a.name,'dat_fn')};
ses    = a.data{strcmp(a.name,'session')};
usable = a.data{strcmp(a.name,'usable')};

bb     = -200;
be     = 1500;
bs     = 100;
bw     = 200;
bins   = [[bb:bs:be-bw]' [bb+bw:bs:be]'];
mbins  = bins(:,2)/1000;

[rM, rSD, rN, rMn] = getML_mRate([Monk 'TRain_MT.txt'], bins, [0 1], [], 0, 0);

m    = rMn(1:7,:,:)-rMn(8:14,:,:);            % select between normalized and unnormalized spikes
sd   = 1./(rN(1:7,:,:)+rN(8:14,:,:));         % weight by # of trials because I don't have sd for normalized spike rate
coh  = [0 0.032 0.064 0.128 0.256 0.512 0.999]';
Lcoh = logical([1 1 1 1 1 1 1]);
               

% plot
lc  = [];
for j = 1:7
    lc = [lc; [0.9 0.9 0.9]-(j-1)/7];
end

hold on
for j = 1:7
    plot(mbins, m(j,:,i), 'o', 'MarkerSize', ms+1, ...
        'MarkerFaceColor', lc(j,:), 'MarkerEdgeColor', 'none')
    mbP = mbins(L);
    t   = [mbP(1):((mbP(end)-mbP(1))/100):mbP(end)]';
    c   = repmat(coh(j),length(t),1);
    A   = [ones(length(t),1), c, t, c.*t];

    plot(t, A*fits(:,i), '-', 'LineWidth', lw, 'Color', lc(j,:))
end
hold off
set(gca, 'xlim', [0 1], 'XTick', [0 0.5 1], 'ylim', [-0.2 1.5], 'YTick', [0 0.5 1], 'YTickLabel', [], ...
          'FontSize', afs, 'LineWidth', alw)
title('Mid', 'FontSize', tfs+1)
text(1, -0.198,  ['Session ' int2str(mses_cy(i))], 'FontSize', tfs-1, ...
     'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom')

 
i=3; 
ah(i+1) = axes;
set(ah(i+1), 'Units', 'Normalized', 'Position', [0.09+i*(asize+0.035) 0.69 asize asize*pw/ph])
hold on
% get data
Monk = 'Cy';
[fits] = getML_MTCTdependence('Cy', 0);
i      = EgInd(4); % example fit
a      = getML_txt([Monk 'TRain_MT.txt']);
fn     = a.data{strcmp(a.name,'dat_fn')};
ses    = a.data{strcmp(a.name,'session')};
usable = a.data{strcmp(a.name,'usable')};

bb     = -200;
be     = 1500;
bs     = 100;
bw     = 200;
bins   = [[bb:bs:be-bw]' [bb+bw:bs:be]'];
mbins  = bins(:,2)/1000;

[rM, rSD, rN, rMn] = getML_mRate([Monk 'TRain_MT.txt'], bins, [0 1], [], 0, 0);

m    = rMn(1:7,:,:)-rMn(8:14,:,:);            % select between normalized and unnormalized spikes
sd   = 1./(rN(1:7,:,:)+rN(8:14,:,:));         % weight by # of trials because I don't have sd for normalized spike rate
coh  = [0 0.032 0.064 0.128 0.256 0.512 0.999]';
Lcoh = logical([1 1 1 1 1 1 1]);
               

% plot
lc  = [];
for j = 1:7
    lc = [lc; [0.9 0.9 0.9]-(j-1)/7];
end

hold on
for j = 1:7
    plot(mbins, m(j,:,i), 'o', 'MarkerSize', ms+1, ...
        'MarkerFaceColor', lc(j,:), 'MarkerEdgeColor', 'none')
    mbP = mbins(L);
    t   = [mbP(1):((mbP(end)-mbP(1))/100):mbP(end)]';
    c   = repmat(coh(j),length(t),1);
    A   = [ones(length(t),1), c, t, c.*t];

    plot(t, A*fits(:,i), '-', 'LineWidth', lw, 'Color', lc(j,:))
end
hold off
set(gca, 'xlim', [0 1], 'XTick', [0 0.5 1], 'ylim', [-0.2 1.5], 'YTick', [0 0.5 1.5], 'YTickLabel', [], ...
          'FontSize', afs, 'LineWidth', alw)
title('Late', 'FontSize', tfs+1)
text(1, -0.198,  ['Session ' int2str(mses_cy(i))], 'FontSize', tfs-1, ...
     'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom')
 
 
 mses_cy(EgInd)
   




%% b: ROC areas during passive fixation task for both monkeys
asize = 0.22;
th = annotation('textbox', [0 0.57 0.05 0.05]);
set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'String', 'b', 'FontSize', ffs, 'FontWeight', 'bold') 
     
fn   = 'CyTRain_MT.txt';
[w_c, r, roc_c, wi_c, roci_c] = getML_tuningProperties(fn, 0);
a             = getML_txt(fn);
fn            = a.data{strcmp(a.name,'dat_fn')};
ses_c         = a.data{find(strcmp(a.name, 'session'))};
use_c         = a.data{find(strcmp(a.name, 'usable'))}; 

fn   = 'ZZTRain_MT.txt';
[w_z, r, roc_z, wi_z, roci_z] = getML_tuningProperties(fn, 0);
a             = getML_txt(fn);
ses_z         = a.data{find(strcmp(a.name, 'session'))};
use_z         = a.data{find(strcmp(a.name, 'usable'))}; 



set(axes, 'Units', 'Normalized', 'Position', [0.1 0.37 0.4 0.2])
hold on
line(repmat(ses_c,1,2)', [roc_c-roci_c roc_c+roci_c]', 'Color', [0.5 0.5 0.5], 'LineWidth', elw)
plot(ses_c, roc_c, 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerSize', ms)
Lgd = ~isnan(roc_c);
[b, bint, rr, rri, stats] = regress(roc_c(Lgd),[ones(size(roc_c(Lgd))) ses_c(Lgd)]);
stats(3)
plot([1:160], b(1)+b(2)*[1:160], '--k', 'LineWidth', lw);
hold off
text(160, 0.5, 'Monkey C', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'FontSize', lfs)
set(gca, 'XLim', [1 160], 'XTick', [1 50 100 150], 'YLim', [0.5 1], 'YTick', [0.5 0.6 0.7 0.8 0.9 1], 'FontSize', afs)
xlabel('Session (d)', 'FontSize', lfs)
ylabel(sprintf('ROC area'), 'FontSize', lfs)
          
set(axes, 'Units', 'Normalized', 'Position', [0.55 0.37 0.4 0.2])
hold on
line(repmat(ses_z,1,2)', [roc_z-roci_z roc_z+roci_z]', 'Color', [0.5 0.5 0.5], 'LineWidth', elw)
plot(ses_z, roc_z, 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerSize', ms)
Lgd = ~isnan(roc_z);
[b, bint, rr, rri, stats] = regress(roc_z(Lgd),[ones(size(roc_z(Lgd))) ses_z(Lgd)]);
stats(3)
plot([1:130], b(1)+b(2)*[1:130], '--k', 'LineWidth', lw);
hold off     
text(130, 0.5, 'Monkey Z', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'FontSize', lfs)
set(gca, 'XLim', [1 130], 'XTick', [1 40 80 120], 'YLim', [0.5 1], 'YTick', [0.5 0.6 0.7 0.8 0.9 1], 'FontSize', afs)
xlabel('Session (d)', 'FontSize', lfs)
%ylabel(sprintf('ROC area'), 'FontSize', lfs)

% th = annotation('textbox', [0.28 0.52 0.5 0.05]);
% set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
%     'String', 'ROC area of pref vs null responses during passive fixation', 'FontSize', lfs) 


     
%% c: tuning width of cell during passive fixation task for both monkeys
th = annotation('textbox', [0 0.28 0.05 0.05]);
set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'String', 'c', 'FontSize', ffs, 'FontWeight', 'bold') 
     
set(axes, 'Units', 'Normalized', 'Position', [0.1 0.08 0.4 0.2])
hold on
line(repmat(ses_c,1,2)', wi_c', 'Color', [0.5 0.5 0.5], 'LineWidth', elw)
plot(ses_c, w_c(:,2), 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerSize', ms)
Lgd = ~isnan(w_c(:,2));
[b, bint, rr, rri, stats] = regress(w_c(Lgd,2),[ones(size(w_c(Lgd,2))) ses_c(Lgd)]);
stats(3)
plot([1:160], b(1)+b(2)*[1:160], '--k', 'LineWidth', lw);
hold off
text(160, 0, 'Monkey C', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'FontSize', lfs)
set(gca, 'XLim', [1 160], 'XTick', [1 50 100 150], 'YLim', [0 60], 'YTick', [0 20 40 60 80], 'FontSize', afs)
xlabel('Session (d)', 'FontSize', lfs)
ylabel(sprintf('Tuning width (degree)'), 'FontSize', lfs)
          
set(axes, 'Units', 'Normalized', 'Position', [0.55 0.08 0.4 0.2])
hold on
line(repmat(ses_z,1,2)', wi_z', 'Color', [0.5 0.5 0.5], 'LineWidth', elw)
plot(ses_z, w_z(:,2), 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerSize', ms)
Lgd = ~isnan(w_z(:,2));
[b, bint, rr, rri, stats] = regress(w_z(Lgd,2),[ones(size(w_z(Lgd,2))) ses_z(Lgd)]);
stats(3)
plot([1:130], b(1)+b(2)*[1:130], '--k', 'LineWidth', lw);
hold off     
text(130, 0, 'Monkey Z', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'FontSize', lfs)
set(gca, 'XLim', [1 130], 'XTick', [1 40 80 120], 'YLim', [0 60], 'YTick', [0 20 40 60 80], 'FontSize', afs)
xlabel('Session (d)', 'FontSize', lfs)
%ylabel(sprintf('ROC area'), 'FontSize', lfs)

% th = annotation('textbox', [0.28 0.52 0.5 0.05]);
% set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
%     'String', 'Tuning width of neural responses during passive fixation', 'FontSize', lfs) 
    
     
     
     

