%% plotML_supp8

%% get performance
tlim = [-inf inf];
[p_c n_c pci_c] = getML_performancePerCoh('CyTRain_psy.txt', tlim, 0);
[p_z n_z pci_z] = getML_performancePerCoh('ZZTRain_psy.txt', tlim, 0);

[pa_c na_c] = getML_performancePerCohAdjusted('CyTRain_psy.txt', tlim, 0);
[pa_z na_z] = getML_performancePerCohAdjusted('ZZTRain_psy.txt', tlim, 0);

%% plot 
% create figure
fh = figure; 
mm2inch = 25.4;
pw = 8;    % paper width, make it nature size of one column, 85mm x 85mm
ph = 10;   % paper height
wysifig(fh, pw, ph) % set it to US letter size

% define figure parameters
alw = 0.5;  % axis line width
afs = 10;    % axis font size
lfs = 10;    % label font size
ffs = 12;   % figure font size
tfs = 10;   % text font size
lw  = 3; % line width
elw = 0.5;
ms  = 4;    % marker size

% define line color
lc1  = [0.3 0.3 0.3];
lc1e = [0.75 0.75 0.75];
lc2  = [1 0.3 0.3];
lc2e = [1 0.75 0.75];

lc  = [];
for i = 1:7
    lc = [lc; [0.9 0.9 0.9]-(i-1)/7];
end


asize = 0.35;
wsize = 20;

%% performance per session (cy)
set(axes, 'Units', 'Normalized', 'Position', [0.1 0.65 asize asize*pw/ph])
hold on
for i = 1:7
   plot(1:length(p_c), p_c(i,:), 'o', ...
        'MarkerSize', ms, 'MarkerFaceColor', lc(i,:), 'MarkerEdgeColor', 'none')
   plot(1:length(p_c), nanrunmean(p_c(i,:),wsize), '-', 'Color', lc(i,:), 'LineWidth', lw) 
end
hold off
set(gca, 'LineWidth', alw, 'FontSize', afs, 'Box', 'off', ...
         'XLim', [1, length(p_c)], 'YLim', [0.4 1], ...
         'XTick', [1 50 100 150], 'YTick', [0.4 0.6 0.8 1])
xlabel('Session (d)', 'FontSize', lfs)
ylabel('Performance (% crt)', 'FontSize', lfs)


%% performance per session (zz)
set(axes, 'Units', 'Normalized', 'Position', [0.6 0.65 asize asize*pw/ph])
hold on
for i = 1:7
   plot(1:length(p_z), p_z(i,:), 'o', ...
        'MarkerSize', ms, 'MarkerFaceColor', lc(i,:), 'MarkerEdgeColor', 'none')
   plot(1:length(p_z), nanrunmean(p_z(i,:),wsize), '-', 'Color', lc(i,:), 'LineWidth', lw)
end
hold off
set(gca, 'LineWidth', alw, 'FontSize', afs, 'Box', 'off', ...
         'XLim', [1, length(p_z)], 'YLim', [0.4 1], ...
         'XTick', [1 50 100 150], 'YTick', [0.4 0.6 0.8 1])
xlabel('Session (d)', 'FontSize', lfs)
ylabel('Performance (% crt)', 'FontSize', lfs)



%% performance per session, adjusted for lapse (cy)
set(axes, 'Units', 'Normalized', 'Position', [0.1 0.3 asize asize*pw/ph])
hold on
for i = 1:7
   plot(1:length(p_c), pa_c(i,:), 'o', ...
        'MarkerSize', ms, 'MarkerFaceColor', lc(i,:), 'MarkerEdgeColor', 'none')
   plot(1:length(p_c), nanrunmean(pa_c(i,:),wsize), '-', 'Color', lc(i,:), 'LineWidth', lw) 
end
hold off
set(gca, 'LineWidth', alw, 'FontSize', afs, 'Box', 'off', ...
         'XLim', [1, length(p_c)], 'YLim', [0.4 1], ...
         'XTick', [1 50 100 150], 'YTick', [0.4 0.6 0.8 1])
xlabel('Session (d)', 'FontSize', lfs)
ylabel('Performance (% crt)', 'FontSize', lfs)



%% performance per session, adjusted for lapse (zz)
set(axes, 'Units', 'Normalized', 'Position', [0.6 0.3 asize asize*pw/ph])
hold on
for i = 1:7
   plot(1:length(p_z), pa_z(i,:), 'o', ...
        'MarkerSize', ms, 'MarkerFaceColor', lc(i,:), 'MarkerEdgeColor', 'none')
   plot(1:length(p_z), nanrunmean(pa_z(i,:),wsize), '-', 'Color', lc(i,:), 'LineWidth', lw)
end
hold off
set(gca, 'LineWidth', alw, 'FontSize', afs, 'Box', 'off', ...
         'XLim', [1, length(p_z )], 'YLim', [0.4 1], ...
         'XTick', [1 50 100 150], 'YTick', [0.4 0.6 0.8 1])
xlabel('Session (d)', 'FontSize', lfs)
ylabel('Performance (% crt)', 'FontSize', lfs)



