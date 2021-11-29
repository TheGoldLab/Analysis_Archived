function plotML_rateByDirCoh(time, rate, ah, lw)
% USAGE:
% time = mean bins
% rate = 14xm matrix, where m = number of bins
%        rows 1-7  are responses for 0-99.9% pref coh,
%             8-14 are responses for 0-99.9% null coh
% ah   = figure handle at while the figure is drawn
% lw   = linewidth

if nargin<2
    return;
end

if nargin < 3
    figure
    ah = axes;
end

if nargin < 4
    lw = 5;
end


% set current axes
set(gcf, 'CurrentAxes', ah)


if size(rate,1)==14
    % prepare plot parameters
    lc  = [];
    for i = 1:7
        lc = [lc; [0.9 0.9 0.9]-(i-1)/7];
    end
    lh  = [line; line; line; line; line; line; line;...
        line; line; line; line; line; line; line];
    set(lh(1:7),  'LineStyle', '-',  'LineWidth', lw)
    set(lh(8:14), 'LineStyle', '--', 'LineWidth', lw)
    for j = 1:7
        set([lh(j) lh(j+7)], 'Color', lc(j,:));
    end


    % plot rate
    for j = 1:7
        set(lh(j),   'XData', time, 'YData', rate(j,:));
        set(lh(j+7), 'XData', time, 'YData', rate(j+7,:));
    end
    
elseif size(rate,1)==7
    % prepare plot parameters
    lc  = [];
    for i = 1:7
        lc = [lc; [0.9 0.9 0.9]-(i-1)/7];
    end
    lh  = [line; line; line; line; line; line; line];
    set(lh(1:7),  'LineStyle', '-',  'LineWidth', lw)
    for j = 1:7
        set([lh(j)], 'Color', lc(j,:));
    end

   % plot rate
    for j = 1:7
        set(lh(j),   'XData', time, 'YData', rate(j,:));
    end
    
end
    
% set axis properties
set(gca, 'Box', 'on',...
          'FontSize', 12)
xlabel('Time (msec)', 'FontSize', 12)
ylabel('Response (spikes/sec)', 'FontSize', 12)

    