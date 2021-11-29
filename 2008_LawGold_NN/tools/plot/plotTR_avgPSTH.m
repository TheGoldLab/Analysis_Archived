function Bses = plotTR_avgPSTH(fn, fn_type, N, binsize, normalized, recompute, savefn)


% check input
if nargin<7    
    % default save path
    savefn = which('plotTR_avgPSTH.m');
    savefn = [savefn(1:end-2) '.mat'];
end

if nargin<6
    recompute = 1;
end

if nargin<5
    normalized = 1;
end

if nargin<4
    binsize = 100;
end

if nargin<3
    N = 3;
end

if nargin<2
    return;
end



% make session array, Bses
a    = getML_txt(fn);
ses  = a.data{strcmp(a.name,'session')};
Ms   = max(ses);
Bses = [1:floor(Ms/N):Ms]';
Bses = [Bses(1:N) [Bses(2:N)-1; Ms]];


if recompute
    % get average psth
    [ym yse yn] = getTR_avgPSTH(fn, fn_type, Bses, binsize, normalized, recompute, savefn);
    save(savefn, 'ym', 'yse', 'yn') 
else
    load(savefn)
end



%% plot average response
% figure parameters:
figsize  = [8 8];
afs      = 35;      % axis font size
tfs      = 25;      % text font size
alw      = 2;       % axis line width
psthlw   = 5;       % psth line width

% plot
ylim_ = [];
fh    = [];

lc  = [];
for i = 1:7
    lc = [lc; [0.9 0.9 0.9]-(i-1)/7];
end 

binsize = 100;
os0     = -200;
maxdd   = 1500;
bins    = [os0:binsize:maxdd-binsize]+binsize/2;

for i = 1:size(Bses,1)
    fh(i) = figure;
    set(fh(i), 'Units', 'inches', 'Position', [0 0 figsize])

    % set default object properties
    lh  = [line; line; line; line; line; line; line;...
           line; line; line; line; line; line; line];
    set(lh(1:7),  'LineStyle', '-',  'LineWidth', psthlw)
    set(lh(8:14), 'LineStyle', '--', 'LineWidth', psthlw)
    for j = 1:7
        set([lh(j) lh(j+7)], 'Color', lc(j,:));
    end
 
    % plot PSTH
    for j = 1:7
        set(lh(j),   'XData', bins, 'YData', ym(:,j,i));
        set(lh(j+7), 'XData', bins, 'YData', ym(:,j+7,i));   
    end
    
    % set axis properties
    set(gca, 'PlotBoxAspectRatioMode', 'manual', 'PlotBoxAspectRatio', [1 1 1],...
             'XLim', [os0+binsize/2 maxdd-1.5*binsize])
    ylim_(i,:) = get(gca, 'YLim');
    xlabel('Time (msec)', 'FontSize', afs)
    ylabel('Response (spikes/sec)', 'FontSize', afs)
    
    % set figure properties
    set(gca, 'LineWidth', alw, 'FontSize', afs, 'Box', 'on')
end

yl = [min(ylim_(:,1)) max(ylim_(:,2))];
for i = 1:size(Bses,1)
    figure(fh(i));
    set(gca, 'YLim', yl);    
    line([0 0], yl, 'LineStyle', '-', 'LineWidth', alw+1, 'Color', 'k');
    text(20, yl(2)-0.1*yl(2),...
          sprintf('Session %d-%d, n=%d', Bses(i,1), Bses(i,2), yn(i)),...
         'FontSize', tfs);
end

%legend(lh([7:-1:1]), '99.9%', '51.2%', '25.6%', '12.8%', '6.4%', '3.2%', '0%');

    