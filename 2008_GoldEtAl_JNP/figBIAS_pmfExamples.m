function fig_ = figBIAS_pmfExamples(num, ex_monk, ex_sss)
% function fig_ = figBIAS_pmfExamples(num, ex_monk, ex_sss)
%
% Plots of two example sessions, each showing:
%   1 % correct pmf vs coh
%   2 bias term(s)

if nargin < 1 || isempty(num)
    num = 8;
end

if nargin < 2 || isempty(ex_monk)
    %ex_monk = 'Ava';
    ex_monk = 'Atticus';
    %ex_monk = 'Cyrus';
    %ex_monk = 'ZZ';
end
monki = find(strcmp(ex_monk, getBIAS_monks));

if nargin < 3 || isempty(ex_sss)
    ex_sss = [20 160];
end

%%
% set up the figure
%%
% units should be in inches, from wysifig
wid  = 4.0; % total width
ht   = 1.2;
cols = {[.4 .6], [.4 .6]};
[axs,fig_] = getBIAS_axes(num, wid, ht, cols);

%%
% Get the data
%%
% data rows are trials, columns are:
%   1 coherence (0...1)
%   2 time (seconds)
%   3 dot dir (-1/0/1)
%   4 choice (-1/1)
%   5 correct (0/1)
%
pmfdat    = FS_loadProjectFile('2008_Bias', 'figBIAS_pmfCompare');
data      = FS_getDotsTrainingData(ex_monk);
sessions  = nonanunique(data(:,1));
os        = ones(1,3);
sp        = {'line_style', '-', 'colors', {0.*os}};
for ss = 1:2

    Lgood = data(:,1) == sessions(ex_sss(ss)) & data(:,2) <= 2 & data(:,3)>=0;
    
    %% PANEL 1: % rt choices PMF
    %
    axes(axs((ss-1)*2+1)); cla reset; hold on;
    Lg = Lgood & data(:,6) > 0.45;
    f = ctPsych_fit(@ddExp3, data(Lg, [5 6 8]), data(Lg, 3));
    ctPsych_plot(@ddExp3,f, data(Lg,[5 6 8 9]), data(Lg,3), 1, [], gca, sp);
    
    ylabel('% right choices')
    xlim([-100 100])

    %% PANEL 2: sample residuals
    %
    axes(axs((ss-1)*2+2)); 
    cla reset; hold on;
    xs = (1:sum(Lgood))';
    ys = pmfdat{monki,4}(Lgood,3,6);  
    plot([1 xs(end)], [0 0], 'k:');
    plot(xs, ys, 'k-');%, 'MarkerSize', 5);
%    plot(xs, nanrunmean(ys,50), 'r-');
    axis([1 1500 -2 2]);
end
