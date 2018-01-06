function fig_ = figBIAS_example(num, ex_monk, ex_ss)
% function fig_ = figBIAS_example(num, ex_monk, ex_ss)
%
% Plots:
%   1 filtered dirs, choices, resids
%   2 autocorrelation dirs, choices, resids
%   3 spectral analysis dirs, choices, resids

if nargin < 1 || isempty(num)
    num = 3;
end
if nargin < 2 || isempty(ex_monk)
    %ex_monk = 'Ava';
    ex_monk = 'Atticus';
    %ex_monk = 'Cyrus';
    %ex_monk = 'ZZ';
end

if nargin < 3 || isempty(ex_ss)
    ex_ss = 185;
end

%%
% set up the figure
%%
% units should be in inches, from wysifig
wid  = 3.0; % total width
hts  = 1.3;
cols = {1,1};
[axs,fig_] = getBIAS_axes(num, wid, hts, cols);

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
data         = FS_getDotsTrainingData(ex_monk);
sessions     = unique(data(:,1));
Lgood        = data(:,1) == sessions(ex_ss) & data(:,2) <= 2 & data(:,3)>=0;
data         = data(Lgood, [5 6 8 9 3]);
fcs          = getBIAS_fcs(ex_monk, Lgood);
fcs(1)       = 0;
ntr          = size(data,1);

Ln0 = data(:,1)>0;
disp(100.*[sum(data(Ln0,4)==1)./sum(Ln0) sum(data(Ln0,3)==1)./sum(Ln0)])

% useful variables
% co    = [1 0 0; (ones(3,1)*linspace(0.7,0,5))'];
co    = [0.7 0.7 0.7; 0 0 0; 1 0 0; 0 1 0; 0 0 1; 1 0 1];
st    = {'-', '-', '.-', '-', '-', '-'};
funs  = {@ddExp3};%; @ddExp4z; @ddExp4fz; @ddExp5fzz};
nfuns = length(funs);
fits  = cell(size(funs));
sdat  = [data(:,3:4) nans(size(data,1), nfuns)];

% do the fits
for ff = 1:nfuns
    if ff < 3
        [fits{ff},s,t,p,r] = ctPsych_fit(funs{ff}, data(:,1:4), data(:,5), [], []);
    else
        [fits{ff},s,t,p,r] = ctPsych_fit(funs{ff}, data(:,1:4), data(:,5), [], [], [], [], fcs);
    end        
    sdat(:,2+ff)    = r(:,4);
end

%% 
%   PANEL 1: PMF
%%
% axes(axs(1)); cla reset; hold on;
% Lgood = data(:,2) > 0.2 & data(:,2) < 0.5;
% 
% % plot ddExp3, ddExp4z
% ctPsych_plot(funs{1}, fits{1}, data(Lgood,1:4), data(Lgood,5), 1, [], axs(1), ...
%     {'line_style', '--'});
% ctPsych_plot(funs{2}, fits{2}, data(Lgood,1:4), data(Lgood,5), 1, [], axs(1), ...
%     {'line_style', '-'});

%% PANEL 1: Filtered dirs (gray), choices (black), 
%       resids 3 (magenta) ... REMOVED:, 4z (blue), 4fz (green), 5fzz (red)
axes(axs(1)); cla reset; hold on;

TAU   = 200;
cfilt = 1./exp((1:TAU*4)./TAU);
cfc   = filter(cfilt, sum(cfilt), sdat(2:end,:));
plot(cfc(:,1), 'k--');
plot(cfc(:,2), 'k-', 'LineWidth', 2);
plot(cfc(:,3), 'k-');

legend({'Stimulus directions', 'Choices', 'Choice residuals'});
plot([0 ntr], [0 0], 'k:');
axis([0 ntr -.2 .35]);
xlabel('Trial number');
ylabel('Filtered choice');


%% PANEL 2: Spectral analysis
%
% axes(axs(2)); cla reset; hold on;
% 
% for ii = 1:3%size(sdat,2)
%     [f, ps] = getBIAS_fft(sdat(2:end,ii));
%     plot(f(2:end), ps(2:end),  st{ii}, 'Color', co(ii,:));
%     plot(.09, ps(1), '.', 'Color', co(ii,:));
% end
% axis([0.09 f(end) 10e-9 1])
% set(gca, 'XScale', 'log', 'YScale', 'log');
% xlabel('Frequency (cycles/1000 trials)');
% ylabel('Power');

%% PANEL 3: Autocorrelation
%
axes(axs(2)); cla reset; hold on;

ACSZ = 50;
plot([0 ACSZ], [0 0], 'k:');

acd = xcorr(sdat(2:end,1), ACSZ, 'coeff');
plot(acd(ACSZ+2:end),  'k--');
acd = xcorr(sdat(2:end,2), ACSZ, 'coeff');
plot(acd(ACSZ+2:end),  'k-', 'LineWidth', 2);
acd = xcorr(sdat(2:end,3), ACSZ, 'coeff');
plot(acd(ACSZ+2:end),  'k-');

axis([0 ACSZ -.08 .15]);
xlabel('Trial lag');
ylabel('Correlation coefficient');
