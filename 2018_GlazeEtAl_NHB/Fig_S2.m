%% Supplementary Figure 2
%
% Adaptivity versus prior precision for the particle-filter and 
%  Bayesian models
%
% Figure-generating code for:
%  Glaze CM, Filipowicz ALS, Kable JW, Balasubramanian V, and Gold JI
%  "A bias-variance trade-off governs individual differences in on-line 
%     learning in an unpredictable environment"
%

%% set up the fig
wid     = 11.6; % total width
ht      = 6;
cols    = {1};
[axs,~] = getPLOT_axes(5, wid, ht, cols, 2, 2, [], 'Glaze et al', true);
set(axs,'Units','normalized');
set(axs, 'FontSize', 12);

%% Load the data
% Where to find the data
[file_list, analysis_data_dir, raw_data_dir] = getDataInfo;

% load it
load(fullfile(analysis_data_dir, 'bayesian_adaptive_sims_logKfixed.mat'));

% set up the axes
axes(axs(1)); cla reset; hold on;

% particle filter simulations
plot(logr0spc(1:2:end),alphamat(:,1:2:end)','linewidth',2)

% Bayesian simulation
plot(logr0spc(1:2:end),alphamatb(1:2:end),'k-','linewidth',2)

% make it pretty
axis([-1.5 8 -.1 1.3])
legend({'sample M=5','10','20','1000','ideal observer'},'fontsize',12)
set(gca,'box','off','ticklength',[0 0],'fontsize',12)

