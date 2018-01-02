%% Supplementary Figure 5
%
% Recover fits using the particle-filter model
%
% Figure-generating code for:
%  Glaze CM, Filipowicz ALS, Kable JW, Balasubramanian V, and Gold JI
%  "A bias-variance trade-off governs individual differences in on-line 
%     learning in an unpredictable environment"
%

%% set up the fig
wid     = 11.6; % total width
ht      = 4;
cols    = {3,3};
[axs,~] = getPLOT_axes(5, wid, ht, cols, 2, 2, [], 'Glaze et al', true);
set(axs,'Units','normalized');
set(axs, 'FontSize', 12);

%% Load the data
% Where to find the data
[file_list, analysis_data_dir, raw_data_dir] = getDataInfo;

% load it
load(fullfile(analysis_data_dir, 'allparam_fixedvarfits_partfilt_sigma_partMfixed_sim.mat'));

% collect some parameters
paramsrl = [psivc(:) whoopsvc(:) logr0vc(:) H0vc(:) logKvc(:)];
params = [paramstrct_sigma_partMfixed.params]';
params = params(:,[2 3 6:8]);

%% report correlations
diag(corr(params,paramsrl(1:size(params,1),:),'type','Spearman'))

%% Sensory noise
axes(axs(1)); cla reset; hold on;
plot(paramsrl(1:size(params,1),1),params(:,1),'ko')
hold on
plot([0 5],[0 5],'k--')
xlabel('simulated sensory noise','fontsize',14)
ylabel('fit sensory noise','fontsize',14)

%% Lapse rate
axes(axs(2)); cla reset; hold on;
plot(paramsrl(1:size(params,1),2),params(:,2),'ko')
hold on
plot([0 .01],[0 .01],'k--')
xlabel('simulated lapse rate','fontsize',14)
ylabel('fit lapse rate','fontsize',14)

%% log prior precision
axes(axs(4)); cla reset; hold on;
plot(paramsrl(1:size(params,1),end-2),params(:,end-2),'ko')
hold on
plot([-1 6],[-1 6],'k--')
xlabel('simulated log(r_0)','fontsize',14)
ylabel('fit log(r_0)','fontsize',14)
xlim([-2 9])
ylim([-2 9])

%% prior mean
axes(axs(5)); cla reset; hold on;
plot(paramsrl(1:size(params,1),end-1),params(:,end-1),'ko')
hold on
plot([0 1],[0 1],'k--')
xlabel('simulated H_0','fontsize',14)
ylabel('fit H_0','fontsize',14)

%% log volatility
axes(axs(6)); cla reset; hold on;
plot(paramsrl(1:size(params,1),end),params(:,end),'ko')
hold on
plot([-9 -2],[-9 -2],'k--')
xlabel('simulated log(K)','fontsize',14)
ylabel('fit log(K)','fontsize',14)
% xlim([-10 -2])
% ylim([-10 -2])

% make 'em pretty
set(axs,'ydir','normal','box','off','ticklength',[0 0],'fontsize',10)

