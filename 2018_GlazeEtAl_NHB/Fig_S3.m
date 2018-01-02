%% Supplementary Figure 3
%
% Properties of the particle-filter and Bayesian models
%
% Figure-generating code for:
%  Glaze CM, Filipowicz ALS, Kable JW, Balasubramanian V, and Gold JI
%  "A bias-variance trade-off governs individual differences in on-line 
%     learning in an unpredictable environment"
%

%% set up the fig
wid     = 17.6; % total width
hts     = [2 1.5 1.5 4];
cols    = {5,3,3,3};
[axs,~] = getPLOT_axes(3, wid, hts, cols, 2, 2, [], 'Glaze et al', true);
set(axs,'Units','normalized');
set(axs, 'FontSize', 12);

%% Panel A: adaptivity vs volatility and prior precisision 
% Load the simulated data
% Where to find the data
[file_list, analysis_data_dir, raw_data_dir] = getDataInfo;

% load the data
load(fullfile(analysis_data_dir, 'bayesian_adaptive_sims1000.mat'));
alphamat1000 = squeeze(alphamat(1,:,:));

load(fullfile(analysis_data_dir, 'bayesian_adaptive_sims3.mat'));
alphamat2 = zeros(4,size(alphamat,2),size(alphamat,3));
alphamat2(1:3,:,:) = alphamat;
alphamat2(4,:,:) = alphamat1000;
alphamat = alphamat2;

% Bayesian model
axes(axs(1));% cla reset; hold on;
imagesc(logr0spc,logKspc(1:end),alphamatb(:,1:end)',[0 .95]);

% Particle-filter models
for ii = 1:4
   axes(axs(1+ii)); %cla reset; hold on;
   imagesc(logr0spc,logKspc(1:end),squeeze(alphamat(5-ii,:,1:end))');
end

% make 'em pretty
set(axs,'ydir','normal','box','off','ticklength',[0 0],'fontsize',12,'clim',[0 .95])
axes(axs(1));
xlabel('Prior precision');
ylabel('Volality');

%% B?D: Posterior distributions for moderate vs low volatility
% load data to use for simulations
[~,analysis_data_dir] = getDataInfo;
load(fullfile(analysis_data_dir, 'model_testbed1.mat'), ...
   'data', 'choice', 'x', 'musgn');

% model parameters, etc
MU_DISTANCE = 30; % distance between triangles
Ks          = [0.01 0.001];
r0          = exp(0.7);
H0          = 0.5;
samples     = [1000 20];
xax         = 1:length(x);

% likelihoods 'n stuff
l1 = normpdf(x, MU_DISTANCE/2,data.sigma);
l2 = normpdf(x,-MU_DISTANCE/2,data.sigma);
F  = data.sigma^2/MU_DISTANCE;

% collect and plot estimates of H
% 1. Ideal-observer model
Hspc = linspace(0.01,0.99,100)';
% for each volatility
for kk = 1:2
   
   % do the simulation
   pH0 = betapdf(Hspc,r0*H0,r0*(1-H0));
   [q1exp,expH_bayes,pHmat_bayes] = ...
      HHMM_mixed_mat_c(l1,l2,Hspc,pH0./sum(pH0),Ks(kk),musgn);
   
   % plot it
   axes(axs(6+(kk-1)*3)); cla reset; hold on;
   imagesc(1:2000,Hspc,pHmat_bayes)
   % plot(xax, expH_bayes, 'k-', 'linewidth', 1.5);
   axis([1 2000 0 1]);
end

%% 2. 20-sample model
H_M         = numel(Hspc);
% for each volatility
for kk = 1:2
   % for each number of samples
   for mm = 1:2
      
      % do the simulation
      H0samps = betarnd(H0*r0,(1-H0)*r0,50000,1);
      [q1samp,Hsamps,Hsampmat] = particle_filter_learnH6_mat( ...
         Ks(kk),H0samps,F,x,musgn,samples(mm),2,zeros(200,1));
      
      % compute prob(H)
      [partM,N] = size(Hsampmat);
      pHmat     = zeros(H_M,N);
      for nn = 1:N
         [N,xout] = hist(Hsampmat(:,nn),Hspc);
         pHmat(:,nn) = N;
      end
      
      % plot it
      axes(axs(7+(kk-1)*3+mm-1)); cla reset; hold on;
      imagesc(1:2000,Hspc,pHmat);
      
      %hs=plot(xax, Hsamps', 'Linewidth', 1.5);
      %set(hs(1), 'Color', 'k');
      %set(hs(2), 'Color', 0.6.*ones(3,1));
      axis([1 2000 0 1]);
   end
end

%% E: Choice variability vs Adaptivity for 20-particle model, varying K
load(fullfile(analysis_data_dir, 'regression_sim_Kvary.mat'), 'paramstrct_regress');

% collect the params
params_regress = [paramstrct_regress.params]';
params_regress = reshape(params_regress,4,41)';

% report correlation coefficient
[r,p] = corr(params_regress(:,2),params_regress(:,3),'type','Spearman')

% plot it
axes(axs(12)); cla reset; hold on;
plot(params_regress(:,2),params_regress(:,3),'ko')
xlabel('Adaptivity');
ylabel('Choice variability');
axis([-0.1 1.8 0 1.5]);

