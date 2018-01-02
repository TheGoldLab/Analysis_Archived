function Fig_5
%% Figure 5
%
% Figure-generating code for
%  Glaze CM, Filipowicz ALS, Kable JW, Balasubramanian V, and Gold JI
%  "A bias-variance trade-off governs individual differences in on-line
%     learning in an unpredictable environment"
%
% Simiulations and adaptivity-choice variability trade off
% for three models:
%   1. ideal-observer model
%   2. 20-sample model
%   3. delta-rule model

%% set up plotz
wid     = 8; % total width
hts     = [0.8.*ones(1,7) 3 3];
cols    = {1,1,1,1,1,1,1,2,2};
[axs,~] = getPLOT_axes(1, wid, hts, cols, 1, 1, [], 'Glaze et al', true);
set(axs,'Units','normalized');

% load data to use for simulations
[~,analysis_data_dir] = getDataInfo;
load(fullfile(analysis_data_dir, 'model_testbed1.mat'), ...
   'data', 'choice', 'x', 'musgn');

% model parameters, etc
MU_DISTANCE = 30; % distance between triangles
K           = 0.01;
r0s         = [2 exp(7)];
H0          = 0.5;
xax         = 1:length(x);

% likelihoods 'n stuff
l1 = normpdf(x, MU_DISTANCE/2,data.sigma);
l2 = normpdf(x,-MU_DISTANCE/2,data.sigma);
F  = data.sigma^2/MU_DISTANCE;

% show the objective hazard rate
axes(axs(1)); cla reset; hold on;
plot(xax, data.H(:,1), 'k--', 'LineWidth', 1.5);

%% collect H estimate
%% 1. Ideal-observer model
Hspc = linspace(0.01,0.99,100)';
for ii = 1:2
   
   % do the simulation
   pH0 = betapdf(Hspc,r0s(ii)*H0,r0s(ii)*(1-H0));
   [q1exp,expH_bayes,pHmat_bayes] = ...
      HHMM_mixed_mat_c(l1,l2,Hspc,pH0./sum(pH0),K,musgn);
   
   % plot it
   axes(axs(1+ii)); cla reset; hold on;
   imagesc(1:2000,Hspc,pHmat_bayes)
   % plot(xax, expH_bayes, 'k-', 'linewidth', 1.5);
   axis([1 2000 0 1]);
end

%% 2. 20-sample model
NUM_SAMPLES = 20;
H_M         = numel(Hspc);
for ii = 1:2
   
   % do the simulation
   H0samps = betarnd(H0*r0s(ii),(1-H0)*r0s(ii),50000,1);
   [q1samp,Hsamps,Hsampmat] = particle_filter_learnH6_mat( ...
      K,H0samps,F,x,musgn,NUM_SAMPLES,2,zeros(200,1));
   
   % compute pH
   [partM,N] = size(Hsampmat);
   pHmat     = zeros(H_M,N);
   for nn = 1:N
      [N,xout] = hist(Hsampmat(:,nn),Hspc);
      pHmat(:,nn) = N;
   end
   
   % plot it
   axes(axs(3+ii)); cla reset; hold on;
   imagesc(1:2000,Hspc,pHmat);
   
   %hs=plot(xax, Hsamps', 'Linewidth', 1.5);
   %set(hs(1), 'Color', 'k');
   %set(hs(2), 'Color', 0.6.*ones(3,1));
   axis([1 2000 0 1]);
end

%% 3. delta-rule model
alphas = exp([-2 -10]);
H0     = 0.5;
for ii = 1:2
   
   % do the simulation
   [Lexp,Lvar,Hexp] = fixedLR_learnH(alphas(ii),H0,F,x,musgn,1,zeros(length(x),1));
   
   % plot it
   axes(axs(5+ii)); cla reset; hold on;
   plot(xax,Hexp', 'k-', 'Linewidth', 1.5);
end

% prettify
set(axs(1:7),'box','off','fontsize',12,'XLim',xax([1 end]),'YLim', [0 1]);
set(axs(1:6),'xtick',[],'ytick',[]);

%% simulated b-v trade-off
sims = {...
   'regression_sim_r0vary_bayesultimate', ...
   'regression_sim_r0vary', ...
   'regression_sim_gammavary_delta_rule_normative'};

for mm = 1:3
   
   % get the best-fitting parameters from the full simulations
   load(fullfile(analysis_data_dir, [sims{mm} '.mat']), 'paramstrct_regress');
   pvals = reshape([paramstrct_regress.params]',4,[])';
   choice_variability = pvals(:,2);
   adaptivity         = pvals(:,3);
   
   % compute correlation coefficient between adaptivity and choice
   %  variability determined from the adaptivity model fit to the
   %  simulations
   [R,P] = corr(choice_variability, adaptivity, 'type', 'Spearman');
   
   % show the scatterplot
   axes(axs(7+mm)); cla reset; hold on;
   plot(choice_variability, adaptivity,'ko');
   set(gca,'ydir','normal','box','off','ticklength',[0 0],'fontsize',12)
   axis([-0.1 1.1 0 0.6])
   title(sprintf('R=%.4f, P=%.4f', R, P));
end
