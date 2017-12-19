function Fig_S6(recompute)
%% Supplementary Figure 6
%
% Figure/stats comparing matched training and test data.
%
% Figure-generating code for:
%  Glaze CM, Filipowicz ALS, Kable JW, Balasubramanian V, and Gold JI
%  "A bias-variance trade-off governs individual differences in on-line 
%     learning in an unpredictable environment"
%

%% set up the fig
wid     = 11.6; % total width
ht      = 4;
cols    = {2,2,2};
[axs,~] = getPLOT_axes(5, wid, ht, cols, 2, 2, [], 'Glaze et al', true);
set(axs,'Units','normalized');
set(axs, 'FontSize', 12);

%% possibly recompute
if nargin>0 && recompute
   fitTrainingData;
end

%% Where the data reside
[~,data_dir] = getDataInfo;

%% load Adaptivity/Particle filter fits to matched training/testing data
% Particle Filter parameters are:
%   1. H0
%   2. logR0
%   3. logK
%   4. vx_scale
%   5. lapse
%
% Adaptivity Model parameters are 
%   1. H_default
%   2. H_slope
%   3. noise in the DV
%   4. lapse rate
load(fullfile(data_dir, 'trainingDataFits.mat'));
fits_adaptivity     = models(1).params;
fits_particleFilter = models(2).params;

%% left column is training data, right is matched testing data
for dd = 1:2
   
   % for each parameter
   for pp = 1:3
      
      % set axis
      axes(axs((pp-1)*2+dd)); cla reset; hold on;
      set(gca, 'FontSize', 12);
      
      if pp == 1
         % b-v trade-off
         xs = fits_adaptivity(:,2,dd);
         ys = fits_adaptivity(:,3,dd);
         xlabel('Adaptivity (mH)')
         ylabel('Choice variability');
         
      elseif pp == 2
         % adaptivity vs logr0
         xs = fits_particleFilter(:,2,dd);
         ys = fits_adaptivity(:,2,dd);
         xlabel('Prior precision');
         ylabel('Adaptivity (mH)')
         
      elseif pp == 3
         % choice variability vs logr0
         xs = fits_particleFilter(:,2,dd);
         ys = fits_adaptivity(:,3,dd);
         xlabel('Prior precision');
         ylabel('Choice variability');
      end
      
      % plotz
      Lg = isfinite(xs) & isfinite(ys);
      [R,P] = corr(xs(Lg), ys(Lg), 'type', 'Spearman');
      plot(xs,ys,'ko');
      title(sprintf('R=%.4f, P=%.4f', R, P))
   end
end
