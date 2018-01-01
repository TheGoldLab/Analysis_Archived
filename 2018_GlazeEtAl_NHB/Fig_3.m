function Fig_3(recompute)
%% Figure 3
%
% Figure-generating code for
%  Glaze CM, Filipowicz ALS, Kable JW, Balasubramanian V, and Gold JI
%  "A bias-variance trade-off governs individual differences in on-line 
%     learning in an unpredictable environment"

% possibly re-do the fits
if nargin==1 && recompute
   fitAdaptivityModel;
end
   
% Where to find the data
[~,data_dir] = getDataInfo;

% Load adaptivity model fits. Parameters are:
%  1: logJ0 = log(Hdefault/(1?Hdefault)
%  2: adaptivity (mH)
%  3: decision noise
%  4: lapse
load(fullfile(data_dir, 'adaptivityModelFits.mat'), 'fits');

%% Set up figure
wid     = 17.6; % total width
ht      = 4;
cols    = {3,3};
[axs,~] = getPLOT_axes(3, wid, ht, cols, 1.5, 1.5, [], 'Glaze et al', true);
set(axs,'Units','normalized','FontSize', 12);

%% B: Scatterplot of adaptivity vs decision variability
axes(axs(2)); cla reset; hold on;
adaptivity  = fits(:,2);
variability = fits(:,3);
plot(adaptivity, variability, 'ko');
[R,P] = corr(adaptivity, variability, 'type', 'Spearman')

