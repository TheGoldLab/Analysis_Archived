function Fig_7
%% Figure 7
%
% Model-free analyses of 20-sample model simulations, using fits
%  to the subjects' data
%
% Figure-generating code for
%  Glaze CM, Filipowicz ALS, Kable JW, Balasubramanian V, and Gold JI
%  "A bias-variance trade-off governs individual differences in on-line 
%     learning in an unpredictable environment"

%% Set up figure
wid     = 17.6; % total width
ht      = 4;
cols    = {3,3};
[axs,~] = getPLOT_axes(2, wid, ht, cols, 1.5, 1.5, [], 'Glaze et al', true);
set(axs,'Units','normalized','FontSize', 12);

%%  just call utility function
plotModelFreeAnalyses('simulations', axs(1:3), axs(4:6));

