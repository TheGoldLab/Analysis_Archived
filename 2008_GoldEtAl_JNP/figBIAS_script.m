%% FigBIAS_script
%
% script for generating figures. duh.
%

nfigs = 15;
figs  = nans(nfigs,1);

%% Figure 2: Residual illustration
%
ff=2; figs(ff) = figBIAS_resids(ff, 'Atticus', 185);

%% Figure 3: Example with
%   run mean, spectral analysis and autocorrelation
ff=3; figs(ff) = figBIAS_example(ff, 'Atticus', 17);

%% Figure 4: Summary of spectral analysis and autocorrelation
%   using no bias, DC bias models
ff=4; figs(ff) = figBIAS_residSummary(ff);

%% Figure 5: Explore different fits to wk
% 
ff=5; figs(ff) = figBIAS_wkFits(ff, 'Atticus', 185);

%% Figure 6: Fits to wk by monkey
%
ff=6; figs(ff) = figBIAS_wkByMonk(ff);

%% Figure 7: resid summary for models with different biases
%
ff=7; figs(ff) = figBIAS_pmfCompare(ff);

%% Figure 8: examples PMF low/high sensitivity
%
ff=8; figs(ff) = figBIAS_pmfExamples(ff);

%% Figure 9: Summary of PMF sensitivity (A) vs session
%
ff=9; figs(ff) = figBIAS_pmfParams(ff);

%% Figure 10: Summary of bias magnitude
%
ff=10; figs(ff) = figBIAS_biasMagnitude(ff);

%% Figure 11: Deviation summary
%
ff=11; figs(ff) = figBIAS_deviationSummary(ff);

%% Figure 12: MT summary
%
ff=12; figs(ff) = figBIAS_MTSummary(ff, 'Cyrus', 82);

%% Figure 13: LIP summary
%
ff=13; figs(ff) = figBIAS_LIPSummary(ff, 'Cyrus', 50);

%% Figure 14: LIP fits
%
ff=14; figs(ff) = figBIAS_LIPROC(ff);

%% Table 1: Saccade summary
%
tableBIAS_saccadeSummary;


for ff = 1:nfigs
    if isfinite(figs(ff))
        print(figs(ff),'-depsc','-loose',sprintf('../figs/raw/%d.eps', ff))
    end
end


%% OLD
%% Figure 12: Deviation fits
%
ff=12; figs(ff) = figBIAS_deviationFits(ff);


% Figure 7: Compare models
%
figs(7) = figBIAS_ddFits(7);

%% Figure 8: Verify local bias
%
figs(8) = figBIAS_verify(8);

%% Figure 12: MT example
%
figs(12) = figBIAS_MTExample(12, 'Cyrus');

%% Figure 14: LIP example
%
figs(14) = figBIAS_LIPExample(14, 'Cyrus');

