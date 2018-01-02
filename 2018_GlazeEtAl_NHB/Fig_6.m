function Fig_6
%% Figure 6
%
% Figure-generating code for
%  Glaze CM, Filipowicz ALS, Kable JW, Balasubramanian V, and Gold JI
%  "A bias-variance trade-off governs individual differences in on-line 
%     learning in an unpredictable environment"

%% set up fig
wid     = 10; % total width
hts     = [5 2 2 2];
cols    = {1,1,1,2};
[axs,~] = getPLOT_axes(6, wid, hts, cols, 2, .5, [], 'Glaze et al', true);
set(axs,'Units','normalized');

%% Get session numbers, sort order
[~, analysis_data_dir, ~, sessions_per_subject] = getDataInfo;

%% A: Compare goodness of fits per model per subject
% celery columns are filename, struct, number of free parameters for each
% model
pnames = { ...
    ... %'allparam_null.mat',                                        'paramstrct_null',                  0; ...
    'allparam_fixedvarfits_bayesultimate_sigmafixed.mat',       'paramstrct_bayesultimate',         5; ...
    'allparam_fixedvarfits_partfilt_sigma_partMfixed5.mat',     'paramstrct_sigma_partMfixed5',     5; ...
    'allparam_fixedvarfits_partfilt_sigma_partMfixed10.mat',    'paramstrct_sigma_partMfixed10',    5; ...
    'allparam_fixedvarfits_partfilt_sigma_partMfixed20.mat',    'paramstrct_sigma_partMfixed',      5; ...
    'allparam_fixedvarfits_fixedLR_sigmafixed.mat',             'paramstrct_fixedLR_sigmafixed',    4};

num_subjects = length(sessions_per_subject);
num_models   = size(pnames, 1);
BICs         = nans(num_subjects, num_models);
for mm = 1:num_models
    load(fullfile(analysis_data_dir, pnames{mm,1}), pnames{mm,2});
    eval(['LLs = [' pnames{mm,2} '.LL]'';']);
    BICs(:,mm) = -2.*LLs+pnames{mm,3}.*log(sessions_per_subject.*2000);
end

% get adaptivity from 20-sample model
params = [paramstrct_sigma_partMfixed.params]';
logr0s = params(:,6);

% subtract 20-sample BIC from each subject, plot
axes(axs(1)); cla reset; hold on;
plot([0 num_models], [0 0], 'k-');
REFi = 4;
[~,I] = sort(logr0s);
for mm = 1:num_models
    if mm ~= REFi        
        BICs(:,mm) = BICs(:,mm) - BICs(:,REFi);
        for ss = 1:num_subjects
            plot(mm, BICs(I(ss),mm), 'ko', 'MarkerSize', 8, ...
                'MarkerFaceColor', ss/num_subjects.*ones(3,1));
        end
        [Rr,Pr] = corr(logr0s, BICs(:,mm), 'type', 'Spearman');
        P = signrank(BICs(:,mm));
        disp(sprintf('%d: r=%.2f, p=%.4f, medp=%.4f', mm, Rr, Pr, P))
        plot(mm+[-0.4 0.4], median(BICs(:,mm)).*[1 1], 'r-', 'LineWidth', 2);
    end
end

%% B,C: Correlation of 20-sample fits with adaptivity, choice variability

% 20-sample params are:
%  1: number of particles ... FIXED (20)
%  2: sensory noise       ... FREE
%  3: lapse               ... FREE
%  4: sigma0              ... UNUSED
%  5: gamma               ... UNUSED
%  6: log prior precision ... FREE
%  7: prior mean          ... FREE
%  8: log volatility (K)  ... FREE
part20fits = [paramstrct_sigma_partMfixed.params]';

% Load adaptivity model fits. Parameters are:
%  1: logJ0 = log(Hdefault/(1?Hdefault)
%  2: adaptivity (mH)
%  3: decision noise
%  4: lapse
load(fullfile(analysis_data_dir, 'adaptivityModelFits.mat'), 'fits');

% get correlations
ainds  = [2 3]; % indices of adaptive parameters to test
nainds = length(ainds);
pinds  = [2 3 6 7 8]; % indices of particle filter parameters to test
npinds = length(pinds);
cdat   = nans(npinds, 2, 2);

% loop through the parameters
for pp = 1:npinds
   for aa = 1:nainds
      [cdat(pp,aa,1), cdat(pp,aa,2)] = corr( ...
         part20fits(:,pinds(pp)), fits(:,ainds(aa)), 'type', 'Spearman');
   end
end

% do the plotz
for aa = 1:2
   axes(axs(aa+1)); cla reset; hold on;
   bar(cdat(:,aa,1));
   title(sprintf('%.4f %.4f %.4f %.4f %.4f', ...
      cdat(1,aa,1), cdat(2,aa,1), cdat(3,aa,1), ...
      cdat(4,aa,1), cdat(5,aa,1)))
   Fp = find(cdat(:,aa,2)<0.05/5);
   if ~isempty(Fp)
      plot(Fp, ones(length(Fp),1), 'r*')
   end
   axis([0 6 -1 1])
end
   
%% D,E Adaptivity (D) and Choice Variability (E) vs prior precision
for aa = 1:2
   axes(axs(aa+3)); cla reset; hold on;
   plot(part20fits(:,6), fits(:,ainds(aa)), 'ko');
   [R,P] = corr(part20fits(:,6), fits(:,ainds(aa)), 'type', 'Spearman');
   title(sprintf('rho=%.4f, p=%.4f', R, P))
   axis([-0.5 8.5 -0.2 1.2])
end
