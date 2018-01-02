function Fig_2(recompute)
%% Figure 2
%
% Figure-generating code for
%  Glaze CM, Filipowicz ALS, Kable JW, Balasubramanian V, and Gold JI
%  "A bias-variance trade-off governs individual differences in on-line 
%     learning in an unpredictable environment"

%% Set up figure
wid     = 17.6; % total width
ht      = 4;
cols    = {3,3,3};
[axs,~] = getPLOT_axes(2, wid, ht, cols, 1.5, 1.5, [], 'Glaze et al', true);
set(axs,'Units','normalized','FontSize', 12);

%% Get adaptivity model fits. 
%  These are generated using:
%  fitAdaptivityModel

% possibly re-do the fits
if nargin==1 && recompute
   fitAdaptivityModel;
end

% Where to find the data
[~, analysis_data_dir] = getDataInfo;

% Parameters are:
%  1: logJ0 = log(Hdefault/(1?Hdefault)
%  2: adaptivity (mH)
%  3: decision noise
%  4: lapse
load(fullfile(analysis_data_dir, 'adaptivityModelFits.mat'), 'fits');

%% B: Histogram of H default
axes(axs(2)); cla reset; hold on;
Hdefault = 1./(1+exp(-fits(:,1)));
[nout,x] = hist(Hdefault,0:.1:1);
bar(x,nout,'facecolor','k')
axis([0 1.0 0 20])
set(gca,'box','off','ticklength',[0 0],'fontsize',12,...
   'XTick',[0 0.5 1], 'YTick', 0:5:20)
xlabel('Hdefault');
ylabel('Count');
disp(sprintf('Hdefault: min=%.2f, median=%.2f, max=%.2f, p=%.10f', ...
   min(Hdefault), median(Hdefault), max(Hdefault), signtest(Hdefault)))

%% C: Histogram of mH
axes(axs(3)); cla reset; hold on;
mH = fits(:,2);
[nout,x] = hist(mH,-.2:.1:1.5);
bar(x,nout,'facecolor','k')
axis([-.3 1.2 0 10])
set(gca,'box','off','ticklength',[0 0],'fontsize',12,...
   'XTick',[0 0.5 1], 'YTick', 0:5:10)
xlabel('mH');
ylabel('Count');
disp(sprintf('mH: min=%.2f, median=%.2f, max=%.2f, p=%.10f', ...
   min(mH), median(mH), max(mH), signtest(mH)))

%% G - I: Pr switch vs trials after CP
%  just call utility function
plotModelFreeAnalyses('data', axs(7:9));

