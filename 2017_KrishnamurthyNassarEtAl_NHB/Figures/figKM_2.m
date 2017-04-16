%% ************************************************************************************
%
%                      FIGURE 2 : ERRORS vs SAC WITH CONTRASTS 
%
% *************************************************************************************
% 
%  THREE COLUMNS:
% 	C1: combine ?initial? and ?subsequent? session data into ?no task? data
% 	C2: prediction data
% 	C3: perceptual report data
% and three rows:
% 	R1: example subject data, scatter ?reported? versus ?actual
% 		- just show raw data points, not medians
% 		- for prediction data (C2), maybe show TACP=1 data in a different color
% 	R2: summary of bias
% 		- option #1: histogram showing average bias per subject. code in terms of ?bias towards straight ahead?
% 		- option #2: scatter plot of reported versus actual, using median values computed per subject for data binned by actual
% 	R3: summary of variance: histogram of per-subject variances re: binned medians



% Josh : the point of R2 is that, on average, subjects were relatively 
% accurate in reporting and predicting sound locations. I?m thinking now a
%line+ribbon plot, where the line is the median of medians, computed per subject,
%of reported/predicted angle per binned true angle, and the ribbon is IQR. 
%So we?d say something like ?across all task conditions, subjects tended to
%make, on average, accurate reports and predictions of the simulated 
%sound-source location, albeit with a slight bias towards frontal space at the lateral extremes?


% Josh : the point of R3 is to quantify the variability in those reports/predictions?
% as I said, I think maybe just as a histogram of per-subject variances re:binned medians. 
%So we?d say something like ?They also exhibited trial-by-trial variability in 
%their perceptual reports and predictions, reflecting inherent uncertainty 
%(inverse of reliability) in both processes. Our goal was to understand how 
%the relative reliabilities of sensory localization cues (C1) and dynamically 
%changing predictions (C2) are combined to govern perceptual reports (C3) 
%in an unpredictable and noisy environment."



getCbColors;

% Use the J. Neurosci figure template
fig_num = 2; % figure number
total_width = 17; % total width
col_height = [4.5, 4.5, 4.5]; % height of each column
cols_per_row = {[0.33, 0.33,0.33], [0.33, 0.33,0.33],[0.33, 0.33,0.33]}; % columns in each row
psh = 1.5; % panel separation height
psw = 2; % panel separation width


[axs,fig_] = getPLOT_axes(fig_num, total_width, col_height, cols_per_row,...
                                 psh, psw, [], '');
set(axs,'Units','normalized');


ax_font_size = 12;
ax_font_size1 = 12;
ms = 8; % marker size;
jit = 0.1; % jitter for smartJitter

markerEdgeCol1 = [0.55 0.55 0.55];
markerEdgeCol2 = [0.1 0.1 0.1];
markerFaceCol1 = [0.75 0.75 0.75];


medianMarkerCol = [0.9 0.2 0.2];
markerFaceAlpha1 = 0.3;
markerEdgeAlpha1 = 0.3;
markerFaceAlpha2 = 0.3;
markerEdgeAlpha2 = 0.5;
markerSize = 16;


diagWidth = 2.5;
%%  --------------------- Row 1 : Example subject data -------------
col1 = [0 161 255]*1.0/255;
col1 = cbColors(4,:);
%col1 = [0.9 0.05 0.05];

example_subj_ID = 4; % pick ID of the example subject
example_idx = training_all_ID == example_subj_ID;
example_training_estimate = training_all_estimate(example_idx);
example_training_outcome = training_all_outcome(example_idx);


selSubj = probeSubjID == example_subj_ID;
sel_CP = probeTAC == 1;
example_pred_CP = probePrediction(selSubj & sel_CP);
example_pred_nonCP = probePrediction(selSubj & ~sel_CP);
example_outcome = probeOutcomes(selSubj);
example_outcome_CP = probeOutcomes(selSubj & sel_CP);
example_outcome_nonCP = probeOutcomes(selSubj & ~sel_CP);
example_mean_nonCP = probeMean(selSubj & ~sel_CP);
example_mean_CP = probeMean(selSubj & sel_CP);
example_estimate = probeEstimate(selSubj);

axes(axs(1)); cla(gca); hold on

% h1 = plot(example_training_outcome,example_training_estimate,'o','markersize',3,'markerfacecolor',cbColors(2,:),...
%              'markeredgecolor',cbColors(2,:));


h2 = plot([-20, 200], [-20, 200], '--','linewidth',diagWidth,'color','k');


% h1 = plot(example_training_outcome,example_training_estimate,'o','markersize',3,'markerfacecolor','None',...
%              'markeredgecolor',markerEdgeCol);

h1 = scatter(example_training_outcome,example_training_estimate,markerSize,'filled',...
    'o','markerfacecolor','None','markeredgecolor',markerEdgeCol2,...
      'markerfacealpha',markerFaceAlpha1,'markeredgealpha',markerEdgeAlpha2);
         
xlim([-20,200])
ylim([-20,200])

%xlabel('True angle')
ylabel('Reported angle ','fontsize',ax_font_size1)
title({'Control task','estimation'},'fontsize',10);
grid on
set(gca,'fontsize',ax_font_size)
axis square
setPLOT_panelLabel(gca, 1,0.18,1.2);
set(gca,'YTick',get(gca,'XTick'))


axes(axs(2)); cla(gca); hold on

% h1 = plot(example_mean_nonCP,example_pred_nonCP,'o','markersize',3,'markerfacecolor',cbColors(2,:),...
%              'markeredgecolor',cbColors(2,:));
% h2 = plot(example_mean_CP,example_pred_CP,'o','markersize',3,'markerfacecolor',col1,...
%              'markeredgecolor',col1);

h3 = plot([-20, 200], [-20, 200], '--','linewidth',diagWidth,'color','k');


% h1 = plot(example_mean_nonCP,example_pred_nonCP,'o','markersize',3,'markerfacecolor','None',...
%              'markeredgecolor',markerEdgeCol);
% h1.LineWidth = 0.2;

h2 = scatter(example_mean_CP,example_pred_CP,0.7*markerSize,'filled',...
    'o','markerfacecolor',markerEdgeCol1,'markeredgecolor',markerEdgeCol1,...
      'markerfacealpha',0.5,'markeredgealpha',0.5,'linewidth',1);

  

h1 = scatter(example_mean_nonCP,example_pred_nonCP,markerSize,'filled',...
    'o','markerfacecolor','None','markeredgecolor',markerEdgeCol2,...
      'markerfacealpha',markerFaceAlpha1,'markeredgealpha',markerEdgeAlpha2);


% h2 = plot(example_mean_CP,example_pred_CP,'o','markersize',3,'markerfacecolor',[0 0 0],...
%              'markeredgecolor',[0 0 0]);


         
xlim([-20,200])
ylim([-20,200])

%xlabel('True angle')
ylabel('Predicted angle ','fontsize',ax_font_size1)
title({'Dynamic task','prediction'},'fontsize',10);
grid on
set(gca,'fontsize',ax_font_size)
axis square
setPLOT_panelLabel(gca, 2,0.18,1.2);
set(gca,'YTick',get(gca,'XTick'))


axes(axs(3)); cla(gca); hold on

% h1 = plot(example_outcome,example_estimate,'o','markersize',3,'markerfacecolor',cbColors(2,:),...
%              'markeredgecolor',cbColors(2,:));


h3 = plot([-20, 200], [-20, 200], '--','linewidth',diagWidth,'color','k');


% h1 = plot(example_outcome,example_estimate,'o','markersize',3,'markerfacecolor','None',...
%              'markeredgecolor',markerEdgeCol);
         
h1 = scatter(example_outcome,example_estimate,markerSize,'filled',...
    'o','markerfacecolor','None','markeredgecolor',markerEdgeCol2,...
      'markerfacealpha',markerFaceAlpha1,'markeredgealpha',markerEdgeAlpha2);
         
xlim([-20,200])
ylim([-20,200])

%xlabel('True angle')
ylabel('Reported angle ','fontsize',ax_font_size1)
title({'Dynamic task','estimation'},'fontsize',10);
grid on
set(gca,'fontsize',ax_font_size)
axis square
setPLOT_panelLabel(gca,3,0.18,1.2);
set(gca,'YTick',get(gca,'XTick'))





%%  --------------------- Row 2 : Binned medians to illustrate bias -------------
col1 = [0 161 255]*1.0/255;
%col1 = [0.9 0.1 0.1];

% ---------------- TRAINING ------------------------
median_subj_reports_training = [];
training_angles = unique(training_all_outcome);

for i0 = 1:length(unique(ID))
    for i1 = 1:length(training_angles)
       
        sel_subj = training_all_ID == i0 & training_all_outcome == training_angles(i1);
        
        aux_reports = training_all_estimate(sel_subj);
        median_subj_reports_training(i0,i1) = median(aux_reports);
    end
end

median_reports_training = median(median_subj_reports_training,1);

median_iqr_training = iqr(median_subj_reports_training);



axes(axs(4)); cla(gca); hold on
h3 = plot([-10, 180], [-10, 180], '--','linewidth',diagWidth,'color','k');

% h2 = plot(training_angles,median_subj_reports_training,'o','markersize',3,'markerfacecolor','None',...
%              'markeredgecolor',markerEdgeCol);
aux1 = repmat(training_angles',29,1);
h2 = scatter(aux1(:),median_subj_reports_training(:),markerSize,'filled',...
    'o','markerfacecolor','None','markeredgecolor',markerEdgeCol2,...
      'markerfacealpha',markerFaceAlpha1,'markeredgealpha',markerEdgeAlpha2);

         
h1 = plot(training_angles,median_reports_training,'o','markersize',7,'markerfacecolor',medianMarkerCol,...
             'markeredgecolor','none');
ylim([-10 180])
xlim([-10 180])
%xlabel('True angle')
ylabel('Median estimate ','fontsize',ax_font_size1)
xlabel('True angle','fontsize',ax_font_size1)

grid on
set(gca,'fontsize',ax_font_size)
axis square
setPLOT_panelLabel(gca, 4,0.18,1.15);
set(gca,'YTick',get(gca,'XTick'))

%---------------------------------------------------





%----------------- PREDICTION ----------------------




median_reports_pred = nanmedian(median_subj_reports_pred);
median_means_pred = nanmedian(median_subj_means_pred);

axes(axs(5)); cla(gca); hold on

h3 = plot([-10, 180], [-10, 180], '--','linewidth',diagWidth,'color','k');

% 
% h1 = plot(median_subj_means_pred,median_subj_reports_pred,'o','markersize',3,'markerfacecolor','None',...
%              'markeredgecolor',markerEdgeCol);

h1 = scatter(median_subj_means_pred(:),median_subj_reports_pred(:),markerSize,'filled',...
    'o','markerfacecolor','None','markeredgecolor',markerEdgeCol2,...
      'markerfacealpha',markerFaceAlpha1,'markeredgealpha',markerEdgeAlpha2);
  
         
h2 = plot(median_means_pred,median_reports_pred,'-','color',medianMarkerCol,...
             'linewidth',3);   

ylim([-10 180])
xlim([-10 180])
%xlabel('True angle')
ylabel('Median prediction ','fontsize',ax_font_size1)
xlabel('True angle','fontsize',ax_font_size1)

grid on
set(gca,'fontsize',ax_font_size)
axis square
setPLOT_panelLabel(gca, 5,0.18,1.15);
set(gca,'YTick',get(gca,'XTick'))


%---------------------------------------------------



% ---------------- ESTIMATION ------------------------



median_reports_est = nanmedian(median_subj_reports_est);
median_outcomes_est = nanmedian(median_subj_outcomes_est);

axes(axs(6)); cla(gca); hold on

h3 = plot([-10, 180], [-10, 180], '--','linewidth',diagWidth,'color','k');


% h1 = plot(median_subj_outcomes_est,median_subj_reports_est,'o','markersize',3,'markerfacecolor','None',...
%              'markeredgecolor',markerEdgeCol);

h1 = scatter(median_subj_outcomes_est(:),median_subj_reports_est(:),markerSize,'filled',...
    'o','markerfacecolor','None','markeredgecolor',markerEdgeCol2,...
      'markerfacealpha',markerFaceAlpha1,'markeredgealpha',markerEdgeAlpha2);
  

h2 = plot(median_outcomes_est,median_reports_est,'-','color',medianMarkerCol,...
             'linewidth',3);   

ylim([-10 180])
xlim([-10 180])
%xlabel('True angle')
ylabel('Median estimate ','fontsize',ax_font_size1)
xlabel('True angle','fontsize',ax_font_size1)


grid on
set(gca,'fontsize',ax_font_size)
axis square
setPLOT_panelLabel(gca, 6,0.18,1.15);
set(gca,'YTick',get(gca,'XTick'))

%---------------------------------------------------








%%  --------------------- Row 3 : Summary of variance -------------

n_bins = 10;

col1 = [0 161 255]*1.0/255;
%col1 = [0.9 0.1 0.1];

% ---------------- TRAINING ------------------------
std_subj_reports_training = [];
training_angles = unique(training_all_outcome);

for i0 = 1:length(unique(ID))
       
        sel_subj = training_all_ID == i0;
        
        aux_reports = training_all_estimate(sel_subj);
        aux_outcomes = training_all_outcome(sel_subj);
        std_subj_reports_training(i0) = std(aux_reports-aux_outcomes);
    
end





%---------------------------------------------------



bin_start = -20;
window_size = 30;
window_step = 15;
bin_end = 200 - window_size;
bin_edges = bin_start:window_step:bin_end;


%----------------- PREDICTION ----------------------

std_subj_reports_pred = [];
median_subj_outcomes_pred = [];

for i0 = 1:length(unique(ID))
        sel_subj = probeSubjID == i0 & probeTAC ~= 1;
       
        aux_pred = probePrediction(sel_subj);
        aux_outcome = probeOutcomes(sel_subj);
        
        std_subj_reports_pred(i0) = std(aux_pred-aux_outcome);
end


% ---------------- ESTIMATION ------------------------

std_subj_reports_est = [];

for i0 = 1:length(unique(ID))
        sel_subj = probeSubjID == i0 & probeTAC ~= 1;
       
        aux_est = probeEstimate(sel_subj);
        aux_outcome = probeOutcomes(sel_subj);
        
        std_subj_reports_est(i0) = std(aux_est-aux_outcome);
end


axes(axs(7)); cla(gca); hold on



%h1 = histogram(std_subj_reports_est,n_bins);

%plot(std_subj_reports_training,std_subj_reports_est, 'o','markersize',6, 'markerfacecolor',[0.8 0.8 0.8], ...
%          'markeredgecolor',cbColors(1,:),'linewidth',1);
      

h1 = scatter(std_subj_reports_training,std_subj_reports_est,36,'filled',...
    'o','markerfacecolor','None','markeredgecolor',markerEdgeCol2,...
      'markerfacealpha',markerFaceAlpha1,'markeredgealpha',1.0);
  
plot([5,35],[5, 35],'--','linewidth',2,'color','k')
xlim([5,35])
ylim([5,35])

xlabel('Training STD','fontsize',12)
ylabel('Estimation STD','fontsize',12)
grid on
set(gca,'fontsize',ax_font_size)
axis square
setPLOT_panelLabel(gca, 7,0.18,1.15);
set(gca,'YTick',get(gca,'XTick'))




axes(axs(8)); cla(gca); hold on


% h0 = plot( std_subj_reports_pred,std_subj_reports_est, 'o','markersize',6, 'markerfacecolor',[0.8 0.8 0.8], ...
%           'markeredgecolor',cbColors(1,:),'linewidth',1);
      
      
h0 = scatter(std_subj_reports_pred,std_subj_reports_est,36,'filled',...
    'o','markerfacecolor','None','markeredgecolor',markerEdgeCol2,...
      'markerfacealpha',markerFaceAlpha1,'markeredgealpha',1.0);
  
plot([10,35],[10, 35],'--','linewidth',2,'color','k')

xlim([8,31])
ylim([8,31])

xlabel('Prediction STD','fontsize',12)
ylabel('Estimation STD','fontsize',12)

%xlim([10 40])
grid on
set(gca,'fontsize',ax_font_size)
axis square

%xlabel('STD(Reported angle)','fontsize',ax_font_size)
%xlim([0 40])

grid on
set(gca,'fontsize',ax_font_size)
setPLOT_panelLabel(gca, 8,0.18,1.15);
set(gca,'YTick',get(gca,'XTick'))


%---------------------------------------------------

axes(axs(9)); cla(gca); hold on


expected_std_subj_reports_est = 1./(sqrt(1./(std_subj_reports_training).^2 + ...
                                        1./(std_subj_reports_pred).^2));
%h1 = histogram(std_subj_reports_est,n_bins);

% a = plot(expected_std_subj_reports_est,std_subj_reports_est,  'o','markersize',6, 'markerfacecolor',[0.8 0.8 0.8], ...
%           'markeredgecolor',cbColors(1,:),'linewidth',1);

      
a = scatter(expected_std_subj_reports_est,std_subj_reports_est,36,'filled',...
    'o','markerfacecolor','None','markeredgecolor',markerEdgeCol2,...
      'markerfacealpha',markerFaceAlpha1,'markeredgealpha',1.0);
plot([8,25],[8, 25],'--','linewidth',2,'color','k')
xlim([8,25])
ylim([8,25])

xlabel({'Expected', 'Estimation STD'},'fontsize',12)
ylabel('Estimation STD','fontsize',12)
grid on
set(gca,'fontsize',ax_font_size)
axis square

setPLOT_panelLabel(gca, 9,0.18,1.15);
set(gca,'YTick',get(gca,'XTick'))

%---------------------------------------------------





kk=annotation('textbox');
set(kk, 'string', 'Krishnamurthy et al, figure 2', 'position', ...
    [0.77 0.9 0.15 0.05], 'EdgeColor', 'none','fontsize',14)


%% OUTPUT STATS FOR CORRELATIONS:

% first do it for the dyamic task
corr_pred_vs_true = [];
corr_est_vs_true = [];
corr_trainingTask = [];
for i1 = 1:length(uniqueIDs)
    sel_subj_noCP = probeSubjID == i1 & probeDay1Sess' == 0 & probeTAC ~= 1;
    sel_subj = probeSubjID == i1 & probeDay1Sess' == 0;
    
    aux_1 = probePrediction(sel_subj_noCP)';
    aux_2 = probeOutcomes(sel_subj_noCP)';
    [c,p] = corr(aux_1,aux_2);
    corr_pred_vs_true(i1) = c;
    
    aux_1 = probeEstimate(sel_subj)';
    aux_2 = probeOutcomes(sel_subj)';
    [c,p] = corr(aux_1,aux_2);
    corr_est_vs_true(i1) = c;
    
    sel_subj_training = training_all_ID == i1;
    aux_train_estimate = training_all_estimate(sel_subj_training);
    aux_train_outcome = training_all_outcome(sel_subj_training);
    [c,p] = corr(aux_train_estimate,aux_train_outcome);
    corr_trainingTask(i1) = c;
    
end


median_corr_pred_vs_true = median(corr_pred_vs_true);
iqr_corr_pred_vs_true = quantile(corr_pred_vs_true,[0.25, 0.75]);

median_corr_est_vs_true = median(corr_est_vs_true);
iqr_corr_est_vs_true = quantile(corr_est_vs_true,[0.25, 0.75]);

median_corr_trainingTask = median(corr_trainingTask);
iqr_corr_trainingTask = quantile(corr_trainingTask,[0.25, 0.75]);

clc
fprintf(' Median corr. Pred vs. True ---> %f  and IQR ---> [%f %f] \n',...
             median_corr_pred_vs_true, iqr_corr_pred_vs_true(1),iqr_corr_pred_vs_true(2));
fprintf(' Median corr. Est vs. True ---> %f  and IQR ---> [%f %f] \n',...
             median_corr_est_vs_true, iqr_corr_est_vs_true(1),iqr_corr_est_vs_true(2));

fprintf(' Median corr. Training ---> %f  and IQR ---> [%f %f] \n',...
             median_corr_trainingTask, iqr_corr_trainingTask(1),iqr_corr_trainingTask(2));


