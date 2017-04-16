%% ************************************************************************************
%
%                      FIGURE 3 : ERRORS vs SAC WITH CONTRASTS 
%
% *************************************************************************************


getCbColors;

% Use the J. Neurosci figure template
fig_num = 2; % figure number
total_width = 17; % total width
col_height = [6, 6, 6]; % height of each column
cols_per_row = {[.6, .35], [.6, .35],[.6, .35]}; % columns in each row
psh = 1; % panel separation height
psw = 1.8; % panel separation width

[axs,fig_] = getPLOT_axes(fig_num, total_width, col_height, cols_per_row,...
                                 psh, psw, [], '');
set(axs,'Units','normalized');


ax_font_size = 14;
ax_font_size1 = 14;
ms = 8; % marker size
jit = 0.1; % jitter for smartJitter

theory_col1 = [0.4 0.4 0.4];
theory_col2 = [0.65 0.65 0.65];
theory_col3 = [0.8 0.8 0.8];

%%  --------------------- Row 1 : Prediction errors & contrasts -------------

axes(axs(1)); cla(gca); hold on

pred_err_lims = [8 37];

auxMu_hi = nanmean(subjPredErrSAC_hi');
auxSE_hi = nanstd(subjPredErrSAC_hi')./sqrt(size(subjPredErrSAC_hi',1));

auxMu_low = nanmean(subjPredErrSAC_low');
auxSE_low = nanstd(subjPredErrSAC_low')./sqrt(size(subjPredErrSAC_low',1));

barInds=[1:6; 1:6];
errBars_hi = [auxMu_hi + auxSE_hi ; auxMu_hi - auxSE_hi];
errBars_low = [auxMu_low + auxSE_low ; auxMu_low - auxSE_low];



auxMu_model_hi = nanmean(modelPredErrSAC_hi');
auxSE_model_hi = nanstd(modelPredErrSAC_hi')./sqrt(size(modelPredErrSAC_hi',1));

auxMu_model_low = nanmean(modelPredErrSAC_low');
auxSE_model_low = nanstd(modelPredErrSAC_low')./sqrt(size(modelPredErrSAC_low',1));

barInds=[1:6; 1:6];
errBars_model_hi = [auxMu_model_hi + auxSE_model_hi ; auxMu_model_hi - auxSE_model_hi];
errBars_model_low = [auxMu_model_low + auxSE_model_low ; auxMu_model_low - auxSE_model_low];

plot([0, 6], [20, 20], '--','linewidth',2,'color',cbColors(2,:))
plot([0, 6], [10, 10], '--','linewidth',2,'color',cbColors(3,:))

plot(barInds,errBars_hi, '-k', 'lineWidth', 1.5)
plot(barInds,errBars_low, '-k', 'lineWidth', 1.5)
h1 = plot(auxMu_hi, 'ob', 'markerSize', ms, 'markerFaceColor', cbColors(2,:), 'markerEdgeColor', 'k', 'lineWidth', 1);
h2 = plot(auxMu_low, 'ob', 'markerSize', ms, 'markerFaceColor', cbColors(3,:), 'markerEdgeColor', 'k', 'lineWidth', 1);

legend([h1,h2],{'high noise','low noise'},'Location','best','Box','off')

% plot(barInds,errBars_model_hi, '-k', 'lineWidth', 1.5)
% plot(barInds,errBars_model_low, '-k', 'lineWidth', 1.5)
plot(auxMu_model_hi, 'db', 'markerSize', ms, 'markerFaceColor', 'None', 'markerEdgeColor',cbColors(2,:) , 'lineWidth', 1)
plot(auxMu_model_low, 'db', 'markerSize', ms, 'markerFaceColor', 'None', 'markerEdgeColor', cbColors(3,:), 'lineWidth', 1)




ylim(pred_err_lims)
xlim([0.5, 6.5])
%xlabel('SAC')
ylabel('Prediction error STD','fontsize',ax_font_size1)
%title('Subject Pred Error','fontsize',14);
grid on
set(gca,'fontsize',ax_font_size)

set(gca,'XTickLabel', '');

setPLOT_panelLabel(gca, 1,0.12,1.05);


% Now plot the CP, TAC and Noise contrasts
axes(axs(2)); cla(gca); hold on


                                       

aux_contrast = [predErr_cp_contrast predErr_sac_contrast_lo predErr_sac_contrast_hi predErr_noise_contrast];

aux_model_contrast = [model_predErr_cp_contrast model_predErr_sac_contrast_low ...
                          model_predErr_sac_contrast_hi model_predErr_noise_contrast];

% Plot this using Matt's smart jitter function
% Note : to compare subject with benchmar, normalisation scale should be
% same

xTickLabels={'CP','SAC(low)','SAC(High)','Noise'}; % 'abs(90-Outcome) 10' ,'abs(90-Outcome) 20'
ll = size(aux_contrast,1);
plot([0 ll*6], [0 0], 'Color',[0.5 0.5 0.5],'linewidth',1)
clear P CI np_P

%which contrasts are significant at 5%?
sig_contrast_idx = [1 2 3 4];
contrast_edge_col_subj = [[0 0 0];cbColors(3,:);cbColors(2,:);[0 0 0]];
contrast_edge_col_theory = [theory_col2;cbColors(3,:);cbColors(2,:);theory_col2];
contrast_line_col = [theory_col3;theory_col3; theory_col3;theory_col3];

for i=1:size(aux_contrast,2)
    xJit = smartJitter(aux_contrast(:,i)./std(aux_contrast(:,i)),.05,5);

    xJit_model = smartJitter(aux_model_contrast(:,i)./std(aux_contrast(:,i)),.05,5);
    
    for i_line = 1:length(uniqueID)
        line([i*6+xJit(i_line), i*6 + 3 + xJit_model(i_line)], ...
             [aux_contrast(i_line,i)./std(aux_contrast(:,i)), aux_model_contrast(i_line,i)./std(aux_contrast(:,i))],...
              'color',contrast_line_col(i,:),'linewidth',0.25);
    end
    
    
    plot(ones(ll, 1).*i.*6+xJit_model + 3 , aux_model_contrast(:,i)./std(aux_contrast(:,i)),...
        'db', 'markerSize', 5, 'markerFaceColor','None',...
        'markerEdgeColor', contrast_edge_col_theory(i,:), 'lineWidth', 0.5);
    
    
    plot(ones(ll, 1).*i.*6+xJit , aux_contrast(:,i)./std(aux_contrast(:,i)),...
        'o', 'markerSize', 5, 'markerFaceColor',[0.96 0.96 0.96],...
        'markerEdgeColor', contrast_edge_col_subj(i,:), 'lineWidth', 1);
    
     xlim([4 6.*(i)+ 4.5])
     
     %y_max = max(aux_y) + 2.0;
    
    if(any(i==sig_contrast_idx))
        plot(i*6,15,'*','markersize',5,'color','k');
    end
end
%set(gca, 'xtick', [2:2:(ll)*2], 'XTickLabel', xTickLabels);
set(gca, 'xtick', [7.5:6:(ll)*6], 'XTickLabel', '');

set(gca,'fontsize',ax_font_size)
ylabel('Contrast','fontsize',ax_font_size1)
grid on


ylim([-5 15.5])
setPLOT_panelLabel(gca, 2,0.18,1.05);




%%  --------------------- Row 2 : Estimation errors & contrasts -------------

axes(axs(3)); cla(gca); hold on

est_err_lims = [7 28]; %26





plot([0 6.5], [nanmean(std_subj_reports_training), nanmean(std_subj_reports_training)],...
         '--k', 'linewidth',2)

auxMu_hi = nanmean(subjEstErrSAC_hi');
auxSE_hi = nanstd(subjEstErrSAC_hi')./sqrt(size(subjEstErrSAC_hi',1));

auxMu_low = nanmean(subjEstErrSAC_low');
auxSE_low = nanstd(subjEstErrSAC_low')./sqrt(size(subjEstErrSAC_low',1));

barInds=[1:6; 1:6];
errBars_hi = [auxMu_hi + auxSE_hi ; auxMu_hi - auxSE_hi];
errBars_low = [auxMu_low + auxSE_low ; auxMu_low - auxSE_low];




auxMu_theory_hi = nanmean(theory_subjective_est_err_hi');
auxSE_theory_hi = nanstd(theory_subjective_est_err_hi')./sqrt(size(theory_subjective_est_err_hi',1));

auxMu_theory_low = nanmean(theory_subjective_est_err_low');
auxSE_theory_low = nanstd(theory_subjective_est_err_low')./sqrt(size(theory_subjective_est_err_low',1));

barInds=[1:6; 1:6];
errBars_theory_hi = [auxMu_theory_hi + auxSE_theory_hi ; auxMu_theory_hi - auxSE_theory_hi];
errBars_theory_low = [auxMu_theory_low + auxSE_theory_low ; auxMu_theory_low - auxSE_theory_low];


%plot(barInds,errBars_theory_hi, '-k', 'lineWidth', 1.5)
%plot(barInds,errBars_theory_low, '-k', 'lineWidth', 1.5)
plot(auxMu_theory_hi, 'db', 'markerSize', ms, 'markerFaceColor', 'None', 'markerEdgeColor',cbColors(2,:) , 'lineWidth', 1)
plot(auxMu_theory_low, 'db', 'markerSize', ms,'markerFaceColor', 'None', 'markerEdgeColor',cbColors(3,:) , 'lineWidth', 1)


plot(barInds,errBars_hi, '-k', 'lineWidth', 1.5)
plot(barInds,errBars_low, '-k', 'lineWidth', 1.5)
plot(auxMu_hi, 'ob', 'markerSize', ms, 'markerFaceColor', cbColors(2,:), 'markerEdgeColor', 'k', 'lineWidth', 1)
plot(auxMu_low, 'ob', 'markerSize', ms, 'markerFaceColor', cbColors(3,:), 'markerEdgeColor', 'k', 'lineWidth', 1)




ylim(est_err_lims)
xlim([0.5, 6.5])
%xlabel('SAC')

%title('Subject Pred Error','fontsize',14);
grid on
set(gca,'fontsize',ax_font_size)

set(gca,'XTickLabel', '');

setPLOT_panelLabel(gca, 3,0.12,1.05);

ylabel('Estimation error STD','fontsize',ax_font_size1)

% Now plot the CP, TAC and Noise contrasts
axes(axs(4)); cla(gca); hold on


                                       

aux_contrast = [estErr_cp_contrast estErr_sac_contrast_lo estErr_sac_contrast_hi estErr_noise_contrast];

aux_model_contrast = [theory_estErr_subjective_cp_contrast ...
                    theory_estErr_subjective_sac_contrast_low ,theory_estErr_subjective_sac_contrast_hi ....
                    theory_estErr_subjective_noise_contrast];


ll = size(aux_contrast,1);
plot([0 ll*6], [0 0], 'Color',[0.5 0.5 0.5],'linewidth',1)
clear P CI np_P

%which contrasts are significant at 5%?
sig_contrast_idx = [1 2 4];
contrast_edge_col_subj = [[0 0 0];cbColors(3,:);cbColors(2,:);[0 0 0]];
contrast_edge_col_theory = [theory_col2;cbColors(3,:);cbColors(2,:);theory_col2];
contrast_line_col = [theory_col3;theory_col3; theory_col3;theory_col3];

for i=1:size(aux_contrast,2)
    xJit = smartJitter(aux_contrast(:,i)./std(aux_contrast(:,i)),.15,.25);

    xJit_model = smartJitter(aux_model_contrast(:,i)./std(aux_contrast(:,i)),.15,.25);
    
    
    for i_line = 1:length(uniqueID)
        line([i*6+xJit(i_line), i*6 + 3 + xJit_model(i_line)], ...
             [aux_contrast(i_line,i)./std(aux_contrast(:,i)), aux_model_contrast(i_line,i)./std(aux_contrast(:,i))],...
              'color',contrast_line_col(i,:),'linewidth',0.25);
    end
    
    
    plot(ones(ll, 1).*i.*6+xJit_model + 3 , aux_model_contrast(:,i)./std(aux_contrast(:,i)),...
        'db', 'markerSize', 5, 'markerFaceColor','None',...
        'markerEdgeColor', contrast_edge_col_theory(i,:), 'lineWidth', 0.5);
    
    
    plot(ones(ll, 1).*i.*6+xJit , aux_contrast(:,i)./std(aux_contrast(:,i)),...
        'o', 'markerSize', 5, 'markerFaceColor',[0.96 0.96 0.96],...
        'markerEdgeColor', contrast_edge_col_subj(i,:), 'lineWidth', 1);
    
     xlim([4 6.*(i)+ 4.5])
     
     %y_max = max(aux_y) + 2.0;
    
    if(any(i==sig_contrast_idx))
        plot(i*6,5.2,'*','markersize',5,'color','k');
    end
end
%set(gca, 'xtick', [2:2:(ll)*2], 'XTickLabel', xTickLabels);
set(gca, 'xtick', [7.5:6:(ll)*6], 'XTickLabel', '');

set(gca,'fontsize',ax_font_size)
ylabel('Contrast','fontsize',ax_font_size1)
grid on
ylim([-6 5.5])

setPLOT_panelLabel(gca, 4,0.18,1.05);






%%  --------------------- Row 3 : Betting & contrasts -------------

axes(axs(5)); cla(gca); hold on

auxMu_hi = nanmean(norm_bet_prop_hi');
auxSE_hi = nanstd(norm_bet_prop_hi')./sqrt(size(norm_bet_prop_hi',1));

auxMu_low = nanmean(norm_bet_prop_low');
auxSE_low = nanstd(norm_bet_prop_low')./sqrt(size(norm_bet_prop_low',1));

barInds=[1:6; 1:6];
errBars_hi = [auxMu_hi + auxSE_hi ; auxMu_hi - auxSE_hi];
errBars_low = [auxMu_low + auxSE_low ; auxMu_low - auxSE_low];



auxMu_theory_hi = nanmean(theory_bet_frequency_hi');
auxSE_theory_hi = nanstd(theory_bet_frequency_hi')./sqrt(size(theory_bet_frequency_hi',1));

auxMu_theory_low = nanmean(theory_bet_frequency_low');
auxSE_theory_low = nanstd(theory_bet_frequency_low')./sqrt(size(theory_bet_frequency_low',1));

barInds=[1:6; 1:6];
errBars_theory_hi = [auxMu_theory_hi + auxSE_theory_hi ; auxMu_theory_hi - auxSE_theory_hi];
errBars_theory_low = [auxMu_theory_low + auxSE_theory_low ; auxMu_theory_low - auxSE_theory_low];




hold on

%plot(barInds,errBars_theory_hi, '-k', 'lineWidth', 1.5)
%plot(barInds,errBars_theory_low, '-k', 'lineWidth', 1.5)
plot(auxMu_theory_hi, 'db', 'markerSize', ms, 'markerFaceColor', 'None', 'markerEdgeColor',cbColors(2,:) , 'lineWidth', 1)
plot(auxMu_theory_low, 'db', 'markerSize', ms,'markerFaceColor', 'None', 'markerEdgeColor',cbColors(3,:) , 'lineWidth', 1)


plot(barInds,errBars_hi, '-k', 'lineWidth', 1)
plot(barInds,errBars_low, '-k', 'lineWidth', 1)
plot(auxMu_hi, 'ob', 'markerSize', ms, 'markerFaceColor', cbColors(2,:), 'markerEdgeColor', 'k', 'lineWidth', 1)
plot(auxMu_low, 'ob', 'markerSize', ms, 'markerFaceColor', cbColors(3,:), 'markerEdgeColor', 'k', 'lineWidth', 1)


xlim([0.5, 6.5])
ylim([-0.25 0.12]);


set(gca,'fontsize',ax_font_size)
xlabel('Sounds after a change-point','fontsize',ax_font_size1)
ylabel({'Normalized bet','frequency'},'fontsize',ax_font_size1)
grid on

set(gca, 'YTick', [-0.2 0 0.1]);

setPLOT_panelLabel(gca, 5,0.12,1.05);


% Now plot the CP, TAC and Noise contrasts
axes(axs(6)); cla(gca); hold on



aux_contrast = [bet_cp_contrast bet_sac_contrast_lo bet_sac_contrast_hi bet_noise_contrast];
aux_model_contrast = [theory_bet_cp_contrast theory_bet_sac_contrast_low(~sel_bet_subj) ...
                 theory_bet_sac_contrast_hi(~sel_bet_subj) theory_bet_noise_contrast];


% Plot this using Matt's smart jitter function

xTickLabels={{'CP'},{'SAC', 'low'},{'SAC', 'high'},{'Noise'}}; % 'abs(90-Outcome) 10' ,'abs(90-Outcome) 20'
ll = size(aux_contrast,1);
plot([0 ll*6], [0 0], 'Color',[0.5 0.5 0.5],'linewidth',1)
clear P CI np_P

%which contrasts are significant at 5%?
sig_contrast_idx = [1 2 4];
contrast_edge_col_subj = [[0 0 0];cbColors(3,:);cbColors(2,:);[0 0 0]];
contrast_edge_col_theory = [theory_col2;cbColors(3,:);cbColors(2,:);theory_col2];
contrast_line_col = [theory_col3;theory_col3; theory_col3;theory_col3];

for i=1:size(aux_contrast,2)
    xJit = smartJitter(aux_contrast(:,i)./std(aux_contrast(:,i)),.15,.25);

    xJit_model = smartJitter(aux_model_contrast(:,i)./std(aux_contrast(:,i)),.15,.25);
    
    
    for i_line = 1:length(uniqueID(~sel_bet_subj))
        line([i*6+xJit(i_line), i*6 + 3 + xJit_model(i_line)], ...
             [aux_contrast(i_line,i)./std(aux_contrast(:,i)), aux_model_contrast(i_line,i)./std(aux_contrast(:,i))],...
              'color',contrast_line_col(i,:),'linewidth',0.25);
    end
    
    
    plot(ones(ll, 1).*i.*6+xJit_model + 3 , aux_model_contrast(:,i)./std(aux_contrast(:,i)),...
        'db', 'markerSize', 5, 'markerFaceColor','None',...
        'markerEdgeColor', contrast_edge_col_theory(i,:), 'lineWidth', 0.5);
    
    
    plot(ones(ll, 1).*i.*6+xJit , aux_contrast(:,i)./std(aux_contrast(:,i)),...
        'o', 'markerSize', 5, 'markerFaceColor',[0.96 0.96 0.96],...
        'markerEdgeColor', contrast_edge_col_subj(i,:), 'lineWidth', 1);
    
     xlim([4 6.*(i)+ 4.5])
     
     %y_max = max(aux_y) + 2.0;
    
    if(any(i==sig_contrast_idx))
        plot(i*6,3.8,'*','markersize',5,'color','k');
    end
end
%set(gca, 'xtick', [2:2:(ll)*2], 'XTickLabel', xTickLabels);
set(gca, 'xtick', [7.5:6:(ll)*6], 'XTickLabel', '');
yTicks = get(gca,'ytick');
xTicks = get(gca, 'xtick');
ax = axis;
vOffset = [0.5 1.0 1.0 0.5];
hOffset = [2.0 1.7 1.5 1.2];

minY = min(yTicks);
for xx = 1:length(xTickLabels)
%Create text box and set appropriate properties
     t = text(xTicks(xx) - hOffset(xx), minY - vOffset(xx),xTickLabels{xx},'fontsize',ax_font_size);
         %'HorizontalAlignment','Right','interpreter', 'latex');
         
end

set(gca,'fontsize',ax_font_size)
ylabel('Contrast','fontsize',ax_font_size1)
grid on

ylim([-4 4])

setPLOT_panelLabel(gca, 6,0.18,1.05);

%Fix the tick labels
%fix_xticklabels;

%%
kk=annotation('textbox');
set(kk, 'string', 'Krishnamurthy et al, figure 3', 'position', ...
    [0.78 0.9 0.15 0.05], 'EdgeColor', 'none','fontsize',14)




figName = '/Users/kamesh/Dropbox/auditoryTaskManuscript/Figure_set_for_writeup/figure_3.pdf';


%%      SPIT OUT THE STATS FOR THE CONTRASTS:

% % SUBJECT PRED ERROR 
% 
% [h_subj_predErr_cp_combined, p_subj_predErr_cp_combined,~,stats_subj_predErr_cp_combined] ...
%         = ttest(predErr_cp_contrast);
% [h_subj_predErr_cp_hi, p_subj_predErr_cp_hi,~,stats_subj_predErr_cp_hi] = ...
%     ttest(predErr_cp_contrast_hi);
% [h_subj_predErr_cp_low, p_subj_predErr_cp_low,~,stats_subj_predErr_cp_low] = ...
%     ttest(predErr_cp_contrast_lo);
% 
% 
% [h_subj_predErr_sac_combined, p_subj_predErr_sac_combined,~,stats_subj_predErr_sac_combined] = ...
%       ttest(predErr_sac_contrast);
% [h_subj_predErr_sac_hi, p_subj_predErr_sac_hi,~,stats_subj_predErr_sac_hi] = ...
%      ttest(predErr_sac_contrast_hi);
% [h_subj_predErr_sac_low, p_subj_predErr_sac_low,~,stats_subj_predErr_sac_low] = ...
%     ttest(predErr_sac_contrast_lo);
% 
% [h_subj_predErr_noise, p_subj_predErr_noise,~,stats_subj_predErr_noise] = ...
%     ttest(predErr_noise_contrast);
% 
% 
% 
% % SUBJECT ESTIMATION ERROR
% 
% [h_subj_estErr_cp_combined, p_subj_estErr_cp_combined, ~, stats_subj_estErr_cp_combined] = ...
%         ttest(estErr_cp_contrast);
% [h_subj_estErr_cp_hi, p_subj_estErr_cp_hi, ~, stats_subj_estErr_cp_hi] = ttest(estErr_cp_contrast_hi);
% [h_subj_estErr_cp_low, p_subj_estErr_cp_low, ~, stats_subj_estErr_cp_low] = ttest(estErr_cp_contrast_lo);
% 
% 
% [h_subj_estErr_sac_combined, p_subj_estErr_sac_combined, ~, stats_subj_estErr_sac_combined] = ...
%     ttest(estErr_sac_contrast);
% [h_subj_estErr_sac_hi, p_subj_estErr_sac_hi, ~, stats_subj_estErr_sac_hi] = ...
%                                                 ttest(estErr_sac_contrast_hi);
% [h_subj_estErr_sac_low, p_subj_estErr_sac_low, ~, stats_subj_estErr_sac_low] = ...
%                                                         ttest(estErr_sac_contrast_lo);
%                                                     
%                                                     
%  [h_subj_estErr_sac_diff, p_subj_estErr_sac_diff, ~, stats_subj_estErr_sac_diff] = ...
%                                          ttest(estErr_sac_contrast_hi - estErr_sac_contrast_lo);              
% 
% [h_subj_estErr_noise, p_subj_estErr_noise, ~, stats_subj_estErr_noise] = ...
%     ttest(estErr_noise_contrast);
% 
% 
% % SUBJECT BETTING
% 
% [h_subj_bet_cp_combined, p_subj_bet_cp_combined, ~, stats_subj_bet_cp_combined] = ...
%                     ttest(bet_cp_contrast);
% [h_subj_bet_cp_hi, p_subj_bet_cp_hi, ~, stats_subj_bet_cp_hi] = ttest(bet_cp_contrast_hi);
% [h_subj_bet_cp_low, p_subj_bet_cp_low, ~, stats_subj_bet_cp_low] = ttest(bet_cp_contrast_lo);
% 
% 
% [h_subj_bet_sac_combined, p_subj_bet_sac_combined,  ~, stats_subj_bet_sac_combined] = ...
%                         ttest(bet_sac_contrast);
% [h_subj_bet_sac_hi, p_subj_bet_sac_hi,  ~, stats_subj_bet_sac_hi] = ttest(bet_sac_contrast_hi);
% [h_subj_bet_sac_low, p_subj_bet_sac_low,  ~, stats_subj_bet_sac_low] = ttest(bet_sac_contrast_lo);
% 
%  [h_subj_bet_sac_diff, p_subj_bet_sac_diff, ~, stats_subj_bet_sac_diff] = ...
%                                          ttest(bet_sac_contrast_lo - bet_sac_contrast_hi);
% 
% [h_subj_bet_noise, p_subj_bet_noise,~, stats_subj_bet_noise] = ttest(bet_noise_contrast);
% 
% % THEORETICAL ESTIMATION ERROR
% 
% [h_theory_estErr_cp_combined, p_theory_estErr_cp_combined] = ttest(theory_estErr_subjective_cp_contrast);
% [h_theory_estErr_cp_hi, p_theory_estErr_cp_hi] = ttest(theory_estErr_subjective_cp_contrast_hi);
% [h_theory_estErr_cp_low, p_theory_estErr_cp_low] = ttest(theory_estErr_subjective_cp_contrast_low);
% 
% 
% [h_theory_estErr_sac_combined, p_theory_estErr_sac_combined] = ttest(theory_estErr_subjective_sac_contrast);
% [h_theory_estErr_sac_hi, p_theory_estErr_sac_hi] = ttest(theory_estErr_subjective_sac_contrast_hi);
% [h_theory_estErr_sac_low, p_theory_estErr_sac_low] = ttest(theory_estErr_subjective_sac_contrast_low);
% 
% [h_theory_estErr_noise, p_theory_estErr_noise] = ttest(theory_estErr_subjective_noise_contrast);
% 
% 
% 
% %THEORY BETTING 
% 
% [h_theory_bet_cp_combined, p_theory_bet_cp_combined] = ttest(theory_bet_cp_contrast);
% [h_theory_bet_sac_combined, p_theory_bet_sac_combined] = ttest(theory_bet_sac_contrast);
% [h_theory_bet_noise, p_theory_bet_noise] = ttest(theory_bet_noise_contrast);


% SUBJECT PRED ERROR 

[p_subj_predErr_cp_combined, h_subj_predErr_cp_combined,stats_subj_predErr_cp_combined] ...
        = signrank(predErr_cp_contrast);
[p_subj_predErr_cp_hi, h_subj_predErr_cp_hi,stats_subj_predErr_cp_hi] = ...
    signrank(predErr_cp_contrast_hi);
[p_subj_predErr_cp_low, h_subj_predErr_cp_low,stats_subj_predErr_cp_low] = ...
    signrank(predErr_cp_contrast_lo);


[p_subj_predErr_sac_combined, h_subj_predErr_sac_combined,stats_subj_predErr_sac_combined] = ...
      signrank(predErr_sac_contrast);
[p_subj_predErr_sac_hi, h_subj_predErr_sac_hi,stats_subj_predErr_sac_hi] = ...
     signrank(predErr_sac_contrast_hi);
[p_subj_predErr_sac_low, h_subj_predErr_sac_low,stats_subj_predErr_sac_low] = ...
    signrank(predErr_sac_contrast_lo);

[p_subj_predErr_noise, h_subj_predErr_noise,stats_subj_predErr_noise] = ...
    signrank(predErr_noise_contrast);



% SUBJECT ESTIMATION ERROR

[p_subj_estErr_cp_combined, h_subj_estErr_cp_combined, stats_subj_estErr_cp_combined] = ...
        signrank(estErr_cp_contrast);
[p_subj_estErr_cp_hi, h_subj_estErr_cp_hi, stats_subj_estErr_cp_hi] = signrank(estErr_cp_contrast_hi);
[p_subj_estErr_cp_low, h_subj_estErr_cp_low, stats_subj_estErr_cp_low] = signrank(estErr_cp_contrast_lo);


[p_subj_estErr_sac_combined, h_subj_estErr_sac_combined, stats_subj_estErr_sac_combined] = ...
    signrank(estErr_sac_contrast);
[p_subj_estErr_sac_hi, h_subj_estErr_sac_hi, stats_subj_estErr_sac_hi] = ...
                                                signrank(estErr_sac_contrast_hi,0,'tail','left');
[p_subj_estErr_sac_low, h_subj_estErr_sac_low, stats_subj_estErr_sac_low] = ...
                                                        signrank(estErr_sac_contrast_lo,0,'tail','left');
                                                    
                                                    
 [p_subj_estErr_sac_diff, h_subj_estErr_sac_diff,stats_subj_estErr_sac_diff] = ...
                                         signrank(estErr_sac_contrast_hi - estErr_sac_contrast_lo);              

[p_subj_estErr_noise, h_subj_estErr_noise, stats_subj_estErr_noise] = ...
    signrank(estErr_noise_contrast);


% SUBJECT BETTING

[p_subj_bet_cp_combined, h_subj_bet_cp_combined, stats_subj_bet_cp_combined] = ...
                    signrank(bet_cp_contrast);
[p_subj_bet_cp_hi, h_subj_bet_cp_hi, stats_subj_bet_cp_hi] = signrank(bet_cp_contrast_hi);
[p_subj_bet_cp_low, h_subj_bet_cp_low, stats_subj_bet_cp_low] = signrank(bet_cp_contrast_lo);


[p_subj_bet_sac_combined, h_subj_bet_sac_combined,  stats_subj_bet_sac_combined] = ...
                        signrank(bet_sac_contrast);
[p_subj_bet_sac_hi, h_subj_bet_sac_hi,  stats_subj_bet_sac_hi] = signrank(bet_sac_contrast_hi,0,'tail','right');
[p_subj_bet_sac_low, h_subj_bet_sac_low,  stats_subj_bet_sac_low] = signrank(bet_sac_contrast_lo,0,'tail','right');

 [p_subj_bet_sac_diff, h_subj_bet_sac_diff, stats_subj_bet_sac_diff] = ...
                                         signrank(bet_sac_contrast_lo - bet_sac_contrast_hi);

[p_subj_bet_noise, h_subj_bet_noise, stats_subj_bet_noise] = signrank(bet_noise_contrast);

% THEORETICAL ESTIMATION ERROR

[p_theory_estErr_cp_combined, h_theory_estErr_cp_combined] = signrank(theory_estErr_subjective_cp_contrast);
[p_theory_estErr_cp_hi, h_theory_estErr_cp_hi] = signrank(theory_estErr_subjective_cp_contrast_hi);
[p_theory_estErr_cp_low, h_theory_estErr_cp_low] = signrank(theory_estErr_subjective_cp_contrast_low);


[p_theory_estErr_sac_combined, h_theory_estErr_sac_combined] = signrank(theory_estErr_subjective_sac_contrast);
[p_theory_estErr_sac_hi, h_theory_estErr_sac_hi] = signrank(theory_estErr_subjective_sac_contrast_hi);
[p_theory_estErr_sac_low, h_theory_estErr_sac_low] = signrank(theory_estErr_subjective_sac_contrast_low);

[p_theory_estErr_noise, h_theory_estErr_noise] = signrank(theory_estErr_subjective_noise_contrast);



%THEORY BETTING 

[p_theory_bet_cp_combined, h_theory_bet_cp_combined] = signrank(theory_bet_cp_contrast);
[p_theory_bet_sac_combined, h_theory_bet_sac_combined] = signrank(theory_bet_sac_contrast);
[p_theory_bet_noise, h_theory_bet_noise] = signrank(theory_bet_noise_contrast);


clc

fprintf(' ---------- SUBJECTS : PREDICTION ERROR CONTRASTS -------------\n\n')


fprintf('Subject predErr CP contrast COMBINED pVal ---> %f  signedrank---> %f\n',...
          p_subj_predErr_cp_combined, stats_subj_predErr_cp_combined.signedrank)
      
fprintf('Subject predErr SAC contrast COMBINED pVal ---> %f  signedrank---> %f\n\n',...
          p_subj_predErr_sac_combined, stats_subj_predErr_sac_combined.signedrank)      
      

fprintf('pVal Subject predErr CP contrast LOW NOISE ---> %f\n',p_subj_predErr_cp_low)
fprintf('pVal Subject predErr SAC contrast LOW NOISE ---> %f\n\n',p_subj_predErr_sac_low)

fprintf('pVal Subject predErr CP contrast HIGH NOISE ---> %f\n',p_subj_predErr_cp_hi)
fprintf('pVal Subject predErr SAC contrast HIGH NOISE ---> %f\n\n',p_subj_predErr_sac_hi)


fprintf('Subject predErr NOISE contrast pVal ---> %f  signedrank---> %f\n\n\n\n',...
          p_subj_predErr_noise, stats_subj_predErr_noise.signedrank)


fprintf(' ---------- SUBJECTS : ESTIMATION ERROR CONTRASTS -------------\n\n')

fprintf('Subject estErr CP contrast COMBINED pVal ---> %f  signedrank---> %f\n',...
          p_subj_estErr_cp_combined, stats_subj_estErr_cp_combined.signedrank)
      
fprintf('Subject estErr SAC contrast COMBINED pVal ---> %f  signedrank---> %f\n\n',...
          p_subj_estErr_sac_combined, stats_subj_estErr_sac_combined.signedrank) 


fprintf('Subject estErr CP contrast LOW pVal ---> %f  signedrank---> %f\n',...
          p_subj_estErr_cp_low, stats_subj_estErr_cp_low.signedrank)
      
fprintf('Subject estErr SAC contrast LOW pVal ---> %f  signedrank---> %f\n\n',...
          p_subj_estErr_sac_low, stats_subj_estErr_sac_low.signedrank) 

fprintf('Subject estErr CP contrast HIGH pVal ---> %f  signedrank---> %f\n',...
          p_subj_estErr_cp_hi, stats_subj_estErr_cp_hi.signedrank)
      
fprintf('Subject estErr SAC contrast HIGH pVal ---> %f  signedrank---> %f\n',...
          p_subj_estErr_sac_hi, stats_subj_estErr_sac_hi.signedrank) 
fprintf('Subject estErr SAC contrast DIFF (LOW - HIGH)  pVal ---> %f  signedrank---> %f\n\n',...
          p_subj_estErr_sac_diff, stats_subj_estErr_sac_diff.signedrank)

fprintf('Subject estErr NOISE contrast pVal ---> %f  signedrank---> %f\n\n\n\n',...
          p_subj_estErr_noise, stats_subj_estErr_noise.signedrank)
      
      



fprintf(' ---------- SUBJECTS :BET CONTRASTS -------------\n\n')


fprintf('Subject Bet CP contrast COMBINED pVal ---> %f  signedrank---> %f\n',...
          p_subj_bet_cp_combined, stats_subj_bet_cp_combined.signedrank)   
fprintf('Subject Bet SAC contrast COMBINED pVal ---> %f  signedrank---> %f\n\n',...
          p_subj_bet_sac_combined, stats_subj_bet_sac_combined.signedrank)   

fprintf('Subject Bet CP contrast LOW pVal ---> %f  signedrank---> %f\n',...
          p_subj_bet_cp_low, stats_subj_bet_cp_low.signedrank)   
fprintf('Subject Bet SAC contrast LOW pVal ---> %f  signedrank---> %f\n',...
          p_subj_bet_sac_low, stats_subj_bet_sac_low.signedrank)
fprintf('Subject Bet SAC contrast DIFF (LOW - HIGH) pVal ---> %f  signedrank---> %f\n\n',...
          p_subj_bet_sac_diff, stats_subj_bet_sac_diff.signedrank)

fprintf('Subject Bet CP contrast HIGH pVal ---> %f  signedrank---> %f\n',...
          p_subj_bet_cp_hi, stats_subj_bet_cp_hi.signedrank)   
fprintf('Subject Bet SAC contrast HIGH pVal ---> %f  signedrank---> %f\n\n',...
          p_subj_bet_sac_hi, stats_subj_bet_sac_hi.signedrank)

fprintf('Subject Bet Noise contrast  pVal ---> %f  signedrank---> %f\n\n\n\n',...
          p_subj_bet_noise, stats_subj_bet_noise.signedrank)   

% --- theory estimation error contrasts

fprintf(' ---------- THEORY : ESTIMATION ERROR CONTRASTS -------------\n\n')

fprintf('pVal Theory Bet CP contrast COMBINED ---> %f\n\n',p_theory_bet_cp_combined)
fprintf('pVal Theory Bet SAC contrast COMBINED ---> %f\n\n',p_theory_bet_sac_combined)

fprintf('pVal Theory Bet Noise contrast ---> %f\n\n\n\n',p_theory_bet_noise)




fprintf(' ---------- SUBJECTS versus MODEL contrasts -------------\n\n')


% tests: ranksum

% Prediction error
% 1) overall difference 
% 2) CP contrast
% 3) EXP contrast
% 4) Noise contrast

% Prediction error

                    %1 -- MEAN
                    %2 -- CP
                    %3 -- EXP
                    %4 -- Noise
                    
 [P_pred_cp,H,STATS]=ranksum(predErr_cp_contrast, model_predErr_cp_contrast)   %2
 [P_pred_sac_lo,H,STATS]=ranksum(predErr_sac_contrast_lo, model_predErr_sac_contrast_low) %3
 [P_pred_sac_hi,H,STATS]=ranksum(predErr_sac_contrast_hi, model_predErr_sac_contrast_hi) %3
 [P_pred_noise,H,STATS]=ranksum(predErr_noise_contrast, model_predErr_noise_contrast) %4

% Estimation error
                                                                      %1
 [P_est_cp,H,STATS]=ranksum(estErr_cp_contrast, theory_estErr_subjective_cp_contrast)   %2
 [P_est_sac_lo,H,STATS]=ranksum(estErr_sac_contrast_lo, theory_estErr_subjective_sac_contrast_low) %3
 [P_est_sac_hi,H,STATS]=ranksum(estErr_sac_contrast_hi, theory_estErr_subjective_sac_contrast_hi) %3
 [P_est_noise,H,STATS]=ranksum(estErr_noise_contrast, theory_estErr_subjective_noise_contrast) %4

% Bet
                                                                      %1
 [P_bet_cp,H,STATS]=ranksum(bet_cp_contrast,  theory_bet_cp_contrast)   %2
  [P_bet_sac_lo,H,STATS]=ranksum(bet_sac_contrast_lo,  theory_bet_sac_contrast_low) %3
 [P_bet_sac_hi,H,STATS]=ranksum(bet_sac_contrast_hi,  theory_bet_sac_contrast_hi) %3
 [P_bet_noise,H,STATS]=ranksum(bet_noise_contrast,  theory_bet_noise_contrast) %4



% [p_theory_estErr_cp_combined, h_theory_estErr_cp_combined] = signrank(theory_estErr_subjective_cp_contrast);
% [p_theory_estErr_cp_hi, h_theory_estErr_cp_hi] = signrank(theory_estErr_subjective_cp_contrast_hi);
% [p_theory_estErr_cp_low, h_theory_estErr_cp_low] = signrank(theory_estErr_subjective_cp_contrast_low);
% 
% 
% [p_theory_estErr_sac_combined, h_theory_estErr_sac_combined] = signrank(theory_estErr_subjective_sac_contrast);
% [p_theory_estErr_sac_hi, h_theory_estErr_sac_hi] = signrank(theory_estErr_subjective_sac_contrast_hi);
% [p_theory_estErr_sac_low, h_theory_estErr_sac_low] = signrank(theory_estErr_subjective_sac_contrast_low);
% 
% [p_theory_estErr_noise, h_theory_estErr_noise] = signrank(theory_estErr_subjective_noise_contrast);



% %THEORY BETTING 
% 
% [p_theory_bet_cp_combined, h_theory_bet_cp_combined] = signrank(theory_bet_cp_contrast);
% [p_theory_bet_sac_combined, h_theory_bet_sac_combined] = signrank(theory_bet_sac_contrast);
% [p_theory_bet_noise, h_theory_bet_noise] = signrank(theory_bet_noise_contrast);

