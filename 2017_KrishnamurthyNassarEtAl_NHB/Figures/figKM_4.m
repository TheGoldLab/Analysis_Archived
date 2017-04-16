%% ************************************************************************************
%
%                      FIGURE 4 : PRIOR WEIGHT FIGURE
%
% *************************************************************************************

getCbColors;

% Use the J. Neurosci figure template
fig_num = 1; % figure number
total_width = 17; % total width
col_height = [7 8]; % height of each column
cols_per_row = {[0.33, 0.33, 0.33], [.65, .4]}; % columns in each row
psh = 1.2; % panel separation height
psw = 2; % panel separation width

[axs,fig_] = getPLOT_axes(fig_num, total_width, col_height, cols_per_row,...
                                 psh, psw, [], '');
set(axs,'Units','normalized');

ax_font_size = 14;
ms = 9; % marker size
jit = 0.05; % jitter for smartJitter
model_contr_marker_col = [0.8 0.1 0];

%%    -----------  Schematic panel 1 ----------------------------
axes(axs(1)); cla(gca); hold on
axis square

subj1_color = [0 0 0];

axis_lims = [-80 80];
example_subjID = 22;

selSubj = probeSubjID == example_subjID & probebCPTAC >= 4 & probeStd == 20;
aux_predErr = probePrediction(selSubj) - probeOutcomes(selSubj);
aux_estErr = probeEstimate(selSubj) - probeOutcomes(selSubj);



ylim(axis_lims)
xlim(axis_lims)

l = lsline;
l.Color = subj1_color;
l.LineWidth = 2;


plot(axis_lims,[0 0],'linestyle','-','color',[0.5 0.5 0.5],'linewidth',1)
plot([0 0],axis_lims,'linestyle','-','color',[0.5 0.5 0.5],'linewidth',1)
plot(axis_lims,axis_lims,'linestyle','-','color',[0.5 0.5 0.5],'linewidth',1)

plot(aux_predErr,aux_estErr,'o','markersize',4,'markeredgecolor',subj1_color,...
          'markerfacecolor','None','linewidth',0.25);


xlabel('Prediction error')
ylabel('Perceptual error')

set(gca,'fontsize',ax_font_size)

setPLOT_panelLabel(gca, 1,0.18,1.2);

title({'Measuring', 'Perceptual Bias'},'fontsize',12)


%%    -----------  Schematic panel 2 ----------------------------
axes(axs(2)); cla(gca); hold on
axis square

%subj1_color = [0.5 0.5 0.5];
subj2_color = [0.2 0.2 0.2];

axis_lims = [-90 90];
example_subjID1 = 23;
example_subjID2 = 29;


selSubj1 = probeSubjID == example_subjID1 & probebCPTAC == 1 & probeStd ~= 0;
aux_predErr1 = probePrediction(selSubj1) - probeOutcomes(selSubj1);
aux_estErr1 = probeEstimate(selSubj1) - probeOutcomes(selSubj1);



Xspace1 = (90 - aux_predErr1)';
X1 = [ones(length(aux_predErr1),1) aux_predErr1' Xspace1];
y1 = aux_estErr1';

[b_1 bint_1] = regress(y1,X1);
l1 = plot(axis_lims,b_1(2)*axis_lims,'linestyle','-','color',subj1_color,'linewidth',2);


% Xspace2 = (90 - aux_predErr2)';
% X2 = [ones(length(aux_predErr2),1) aux_predErr2' Xspace2];
% y2 = aux_estErr2';
% 
% [b_2 bint_2] = regress(y2,X2);
% l2 = plot(axis_lims,b_2(2)*axis_lims,'linestyle','-','color',subj2_color,'linewidth',3);




ylim(axis_lims)
xlim(axis_lims)

%l = lsline;
%l.Color = [0.4 0.4 0.4];
%l.LineWidth = 2;


plot(axis_lims,[0 0],'linestyle','-','color',[0.5 0.5 0.5],'linewidth',1)
plot([0 0],axis_lims,'linestyle','-','color',[0.5 0.5 0.5],'linewidth',1)
plot(axis_lims,axis_lims,'linestyle','-','color',[0.5 0.5 0.5],'linewidth',1)



plot(aux_predErr1,aux_estErr1,'o','markersize',4,'markeredgecolor',subj1_color, ...
           'markerfacecolor','None','linewidth',0.25);

xlabel('Prediction error')

set(gca,'fontsize',ax_font_size)


%legend([l1],'S1','S2','fontsize',12,'location','southeast')
setPLOT_panelLabel(gca, 2,0.18,1.2);

title('Change-point trials','fontsize',12)

%%    -----------  Schematic panel 3 ----------------------------
axes(axs(3)); cla(gca); hold on
axis square

%subj1_color = [0.5 0.5 0.5];
subj2_color = [0.2 0.2 0.2];

axis_lims = [-90 90];



selSubj1 = probeSubjID == example_subjID1 & probebCPTAC >= 2 & probeStd ~= 0;
aux_predErr1 = probePrediction(selSubj1) - probeOutcomes(selSubj1);
aux_estErr1 = probeEstimate(selSubj1) - probeOutcomes(selSubj1);

% selSubj2 = probeSubjID == example_subjID2 & probebCPTAC >= 5 & probeStd ~= 0;
% aux_predErr2 = probePrediction(selSubj2) - probeOutcomes(selSubj2);
% aux_estErr2 = probeEstimate(selSubj2) - probeOutcomes(selSubj2);





% plot(aux_predErr2,aux_estErr2,'x','markersize',6,'markeredgecolor',[0.3 0.3 0.3], ...
%            'markerfacecolor',[0.8 0.8 0.8],'linewidth',2);


       


Xspace1 = (90 - aux_predErr1)';
X1 = [ones(length(aux_predErr1),1) aux_predErr1' Xspace1];
y1 = aux_estErr1';

[b_1 bint_1] = regress(y1,X1);
l1 = plot(axis_lims,b_1(2)*axis_lims,'linestyle','-','color',subj1_color,'linewidth',2);


% Xspace2 = (90 - aux_predErr2)';
% X2 = [ones(length(aux_predErr2),1) aux_predErr2' Xspace2];
% y2 = aux_estErr2';
% 
% [b_2 bint_2] = regress(y2,X2);
% l2 = plot(axis_lims,b_2(2)*axis_lims,'linestyle','-','color',subj2_color,'linewidth',3);
% 



ylim(axis_lims)
xlim(axis_lims)

%l = lsline;
%l.Color = [0.4 0.4 0.4];
%l.LineWidth = 2;


plot(axis_lims,[0 0],'linestyle','-','color',[0.5 0.5 0.5],'linewidth',1)
plot([0 0],axis_lims,'linestyle','-','color',[0.5 0.5 0.5],'linewidth',1)
plot(axis_lims,axis_lims,'linestyle','-','color',[0.5 0.5 0.5],'linewidth',1)


plot(aux_predErr1,aux_estErr1,'o','markersize',4,'markeredgecolor',subj1_color, ...
           'markerfacecolor','None','linewidth',0.25);

%legend([l1],'S1','S2','fontsize',12,'location','southeast')
xlabel('Prediction error')


set(gca,'fontsize',ax_font_size)
setPLOT_panelLabel(gca, 3,0.18,1.2);

title('Nonchange-point trials','fontsize',12)


%%    -----------  SUBJECT AND MODEL PRIOR WEIGHTS WITH CONTRASTS

axes(axs(4)); cla(gca); hold on
pltIdx = 1:6; %indices corresponding to SAC 1-6


% NOW PLOT THE SUBJECT PWt vs. SAC


count = 1;
pointHandles = [];
delta = [0 0.1];    %spacing between points for different stds

% plot(repmat((1:numTAC)+delta(count), 2, 1)-jit,...
%     [theoryBias_subjective_avg_low(pltIdx) - theoryBias_subjective_Int_low(pltIdx) ...
%     theoryBias_subjective_avg_low(pltIdx) + theoryBias_subjective_Int_low(pltIdx)]',...
%     cbColors(2,:),'linewidth',1);
pointHandles(end+1) = plot((1:numTAC)+delta(count)-jit,...
    theoryBias_subjective_avg_low(pltIdx), 'dk', 'markerEdgeColor', cbColors(3,:), ...
    'markerFaceColor', 'None',...
    'lineWidth', 1, 'markerSize', ms);
% 
% plot(repmat((1:numTAC)+delta(count), 2, 1)-jit,...
%     [theoryBias_subjective_avg_hi(pltIdx) - theoryBias_subjective_Int_hi(pltIdx) ...
%     theoryBias_subjective_avg_hi(pltIdx) + theoryBias_subjective_Int_hi(pltIdx)]',...
%     cbColors(3,:),'linewidth',1);
pointHandles(end+1) = plot((1:numTAC)+delta(count)-jit,...
    theoryBias_subjective_avg_hi(pltIdx), 'dk', 'markerEdgeColor', cbColors(2,:), ...
    'markerFaceColor', 'None',...
    'lineWidth', 1, 'markerSize', ms);

plot(jit+repmat((1:numTAC)+delta(count), 2, 1),...
    [subjAvgBetas_lo(pltIdx)-subjBetaInt_lo(pltIdx) subjAvgBetas_lo(pltIdx)+subjBetaInt_lo(pltIdx)]',...
    'k','linewidth',1.5);
h1 =  plot(jit+(1:numTAC)+delta(count),...
    subjAvgBetas_lo(pltIdx), 'ok', 'markerEdgeColor', 'k', ...
    'markerFaceColor', cbColors(3,:),...
    'lineWidth', 1, 'markerSize', ms);

plot(repmat((1:numTAC)+delta(count), 2, 1)-jit,...
    [subjAvgBetas_hi(pltIdx)-subjBetaInt_hi(pltIdx) subjAvgBetas_hi(pltIdx)+subjBetaInt_hi(pltIdx)]',...
    'k','linewidth',1.5);
h2 = plot((1:numTAC)+delta(count)-jit,...
    subjAvgBetas_hi(pltIdx), 'ok', 'markerEdgeColor', 'k', ...
    'markerFaceColor', cbColors(2,:),...
    'lineWidth', 1, 'markerSize', ms);

legend([h2,h1],{'high noise','low noise'},'Location','best','Box','off')


count = count + 1;

%legend(pointHandles,'LOW', 'HIGH','location','NorthWest');

grid on
ylabel('Perceptual Bias','fontsize',ax_font_size);
xlabel('Sounds after a change-point');
%title('SUBJECT PRIOR WEIGHT','fontsize',ax_font_size)
xlim([0.5 6.5])
set(gca, 'fontsize',ax_font_size)
set(gca,'YTick',[0.1 0.2 0.3 0.4 0.5 0.6 0.7])


setPLOT_panelLabel(gca, 4,0.12,1.1);




% PLOT THE CONTRASTS
% Now plot the CP, TAC and Noise contrasts
axes(axs(5)); cla(gca); hold on



aux_contrast = [pw_cp_contrast_combined pw_sac_contrast_lo pw_sac_contrast_hi pw_noise_contrast];
aux_theory_subjective_contrast = [theory_subjective_pw_cp_contrast theory_subjective_pw_sac_contrast_low ...
                                    theory_subjective_pw_sac_contrast_hi theory_subjective_pw_noise_contrast];
%aux_theory_subjective_contrast = bsxfun(@times,aux_theory_subjective_contrast,...
%                                        1./std(aux_theory_subjective_contrast));

aux_model_contrast = aux_theory_subjective_contrast;
                                    


% 
% xTickLabels={'CP','SAC(low)','SAC(High)','Noise'}; % 'abs(90-Outcome) 10' ,'abs(90-Outcome) 20'
% ll = size(aux_contrast,1);
% plot([0 ll*6], [0 0], 'k','linewidth',2)
% clear P CI np_P
% 
% %which contrasts are significant at 5%?
% sig_contrast_idx = [1 3 4];
% contrast_edge_col_subj = [[0 0 0];cbColors(3,:);cbColors(2,:);[0 0 0]];
% contrast_edge_col_theory = [theory_col2;cbColors(3,:);cbColors(2,:);theory_col2];
% contrast_line_col = [theory_col3;theory_col3; theory_col3;theory_col3];
% 
% for i=1:size(aux_contrast,2)
%     xJit = smartJitter(aux_contrast(:,i)./std(aux_contrast(:,i)),.15,.25);
% 
%     xJit_model = smartJitter(aux_theory_subjective_contrast(:,i)./std(aux_contrast(:,i)),.15,.25);
%     
% % %   
% %     aux_x = zeros(58,1);
% %     aux_x(1:2:end) = (ones(ll, 1).*);
% %     aux_x(2:2:end) = (ones(ll, 1).*i.*6+xJit_model + 3);
% %     aux_y = zeros(58,1);
% %     aux_y(1:2:end) = (aux_contrast(:,i)./std(aux_contrast(:,i)));
% %     aux_y(2:2:end) = (aux_model_contrast(:,i)./std(aux_contrast(:,i)));
% %     
%     for i_line = 1:length(uniqueID)
%         line([i*6+xJit(i_line), i*6 + 3 + xJit_model(i_line)], ...
%              [aux_contrast(i_line,i)./std(aux_contrast(:,i)), aux_model_contrast(i_line,i)./std(aux_contrast(:,i))],...
%               'color',contrast_line_col(i,:),'linewidth',0.25);
%     end
%     
%     
%     plot(ones(ll, 1).*i.*6+xJit_model + 3 , aux_model_contrast(:,i)./std(aux_contrast(:,i)),...
%         'db', 'markerSize', 5, 'markerFaceColor','None',...
%         'markerEdgeColor', contrast_edge_col_theory(i,:), 'lineWidth', 0.5);
%     
%     
%     plot(ones(ll, 1).*i.*6+xJit , aux_contrast(:,i)./std(aux_contrast(:,i)),...
%         'o', 'markerSize', 5, 'markerFaceColor',[0.96 0.96 0.96],...
%         'markerEdgeColor', contrast_edge_col_subj(i,:), 'lineWidth', 1);
%     
%      xlim([4 6.*(i)+ 4.5])
%      
%      y_max = max(aux_y) + 2.0;
%     
%     if(any(i==sig_contrast_idx))
%         plot(i*6,3,'*','markersize',5,'color','k');
%     end
% end
% %set(gca, 'xtick', [2:2:(ll)*2], 'XTickLabel', xTickLabels);
% set(gca, 'xtick', [7.5:6:(ll)*6], 'XTickLabel', xTickLabels);
% 
% set(gca,'fontsize',ax_font_size)
% ylabel('Contrast','fontsize',ax_font_size1)
% grid on



xTickLabels={{'CP'},{'SAC', 'low'},{'SAC', 'high'},{'Noise'}}; % 'abs(90-Outcome) 10' ,'abs(90-Outcome) 20'
ll = size(aux_contrast,1);
plot([0 ll*6], [0 0], 'Color',[0.5 0.5 0.5],'linewidth',1)
clear P CI np_P

%which contrasts are significant at 5%?
sig_contrast_idx = [1 3 4];
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
         plot(i*6,3.5,'*','markersize',5,'color','k');
     end
end
%set(gca, 'xtick', [2:2:(ll)*2], 'XTickLabel', xTickLabels);
set(gca, 'xtick', [7.5:6:(ll)*6], 'XTickLabel', '');
yTicks = get(gca,'ytick');
xTicks = get(gca, 'xtick');
ax = axis;
vOffset = [0.4 0.8 0.8 0.4];
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

ylim([-5.5 4.5])

setPLOT_panelLabel(gca, 5,0.18,1.1);


%%

kk=annotation('textbox');
set(kk, 'string', 'Krishnamurthy et al, figure 4', 'position', ...
    [0.78 0.85 0.15 0.05], 'EdgeColor', 'none','fontsize',14)

figName = '/Users/kamesh/Dropbox/auditoryTaskManuscript/Figure_set_for_writeup/figure_4.pdf';


%%      SPIT OUT THE STATS FOR THE CONTRASTS:

% [h_subj_pw_cp_combined, p_subj_pw_cp_combined, ~, stats_subj_pw_cp_combined] = ttest(pw_cp_contrast_combined);
% [h_subj_pw_cp_hi, p_subj_pw_cp_hi, ~, stats_subj_pw_cp_hi] = ttest(pw_cp_contrast_hi);
% [h_subj_pw_cp_hi, p_subj_pw_cp_low, ~,stats_subj_pw_cp_low] = ttest(pw_cp_contrast_lo);
% 
% 
% [h_subj_pw_sac_combined, p_subj_pw_sac_combined, ~, stats_subj_pw_sac_combined] = ttest(pw_sac_contrast_combined);
% [h_subj_pw_sac_hi, p_subj_pw_sac_hi, ~, stats_subj_pw_sac_hi] = ttest(pw_sac_contrast_hi);
% [h_subj_pw_sac_hi,  p_subj_pw_sac_low, ~, stats_subj_pw_sac_low] = ttest(pw_sac_contrast_lo);
% 
% [h_subj_pw_noise, p_subj_pw_noise, ~, stats_subj_pw_noise] = ttest(pw_noise_contrast);


[p_subj_pw_cp_combined, h_subj_pw_cp_combined, stats_subj_pw_cp_combined] = signrank(pw_cp_contrast_combined,0);
[p_subj_pw_cp_hi, h_subj_pw_cp_hi, stats_subj_pw_cp_hi] = signrank(pw_cp_contrast_hi,0);
[p_subj_pw_cp_low, h_subj_pw_cp_low,stats_subj_pw_cp_low] = signrank(pw_cp_contrast_lo,0);


[p_subj_pw_sac_combined, h_subj_pw_sac_combined,stats_subj_pw_sac_combined] = signrank(pw_sac_contrast_combined,0,'tail','right');
[p_subj_pw_sac_hi, h_subj_pw_sac_hi, stats_subj_pw_sac_hi] = signrank(pw_sac_contrast_hi,0);
[p_subj_pw_sac_low,  h_subj_pw_sac_low, stats_subj_pw_sac_low] = signrank(pw_sac_contrast_lo,0);

[p_subj_pw_noise, h_subj_pw_noise, stats_subj_pw_noise] = signrank(pw_noise_contrast,0);


clc


fprintf(' ---------- SUBJECTS : PRIOR BIAS CONTRASTS -------------\n\n')

fprintf('Subject prior bias CP contrast COMBINED pVal ---> %f  signedrank---> %f\n',...
          p_subj_pw_cp_combined, stats_subj_pw_cp_combined.signedrank)
      
fprintf('Subject prior bias SAC contrast COMBINED pVal ---> %f  signedrank---> %f\n\n',...
          p_subj_pw_sac_combined, stats_subj_pw_sac_combined.signedrank) 


fprintf('Subject prior bias CP contrast LOW pVal ---> %f  signedrank---> %f\n',...
          p_subj_pw_cp_low, stats_subj_pw_cp_low.signedrank)
      
fprintf('Subject prior bias SAC contrast LOW pVal ---> %f  signedrank---> %f\n\n',...
          p_subj_pw_sac_low, stats_subj_pw_sac_low.signedrank) 

fprintf('Subject prior bias CP contrast HIGH pVal ---> %f  signedrank---> %f\n',...
          p_subj_pw_cp_hi, stats_subj_pw_cp_hi.signedrank)
      
fprintf('Subject prior bias SAC contrast HIGH pVal ---> %f  signedrank---> %f\n',...
          p_subj_pw_sac_hi, stats_subj_pw_sac_hi.signedrank) 


fprintf('Subject prior bias NOISE contrast pVal ---> %f  signedrank---> %f\n\n\n\n',...
          p_subj_pw_noise, stats_subj_pw_noise.signedrank)
      
      

      
      
      
 [P_pw_cp,H,STATS]=ranksum(theory_subjective_pw_cp_contrast, pw_cp_contrast_combined)   %2
 [P_pw_sac_lo,H,STATS]=ranksum(theory_subjective_pw_sac_contrast_low, pw_sac_contrast_lo) %3
 [P_pw_sac_hi,H,STATS]=ranksum(theory_subjective_pw_sac_contrast_hi, pw_sac_contrast_hi) %3
 [P_pw_noise,H,STATS]=ranksum(theory_subjective_pw_noise_contrast, pw_noise_contrast) %4     
      