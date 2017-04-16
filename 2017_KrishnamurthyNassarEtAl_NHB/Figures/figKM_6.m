%% ************************************************************************************
%
%                      FIGURE 6 : THE LINEAR MODEL
%
% *************************************************************************************

getCbColors;

% Use the J. Neurosci figure template
fig_num = 6; % figure number
total_width = 17; % total width
col_height = [5, 8]; % height of each column
cols_per_row = {[0.3 0.3 0.4], [0.9]}; % columns in each row

col_height = [5, 5, 6]; % height of each column
cols_per_row = {[.45, .45], [.45, .45], [0.9]}; % columns in each row

psh = 2; % panel separation height
psw = 2; % panel separation width


[axs,fig_] = getPLOT_axes(fig_num, total_width, col_height, cols_per_row,...
                                 psh, psw, [], '');
set(axs,'Units','normalized');


ax_font_size = 12;
ms = 8; % marker size
jit = 0.1; % jitter for smartJitter


%%  --------------------- Stability and Uncertainty terms vs. SAC -------------

axes(axs(3)); cla(gca); hold on


stab_plot_lims = [0.4 1.0];

auxMu_hi = nanmean(stabTerm_SAC_hi');
auxSE_hi = nanstd(stabTerm_SAC_hi')./sqrt(size(stabTerm_SAC_hi',1));

auxMu_low = nanmean(stabTerm_SAC_low');
auxSE_low = nanstd(stabTerm_SAC_low')./sqrt(size(stabTerm_SAC_low',1));

barInds=[1:6; 1:6];
errBars_hi = [auxMu_hi + auxSE_hi ; auxMu_hi - auxSE_hi];
errBars_low = [auxMu_low + auxSE_low ; auxMu_low - auxSE_low];

plot(barInds,errBars_hi, '-k', 'lineWidth', 1.5)
plot(barInds,errBars_low, '-k', 'lineWidth', 1.5)
h1 = plot(auxMu_hi, 'ob', 'markerSize', ms, 'markerFaceColor', cbColors(2,:), 'markerEdgeColor', 'k', 'lineWidth', 1);
h2 = plot(auxMu_low, 'ob', 'markerSize', ms, 'markerFaceColor', cbColors(3,:), 'markerEdgeColor', 'k', 'lineWidth', 1);
%legend([h1,h2],{'high noise','low noise'},'fontsize',ax_font_size,'Location','best')

ylim(stab_plot_lims)
xlim([0.5, 6.5])
xlabel('SAC')
ylabel('Prior relevance','fontsize',ax_font_size)
%title('Subject Pred Error','fontsize',14);
grid on
set(gca,'fontsize',ax_font_size)
setPLOT_panelLabel(gca, 3,0.18,1.2);




% now plot the uncertainty term
axes(axs(4)); cla(gca); hold on

stab_plot_lims = [-0.15 0.15];
delta = 0.2;

auxMu_hi = nanmean(uncTerm_SAC_hi');
auxSE_hi = nanstd(uncTerm_SAC_hi')./sqrt(size(uncTerm_SAC_hi',1));

auxMu_low = nanmean(uncTerm_SAC_low');
auxSE_low = nanstd(uncTerm_SAC_low')./sqrt(size(uncTerm_SAC_low',1));

barInds=[1:6; 1:6];
errBars_hi = [auxMu_hi + auxSE_hi ; auxMu_hi - auxSE_hi];
errBars_low = [auxMu_low + auxSE_low ; auxMu_low - auxSE_low];

plot(barInds,errBars_hi, '-k', 'lineWidth', 1.5)
plot(barInds+delta,errBars_low, '-k', 'lineWidth', 1.5)
plot(barInds,auxMu_hi, 'ob', 'markerSize', ms, 'markerFaceColor', cbColors(2,:), 'markerEdgeColor', 'k', 'lineWidth', 1)
plot(barInds+delta,auxMu_low, 'ob', 'markerSize', ms, 'markerFaceColor', cbColors(3,:), 'markerEdgeColor', 'k', 'lineWidth', 1)

%ylim(stab_plot_lims)
xlim([0.5, 6.5])
xlabel('SAC')
ylabel('Prior reliability','fontsize',ax_font_size)
%title('Subject Pred Error','fontsize',14);
grid on
set(gca,'fontsize',ax_font_size)
setPLOT_panelLabel(gca, 4,0.18,1.2);





%% --------------- PW comparison of subjects with the generative model 
axes(axs(1)); cla(gca); hold on
count = 1;


pltIdx = 1:6; %indices corresponding to SAC 1-6


% PLOT THE MODEL PRIOR WT vs. SAC
tim = (1:6)';
mns = genModAvgBetas_split_lo(pltIdx);
ses = genModBetaInt_split_lo(pltIdx);
hP = patch([tim; flipdim(tim,1)], [mns-ses; flipdim(mns+ses,1)], ...
            cbColors(min(count,5),:));
set(hP, 'LineStyle', 'none');
set(hP,'FaceColor', cbColors(3,:));
set(hP,'FaceAlpha',0.5);

mns = genModAvgBetas_split_hi(pltIdx);
ses = genModBetaInt_split_hi(pltIdx);
hP = patch([tim; flipdim(tim,1)], [mns-ses; flipdim(mns+ses,1)], ...
            cbColors(min(count,5),:));
set(hP, 'LineStyle', 'none');
set(hP,'FaceColor', cbColors(2,:));
set(hP,'FaceAlpha',0.5);


count = 1;
pointHandles = [];
delta = [0 0.1];    %spacing between points for different stds


plot(jit+repmat((1:numTAC)+delta(count), 2, 1),...
    [subjAvgBetas_allTrials_lo(pltIdx)-subjBetaInt_allTrials_lo(pltIdx) subjAvgBetas_allTrials_lo(pltIdx)+subjBetaInt_allTrials_lo(pltIdx)]',...
    'k','linewidth',1.5);
h1 = plot(jit+(1:numTAC)+delta(count),...
    subjAvgBetas_allTrials_lo(pltIdx), 'ok', 'markerEdgeColor', 'k', ...
    'markerFaceColor', cbColors(3,:),...
    'lineWidth', 1, 'markerSize', ms);

plot(repmat((1:numTAC)+delta(count), 2, 1)-jit,...
    [subjAvgBetas_allTrials_hi(pltIdx)-subjBetaInt_allTrials_hi(pltIdx) subjAvgBetas_allTrials_hi(pltIdx)+subjBetaInt_allTrials_hi(pltIdx)]',...
    'k','linewidth',1.5);
h2 = plot((1:numTAC)+delta(count)-jit,...
    subjAvgBetas_allTrials_hi(pltIdx), 'ok', 'markerEdgeColor', 'k', ...
    'markerFaceColor', cbColors(2,:),...
    'lineWidth', 1, 'markerSize', ms);

legend([h2,h1],{'high noise','low noise'},'fontsize',ax_font_size,'Location','best','Box','off')


count = count + 1;

%legend(pointHandles,'LOW', 'HIGH','location','SouthEast');

grid on
ylabel('Perceptual bias');

xlim([0 6.5])
set(gca, 'fontsize',ax_font_size)

setPLOT_panelLabel(gca, 1,0.18,1.2);



axes(axs(2)); cla(gca); hold on
count = 1;

pltIdx = 1:6; %indices corresponding to SAC 1-6


% PLOT THE MODEL PRIOR WT vs. SAC
tim = (1:6)';
mns = linearApprxAvgBetas_split_lo(pltIdx);
ses = linearApprxBetaInt_split_lo(pltIdx);
hP = patch([tim; flipdim(tim,1)], [mns-ses; flipdim(mns+ses,1)], ...
            cbColors(min(count,5),:));
set(hP, 'LineStyle', 'none');
set(hP,'FaceColor', cbColors(3,:));
set(hP,'FaceAlpha',0.5);

mns = linearApprxAvgBetas_split_hi(pltIdx);
ses = linearApprxBetaInt_split_hi(pltIdx);
hP = patch([tim; flipdim(tim,1)], [mns-ses; flipdim(mns+ses,1)], ...
            cbColors(min(count,5),:));
set(hP, 'LineStyle', 'none');
set(hP,'FaceColor', cbColors(2,:));
set(hP,'FaceAlpha',0.5);


count = 1;
pointHandles = [];
delta = [0 0.1];    %spacing between points for different stds


plot(jit+repmat((1:numTAC)+delta(count), 2, 1),...
    [subjAvgBetas_allTrials_lo(pltIdx)-subjBetaInt_allTrials_lo(pltIdx) subjAvgBetas_allTrials_lo(pltIdx)+subjBetaInt_allTrials_lo(pltIdx)]',...
    'k','linewidth',1.5);
pointHandles(end+1) = plot(jit+(1:numTAC)+delta(count),...
    subjAvgBetas_allTrials_lo(pltIdx), 'ok', 'markerEdgeColor', 'k', ...
    'markerFaceColor', cbColors(3,:),...
    'lineWidth', 1, 'markerSize', ms);

plot(repmat((1:numTAC)+delta(count), 2, 1)-jit,...
    [subjAvgBetas_allTrials_hi(pltIdx)-subjBetaInt_allTrials_hi(pltIdx) subjAvgBetas_allTrials_hi(pltIdx)+subjBetaInt_allTrials_hi(pltIdx)]',...
    'k','linewidth',1.5);
pointHandles(end+1) = plot((1:numTAC)+delta(count)-jit,...
    subjAvgBetas_allTrials_hi(pltIdx), 'ok', 'markerEdgeColor', 'k', ...
    'markerFaceColor', cbColors(2,:),...
    'lineWidth', 1, 'markerSize', ms);

count = count + 1;

%legend(pointHandles,'LOW', 'HIGH','location','SouthEast');

grid on

% xlabel('SAC');
xlim([0 6.5])
set(gca, 'fontsize',ax_font_size)

setPLOT_panelLabel(gca, 2,0.18,1.2);

%% ----------------- LINEAR MODEL V1a REGRESSION PLOT-----------------------


axes(axs(5)); cla(gca); hold on
getCbColors
% Plot this using Matt's smart jitter function

% LABELS FOR MODEL 1B:
% xTickLabels={ 'intercept', 'PE','PE*relevance', 'PE*reliability','PE*Bet', 'PE*noise',...
%                  'PE*CentDist' , 'Center bias', 'Subjective bias'}; 
 % LABELS FOR MODEL 1:                
 xTickLabels={ 'intercept', 'PE','PE*relevance', 'PE*reliability','PE*Bet', ...
                  'Center bias', 'Subjective bias'}; 

ll=size(beta_v1,1);

%which ones are significant
[is_beta_v1_sig, p_beta_v1] = ttest(beta_v1);

% don't include subjects with all or none bets for the bet term:
[is_beta_v1_sig(5), p_beta_v1(5)] = ttest(beta_v1(subjNormalBets,5));




clear P CI np_P
hold on
wagerIdx = 5;  % Column of wager

aux_col1 = [139,0,139]/256;
aux_col2 = [0.7 0.6 0.2]*0.85;

fig6_colors = zeros(10,3);
fig6_colors(2,:) = [0.65 0.65 0.65];
fig6_colors(3,:) = [0.7 0.2 0.3];
fig6_colors(4,:) = [0.1 0.4 0.05];
fig6_colors(5:end,:) = cbColors(3:end,:);
%fig6_colors(6,:) = fig6_colors(7,:)*0.97;
fig6_colors(5,:) = aux_col2;
fig6_colors(7,:) = aux_col1;
fig6_colors(6,:) = cbColors(6,:)*0.8;
fig6_colors(9,:) = [0.4 0.35 0.55];
fig6_colors(8,:) = [0.55 0.35 0.4];

%fig6_colors(5,:) = 
plot([0 ll*4], [0 0], '-','Color',[0.5 0.5 0.5])
for i=1:size(beta_v1,2)
    if(i ~= wagerIdx)
    auxBeta = beta_v1(:,i);
    
    xJit = smartJitter(auxBeta./std(auxBeta),.2,.2);
    plot(ones(length(auxBeta), 1).*i.*4+xJit , (auxBeta)./std(auxBeta),...
        'o', 'markerSize', 7, 'markerFaceColor', fig6_colors(mod(i,9)+1,:),...
        'markerEdgeColor', 'k', 'lineWidth', 1);
    xlim([0 4.*i+1.75])
    
    
    
    modAvBeta=nanmean(avBeta_v1_sim(:,i)./std(auxBeta));
    
    % OK, instead of 95% confidence over mean, lets do 95% confidence over
    % individual subject estimates
    allModCoeffs=squeeze((beta_v1_sim(:,i,:)))./std(auxBeta);
    else
       auxBeta = beta_v1(subjNormalBets,i);
       xJit = smartJitter(auxBeta./std(auxBeta),.2,.2);
    plot(ones(length(auxBeta), 1).*i.*4+xJit , (auxBeta)./std(auxBeta),...
        'o', 'markerSize', 7, 'markerFaceColor', fig6_colors(mod(i,9)+1,:),...
        'markerEdgeColor', 'k', 'lineWidth', 1);
    xlim([0 4.*i+1.75])
    
    modAvBeta=nanmean(avBeta_v1_sim(subjNormalBets,i)./std(auxBeta));
    
    
    allModCoeffs=squeeze((beta_v1_sim(subjNormalBets,i,:)))./std(auxBeta);
    
    
    end
    
    CI_mean = prctile(nanmean(allModCoeffs, 1),[2.5, 97.5]);
    allModCoeffs=allModCoeffs(:);
    CI = prctile(allModCoeffs,[2.5, 97.5]);
    
   
     
    plot([i.*4+1.5, i.*4+1.5] , CI, ...
        '-', 'markerSize', 8, 'color', fig6_colors(mod(i,9)+1,:), 'lineWidth', 1);
    
       plot([i.*4+1.5, i.*4+1.5] , CI_mean, ...
        '-', 'markerSize', 8, 'color', fig6_colors(mod(i,9)+1,:), 'lineWidth', 10);
 
%     plot(i.*4+1.5 , modAvBeta,...
%         'd', 'markerSize', 8, 'markerFaceColor', [1 1 1],...
%         'markerEdgeColor', fig6_colors(mod(i,9)+1,:), 'lineWidth', 1.75);
    
    if(is_beta_v1_sig(i))
       plot(i*4,6.1,'*','markersize',5,'color','k'); 
    end
    
end

ylim([-4.1 6.5])
set(gca, 'xtick', [4:4:ll*4], 'XTickLabel', xTickLabels);
set(gca,'fontsize',14)
%title(['Relation of model variables with Subject Behaviour'])
grid on

ylabel('$\beta$','fontsize',18,'interpreter','latex')
set(axs(5),'XTickLabelRotation',90);

setPLOT_panelLabel(gca, 5,0.07,1.07);

%%
figName = '/Users/kamesh/Dropbox/auditoryTaskManuscript/Figure_set_for_writeup/figure_6.pdf';


kk=annotation('textbox');
set(kk, 'string', 'Krishnamurthy et al, figure 6', 'position', ...
    [0.77 0.88 0.15 0.05], 'EdgeColor', 'none','fontsize',14)
