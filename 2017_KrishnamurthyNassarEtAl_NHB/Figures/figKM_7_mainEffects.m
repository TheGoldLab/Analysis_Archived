%% ************************************************************************************
%
%                      FIGURE 7 : MAIN PUPIL REGRESSION
%
% *************************************************************************************

getCbColors;

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


fig7_colors = zeros(6,3);
fig7_colors(1,:) = fig6_colors(4,:);
fig7_colors(2,:) = fig6_colors(5,:);
fig7_colors(3,:) = fig6_colors(6,:);
fig7_colors(4,:) = [220 70 5]*1.0/256;


resIdx=5;

% Use the J. Neurosci figure template
fig_num = 7; % figure number
total_width = 15; % total width
col_height = [7, 2, 2, 7]; % height of each column
cols_per_row = {[0.15 0.85],[0.18 0.32 0.32 0.18], [0.18 0.32 0.32 0.18], [0.15 0.85]}; % columns in each row
psh = 2.3; % panel separation height
psw = 1.5; % panel separation width


[axs,fig_] = getPLOT_axes(fig_num, total_width, col_height, cols_per_row,...
    psh, psw, [], '');
set(axs,'Units','normalized');


ax_font_size = 14;
ms = 8; % marker size
jit = 0.1; % jitter for smartJitter


%time stamps for pupil trace
tim = linspace(0,size(avg_subj_pup_trace_KK,2)/60.0,size(avg_subj_pup_trace_KK,2))';






%%  ----------------------- HOW DOES THE AVERAGE PUPIL LOOK ------------------

axes(axs(2)); cla(gca); hold on

pupColor = [0.05 0.05 0.05];

mns = nanmean(avg_subj_pup_trace_KK,1)';
ses = nanstd(avg_subj_pup_trace_KK)'/sqrt(size(avg_subj_pup_trace_KK,1));

l = line([0 2.5],[0,0],'color',[0.5 0.5 0.5],'linestyle','-');


hP = patch([tim; flipdim(tim,1)], [mns-ses; flipdim(mns+ses,1)], ...
    cbColors(7,:));
set(hP, 'LineStyle', 'none');
set(hP,'FaceColor', pupColor);
set(hP,'FaceAlpha',0.2);
aux = plot(tim, mns, 'color',pupColor,'linewidth',3);
hP(end+1) = aux;
set(aux,'LineWidth',3)



plot([tim(respTime), tim(respTime)], [-.05, .15], '--r')


set(axs(2),'XTick',[0 2.5])
set(axs(2),'XTickLabel',{'0','2.5'})

ylabel('Evoked pupil response (mm)')
xl_txt1 = xlabel(gca,'Time re: stimulus onset (s)');
aux_p1 = xl_txt1.Position;
aux_p1(2) = aux_p1(2) + 0.01;
xl_txt1.Position = aux_p1;
set(gca,'fontsize',ax_font_size)
setPLOT_panelLabel(gca, 1,0.2,1.1);


set(axs(1),'visible','off')








%% ------------------------ PUPIL BINNED BY BEHAVIOURAL VARIABLES --------------------------

marker_face_col = [0.9 0.9 0.9];
ms = 6; % marker size
ma = 60; % marker area
y_min = -0.03;
y_max = 0.03;



% --------- %
axes(axs(4)); cla(gca); hold on
hold on
[rho, pVals] = corr(nanmean(MC_pup(:,1:3),2),Stab,'type','spearman');
aux1 = nanmean(medStab);
aux2 = nanmean(stabBinBaseline);
%plot(aux1, aux2, 'o','markersize',ms, 'markerfacecolor',marker_face_col,...
%    'markeredgecolor',[0 0 0],'linewidth',1);
scatter(aux1, aux2,ma,'filled','markerfacecolor',fig7_colors(1,:), ...
    'markeredgecolor',fig7_colors(1,:),'markerfacealpha',0.5,...
     'markeredgealpha',1.0,'linewidth',1);
ylim([y_min y_max])

plot([nanmean(medStab); nanmean(medStab)], ...
    [nanmean(stabBinBaseline)+ste_stabBaseline; nanmean(stabBinBaseline)-ste_stabBaseline],  '-k')
ylabel({'Baseline'; '(z-scored)'})

setPLOT_panelLabel(gca, 2);
set(gca,'fontsize',ax_font_size)
titlestr = ['$p: ' num2str(pVals,'%1.2f') '$'];
t1 =title(titlestr,'fontsize',12,'interpreter','latex');


% --------- %
axes(axs(8)); cla(gca); hold on
hold on
aux1 = nanmean(medStab);
aux2 = nanmean(stabBinResponse);
[rho, pVals] = corr(MC_respResidual,Stab,'type','spearman');
%plot(aux1, aux2, 'o','markersize',ms, 'markerfacecolor',marker_face_col,...
%    'markeredgecolor',[0 0 0],'linewidth',1);
scatter(aux1, aux2,ma,'filled','markerfacecolor',fig7_colors(1,:), ...
    'markeredgecolor',fig7_colors(1,:),'markerfacealpha',0.5,...
     'markeredgealpha',1.0,'linewidth',1);
ylim([y_min y_max])
h = lsline;
h.LineWidth = 1;
h.LineStyle = '-';
h.Color = fig7_colors(1,:);
plot([nanmean(medStab); nanmean(medStab)], ...
    [nanmean(stabBinResponse)+ste_stabResponse; nanmean(stabBinResponse)-ste_stabResponse],  '-k')
ylabel({'Evoked'; '(z-scored)'})
xlabel('Relevance','color',fig7_colors(1,:),'fontweight','bold')
setPLOT_panelLabel(gca, 5);
set(gca,'fontsize',ax_font_size)
%titlestr = ['$\mathbf{\beta}$ : ' num2str(rho,'%1.3f') ' p: ' num2str(pVals,'%1.3f')];
%titlestr = [' p: ' num2str(pVals,'%1.2f')];
titlestr = ['$p<0.01$'];
t1 =title(titlestr,'fontsize',12,'interpreter','latex');


% --------- %
axes(axs(5)); cla(gca); hold on
hold on
aux1 = nanmean(medRel);
aux2 = nanmean(relBinBaseline);
[rho, pVals] = corr(nanmean(MC_pup(:,1:3),2),Rel,'type','spearman');
%plot(aux1, aux2, 'o','markersize',ms, 'markerfacecolor',marker_face_col,...
%    'markeredgecolor',[0 0 0],'linewidth',1);
scatter(aux1, aux2,ma,'filled','markerfacecolor',fig7_colors(2,:), ...
    'markeredgecolor',fig7_colors(2,:),'markerfacealpha',0.5,...
     'markeredgealpha',1.0,'linewidth',1);
h = lsline;
h.LineWidth = 1;
h.LineStyle = '-';
h.Color = fig7_colors(2,:);

plot([nanmean(medRel); nanmean(medRel)], ...
    [nanmean(relBinBaseline)+ste_relBaseline; nanmean(relBinBaseline)-ste_relBaseline],  '-k')
ylim([y_min y_max])
xlim([0.3 0.75])
set(gca,'yticklabel','');
setPLOT_panelLabel(gca, 3);
set(gca,'fontsize',ax_font_size)
%titlestr = ['$\mathbf{\beta}$ : ' num2str(rho,'%1.3f') ' p: ' num2str(pVals,'%1.3f')];
titlestr = ['$p<0.01$'];
t1 =title(titlestr,'fontsize',12,'interpreter','latex');


% --------- %
axes(axs(9)); cla(gca); hold on
hold on
[rho, pVals] = corr(MC_respResidual,Rel,'type','spearman');
aux1 = nanmean(medRel);
aux2 =  nanmean(relBinResponse);
%plot(aux1, aux2, 'o','markersize',ms, 'markerfacecolor',marker_face_col,...
%    'markeredgecolor',[0 0 0],'linewidth',1);
scatter(aux1, aux2,ma,'filled','markerfacecolor',fig7_colors(2,:), ...
    'markeredgecolor',fig7_colors(2,:),'markerfacealpha',0.5,...
     'markeredgealpha',1.0,'linewidth',1);
plot([nanmean(medRel); nanmean(medRel)], ...
    [nanmean(relBinResponse)+ste_relResponse; nanmean(relBinResponse)-ste_relResponse],  '-k')

ylim([y_min y_max])
xlim([0.3 0.75])
xlabel('Reliability','color',fig7_colors(2,:),'fontweight','bold')
set(gca,'yticklabel','');
setPLOT_panelLabel(gca, 6);
set(gca,'fontsize',ax_font_size)
titlestr = ['$p: ' num2str(pVals,'%1.2f') '$'];
t1 =title(titlestr,'fontsize',12,'interpreter','latex');


% --------- %
axes(axs(6)); cla(gca); hold on
hold on
[rho, pVals] = corr(nanmean(MC_pup(:,1:3),2),Bet,'type','spearman');
%plot([0, 1], nanmean(betBinBaseline(subjNormalBets,:)), 'o','markersize',ms, 'markerfacecolor',marker_face_col,...
%    'markeredgecolor',[0 0 0],'linewidth',1)
scatter([0 1], nanmean(betBinBaseline(subjNormalBets,:)),ma,'filled','markerfacecolor',fig7_colors(3,:), ...
    'markeredgecolor',fig7_colors(3,:),'markerfacealpha',0.5,...
     'markeredgealpha',1.0,'linewidth',1);
ylim([y_min y_max])
plot([0, 1; 0, 1], ...
    [nanmean(betBinBaseline(subjNormalBets,:))+ste_betBaseline; nanmean(betBinBaseline(subjNormalBets,:))-ste_betBaseline],  '-k')
xlim([-.3, 1.1])
set(gca,'yticklabel','');

setPLOT_panelLabel(gca, 4);
set(gca,'fontsize',ax_font_size)
titlestr = ['$p: ' num2str(pVals,'%1.2f') '$'];
t1 =title(titlestr,'fontsize',12,'interpreter','latex');

% --------- %
axes(axs(10)); cla(gca); hold on
hold on
[rho, pVals] = corr(MC_respResidual,Bet,'type','spearman');
%plot([0, 1], nanmean(betBinResponse(subjNormalBets,:)), 'o','markersize',ms, 'markerfacecolor',marker_face_col,...
%    'markeredgecolor',[0 0 0],'linewidth',1)
scatter([0 1], nanmean(betBinResponse(subjNormalBets,:)),ma,'filled','markerfacecolor',fig7_colors(3,:), ...
    'markeredgecolor',fig7_colors(3,:),'markerfacealpha',0.5,...
     'markeredgealpha',1.0,'linewidth',1);
ylim([y_min y_max])
plot([0, 1; 0, 1], ...
    [nanmean(betBinResponse(subjNormalBets,:))+ste_betResponse; nanmean(betBinResponse(subjNormalBets,:))-ste_betResponse],  '-k')
xlabel('Confidence','color',fig7_colors(3,:),'fontweight','bold')
xlim([-.3, 1.1])
set(gca,'yticklabel','');
setPLOT_panelLabel(gca, 7);
set(gca,'fontsize',ax_font_size)
%titlestr = [' p: ' num2str(pVals,'%1.3f')];
titlestr = [' $p<0.01$'];
t1 =title(titlestr,'fontsize',12,'interpreter','latex');

% --------- %
set(axs(3),'visible','off')
set(axs(7),'visible','off')

%% ----------------------- HOW DO THE REGRESSION COEFFICIENTS LOOK ----------------
%

axes(axs(11)); cla(gca); hold on

cbColors(1,:) = cbColors(7,:);

% 1 - Relevance
% 2 - Reliability
% 3 - Bet
% 4 - Residual

% avg pupil betas for plotting
avg_pupilBetas_base = mean(baseCoeffs(:,1:4));
std_pupilBetas_base = std(baseCoeffs(:,1:4))/sqrt(size(baseCoeffs,1));


h = [];
h_B = [];
betaOffset = -0.4;
delta = 0.09;
count = 1;
plotLegend = {};
avgPupilCoeffs = [];
%significantEpochs = [];


baseline_width = 0.4*ones(4,1); %[0.4 2 0.5 0.5 2 0.4];
baseline_width(P_base(1:4) < 0.05) = 2;

baseline_markersize = 4*ones(4,1);   %[5 5 9 5 5 9];
baseline_markersize(P_base(1:4) < 0.05) = 8;

pupilRegIdx = [1:4];
labels={ 'Relevance','Reliability','Bet','Residual'};

fig7_colors = zeros(6,3);
fig7_colors(1,:) = fig6_colors(4,:);
fig7_colors(2,:) = fig6_colors(5,:);
fig7_colors(3,:) = fig6_colors(6,:);
fig7_colors(4,:) = [220 70 5]*1.0/256;


l = line([betaOffset-0.1 0],[0,0],'color',[0.5 0.5 0.5]);
l.LineWidth = 1;
l.LineStyle = '-';

for i1 = pupilRegIdx
    
    
    yB = [avg_pupilBetas_base(i1) - std_pupilBetas_base(i1) ...
        avg_pupilBetas_base(i1) + std_pupilBetas_base(i1)];
    xB = betaOffset + count*delta;
    
    aux = plot(xB,avg_pupilBetas_base(i1),'o','MarkerSize', baseline_markersize(i1),...
        'MarkerFaceColor',fig7_colors(i1,:),...
        'MarkerEdgeColor',fig7_colors(i1,:));
    
    plot([xB xB],yB,'color',fig7_colors(i1,:),'linewidth',baseline_width(i1));
    
    h_B(end+1) = aux;
    
    count = count + 1;
    
end

ylim([-0.08, 0.024])



ylabel( '\beta : baseline','fontsize',16,'fontname','Arial')
set(axs(11),'YTick',[-0.06,0,0.02])
set(axs(11),'XTickLabel','')

set(axs(11),'FontSize',16)
setPLOT_panelLabel(gca, 8,0.22,1.06);




%% Now plot the coefficients in the evoked trace


axes(axs(12)); cla(gca); hold on
set(axs(12),'YAxisLocation','right')

%which traces survive multiple comparisons?
is_trace_sig = any(significantEpochs);

alpha_trace = [0.4 0.1 0.2 0.4];
width_trace = [2.5 0.75 1.5 2.5];
sig_width = 0.0003;

sigLineLevel = [-0.01 0 -0.01 -0.01 -0.01];
sigStep = [1e-4 5e-4 10e-4 20e-4];

avg_pupilBetas = squeeze(mean(respCoeffs,1));
std_pupilBetas = squeeze(std(respCoeffs)/sqrt(size(respCoeffs,1)));

l = line([0 2.5],[0,0],'color',[0.5 0.5 0.5]);
l.LineWidth = 1;
l.LineStyle = '-';

count = 1;
pupilRegIdx = [2 3 4 1];
for i1 = pupilRegIdx
    
    mns = avg_pupilBetas(:,i1);
    ses = std_pupilBetas(:,i1);
    
    hP = patch([tim; flipdim(tim,1)], [mns-ses; flipdim(mns+ses,1)], ...
        cbColors(min(count,5),:));
    set(hP, 'LineStyle', 'none');
    set(hP,'FaceColor', fig7_colors(i1,:));
    set(hP,'FaceAlpha',alpha_trace(i1));
    
    aux = plot(tim, mns, 'color',fig7_colors(i1,:),'linewidth',width_trace(i1));
    h(end+1) = aux;
    
    
    
    %     sigIdx = nan(length(tim),1);
    %     sigIdx(significantEpochs(:,i1)) = sigLineLevel(i1) + sigStep(i1);
    %     l = length(sigIdx);
    %     sigIdx1 = nan(length(tim),1);
    %     sigIdx1(1:4:l) = sigIdx(1:4:l);
    
    if(is_trace_sig(i1))
        sig_time = tim(significantEpochs(:,i1));
        sig_level = sigLineLevel(i1) - sigStep(i1);
        aux_sig = sig_level*ones(length(sig_time),1);
        hP_sig = patch([sig_time; flipdim(sig_time,1)], [aux_sig; flipdim(aux_sig+sig_width,1)], ...
            fig7_colors(i1,:));
        
        set(hP_sig, 'LineStyle', 'none');
        set(hP_sig,'FaceColor', fig7_colors(i1,:));
        set(hP_sig,'FaceAlpha',0.7);
    end
    %aux_1 = plot(tim,sigIdx1,'*','color',fig7_colors(i1,:),'linewidth',2);
    %aux_1.Parent = ax_2;
    %aux_1.MarkerSize = 4;
    %aux_1.MarkerFaceColor  = fig7_colors(i1,:);
    
    
    plotLegend(count) = labels(i1);
    
    
    
    count = count + 1;
    
end

ylim([-0.014, 0.025])
xlim([0, 2.5])

set(axs(12),'YTick',[-0.01,0,0.01,0.02])

set(axs(12),'XTick',[0 2.5])

set(axs(12),'FontSize',16)

leg_text1 = text(0.01,0.92,'Residual','fontsize',13,...
     'color',fig7_colors(4,:),'Units','normalized',...
     'fontweight','bold');
% leg_text1.FontName = 'Arial Narrow';
leg_text1.FontSize = 16;

%leg = legend(h,plotLegend,'Location','NorthWest','fontsize',12);
% leg_text1 = text(0.01,0.96,'Relevance','fontsize',13,...
%     'color',fig7_colors(1,:),'Units','normalized',...
%     'fontweight','bold');
% leg_text1.FontName = 'Arial Narrow';
% leg_text1.FontSize = 14;
% leg_text2 = text(0.01,0.89,'Reliability','fontsize',13,...
%     'color',fig7_colors(2,:),'Units','normalized',...
%     'fontweight','bold');
% leg_text2.FontName = 'Arial Narrow';
% leg_text2.FontSize = 14;
% leg_text3 = text(0.01,0.82,'Confidence','fontsize',13,...
%     'color',fig7_colors(3,:),'Units','normalized',...
%     'fontweight','bold');
% leg_text3.FontName = 'Arial Narrow';
% leg_text3.FontSize = 14;
% leg_text4 = text(0.01,0.75,'Residual','fontsize',13,...
%     'color',fig7_colors(4,:),'Units','normalized',...
%     'fontweight','bold');
% leg_text4.FontName = 'Arial Narrow';
% leg_text4.FontSize = 14;

ylabel({'\beta : evoked'},'fontsize',16,'fontname','Arial')
xl_txt2 = xlabel(gca,'Time re: stimulus onset (s)');
aux_p2 = xl_txt2.Position;
aux_p2(2) = aux_p2(2) + 0.004;
xl_txt2.Position = aux_p2

%xlabel(gca,'Time(s)')
%xlim([-0.4 tim(end)]);
%ylim([-0.42 0.42]);
line([0 0],ylim,'color','k','linewidth',1)

%title('Within-subject regression','fontsize',12)

set(gca,'fontsize',ax_font_size)
setPLOT_panelLabel(gca, 9,0.06,1.06);


kk=annotation('textbox');
set(kk, 'string', 'Krishnamurthy, Nassar et al. figure 7', 'position', ...
    [0.61 0.9 0.10 0.1], 'EdgeColor', 'none','fontsize',14)

% fig_title=annotation('textbox');
% set(fig_title, 'string', 'Within subject pupil effects', 'position', ...
%     [0.05 0.88 0.10 0.1], 'EdgeColor', 'none','fontsize',14)

figName = '~/Dropbox/auditoryTaskManuscript/Figure_set_for_writeup/figure_7_withinSubj.pdf';

%print(gcf,figName,'-dpdf','-fillpage','-painters')

