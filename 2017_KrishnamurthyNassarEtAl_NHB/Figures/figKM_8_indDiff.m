%% ************************************************************************************
%
%                      FIGURE 8 : Pupil individual differences
%
% *************************************************************************************

getCbColors;

resIdx=5;

% Use the J. Neurosci figure template
fig_num = 8; % figure number
total_width = 15; % total width
col_height = [4 4 7]; % height of each column
cols_per_row = {[ 0.45 0.45], [0.45  0.45], [.15 0.75]}; % columns in each row
psh = 2.8; % panel separation height
psw = 1.8; % panel separation width


[axs,fig_] = getPLOT_axes(fig_num, total_width, col_height, cols_per_row,...
    psh, psw, [], '');
set(axs,'Units','normalized');


ax_font_size = 14;
ms = 8; % marker size
jit = 0.1; % jitter for smartJitter


%time stamps for pupil trace
tim = linspace(0,size(avg_subj_pup_trace_KK,2)/60.0,size(avg_subj_pup_trace_KK,2))';



%%  ----------------------- HOW DOES THE AVERAGE PUPIL LOOK ------------------

marker_face_col = [0.9 0.9 0.9];
ms = 6; % marker size
y_min = -0.03;
y_max = 0.03;

axes(axs(1)); cla(gca); hold on
aux1 = beta_v1(:,2);
aux2 = subMeanBaseline';
[rho, pVals] = corr(aux1,aux2,'type','spearman');
[RHO1,PVAL1]=corr(beta_v1(:,2), subMeanBaseline');
plot(aux1, aux2, 'o','markersize',ms, 'markerfacecolor',marker_face_col,...
    'markeredgecolor',[0 0 0],'linewidth',1);
%ylim([y_min y_max])
ylabel({'Raw baseline'; '(mm)'})
setPLOT_panelLabel(gca, 1);
set(gca,'fontsize',ax_font_size)
titlestr = ['$p: ' num2str(pVals,'%1.2f') '$'];
t1 =title(titlestr,'fontsize',12,'interpreter','latex');


axes(axs(2)); cla(gca); hold on
aux1 = beta_v1(:,3);
aux2 = subMeanBaseline';
[rho, pVals] = corr(aux1,aux2,'type','spearman');
[RHO2,PVAL2]=corr(beta_v1(:,3), subMeanBaseline');
plot(aux1, aux2, 'o','markersize',ms, 'markerfacecolor',marker_face_col,...
    'markeredgecolor',[0 0 0],'linewidth',1);
setPLOT_panelLabel(gca, 2);
set(gca,'yticklabel','');
set(gca,'fontsize',ax_font_size)
titlestr = ['$p: ' num2str(pVals,'%1.2f') '$'];
t1 =title(titlestr,'fontsize',12,'interpreter','latex');


axes(axs(3)); cla(gca); hold on
aux1 = beta_v1(:,2);
aux2 = subjMeanRespResid;
[rho, pVals] = corr(aux1,aux2,'type','spearman');
[RHO3,PVAL3]=corr(beta_v1(:,2), subjMeanRespResid);
plot(aux1, aux2, 'o','markersize',ms, 'markerfacecolor',marker_face_col,...
    'markeredgecolor',[0 0 0],'linewidth',1);
h = lsline;
h.LineWidth = 1;
h.LineStyle = '-';
h.Color = [0 0 0];
setPLOT_panelLabel(gca, 3);
ylabel({'Evoked response'; '(mm)'})
xlabel('Overall bias')
set(gca,'fontsize',ax_font_size)
titlestr = ['$p: ' num2str(pVals,'%1.2f') '$'];
t1 =title(titlestr,'fontsize',12,'interpreter','latex');


axes(axs(4)); cla(gca); hold on
aux1 = beta_v1(:,3);
aux2 = subjMeanRespResid;
[rho, pVals] = corr(aux1,aux2,'type','spearman');
[RHO4,PVAL4]=corr(beta_v1(:,3), subjMeanRespResid);
plot(aux1, aux2, 'o','markersize',ms, 'markerfacecolor',marker_face_col,...
    'markeredgecolor',[0 0 0],'linewidth',1);
h = lsline;
h.LineWidth = 1;
h.LineStyle = '-';
h.Color = [0 0 0];
setPLOT_panelLabel(gca, 4);
xlabel('Relevance-dependent bias')
set(gca,'yticklabel','');
set(gca,'fontsize',ax_font_size)
titlestr = ['$p<0.01$'];
t1 =title(titlestr,'fontsize',12,'interpreter','latex');




disp(sprintf('perceptual bias -- pupil baseline:    R = %g, P = %g', RHO1, PVAL1));
disp(sprintf('relevance dep bias -- pupil baseline: R = %g, P = %g', RHO2, PVAL2));
disp(sprintf('perceptual bias -- pupil response:    R = %g, P = %g', RHO3, PVAL3));
disp(sprintf('relevance dep bias -- pupil response: R = %g, P = %g', RHO4, PVAL4));

% --------- %
%set(axs(2),'visible','off')
%set(axs(5),'visible','off')


%% ----------------------- HOW DO THE REGRESSION COEFFICIENTS LOOK ----------------
%
axes(axs(5)); cla(gca); hold on

aux = jet(10);
corrCol(1,:) = aux(3,:);
corrCol(2,:) = aux(7,:)*0.8;

corrCol(1,:) = [0.6 0.3 0.3];
corrCol(2,:) = [0.3 0.3 0.6];
xB=0:1;

count = 1;

l = line(xB+[-.5, .5],[0,0],'color',[0.5 0.5 0.5]);
l.LineWidth = 1;
l.LineStyle = '-';

% SHARED
plot([xB(1); xB(1)],base_pupCorrInt(2,:)','color',corrCol(1,:),'linewidth',2);

aux_1 = plot(xB(1), baseIndPupCorr(2,:),'o','MarkerSize', 6,...
    'MarkerFaceColor',corrCol(1,:),...
    'MarkerEdgeColor',corrCol(1,:));

% UNIQUE
plot([xB(2); xB(2)],base_pupCorrInt(3,:)','color',corrCol(2,:),'linewidth',2);

aux_2 = plot(xB(2), baseIndPupCorr(3,:),'o','MarkerSize', 6,...
    'MarkerFaceColor',corrCol(2,:),...
    'MarkerEdgeColor',corrCol(2,:));


ylim([-0.6 0.6]);
ylabel( '\beta : baseline ','fontsize',16)

xlim([-0.5 1.5]);

set(axs(5),'YTick',[-0.5,0,0.5])
set(axs(5),'XTickLabel','')
set(axs(5),'FontSize',16)
setPLOT_panelLabel(gca, 5,0.22,1.06);




axes(axs(6)); cla(gca); hold on
set(axs(6),'YAxisLocation','right')



h = [];

l  = line([0 2.5],[0,0],'color',[0.5 0.5 0.5],'linestyle','-');
l.LineWidth = 1.5;




%which traces survive multiple comparisons?
sig_width = 0.004;

sigLineLevel = [-0.20 -0.2];
sigStep = [1.5e-2 3e-2];


if(any( sharedSigTimes))
    sig_time = tim(sharedSigTimes>0);
    sig_level = sigLineLevel(1) - sigStep(1);
    aux_sig = sig_level*ones(length(sig_time),1);
    hP_sig = patch([sig_time; flipdim(sig_time,1)], [aux_sig; flipdim(aux_sig+sig_width,1)], ...
        fig7_colors(i1,:));
    
    set(hP_sig, 'LineStyle', 'none');
    set(hP_sig,'FaceColor', corrCol(1,:));
    set(hP_sig,'FaceAlpha',0.7);
    
end


if(any( uniqueSigTimes))
    sig_time = tim(uniqueSigTimes>0);
    sig_level = sigLineLevel(2) - sigStep(2);
    aux_sig = sig_level*ones(length(sig_time),1);
    hP_sig = patch([sig_time; flipdim(sig_time,1)], [aux_sig; flipdim(aux_sig+sig_width,1)], ...
        fig7_colors(i1,:));
    
    set(hP_sig, 'LineStyle', 'none');
    set(hP_sig,'FaceColor', corrCol(2,:));
    set(hP_sig,'FaceAlpha',0.7);
    
end

% Plot the shared correlation
mns = indPupCorr(:,2)';
ses = squeeze(b_int(:,2,:));

hP = patch([tim; flipdim(tim,1)], [ses(:,1); flipdim(ses(:,2),1)], ...
    cbColors(2,:));
set(hP, 'LineStyle', 'none');
set(hP,'FaceColor', corrCol(1,:));
set(hP,'FaceAlpha',0.2);
aux = plot(tim, mns, 'color',corrCol(1,:),'linewidth',3);
h(end+1) = aux;
set(aux,'LineWidth',3)


% Plot the unique correlation
mns = indPupCorr(:,3)';
ses = squeeze(b_int(:,3,:));
hP = patch([tim; flipdim(tim,1)], [ses(:,1); flipdim(ses(:,2),1)], ...
    cbColors(2,:));
set(hP, 'LineStyle', 'none');
set(hP,'FaceColor', corrCol(2,:));
set(hP,'FaceAlpha',0.2);
aux = plot(tim, mns, 'color',corrCol(2,:),'linewidth',3);
h(end+1) = aux;
set(aux,'LineWidth',3)

%grid on
ylabel({'\beta : evoked'},'fontsize',16)

xlim([0,2.5]);
ylim([-0.23 0.15])
line([0 0],ylim,'color','k','linewidth',1)

%
%
%set(axs(6),'YTick',[-1.0,0,1.0])


% line indicating peak response (for scatter)
plot([tim(respTime), tim(respTime)], ylim, '--r')

set(axs(6),'FontSize',16)
xlabel(gca,'Time re: stimulus onset (s)');
%title('Individual-differences regression','fontsize',12)
%legend(h,'Shared','Unique','Location','best','fontsize',12)

leg_text5 = text(0.01,0.94,'Shared','fontsize',13,...
    'color',corrCol(1,:),'Units','normalized',...
    'fontweight','bold');
leg_text6 = text(0.01,0.86,'Unique','fontsize',13,...
    'color',corrCol(2,:),'Units','normalized',...
    'fontweight','bold');

line([0 0],ylim,'color','k','linewidth',1)

setPLOT_panelLabel(gca, 6,0.06,1.06);

kk=annotation('textbox');
set(kk, 'string', 'Krishnamurthy, Nassar et al. figure 8', 'position', ...
    [0.78 0.92 0.15 0.05], 'EdgeColor', 'none','fontsize',14)




figName = '~/Dropbox/auditoryTaskManuscript/Figure_set_for_writeup/figure_8.pdf';