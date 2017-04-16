%% ************************************************************************************
%
%                      FIGURE 5 : BEHAVIOUR INDIVIDUAL DIFFERENCES
%
% *************************************************************************************


getCbColors;

% Use the J. Neurosci figure template
fig_num = 5; % figure number
total_width = 17; % total width
col_height = [5 5 5 ]; % height of each column
cols_per_row = { [0.45, 0.45],[0.45, 0.45],[0.45, 0.45]}; % columns in each row
psh = 3; % panel separation height
psw = 3.5; % panel separation width


[axs,fig_] = getPLOT_axes(fig_num, total_width, col_height, cols_per_row,...
                                 psh, psw, [], '');
set(axs,'Units','normalized');


ax_font_size = 12;
ax_font_size1 = 12;
ax_font_size2 = 10;
ms = 9; % marker size
jit = 0.05; % jitter for smartJitter

marker_face_col = [0.9 0.9 0.9];


figName = '/Users/kamesh/Dropbox/auditoryTaskManuscript/Figure_set_for_writeup/figure_5.pdf';

%% -------- Residula prior weight vs. pred err std
% 

% MATT creates the most basic measures of:

% key variables: 
% average non-changepoint prediction error dispesion:  predErr_avNcp_contrast ;
% average non-changepoint prior weight              :  pw_avNcp_contrast
% likelihood dispersion:                            :  train_PE_subject

% train_PE_subject ; % training likelihood width

%pw_avNcp_contrast = zscore(pw_avNcp_contrast);
aux_train_PE_subject = (train_PE_subject);
aux_predErr_avNcp_contrast = (predErr_avNcp_contrast);
xMat=[ones(size(train_PE_subject)), aux_train_PE_subject, aux_predErr_avNcp_contrast];
STATS=regstats( pw_avNcp_contrast, xMat(:,2:end)) ;

% this is the p values for t-tests on beta ~=0
% OK: statistics check out. Both factors seem to contribute to overall
% prior weight:

clc
disp(sprintf('stats for paper: fstat and associated p:  \n %s \n %s', num2str(STATS.fstat.f), num2str(STATS.fstat.pval)))
disp(sprintf('perceptual reliability stats for paper: beta/t/p:  \n %s \n %s \n %s', num2str(STATS.tstat.beta(2)), num2str(STATS.tstat.t(2)),  num2str(STATS.tstat.pval(2))))
disp(sprintf('prediction error  stats for paper: beta/t/p:  \n %s \n %s \n %s', num2str(STATS.tstat.beta(3)), num2str(STATS.tstat.t(3)),  num2str(STATS.tstat.pval(3))))




disp(sprintf('prior weight effects: \n training error width p = %0.2g \n prediction error width p = %0.2g', STATS.tstat.pval (2), STATS.tstat.pval (3) ))

par_rho_1 = STATS.tstat.beta(2);
pVals_1 = STATS.tstat.pval(2);
par_rho_2 = STATS.tstat.beta(3);
pVals_2 = STATS.tstat.pval(3);
% OK, now lets get residuals for plots:


[B1,BINT,R_wTraining,RINT,STATS] = regress( pw_avNcp_contrast, xMat(:, [1, 2])) ;
[B2,BINT,R_wPredErrs,RINT,STATS] = regress( pw_avNcp_contrast, xMat(:, [1, 3])) ;




axes(axs(1)); cla(gca); hold on

aux1=train_PE_subject;
aux2=R_wPredErrs;
[rho pVals] = corr(aux1, aux2,'type','spearman');
plot(aux1, aux2, 'o','markersize',8, 'markerfacecolor',marker_face_col,...
     'markeredgecolor',[0 0 0],'linewidth',1);
h = lsline;
h.LineWidth = 2;
h.LineStyle = '-';
h.Color = [0 0 0];

titlestr = ['$\mathbf{\beta}$ : ' num2str(par_rho_1,'%1.3f') ' p: ' num2str(pVals_1,'%1.3f')];
t1 =title(titlestr,'fontsize',14,'interpreter','latex');

%title([' RHO ' num2str(rho) '  pVal ' num2str(pVals)],'fontsize',ax_font_size1);
%ylabel({'Perceptual bias residual',' for prediction accuracy'},'fontsize',ax_font_size1)
ylabel('Perceptual bias residual','fontsize',ax_font_size1)
xlabel('Control error STD (deg)','fontsize',ax_font_size1)
grid on
set(gca,'fontsize',ax_font_size)


setPLOT_panelLabel(gca, 1,0.15,1.15);


%% -------- Average prior weight vs. subject estimation error

axes(axs(2)); cla(gca); hold on
%aux_avg_pw = mean(subj_pw(2:7,:) + subj_pw(8:13,:));

aux2 = R_wTraining;
aux1 = predErr_avNcp_contrast; % training likelihood width
%aux1 = subj_estErr_std';


[rho pVals] = corr(aux1, aux2,'type','spearman');
plot(aux1, aux2, 'o','markersize',8, 'markerfacecolor',marker_face_col,...
     'markeredgecolor',[0 0 0],'linewidth',1);
h = lsline;
h.LineWidth = 2;
h.LineStyle = '-';
h.Color = [0 0 0];

titlestr = ['$\mathbf{\beta}$ : ' num2str(par_rho_2,'%1.3f') ' p: ' num2str(pVals_2,'%1.3f')];
t1 =title(titlestr,'fontsize',14,'interpreter','latex');

%title([' RHO ' num2str(rho) '  pVal ' num2str(pVals)],'fontsize',ax_font_size1);
%ylabel({'Perceptual bias residual','for control task accuracy'},'fontsize',ax_font_size1)
ylabel('Perceptual bias residual','fontsize',ax_font_size1)
xlabel('Prediction error STD (deg)','fontsize',ax_font_size1)
grid on
set(gca,'fontsize',ax_font_size)


% 
setPLOT_panelLabel(gca, 2,0.15,1.15);




%% -------- Residula prior weight vs. pred err std for high & LOW wager

uniqueStds = [10,20];
uniqueIDs = unique(ID);
% calculate prior weights vs TAC for the model
clear  beta_subj_highWager_pw beta_subj_lowWager_pw

for i0 = find(~sel_bet_subj)%1:length(uniqueIDs)
    
    if(use_day1)
        selID = probeSubjID == uniqueIDs(i0);
    else
        selID_lowWager = probeSubjID == uniqueIDs(i0) & probeDay1Sess' == 0 & probeWager == 0;
        selID_highWager = probeSubjID == uniqueIDs(i0) & probeDay1Sess' == 0 & probeWager == 1;
    end
    
    selSTD_lowWager = probeStd(selID_lowWager)';
    selSTD_highWager = probeStd(selID_highWager)';
    
    clear Bsubj stats
    subjPredErr = probeOutcomes' - probePrediction';
    subjPercErr = probeOutcomes' - probeEstimate';
    
    perceptualErr_lowWager = subjPercErr(selID_lowWager)-nanmean(subjPercErr(selID_lowWager));
    predictionErr_lowWager = subjPredErr(selID_lowWager) -nanmean(subjPredErr(selID_lowWager));
    
    perceptualErr_highWager = subjPercErr(selID_highWager)-nanmean(subjPercErr(selID_highWager));
    predictionErr_highWager = subjPredErr(selID_highWager) -nanmean(subjPredErr(selID_highWager));
    
    %create the decision matrix (high wager)
    X0_highWager = ones(length(predictionErr_highWager),1);
    
    % columns corresponding to b1 - b6
    Xpred_highWager = zeros(length(predictionErr_highWager),2*numTAC);
    for i4 = 1:numTAC
        Xpred_highWager(probeTAC(selID_highWager) == i4, i4) = 1;
        Xpred_highWager(probeTAC(selID_highWager) == i4, i4+6) = 1;
    end
    Xpred_highWager(:,1:6) = bsxfun(@times,Xpred_highWager(:,1:6),...
                               predictionErr_highWager.*(selSTD_highWager==10));
    Xpred_highWager(:,7:12) = bsxfun(@times,Xpred_highWager(:,7:12),...
                                predictionErr_highWager.*(selSTD_highWager==20));
    Xpred_highWager = bsxfun(@minus,Xpred_highWager,nanmean(Xpred_highWager,1));
    
    
    rawPred_highWager = probePrediction(selID_highWager)';
    
    % how far left of center are we?
    XSpace_highWager = (90 - rawPred_highWager);
    XSpace_highWager = (XSpace_highWager - nanmean(XSpace_highWager))./std(XSpace_highWager);
    
    
    Xsubj_highWager = [ones(size(Xpred_highWager,1),1) Xpred_highWager XSpace_highWager];
    
    %[Bsubj, stats] = robustfit(Xsubj,perceptualErr,[],[], 'off');
    
    [Bsubj_highWager, ~,~,LLsubj_highWager] = regressW_mike(perceptualErr_highWager, ...
        trialExpStd(selID_highWager),Xsubj_highWager);
    
    beta_subj_highWager_pw(:,i0) = Bsubj_highWager;
    
    
    
    %create the decision matrix (low wager)
    X0_lowWager = ones(length(predictionErr_lowWager),1);
    
    % columns corresponding to b1 - b6
    Xpred_lowWager = zeros(length(predictionErr_lowWager),2*numTAC);
    for i4 = 1:numTAC
        Xpred_lowWager(probeTAC(selID_lowWager) == i4, i4) = 1;
        Xpred_lowWager(probeTAC(selID_lowWager) == i4, i4+6) = 1;
    end
    Xpred_lowWager(:,1:6) = bsxfun(@times,Xpred_lowWager(:,1:6),...
                               predictionErr_lowWager.*(selSTD_lowWager==10));
    Xpred_lowWager(:,7:12) = bsxfun(@times,Xpred_lowWager(:,7:12),...
                                predictionErr_lowWager.*(selSTD_lowWager==20));
    Xpred_lowWager = bsxfun(@minus,Xpred_lowWager,nanmean(Xpred_lowWager,1));
    
    
    rawPred_lowWager = probePrediction(selID_lowWager)';
    
    % how far left of center are we?
    XSpace_lowWager = (90 - rawPred_lowWager);
    XSpace_lowWager = (XSpace_lowWager - nanmean(XSpace_lowWager))./std(XSpace_lowWager);
    
    
    Xsubj_lowWager = [ones(size(Xpred_lowWager,1),1) Xpred_lowWager XSpace_lowWager];
    
    %[Bsubj, stats] = robustfit(Xsubj,perceptualErr,[],[], 'off');
    
    [Bsubj_lowWager, ~,~,LLsubj_lowWager] = regressW_mike(perceptualErr_lowWager, ...
        trialExpStd(selID_lowWager),Xsubj_lowWager);
    
    beta_subj_lowWager_pw(:,i0) = Bsubj_lowWager;
    
    
end



subjPredErrSAC_highWager_hi = [];
subjPredErrSAC_highWager_low = [];
subjEstErrSAC_highWager_hi = [];
subjEstErrSAC_highWager_low = [];

subjPredErrSAC_lowWager_hi = [];
subjPredErrSAC_lowWager_low = [];
subjEstErrSAC_lowWager_hi = [];
subjEstErrSAC_lowWager_low = [];



uniqueStds = unique(fixStd);
numTAC = 6;




clear subj_est_err_disp subj_pred_err_disp


for i1 = find(~sel_bet_subj)
    
    subj_bet_freq(i1) = nanmean(probeWager(probeSubjID == i1));
    
    
    
    for i2 = 1:numTAC
        
        if(use_day1)
            selSAC_hi = probebCPTAC == i2 & probeSubjID == i1 & probeStd == 20;
            selSAC_low = probebCPTAC == i2 & probeSubjID == i1 & probeStd == 10;
        else
            selSAC_highWager_hi = probebCPTAC == i2 & probeWager == 1 &...
                 probeSubjID == i1 & probeStd == 20 & probeDay1Sess' == 0 ;
            selSAC_highWager_low = probebCPTAC == i2 & probeWager == 1 &... 
                probeSubjID == i1 & probeStd == 10 & probeDay1Sess' == 0 ;
            
             selSAC_lowWager_hi = probebCPTAC == i2 & probeWager == 0 &...
                 probeSubjID == i1 & probeStd == 20 & probeDay1Sess' == 0 ;
            selSAC_lowWager_low = probebCPTAC == i2 & probeWager == 0 &... 
                probeSubjID == i1 & probeStd == 10 & probeDay1Sess' == 0 ;
        end
        
        
        
        subjPredErrSAC_highWager_hi(i2,i1) = nanstd(subjPredErr(selSAC_highWager_hi));
        subjPredErrSAC_highWager_low(i2,i1) = nanstd(subjPredErr(selSAC_highWager_low));
        subjPredErrSAC_lowWager_hi(i2,i1) = nanstd(subjPredErr(selSAC_lowWager_hi));
        subjPredErrSAC_lowWager_low(i2,i1) = nanstd(subjPredErr(selSAC_lowWager_low));
        
        subjEstErrSAC_highWager_hi(i2,i1) = nanstd(subjEstErr(selSAC_highWager_hi));
        subjEstErrSAC_highWager_low(i2,i1) = nanstd(subjEstErr(selSAC_highWager_low));
        subjEstErrSAC_lowWager_hi(i2,i1) = nanstd(subjEstErr(selSAC_lowWager_hi));
        subjEstErrSAC_lowWager_low(i2,i1) = nanstd(subjEstErr(selSAC_lowWager_low));
        
    end
end










%% -------- PW(2-5) vs. Pred Err disp (2-5) (Low NOISE)

axes(axs(3)); cla(gca); hold on


aux2 = pw_sac_contrast_lo;
%aux2 = pw_noise_contrast;
aux1 = predErr_sac_contrast_lo;

%aux1 = predErrRel_sac_contrast;

[rho pVals] = corr(aux1, aux2,'type','spearman');
plot(aux1, aux2, 'o','markersize',8, 'markerfacecolor',cbColors(3,:),...
     'markeredgecolor',[0 0 0],'linewidth',1);

xlim([-10 5])
ylim([-0.16 0.25])


h = lsline;
h.LineWidth = 2;
h.LineStyle = '-';
h.Color = [0 0 0];

rho_str = ['RHO : ' num2str(rho,'%1.2f')];
p_str = ['p : ' num2str(pVals,'%1.3f')];
grid on

ylabel({'Perceptual bias Exp', 'contrast (low noise)'},'fontsize',ax_font_size1)
xlabel({'Prediction error' 'Exp contrast (low noise)'},'fontsize',ax_font_size1)
set(gca,'fontsize',ax_font_size1)
%legend('LOW NOISE')
%t2 = title(['RHO : ' rho_str '  p : ' p_str],'fontsize',ax_font_size2);

titlestr = ['$\mathbf{\rho}$ : ' num2str(rho,'%1.2f') ' p: ' num2str(pVals,'%1.3f')];
t2 =title(titlestr,'fontsize',14,'interpreter','latex');

%legend('low noise','Location','best')
set(gca,'YTick',[-0.1 0 0.2])





setPLOT_panelLabel(gca, 3,0.15,1.1);



%% -------- PW(2-5) vs. Pred Err disp (2-5) (HIGH NOISE)

axes(axs(4)); cla(gca); hold on


aux2 = pw_sac_contrast_hi;
%aux2 = pw_noise_contrast;
aux1 = predErr_sac_contrast_hi;

%aux1 = predErrRel_sac_contrast;

[rho pVals] = corr(aux1, aux2,'type','spearman');
plot(aux1, aux2, 'o','markersize',8, 'markerfacecolor',cbColors(2,:),...
     'markeredgecolor',[0 0 0],'linewidth',1);
xlim([-10 5])
ylim([-0.15 0.25])

% h = lsline;
% h.LineWidth = 2;
% h.LineStyle = '-';
% h.Color = [0 0 0];

rho_str = ['RHO : ' num2str(rho,'%1.2f')];
p_str = ['p : ' num2str(pVals,'%1.3f')];
grid on

ylabel({'Perceptual bias Exp', 'contrast (high noise)'},'fontsize',ax_font_size1)
xlabel({'Prediction error' 'Exp contrast (high noise)'},'fontsize',ax_font_size1)
set(gca,'fontsize',ax_font_size1)
%legend('LOW NOISE')
%t2 = title(['RHO : ' rho_str '  p : ' p_str],'fontsize',ax_font_size2);

titlestr = ['$\mathbf{\rho}$ : ' num2str(rho,'%1.2f') ' p: ' num2str(pVals,'%1.3f')];
t2 =title(titlestr,'fontsize',14,'interpreter','latex');

%legend('high noise','Location','best')
set(gca,'YTick',[-0.1 0 0.2])





setPLOT_panelLabel(gca, 4,0.15,1.1);



%% -------- PW CP contrast vs. Pred Err disp on SAC 2 

axes(axs(5)); cla(gca); hold on


aux2 = (pw_cp_contrast_hi + pw_cp_contrast_lo)*0.5;

aux1 = (predErr_cp_contrast_hi + predErr_cp_contrast_lo)*0.5;

aux1 = 1*mean(subjPredErrSAC_hi(2:end,:)' + subjPredErrSAC_low(2:end,:)',2);

%aux1 = predErrRel_cp_contrast;


%aux2 = pw_cp_contrast_combined;
%aux1 = predErr_cp_contrast;

[rho pVals] = corr(aux1, aux2,'type','spearman');
plot(aux1, aux2, 'o','markersize',8, 'markerfacecolor',marker_face_col,...
     'markeredgecolor',[0 0 0],'linewidth',1);
 
set(gca,'YTick',[-0.1 -0.05 0])

xlim([30 60])
ylim([-0.11 0])

h = lsline;
h.LineWidth = 2;
h.LineStyle = '-';
h.Color = [0 0 0];


rho_str = ['\rho : ' num2str(rho,'%1.2f')];
p_str = ['p : ' num2str(pVals,'%1.3f')];

grid on

ylabel({'Perceptual bias CP contrast'},'fontsize',ax_font_size1)
xlabel({'Prediction', 'error STD(SAC 2) (deg)'},'fontsize',ax_font_size1)
set(gca,'fontsize',ax_font_size1)
%legend('LOW NOISE')
titlestr = ['$\mathbf{\rho}$ : ' num2str(rho,'%1.2f') ' p: ' num2str(pVals,'%1.3f')];
t1 =title(titlestr,'fontsize',14,'interpreter','latex');






setPLOT_panelLabel(gca, 5,0.15,1.1);

%% -------- PW Noise contrast vs. Pred Err disp noise Contrast

axes(axs(6)); cla(gca); hold on


aux2 = pw_noise_contrast;
%aux2 = pw_noise_contrast;
aux1 = predErr_noise_contrast;


%aux1 = predErrRel_noise_contrast;

[rho pVals] = corr(aux1, aux2,'type','spearman');
plot(aux1, aux2, 'o','markersize',8, 'markerfacecolor',marker_face_col,...
     'markeredgecolor',[0 0 0],'linewidth',1);
 
ylim([-0.31 0.2])
xlim([4 20])


h = lsline;
h.LineWidth = 2;
h.LineStyle = '-';
h.Color = [0 0 0];


rho_str = ['RHO : ' num2str(rho,'%1.2f')];
p_str = ['p : ' num2str(pVals,'%1.3f')];
grid on

ylabel({'Perceptual bias', 'Noise contrast'},'fontsize',ax_font_size1)
xlabel({'Prediction error' 'noise contrast'},'fontsize',ax_font_size1)
set(gca,'fontsize',ax_font_size1)
%t3 = title(['RHO : ' rho_str '  p : ' p_str],'fontsize',ax_font_size2);

titlestr = ['$\mathbf{\rho}$ : ' num2str(rho,'%1.2f') ' p: ' num2str(pVals,'%1.3f')];
t3 =title(titlestr,'fontsize',14,'interpreter','latex');

set(gca,'YTick',[-0.3 -0.15  0 0.15])




setPLOT_panelLabel(gca, 6,0.15,1.1);








%%

kk=annotation('textbox');
set(kk, 'string', 'Krishnamurthy et al, figure 5', 'position', ...
    [0.75 0.92 0.15 0.05], 'EdgeColor', 'none','fontsize',14)


figName = '/Users/kamesh/Dropbox/auditoryTaskManuscript/Figure_set_for_writeup/figure_5.pdf';




