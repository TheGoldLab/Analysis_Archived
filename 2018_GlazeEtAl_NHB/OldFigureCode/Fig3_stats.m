clear
clc
close all

dirnm = pwd;
data_dir = [dirnm(1:strfind(dirnm,'Dropbox')+19) filesep 'fit_data' filesep];

load([data_dir 'allparam_fixedvarfits_regression_crossval.mat'],'alphamat','H0mat','psimat')
load([data_dir 'allparam_fixedvarfits_regression.mat'])

[correctvc,RTvc,subjvc,Hvc,rvc,musgnvc,evc,initialevc,blockindvc,isprevc,sessionindsvc,paramat,alphavc,datacorrectvc] = get_modelfree_data;

pcorrectnorm = grpstats(correctvc-datacorrectvc,subjvc);

paramat_regress = [paramstrct_regress.params]';
paramat_regress = reshape(paramat_regress,4,41)';

alpha = paramat_regress(:,2);
psi = paramat_regress(:,3);
H0 = paramat_regress(:,1);
H0 = 1./(1+exp(-H0));

[rnoise,pnoise] = corr(alpha,psi,'type','Spearman');
[rstab,pstab] = corr(psimat(:,1),psimat(:,2),'type','Spearman');

psistats = sprintf('noise corr= %0.3g, p=%0.4g', rnoise, pnoise)
psistats2 = sprintf('noise stab corr= %0.3g, p=%0.4g', rstab, pstab)


figure
plot(alpha,psi,'ko')
set(gca,'box','off','ticklength',[0 0],'xtick',[],'ytick',[])
xlim([-.1 1.1])
ylim([0 1.5])
set(gcf,'Position',[440   702   128    96])

figure
plot(psimat(:,1),psimat(:,2),'ko')
set(gca,'box','off','ticklength',[0 0],'xtick',[],'ytick',[])
xlim([0 1.5])
ylim([0 1.5])
set(gcf,'Position',[440   702   128    96])

