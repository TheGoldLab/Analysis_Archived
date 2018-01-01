clear
clc
close all

dirnm = pwd;
data_dir = [dirnm(1:strfind(dirnm,'Dropbox')+19) filesep 'fit_data' filesep];


load([data_dir 'allparam_fixedvarfits_fixedLRapprox_sigmafixed.mat'],'paramstrct_fixedLRapprox_sigmafixed');

load([data_dir 'allparam_fixedvarfits_fixedLR_sigmafixed.mat'],'paramstrct_fixedLR_sigmafixed');
load([data_dir 'allparam_fixedvarfits_partfilt_sigma_partMfixed20.mat'],'paramstrct_sigma_partMfixed');
load([data_dir 'allparam_fixedvarfits_partfilt_sigma_partMfixed5.mat'],'paramstrct_sigma_partMfixed5');
load([data_dir 'allparam_fixedvarfits_partfilt_sigma_partMfixed10.mat'],'paramstrct_sigma_partMfixed10');

load([data_dir 'allparam_fixedvarfits_bayesultimate_sigmafixed.mat'],'paramstrct_bayesultimate');

paramstrct_bayesultimate_sigmafixed = paramstrct_bayesultimate;

load([data_dir 'allparam_null.mat'],'paramstrct_null')

LLnull = [paramstrct_null.LL]';
LLfixedLRapprox_sigmafixed = [paramstrct_fixedLRapprox_sigmafixed.LL]';
LLfixedLR_sigmafixed = [paramstrct_fixedLR_sigmafixed.LL]';
LLbayes_sigmafixed = [paramstrct_bayesultimate_sigmafixed.LL]';
LL20_sigma_partMfixed = [paramstrct_sigma_partMfixed.LLrl]';
LL10_sigma_partMfixed = [paramstrct_sigma_partMfixed10.LLrl]';
LL5_sigma_partMfixed = [paramstrct_sigma_partMfixed5.LLrl]';

[substrct,sessionNvc,sigmamedvc,Hlowvc,Hhighvc,signaledmedvc,signaledfirstvc] = get_subjstats;
indsval = (sessionNvc>1 & signaledfirstvc>0);
sessionN = sessionNvc(indsval);

BICnull = -2*LLnull;
BIC5_sigmafixed = -2*LL5_sigma_partMfixed+5*log(sessionN*2000);
BIC10_sigmafixed = -2*LL10_sigma_partMfixed+5*log(sessionN*2000);
BIC20_sigmafixed = -2*LL20_sigma_partMfixed+5*log(sessionN*2000);
BICbayes_sigmafixed = -2*LLbayes_sigmafixed+5*log(sessionN*2000);
BICfixedLR_sigmafixed = -2*LLfixedLR_sigmafixed+4*log(sessionN*2000);
BICfixedLRapprox_sigmafixed = -2*LLfixedLRapprox_sigmafixed+4*log(sessionN*2000);

BICmat = [BIC5_sigmafixed BIC10_sigmafixed BIC20_sigmafixed BICbayes_sigmafixed BICfixedLR_sigmafixed BICfixedLRapprox_sigmafixed BICnull];

BICdiff = BICmat - repmat(BICmat(:,3),1,size(BICmat,2));
BICdiffmed = mean(BICdiff);
bootsamps =  bootstrp(1000,@mean,BICdiff);
BICdiffsem1 = prctile(bootsamps,25);
BICdiffsem2 = prctile(bootsamps,75);
ub = BICdiffsem2-BICdiffmed;
lb = BICdiffmed-BICdiffsem1;

BICdiffsem = std(bootsamps);

pvc = zeros(7,1);
for n = 1:7
   pvc(n) = signrank(BICdiff(:,n)); 
end

modN = 5;

figure
bar(BICdiffmed(:,1:modN),'facecolor','w')
hold on

% ISSUE HERE IS THAT SAMPLED DISTRIBUTIONS OF MEDIANS ARE VERY SKEWED POSITIVELY, SO
% PLOTTING SEM EXAGGERATES ERROR ON THE DIFFERENCE SIGN
% errorbar(1:modN,BICdiffmed(:,1:modN),BICdiffsem(:,1:modN),'k.','linewidth',2)
errorbar(1:modN,BICdiffmed(:,1:modN),BICdiffsem(:,1:modN),'k.','linewidth',2)
plot([1 modN],[0 0],'k--')

set(gca,'box','off','ticklength',[0 0],'fontsize',12,'xtick',[])
xlim([0 modN+1])
set(gcf,'Position',[440   640   262   158])

BICstats1 = sprintf('M=5 diff median= %0.3g +/- %0.3g, p= %0.4g \nM=10 diff median= %0.3g +/- %0.3g, p= %0.4g \nideal diff median= %0.3g +/- %0.3g, p= %0.4g\nfixedLR diff median= %0.3g +/- %0.3g, p= %0.4g\nnull diff median= %0.3g +/- %0.3g, p= %0.4g',...
    BICdiffmed(1), BICdiffsem(1), pvc(1),...
    BICdiffmed(2), BICdiffsem(2), pvc(2),...
    BICdiffmed(4), BICdiffsem(4), pvc(4),...
    BICdiffmed(5), BICdiffsem(5), pvc(5),...
    BICdiffmed(6), BICdiffsem(6), pvc(6),...
    BICdiffmed(7), BICdiffsem(7), pvc(7))
