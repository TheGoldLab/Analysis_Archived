clear
clc
close all


load ../fit_data/allparam_fixedvarfits_regression_crossval.mat alphamat H0mat psimat
load ../fit_data/allparam_fixedvarfits_regression.mat
load ../fit_data/allparam_fixedvarfits_partfilt_sigma_partMfixed20.mat 

% 
% Fig 7 stats:
% -H0 correlation between regression and sampling model
% -correlation b/w r0 and adaptability, choice noise 
% -correlation b/w r0 across partitions (and do for K?)
% 




load ../fit_data/allparam_fixedvarfits_partfilt_sigma_partMfixed_crossval.mat
params11 = [paramstrct_sigma_partMfixed_crossval(:,1,1).params]';
params12 = [paramstrct_sigma_partMfixed_crossval(:,1,2).params]';
params21 = [paramstrct_sigma_partMfixed_crossval(:,2,1).params]';
params22 = [paramstrct_sigma_partMfixed_crossval(:,2,2).params]';
[rKcross,pKcross] = corr(params22(:,end),params21(:,end),'type','Spearman','tail','right');
[rr0cross,pr0cross] = corr(params12(:,end-2),params11(:,end-2),'type','Spearman','tail','right');


paramat = [paramstrct_sigma_partMfixed.params]';

paramat_regress = [paramstrct_regress.params]';
paramat_regress = reshape(paramat_regress,4,41)';

alpha = paramat_regress(:,2);
psi = paramat_regress(:,3);
H0r = paramat_regress(:,1);
H0r = 1./(1+exp(-H0r));

[rH0H0,pH0H0] = corr(paramat(:,7),H0r,'type','Spearman');

[rnoiser0,pnoiser0] = corr(paramat(:,6),psi,'type','Spearman');
[ralphar0,palphar0] = corr(paramat(:,6),alpha,'type','Spearman');

[rnoiseK,pnoiseK] = corr(paramat(:,end),psi,'type','Spearman');
[ralphaK,palphaK] = corr(paramat(:,end),alpha,'type','Spearman');

H0stats1 = sprintf('H0-H0 corr= %0.3g, p=%0.4g', rH0H0, pH0H0)

r0stats1 = sprintf('r0-alpha corr= %0.3g, p=%0.4g \nr0-psi corr= %0.3g, p=%0.4g', ralphar0, palphar0, rnoiser0, pnoiser0)
r0stats2 = sprintf('r0 stab corr= %0.3g, p=%0.4g', rr0cross, pr0cross)

Kstats1 = sprintf('K-alpha corr= %0.3g, p=%0.4g \nK-psi corr= %0.3g, p=%0.4g', ralphaK, palphaK, rnoiseK, pnoiseK)
Kstats2 = sprintf('K stab corr= %0.3g, p=%0.4g', rKcross, pKcross)

[rnoise_xnoise,pnoise_xnoise] = corr(paramat(:,2),psi,'type','Spearman');
xnoisestats = sprintf('noise-xnoise corr= %0.3g, p=%0.4g', rnoise_xnoise,pnoise_xnoise)
