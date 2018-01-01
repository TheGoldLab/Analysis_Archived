clear
clc
close all

dirnm = pwd;
data_dir = [dirnm(1:strfind(dirnm,'Dropbox')+19) filesep 'fit_data' filesep];

load([data_dir 'allparam_fixedvarfits_regression_crossval.mat'],'alphamat','H0mat','psimat')
load([data_dir 'allparam_fixedvarfits_regression.mat'],'paramstrct_regress')

[correctvc,RTvc,subjvc,Hvc,rvc,musgnvc,evc,initialevc,blockindvc,isprevc,sessionindsvc,paramat,alphavc,datacorrectvc] = get_modelfree_data;

pcorrectnorm = grpstats(correctvc-datacorrectvc,subjvc);

paramat_regress = [paramstrct_regress.params]';
paramat_regress = reshape(paramat_regress,4,41)';

alpha = paramat_regress(:,2);
psi = paramat_regress(:,3);
H0 = paramat_regress(:,1);
H0 = 1./(1+exp(-H0));

[rcorrect,pcorrect] = corr(alpha,pcorrectnorm,'type','Spearman');
[rcorrectH0,pcorrectH0] = corr(H0,pcorrectnorm,'type','Spearman');

alphamed = median(alpha);
alphasem = std(bootstrp(500,@median,alpha));
alphamin = min(alpha);
alphamax = max(alpha);
alphap = signrank(alpha);

H0med = median(H0);
H0sem = std(bootstrp(500,@median,H0));
H0min = min(H0);
H0max = max(H0);

[rstab,pstab] = corr(alphamat(:,1),alphamat(:,2),'type','Spearman');
[rstabH0,pstabH0] = corr(H0mat(:,1),H0mat(:,2),'type','Spearman');

alphaminsem = std(bootstrp(500,@min,alpha));
alphamaxsem = std(bootstrp(500,@max,alpha));

alphastats1 = sprintf('alpha median= %0.3g +/- %0.3g \nalphap= %0.4g \nalphamin=%0.3g, alphamax=%0.3g', alphamed, alphasem, alphap, alphamin, alphamax)
alphastats2 = sprintf('alpha correct corr= %0.3g, p= %0.4g', rcorrect, pcorrect)
alphastats3 = sprintf('alpha stability corr= %0.3g, p= %0.4g', rstab, pstab)

H0stats1 = sprintf('H0 median= %0.3g +/- %0.3g \nH0min=%0.3g, H0max=%0.3g', H0med, H0sem, H0min, H0max)
H0stats2 = sprintf('H0 correct corr= %0.3g, p= %0.4g', rcorrectH0, pcorrectH0)
H0stats3 = sprintf('H0 stability corr= %0.3g, p= %0.4g', rstabH0, pstabH0)


% PARAMETER DISTRIBUTION PLOTS
figure
[nout,x] = hist(alpha,-.2:.1:1.5);
bar(x,nout,'facecolor','k')
set(gcf,'Position',[440   672   271   126])
xlim([-.3 1.2])
set(gca,'box','off','ticklength',[0 0],'fontsize',12)

figure
[nout,x] = hist(H0,0:.1:1);
bar(x,nout,'facecolor','k')
set(gcf,'Position',[440   672   271   126])
xlim([-.1 1.1])
set(gca,'box','off','ticklength',[0 0],'fontsize',12)


% STABILITY PLOTS
figure
plot(alphamat(:,1),alphamat(:,2),'ko')
set(gca,'box','off','ticklength',[0 0],'xtick',[],'ytick',[])
xlim([-.1 1.1])
ylim([-.1 1.1])
set(gcf,'Position',[440   702   128    96])

figure
plot(H0mat(:,1),H0mat(:,2),'ko')
set(gca,'box','off','ticklength',[0 0],'xtick',[],'ytick',[])
xlim([0 1])
ylim([0 1])
set(gcf,'Position',[440   702   128    96])


load([data_dir 'allparam_fixedvarfits_regression_crossval_nonrand.mat'],'alphamat','H0mat','psimat')

alphachangemed = median(alphamat(:,2)-alphamat(:,1));
alphachangesem = std(bootstrp(500,@median,alphamat(:,2)-alphamat(:,1)));

alphachangep = signrank(alphamat(:,2)-alphamat(:,1));

H0changemed = median(H0mat(:,2)-H0mat(:,1));
H0changesem = std(bootstrp(500,@median,H0mat(:,2)-H0mat(:,1)));

H0changep = signrank(H0mat(:,2)-H0mat(:,1));

alphachangestats1 = sprintf('alpha change median= %0.3g +/- %0.3g \nalphap= %0.4g', alphachangemed, alphachangesem, alphachangep)
H0changestats1 = sprintf('H0 change median= %0.3g +/- %0.3g \nalphap= %0.4g', H0changemed, H0changesem, H0changep)

