clear
clc
close all

load bayesian_adaptive_sims1000.mat
alphamat1000 = squeeze(alphamat(1,:,:));

load bayesian_adaptive_sims3.mat

alphamat2 = zeros(4,size(alphamat,2),size(alphamat,3));
alphamat2(1:3,:,:) = alphamat;
alphamat2(4,:,:) = alphamat1000;
alphamat = alphamat2;

figure;
subplot(1,5,1)
imagesc(logr0spc,logKspc(1:end),alphamatb(:,1:end)',[0 .95]);

for i = 1:4
    subplot(1,5,i+1)
    imagesc(logr0spc,logKspc(1:end),squeeze(alphamat(5-i,:,1:end))');
end

for i = 1:5
    subplot(1,5,i)
    set(gca,'ydir','normal','box','off','ticklength',[0 0],'fontsize',12,'clim',[0 .95])
end

set(gcf,'Position',[324   683   750   115])

load bayesian_adaptive_sims_logKfixed.mat

figure; 
plot(logr0spc(1:2:end),alphamat(:,1:2:end)','linewidth',2)
hold on
plot(logr0spc(1:2:end),alphamatb(1:2:end),'k-','linewidth',2)

xlim([-1.5 8])
ylim([-.1 1.3])

% legend({'sample M=5','10','20','1000','ideal observer'},'fontsize',12)

set(gca,'box','off','ticklength',[0 0],'fontsize',12)

set(gcf,'Position',[440   660   221   138])
