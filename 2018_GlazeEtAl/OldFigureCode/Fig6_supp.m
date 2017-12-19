clear
clc
close all


load ../fit_data/allparam_fixedvarfits_partfilt_sigma_partMfixed_sim.mat
paramsrl = [psivc(:) whoopsvc(:) logr0vc(:) H0vc(:) logKvc(:)];
params = [paramstrct_sigma_partMfixed.params]';
params = params(:,[2 3 6:8]);

diag(corr(params,paramsrl(1:size(params,1),:),'type','Spearman'))

figure;
subplot(1,5,1)
plot(paramsrl(1:size(params,1),1),params(:,1),'ko')
hold on
plot([0 5],[0 5],'k--')
xlabel('simulated sensory noise','fontsize',14)
ylabel('fit sensory noise','fontsize',14)

subplot(1,5,2)
plot(paramsrl(1:size(params,1),2),params(:,2),'ko')
hold on
plot([0 .01],[0 .01],'k--')
xlabel('simulated lapse rate','fontsize',14)
ylabel('fit lapse rate','fontsize',14)

subplot(1,5,3)
plot(paramsrl(1:size(params,1),end-2),params(:,end-2),'ko')
hold on
plot([-1 6],[-1 6],'k--')
xlabel('simulated log(r_0)','fontsize',14)
ylabel('fit log(r_0)','fontsize',14)
xlim([-2 9])
ylim([-2 9])

subplot(1,5,4)
plot(paramsrl(1:size(params,1),end-1),params(:,end-1),'ko')
hold on
plot([0 1],[0 1],'k--')
xlabel('simulated H_0','fontsize',14)
ylabel('fit H_0','fontsize',14)


subplot(1,5,5)
plot(paramsrl(1:size(params,1),end),params(:,end),'ko')
hold on
plot([-9 -2],[-9 -2],'k--')
xlabel('simulated log(K)','fontsize',14)
ylabel('fit log(K)','fontsize',14)
% xlim([-10 -2])
% ylim([-10 -2])

for i = 1:5
    subplot(1,5,i)
    set(gca,'ydir','normal','box','off','ticklength',[0 0],'fontsize',10)
end

set(gcf,'Position',[434   694   836   111])