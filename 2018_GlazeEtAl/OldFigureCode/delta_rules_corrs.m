clear
clc
close all


load regression_sim_r0vary.mat paramstrct_regress

params_regress = [paramstrct_regress.params]';
params_regress = reshape(params_regress,4,41)';

figure
plot(params_regress(:,2),params_regress(:,3),'o','markersize',8)
hold on

[r,p] = corr(params_regress(:,2:3),'type','Spearman')



load regression_sim_gammavary_delta_rule_normative.mat paramstrct_regress
params_regress = [paramstrct_regress.params]';
params_regress = reshape(params_regress,4,41)';


[r,p] = corr(params_regress(:,2:3),'type','Spearman')

plot(params_regress(:,2),params_regress(:,3),'ko','markersize',8)

% load regression_sim_gammavary_delta_rule_leakyint.mat paramstrct_regress
% params_regress = [paramstrct_regress.params]';
% params_regress = reshape(params_regress,4,41)';
% 
% [r,p] = corr(params_regress(:,2:3),'type','Spearman')
% 
% plot(params_regress(:,2),params_regress(:,3),'ro')
% 

set(gca,'ydir','normal','box','off','ticklength',[0 0],'fontsize',12)

set(gcf,'Position',[440   615   261   183])

% xlim([-.1 1.1])
% ylim([0 .6])



% 
% figure
% subplot(3,1,1)
% plot(data.H(:,1),'k-','linewidth',2)
% 
% subplot(3,1,2)
% plot(data.muinds,'k-','linewidth',1)
% 
% subplot(3,1,3)
% plot(x,'k-','linewidth',1)
% 
% for i = 1:3
%    subplot(3,1,i)
%    set(gca,'box','off','ticklength',[0 0],'xtick',[],'ytick',[])
% end
% 
% 
% 
% set(gcf,'Position',[440   685   900   113])


