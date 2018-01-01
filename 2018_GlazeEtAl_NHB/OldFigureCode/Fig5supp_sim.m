clear
clc
close all

load regression_sim_Kvary.mat
params_regress = [paramstrct_regress.params]';
params_regress = reshape(params_regress,4,41)';

figure
plot(params_regress(:,2),params_regress(:,3),'ko')

set(gca,'ydir','normal','box','off','ticklength',[0 0],'fontsize',12)
% xlim([0 .7])
% ylim([-.1 1.1])
