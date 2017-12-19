clear
clc
close all

spcM = 30;

load model_testbed1.mat data choice x musgn

Hspc = linspace(.01,.99,spcM);

Mvc = [1000 20];

% FOR GENERATING DIFFERENT K
K1 = .01;
K2 = .001;
r01 = 2;
r02 = 2;
H0 = 0.5;



% FOR GENERATING DIFFERENT r0
K1 = .01;
K2 = .01;
r01 = 2;
r02 = exp(7);
H0 = 0.5;


l1 = normpdf(x,15,data.sigma);
l2 = normpdf(x,-15,data.sigma);

Hspc2 = linspace(.01,.99,100);
Hspc2 = Hspc2(:);
pH01 = betapdf(Hspc2,r01*H0,r01*(1-H0));
pH01 = pH01 / sum(pH01);

pH02 = betapdf(Hspc2,r02*H0,r02*(1-H0));
pH02 = pH02 / sum(pH02);

[q1exp1,expH_bayes1,pHmat_bayes1] = ...
    HHMM_mixed_mat_c(l1,l2,Hspc2,pH01,K1,musgn);

[q1exp2,expH_bayes2,pHmat_bayes2] = ...
    HHMM_mixed_mat_c(l1,l2,Hspc2,pH02,K2,musgn);

figure
subplot(2,3,1)
plot(1:2000,expH_bayes1,'k-','linewidth',2)

subplot(2,3,4)
plot(1:2000,expH_bayes2,'k-','linewidth',2)

H0samps1 = betarnd(H0*r01,(1-H0)*r01,50000,1);
H0samps2 = betarnd(H0*r02,(1-H0)*r02,50000,1);
F = data.sigma^2/30;

[q1samps,Hsamps1] = particle_filter_learnH(K1,H0samps1,F,x,musgn,Mvc(1),3,zeros(200,1));
[q1samps,Hsamps2] = particle_filter_learnH(K2,H0samps2,F,x,musgn,Mvc(1),3,zeros(200,1));

subplot(2,3,2)
plot(1:2000,Hsamps1')

subplot(2,3,5)
plot(1:2000,Hsamps2')


[q1samps,Hsamps1] = particle_filter_learnH(K1,H0samps1,F,x,musgn,Mvc(2),3,zeros(200,1));
[q1samps,Hsamps2] = particle_filter_learnH(K2,H0samps2,F,x,musgn,Mvc(2),3,zeros(200,1));

subplot(2,3,3)
plot(1:2000,Hsamps1')

subplot(2,3,6)
plot(1:2000,Hsamps2')

for i = 1:6
   subplot(2,3,i)
   set(gca,'ydir','normal','box','off','ticklength',[0 0],'fontsize',12,'xtick',[],'ytick',[])
   ylim([0 1])
%    hold on
%    plot(data.H(:,1),'k--')
    
end

set(gcf,'Position',[440   685   801   113])




