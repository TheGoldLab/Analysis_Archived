clear
clc
close all

spcM = 30;

load model_testbed1.mat data choice x musgn

Hspc = linspace(.01,.99,spcM);

Mvc = [1000 20];

K1 = .005;
K2 = .005;
r01 = 2;
r02 = 2;
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
imagesc(1:2000,Hspc,pHmat_bayes1)

subplot(2,3,4)
imagesc(1:2000,Hspc,pHmat_bayes2)

H0samps1 = betarnd(H0*r01,(1-H0)*r01,50000,1);
H0samps2 = betarnd(H0*r02,(1-H0)*r02,50000,1);
F = data.sigma^2/30;

[Lexp,Lvar,Hsampmat] = particle_filter_learnH6_mat(K1,H0samps1,F,x,musgn,Mvc(1),1,zeros(200,1));
pHmat11 = calc_pH_partfilt(Hsampmat,Hspc);

[Lexp,Lvar,Hsampmat] = particle_filter_learnH6_mat(K2,H0samps2,F,x,musgn,Mvc(1),1,zeros(200,1));
pHmat21 = calc_pH_partfilt(Hsampmat,Hspc);

subplot(2,3,2)
imagesc(1:2000,Hspc,pHmat11)

subplot(2,3,5)
imagesc(1:2000,Hspc,pHmat21)


[Lexp,Lvar,Hsampmat] = particle_filter_learnH6_mat(K1,H0samps1,F,x,musgn,Mvc(2),1,zeros(200,1));
pHmat12 = calc_pH_partfilt(Hsampmat,Hspc);

[Lexp,Lvar,Hsampmat] = particle_filter_learnH6_mat(K2,H0samps2,F,x,musgn,Mvc(2),1,zeros(200,1));
pHmat22 = calc_pH_partfilt(Hsampmat,Hspc);

subplot(2,3,3)
imagesc(1:2000,Hspc,pHmat12)

subplot(2,3,6)
imagesc(1:2000,Hspc,pHmat22)

for i = 1:6
   subplot(2,3,i)
   set(gca,'ydir','normal','box','off','ticklength',[0 0],'fontsize',12,'xtick',[],'ytick',[])
    
end

set(gcf,'Position',[440   685   801   113])



figure
subplot(3,1,1)
plot(data.H(:,1),'k-','linewidth',2)

subplot(3,1,2)
plot(data.muinds,'k-','linewidth',1)

subplot(3,1,3)
plot(x,'k-','linewidth',1)

for i = 1:3
   subplot(3,1,i)
   set(gca,'box','off','ticklength',[0 0],'xtick',[],'ytick',[])
end



set(gcf,'Position',[ 1099         685         241         113])


