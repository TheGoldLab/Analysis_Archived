clear
clc
close all

K1 = .01;
K2 = .001;
N = 2;

figure
subplot(2,1,1)
[muvc1,rvc1] = generate_cp_cont(K1,2000,.5,.5,.01,.99);
[muvc2,rvc2] = generate_cp_cont(K1,2000,.5,.5,.01,.99);
plot(muvc1,'linewidth',2)
hold on
plot(muvc2,'linewidth',2);

subplot(2,1,2)
[muvc1,rvc1] = generate_cp_cont(K2,2000,.5,.5,.01,.99);
[muvc2,rvc2] = generate_cp_cont(K2,2000,.5,.5,.01,.99);
plot(muvc1,'linewidth',2)
hold on
plot(muvc2,'linewidth',2);

for i = 1:2
    subplot(2,1,i)
    set(gca,'box','off','ticklength',[0 0],'fontsize',12,'xtick',[],'ytick',[])
end
set(gcf,'Position',[440   571   331   227])

