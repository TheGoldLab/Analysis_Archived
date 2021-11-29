clear
clc
close all

% CODE USED TO GENERATE STIMULUS SEQUENCES

Hlen = [400 200 400];
Hid = [.01 .3 .01];
lid = [.8 1 .8];


Hlen = [200 500 100 200];
Hid = [.99 .01 .99 .01];
lid = [.8 .8 .8 .8];

N = sum(Hlen);
blockN = numel(Hlen);
Hvc = zeros(N,1);
lvc = Hvc;

k = 1;

for i = 1:blockN
    kinds = k:k+Hlen(i)-1;
    Hvc(kinds) = Hid(i);
    lvc(kinds) = lid(i);
    k = k + Hlen(i);
end

z = generate_markovchain3(Hvc,N);
z = z + 1;

stimvc = zeros(N,1);
stimvc(z==1) = rand(sum(z==1),1)>lvc(z==1);
stimvc(z==2) = rand(sum(z==2),1)<lvc(z==2);

stimvc = stimvc+1;

figure
subplot(3,1,1)
plot(Hvc)

subplot(3,1,2)
plot(z)

subplot(3,1,3)
plot(stimvc)


stimvc1 = stimvc(1:Hlen(1));
mean(stimvc1(2:end)~=stimvc1(1:end-1))
stimvc2 = stimvc(Hlen(1)+1:Hlen(1)+Hlen(2));
mean(stimvc2(2:end)~=stimvc2(1:end-1))
stimvc3 = stimvc(Hlen(1)+Hlen(2)+1:end);
mean(stimvc3(2:end)~=stimvc3(1:end-1))
