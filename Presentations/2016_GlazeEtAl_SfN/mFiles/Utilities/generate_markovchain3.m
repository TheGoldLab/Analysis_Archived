function z = generate_markovchain3(Hvc,N)

% x = zeros(N,1);
% z = double(rand>0);
% x(1) = z;
% 
% for n = 1:N-1
%     ch = rand<H(n);
%     
%     if ch
%         z = 1-z;
%     end
%     x(n+1) = z;
% end
% 

ch1 = rand(N,1)<Hvc;
ch1(1) = 1;

znum = sum(ch1);
zid = cumsum(ch1);
zunique = mod(1:znum,2);
if rand
    zunique = 1-zunique;
end

% cp1 = [1;find(ch1)];
z = zunique(zid);