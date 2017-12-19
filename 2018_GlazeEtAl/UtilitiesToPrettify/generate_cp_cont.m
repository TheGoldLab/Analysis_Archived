function [muvc,rvc] = generate_cp_cont(lambda,vclen,mu0,sigma0,mumin,mumax)

cp = double(rand(vclen,1)<lambda);
cploc = find(cp);
cploc(1) = 1;
cploc = [cploc;vclen+1];

muvc = zeros(vclen,1);

rvc = muvc;

if nargin==4
    muspc = normrnd(mu0,sigma0,length(cploc)-1,1);
else
    muspc = mumin + (mumax-mumin)*rand(length(cploc)-1,1);
end

for cpind = 1:length(cploc)-1
    mutmp = muspc(cpind);
    muvc(cploc(cpind):cploc(cpind+1)-1) = mutmp;
    
    rvc(cploc(cpind):cploc(cpind+1)-1) = 1:(cploc(cpind+1)-cploc(cpind));
end