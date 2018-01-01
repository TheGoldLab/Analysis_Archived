function [signalvc,rvc] = generate_cp(lambda,vclen,mu1,mu2)

cp = double(rand(vclen,1)<lambda);
cploc = find(cp);
cploc(1) = 1;
cploc = [cploc;vclen+1];

muvc = zeros(vclen,1);

mutmp = rand>.5;

rvc = muvc;

for cpind = 1:length(cploc)-1
    mutmp = 1 - mutmp;
    muvc(cploc(cpind):cploc(cpind+1)-1) = mutmp;
    
    rvc(cploc(cpind):cploc(cpind+1)-1) = 1:(cploc(cpind+1)-cploc(cpind));
end

signalvc = zeros(vclen,1);
signalvc(muvc==1)=mu1;
signalvc(muvc==0)=mu2;