function [pHmat,varH,expH] = calc_pH_bayesian(prmat,Hmat,r0,Hspc)

M = numel(Hspc);
Hspc = Hspc(:);
Hspc2 = Hspc(:);
[R,N] = size(prmat);
pHmat = zeros(M,N);

Hmatorg = repmat(Hspc,1,N);
varH = zeros(N,1);
expH = varH;

for n = 1:N
 %   Hmatmp = repmat(Hspc,1,n);
 
    indsval = find(prmat(:,n)>eps);
    
    Hmatmp = Hmatorg(:,indsval);
    
    rmatmp = repmat(indsval',M,1)-1;
    Hexptmp = repmat(Hmat(indsval,n)',M,1);
    
    alphamat = (r0+rmatmp).*Hexptmp;
    betamat = (r0+rmatmp).*(1-Hexptmp);
    
    pH = betapdf(Hmatmp,alphamat,betamat);
    
    pHexp = pH*prmat(indsval,n);
    
    pHmat(:,n) = pHexp;
    
    pHexp = pHexp/sum(pHexp);
    
    expH(n) = sum(pHexp.*Hspc);
    
    varH(n) = sum(pHexp.*Hspc2)-expH(n)^2;
    
    
    
end