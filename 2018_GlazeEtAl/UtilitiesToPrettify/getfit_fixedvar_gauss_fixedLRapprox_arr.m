function [LL,LLvc,LP] = getfit_fixedvar_gauss_fixedLRapprox_arr(xarr,choicearr,musgnarr,sigmavc,mu,params,eta)

sigma0 = params(5);
gamma = params(6);
psi = params(3);
whoops = params(4);
logalpha = params(1);
H0 = params(2);

alpha = exp(logalpha);

arrN = numel(xarr);
LLvc = zeros(arrN,1);

parfor arrind = 1:arrN
    
    sigma = sigmavc(arrind);
    x = xarr{arrind};
    choice = choicearr{arrind};
    musgn = musgnarr{arrind};
    
    sigmasubj = sigma0 + gamma*(sigma-sigma0);
    F = sigmasubj^2 / mu;
    
    psirl = psi/F;
   
    Lexp = fixedLR_learnH6_approx(alpha,H0,F,x,musgn,1,eta);
  
    choice1 = 1./(1+exp(-Lexp/psirl));
    choice1 = (1-whoops)*choice1 + whoops*(1-choice1);
    LL = sum(choice.*log(choice1)+(1-choice).*log(1-choice1));
    
    LLvc(arrind) = LL;
    
end

LL = sum(LLvc);
LP = 0;
% 
% LP = lognormpdf(log(r0_final),2,1.5) + ...
%     logbetapdf(K_final,1,10) + ...
%     logbetapdf(H0_final,2,2.5);

% LP = lognormpdf(logr0_final,log(2),1) + ...
%     logbetapdf(K_final,5/4+1,500/4) + ...
%     logbetapdf(H0_final,10,15) + ...
%     logbetapdf(gamma_final,20,20) + ...
%     loggampdf(sigma0_final,21,10/20) + ...
%     logbetapdf(whoops,1,100);

% LP = lognormpdf(logr0_final,log(2),2) + ...
%     logbetapdf(K_final,5/4+1,100/4) + ...
%     logbetapdf(H0_final,2,2.5) + ...
%     logbetapdf(gamma_final,2,3) + ...
%     loggampdf(sigma0_final,11,1);

LL = LL + LP;
