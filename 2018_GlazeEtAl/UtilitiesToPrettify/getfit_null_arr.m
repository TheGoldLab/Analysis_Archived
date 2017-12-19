function [LL,LLvc,LP] = getfit_null_arr(xarr,choicearr,musgnarr,sigmavc,mu,params)

psi = params(1);
whoops = params(2);

arrN = numel(xarr);
LLvc = zeros(arrN,1);

parfor arrind = 1:arrN
    
    sigma = sigmavc(arrind);
    x = xarr{arrind};
    choice = choicearr{arrind};

    choice1 = .5 + .5*erf(x/(sqrt(2)*psi));
    
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
