function [LL,LLvc,LP] = getfit_fixedvar_gauss_partfilt_arr(xarr,choicearr,musgnarr,sigmavc,mu,simN,params,eta0)

simN = ceil(simN);

H0_final = params(7);
r0_final = exp(params(6));
K_final = exp(params(8));
sigma0_final = params(4);
gamma_final = params(5);
psi_final = params(2);
partM = ceil(params(1));
whoops = params(3);
% 
% logr0_final = log(r0_final);

H0samps = betarnd(r0_final*H0_final,r0_final*(1-H0_final),10000,1);
% H0samps = .001 + .998*H0samps;
H0samps = max(min(H0samps,.999),.001);

% logJ0_final = log(H0_final) - log(1-H0_final);
% logJ0samps = normrnd(logJ0_final,1/r0_final,10000,1);
% H0samps = 1./(1+exp(-logJ0samps));

eta = psi_final*eta0;

arrN = numel(xarr);
LLvc = zeros(arrN,1);

parfor arrind = 1:arrN
    
    sigma = sigmavc(arrind);
    x = xarr{arrind};
    choice = choicearr{arrind};
    musgn = musgnarr{arrind};
    
    sigmasubj = sigma0_final + gamma_final*(sigma-sigma0_final);
    F = sigmasubj^2 / mu;
   
    [Lexp,Lvar] = particle_filter_learnH6(K_final,H0samps,F,x,musgn,partM,simN,eta);
    
    Lvar = max(Lvar,.01);

    choice1 = .5 + .5*erf(Lexp./sqrt(2*Lvar));
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
