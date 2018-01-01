function params = fit_dsprt_blocked(x,choice,musgn,blockinds,blockinds_sigma)

whoops = .002;

logdiffSq1 = -1*(x-15).^2;
logdiffSq2 = -1*(x+15).^2;

blockU = unique(blockinds);
blockUsigma = unique(blockinds_sigma);
blockM = numel(blockU);
blockMsigma = numel(blockUsigma);

ub = [ones(1,blockM) 20*ones(1,blockMsigma) 10];
params0 = [.4*ones(1,blockM) 10*ones(1,blockMsigma) .5];
lb = [zeros(1,blockM) 2*ones(1,blockMsigma) 0];

model = @fit_dsprt_min;
%
% if matlabpool('size') == 0
%     matlabpool open
% end

% opts = optimoptions(@fmincon,'Algorithm','interior-point');
% ms = MultiStart('UseParallel','always');
% 
% problem = createOptimProblem('fmincon','objective',...
%     model,'x0',params0,'lb',lb,'ub',ub,'options',opts);
% 
% [params,fval] = run(ms,problem,16);


% options = psoptimset('maxiter',5000,'Display','none','UseParallel',true,'CompleteSearch', 'on','CompletePoll', 'on');

options = psoptimset('maxiter',5000,'Display','none','UseParallel',false,'CompleteSearch', 'off','CompletePoll', 'off');

% options = optimoptions('patternsearch','UseParallel', true, 'UseCompletePoll', false, 'UseVectorized', false);

[params] = patternsearch(model,params0,[],[],[],[],lb,ub,[],options);

    function e = fit_dsprt_min(params)
        
        Hset = params(1:blockM);
        sigmaset = params((blockM+1):(blockM+blockMsigma));
        psi = params(end);
        
        sigmavc = sigmaset(blockinds_sigma)';
        
%         l1 = normpdf(x,15,sigmavc);
%         l2 = normpdf(x,-15,sigmavc);
        
        l1 = exp(logdiffSq1./(2*sigmavc.^2));
        l2 = exp(logdiffSq2./(2*sigmavc.^2));
        
        Hvc = Hset(blockinds)';
        
        [q1,q2] = dsprt_mixed_Hvc2(l1,l2,Hvc,musgn);
        L = log(q1)-log(q2);
        
    %    choicehat = .5 + .5*erf(L/(sqrt(2)*psi));
        choicehat = 1./(1+exp(-L/psi));
        choicehat = (1-whoops)*choicehat+whoops*(1-choicehat);
 
        e = -sum((1-choice).*log(1-choicehat)) - sum(choice.*log(choicehat));
        
%         prec = 1/psi^2;
        
%         logp_psi = sum(loggampdf(prec,1.05,20));
%         e = e - logp_psi;

        logJ = log((1-Hset)./Hset);
        dlogJ = diff(logJ);
        
        logp_dJ = sum(lognormpdf(dlogJ,0,2));
        
        e = e - logp_dJ;
       
        dsigma = diff(sigmaset);
        
        logp_dsigma = sum(lognormpdf(dsigma,0,2));
        
        e = e - logp_dsigma;
        
        
    end

end