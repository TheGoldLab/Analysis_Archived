function paramstrct = fit_null_multsession(xarr,choicearr,musgnarr,sigmavc,mu)

params0 = [1 .002];
ub = [10 .2];
lb = [0 0];

model = @fit_dsprt_min;


opts = optimoptions(@fmincon,'Algorithm','interior-point');

% Multistart optimization object runs the fitting from multiple starting
% points (params0 varied)
ms = MultiStart('UseParallel','always','Display','iter','StartPointsToRun','bounds');

problem = createOptimProblem('fmincon','objective',...
    model,'x0',params0,'lb',lb,'ub',ub,'options',opts);

[params,fval] = run(ms,problem,4);
LLfinal = -fval;

paramstrct.params = params;
paramstrct.LL = LLfinal;

    function e = fit_dsprt_min(params)
 
        LL = getfit_null_arr(xarr,choicearr,musgnarr,sigmavc,mu,params);
        e = -1*LL;
        
    end

end