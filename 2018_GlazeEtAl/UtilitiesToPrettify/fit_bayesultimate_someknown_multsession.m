function [paramstrct] = fit_bayesultimate_someknown_multsession(xarr,choicearr,musgnarr,sigmavc,mu,fixedparamstrct,Hspc)


ub = fixedparamstrct.ub(fixedparamstrct.free);
lb = fixedparamstrct.lb(fixedparamstrct.free);
params0 = fixedparamstrct.params0(fixedparamstrct.free);

params20 = zeros(numel(fixedparamstrct.free),1);
params20(~fixedparamstrct.free) = fixedparamstrct.params0(~fixedparamstrct.free);

model = @fit_dsprt_min;

% options = psoptimset('maxiter',5000,'Display','iter','CompletePoll','on','TolMesh',1e-4,'CompleteSearch','on','UseParallel',true,'ScaleMesh','on','InitialMeshSize',1,'SearchMethod',@searchlhs);

% options = psoptimset('maxiter',5000,'Display','iter','CompletePoll','on','TolMesh',1e-6,'CompleteSearch','on','UseParallel',true,'ScaleMesh','on','InitialMeshSize',1);
% [params,fval] = patternsearch(model,params0,[],[],[],[],lb,ub,[],options);

% options = saoptimset('maxiter',5000,'Display','iter','StallIterLimit',500,'AnnealingFcn',@annealingboltz,'TemperatureFcn',@temperatureboltz,'ReannealInterval',100);
% [params,fval] = simulannealbnd(model,params0,lb,ub,options);
% 
% 
% startpts = [params0;ub-.01;lb+.01;.5*(ub+lb)];
% startpts = CustomStartPointSet(startpts);

opts = optimoptions(@fmincon,'Algorithm','interior-point');

% Multistart optimization object runs the fitting from multiple starting
% points (params0 varied)
ms = MultiStart('UseParallel','always','Display','iter','StartPointsToRun','bounds');

problem = createOptimProblem('fmincon','objective',...
    model,'x0',params0,'lb',lb,'ub',ub,'options',opts);

[params,fval] = run(ms,problem,4);

paramsfinal = params20;
paramsfinal(fixedparamstrct.free) = params;

paramstrct = fixedparamstrct;

paramstrct.LL = -1*fval;
paramstrct.params = paramsfinal;

    function e = fit_dsprt_min(params)
        
        params2 = params20;
        params2(fixedparamstrct.free) = params;

        LL = getfit_fixedvar_gauss_bayesultimate_arr(xarr,choicearr,musgnarr,sigmavc,mu,params2,Hspc);
        e = -1*LL;
        
    end

end