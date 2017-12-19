function paramstrct = fit_Hregress_someknown_multsession(xarr,choicearr,musgnarr,logJarr,sigmavc,mudist,paramstrct)

arrN = numel(xarr);
l1arr = cell(arrN,1);
l2arr = l1arr;

for arrind2 = 1:arrN
    l1arr{arrind2} = normpdf(xarr{arrind2},mudist/2,sigmavc(arrind2));
    l2arr{arrind2} = normpdf(xarr{arrind2},-mudist/2,sigmavc(arrind2));
end

ub = paramstrct.ub(paramstrct.free);
params0 = paramstrct.params0(paramstrct.free);
lb = paramstrct.lb(paramstrct.free);

paramsall0 = paramstrct.params0;


model = @fit_dsprt_min;

opts = optimoptions(@fmincon,'Algorithm','interior-point');
ms = MultiStart('UseParallel','always');
problem = createOptimProblem('fmincon','objective',...
    model,'x0',params0,'lb',lb,'ub',ub,'options',opts);

[params,fval] = run(ms,problem,16);

paramstrct.params = paramsall0;
paramstrct.params(paramstrct.free) = params;
paramstrct.LL = -fval;

   function e = fit_dsprt_min(params)

        paramsall = paramsall0;
        paramsall(paramstrct.free) = params;
        paramsall(~paramstrct.free) = paramstrct.params0(~paramstrct.free);
        e = getfit_Hregress(l1arr,l2arr,choicearr,musgnarr,logJarr,paramsall);

    end

end