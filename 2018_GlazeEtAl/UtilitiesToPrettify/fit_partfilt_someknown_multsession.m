function [paramstrct] = fit_partfilt_someknown_multsession(xarr,choicearr,musgnarr,sigmavc,mu,fixedparamstrct)

eta0 = randn(20000,1);

ub = fixedparamstrct.ub(fixedparamstrct.free);
lb = fixedparamstrct.lb(fixedparamstrct.free);
params0 = fixedparamstrct.params0(fixedparamstrct.free);

params20 = zeros(numel(fixedparamstrct.free),1);
params20(~fixedparamstrct.free) = fixedparamstrct.params0(~fixedparamstrct.free);

model = @fit_dsprt_min;

% options = psoptimset('maxiter',5000,'Display','iter','CompletePoll','on','TolMesh',1e-4,'CompleteSearch','on','UseParallel',true,'ScaleMesh','on','InitialMeshSize',1,'SearchMethod',@searchlhs);
% % 
options = psoptimset('MaxIter',5000,'Display','iter','CompletePoll','on','TolMesh',1e-6,'UseParallel',true,'ScaleMesh','on','InitialMeshSize',1,'TimeLimit',900,'CompleteSearch','off');
[params,fval] = patternsearch(model,params0,[],[],[],[],lb,ub,[],options);

% options = saoptimset('maxiter',1000,'Display','iter','InitialTemperature',100);
% [params,fval] = simulannealbnd(model,params0,lb,ub,options);

% opts = optimoptions('ga','MaxTime',3600,'Display','iter','UseParallel',true,'MaxGenerations',20,'PopulationSize',50);

% partM0 = [5*ones(12,1);10*ones(12,1);20*ones(13,1);40*ones(12,1)];
% psi0 = linspace(1,3,49);
% psi0 = psi0(randperm(49));
% whoops0 = linspace(.001,.015,49);
% whoops0 = whoops0(randperm(49));
% sigma0 = linspace(7,11,49);
% sigma0 = sigma0(randperm(49));
% gamma = .7*ones(49,1);
% logr0 = [linspace(.5,6,12)';linspace(1,7,12)';linspace(1,10,13)';linspace(1,10,12)'];
% H0 = linspace(.1,.7,49);
% H0 = H0(randperm(49));
% logK = linspace(log(.0001),log(.1),49);
% logK = logK(randperm(49));

% params0 = [params0;[partM0 psi0(:) whoops0(:) sigma0(:) gamma(:) logr0(:) H0(:) logK(:)]];
% 
% if nargin>=7
%     params0 = paramat0;
% else
%     params0 = fixedparamstrct.params0;
% end
% 
% opts.InitialPopulationMatrix = params0;
% 
% [params,fval] = ga(model,size(params0,2),[],[],[],[],lb,ub,[],1,opts);

paramsfinal = params20;
paramsfinal(fixedparamstrct.free) = params;

paramstrct = fixedparamstrct;

paramstrct.LL = -1*fval;
paramstrct.params = paramsfinal;

paramstrct.LLrl = getfit_fixedvar_gauss_partfilt_arr(xarr,choicearr,musgnarr,sigmavc,mu,10000,paramstrct.params,eta0);

    function e = fit_dsprt_min(params,simN)
        
        if nargin==1
            simN = 10000;
        end
        
        params2 = params20;
        params2(fixedparamstrct.free) = params;
        
%         simN = ceil(params2(2)*100 + 500/params2(6));
%         simN = max(min(simN,1500),10);
 
        LL = getfit_fixedvar_gauss_partfilt_arr(xarr,choicearr,musgnarr,sigmavc,mu,simN,params2,eta0);
        e = -1*LL;
        
    end

end