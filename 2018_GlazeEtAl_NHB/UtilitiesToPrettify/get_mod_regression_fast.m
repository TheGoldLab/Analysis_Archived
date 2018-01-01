function alphavc = get_mod_regression_fast(params,genmod,logJarr,xarr,l1arr,l2arr,musgnarr,sigmavc,subjvc)
subjU = unique(subjvc);
subjN = numel(subjU);

subjN = 1;

alphavc = zeros(subjN,1);

blockinds = sort(repmat([1:10]',200,1));

switch genmod
    case 'partfilt'
        H0samps = betarnd(params(2)*exp(params(1)),(1-params(2))*exp(params(1)),20000,1);
    case 'bayesultimate'
        Hspc = linspace(.01,.99,100);
        Hspc = Hspc(:);
        pH0 = betapdf(Hspc,params(2)*exp(params(1)),(1-params(2))*exp(params(1)));
        pH0 = pH0 / sum(pH0);
    case 'delta_rule_norm'
        alpha = exp(params(1));
        H0 = params(2);   
    case 'delta_rule_leakint'
        alpha = exp(params(1));
        H0 = params(2);
end


for subjind = 1:subjN
    sessioninds = find(subjvc==subjU(subjind));
    sessionN = numel(sessioninds);
    logJall = zeros(sessionN*10,1);
    logJexpall = logJall;
    
    k = 1;
    kindsarr = cell(sessionN,1);
    
    for sessionind = 1:sessionN
        kindsarr{sessionind} = [k:k+9];
        k = k + 10;
    end
    
    for sessionind = 1:sessionN
        logJ = logJarr{sessioninds(sessionind)};
        musgn = musgnarr{sessioninds(sessionind)};
        l1 = l1arr{sessioninds(sessionind)};
        l2 = l2arr{sessioninds(sessionind)};
        x = xarr{sessioninds(sessionind)};
        
        switch genmod
            case 'bayesultimate'
                [q1samp,Hexp] = HHMM_mixed_c(l1,l2,Hspc,pH0,exp(params(3)),musgn);

            case 'partfilt'
                F = sigmavc(sessioninds(sessionind))^2/30;
                [Lexp,Lvar,Hexp] = particle_filter_learnH7(exp(params(3)),H0samps,F,x,musgn,params(4),1,zeros(10,1));
                                   
            case 'delta_rule_norm'
                F = sigmavc(sessioninds(sessionind))^2/30;
                [Lexp,Lvar,Hexp] = fixedLR_learnH(alpha,H0,F,x,musgn,1,zeros(length(x),1));

            case 'delta_rule_leakint'
                F = sigmavc(sessioninds(sessionind))^2/30;
                [Lexp,Lvar,Hexp] = fixedLR_learnH_approx(alpha,1+30/sigmavc(sessioninds(sessionind)),H0,F,x,musgn,1,zeros(length(x),1));

        end

        logJexp = log(Hexp) - log(1-Hexp);
        kindstmp = kindsarr{sessionind};
        logJall(kindstmp) = logJ;
        logJexpall(kindstmp) = grpstats(logJexp,blockinds);
        
    end
    
    beta = regress(logJexpall,[ones(k-1,1) logJall]);
    alphavc(subjind) = beta(2);
    
end
