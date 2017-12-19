function paramstrct_regress = get_mod_regression(params,genmod)

dirnm = pwd;
data_dir = [dirnm(1:strfind(dirnm,'Dropbox')+19) filesep 'data' filesep];

load([data_dir 'subj_info.mat'])
load([data_dir 'task_params.mat'])

subjN = 56;


[substrct,sessionNvc,sigmamedvc,Hlowvc,Hhighvc,signaledmedvc,signaledfirstvc] = get_subjstats;
indsval = (sessionNvc>1 & signaledfirstvc>0);



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


fixedparamstrct.free = logical([1 1 1 0]);
fixedparamstrct.ub = [10 2 10 .02];
fixedparamstrct.lb = [-10 -1 0 .0001];
fixedparamstrct.params0 = [0 .5 1 0];

mudist = 30;

for subjind = 1:length(subjids)
    
    subjind
    subjid = subjids{subjind};
    flnm = [data_dir 'data_' subjid '.mat'];
    
    choicearr = cell(10,1);
    xarr = choicearr;
    musgnarr = choicearr;
    logJarr = choicearr;
    sigmavc = zeros(10,1);
    
    session_k = 0;
    
    if exist(flnm)
        load(flnm)
        
        if session_ind>1
            
            for condind = 1:session_ind-1
                
                eval(['data=data' num2str(condind) ';'])
                
                midpt = mean(data.muall(:,1));
                x = data.X(:,1)-midpt;
                F = data.sigma^2 / 30;
                
                LLR = x / F;
                
                if data.H(1)>1
                    data.H = 1./data.H;
                end
                
                musgn = sign(data.muall(data.muinds,1)-midpt);
                musgn(data.signaled==0) = 0;
                
                musgn = [0;musgn(1:end-1)];
                
                cpvc = double(data.r==1);
                
                if data.Hset(1)>1
                    data.Hset = 1./data.Hset;
                end
                
                
                switch genmod
                    case 'bayesultimate'
                        sigma0 = data.sigma;
                        l1 = normpdf(x,15,sigma0);
                        l2 = normpdf(x,-15,sigma0);
                        [q1samp,Hexp] = HHMM_mixed_c(l1,l2,Hspc,pH0,exp(params(3)),musgn);
                        choice = double(q1samp>.5);
                        
                    case 'partfilt'
                        sigma0 = data.sigma;
                        F = sigma0^2/30;
                        Lexp = particle_filter_learnH6(exp(params(3)),H0samps,F,x,musgn,params(4),1,zeros(10,1));
                        choice = double(Lexp>0);
                        
                    case 'delta_rule_norm'
                        sigma0 = data.sigma;
                        F = sigma0^2/30;
                        Lexp = fixedLR_learnH6(alpha,H0,F,x,musgn,1,zeros(length(x),1));
                        
                        choice = double(Lexp>0);
                        
                    case 'delta_rule_leakint'
                        sigma0 = data.sigma;
                        F = sigma0^2/30;
                        Lexp = fixedLR_learnH6_approx(alpha,H0,F,x,musgn,1,zeros(length(x),1));
                        
                        choice = double(Lexp>0);

                end
                
                session_k = session_k + 1;
                
                choicearr{condind} = choice;
                xarr{condind} = x;
                
                musgnarr{condind} = musgn;
                logJarr{condind} = log(data.H(:,1)./(1-data.H(:,1)));
                sigmavc(condind) = data.sigma;
                
            end
            
        end
        
    end
    
    choicearr = choicearr(1:session_k);
    musgnarr = musgnarr(1:session_k);
    xarr = xarr(1:session_k);
    logJarr = logJarr(1:session_k);
    sigmavc = sigmavc(1:session_k);
    
    if indsval(subjind)
        
        subjind
        
        paramstrct_regress(subjind) = fit_Hregress_someknown_multsession(xarr,choicearr,musgnarr,logJarr,sigmavc,mudist,fixedparamstrct);
        
    end
    
end


