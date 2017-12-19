clear
clc
close all

dirnm = pwd;
data_dir = [dirnm(1:strfind(dirnm,'Dropbox')+19) filesep 'data' filesep];
load([data_dir 'subj_info.mat'])
load([data_dir 'task_params.mat'])

all_k = 0;

subjN = 56;

genmod = 'data';
genmod = 'partfilt';
genmod = 'ideal';
genmod = 'delta_rule_leakyint';

if ~strcmp(genmod,'data')
%     % for Fig 5 main, this set of parameters for simulations gives us
%     % changes in adaptability by prior precision
%     logKvc = log(.005)*ones(subjN,1);
%     logr0vc = linspace(0,10,subjN);
%     logr0vc = logr0vc(randperm(subjN));
    
%     for Fig 5 supplemental, this set of parameters for simulations gives us
%     changes in adaptability by K
    logKvc = linspace(-10,log(.05),subjN);
    logKvc = logKvc(randperm(subjN));
    logr0vc = 2*ones(subjN,1);
    
    logalphavc = linspace(-10,log(.1),subjN);
    logalphavc = logalphavc(randperm(subjN));
    
end

mudist = 30;
sigmamn = zeros(subjN,1);

[substrct,sessionNvc,sigmamedvc,Hlowvc,Hhighvc,signaledmedvc,signaledfirstvc] = get_subjstats;
indsval = (sessionNvc>1 & signaledfirstvc>0);

fixedparamstrct.free = logical([1 1 1 1]);
fixedparamstrct.ub = [10 2 10 .02];
fixedparamstrct.lb = [-10 -1 0 .0001];
fixedparamstrct.params0 = [0 .5 1 .0025];

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
                    case 'data'
                        choice = double(data.pred==2);
                    case 'ideal'

                        r0 = exp(logr0vc(subjind));
                        H0 = .4;
                        K = exp(logKvc(subjind));
                        sigma0 = data.sigma;
                        l1 = normpdf(x,15,sigma0);
                        l2 = normpdf(x,-15,sigma0);
      
                        Hspc = linspace(.01,.99,100);
                        Hspc = Hspc(:);
                        pH0 = betapdf(Hspc,r0*H0,r0*(1-H0));
                        pH0 = pH0 / sum(pH0);
          
                        q1exp = ...
                            HHMM_mixed_c(l1,l2,Hspc,pH0,K,musgn);
                        
                        choice = double(q1exp>.5);
                        choice = choice(:);
                        
                    case 'partfilt'
                        partM = 20;
                        r0 = exp(logr0vc(subjind));
                        H0 = .4;
                        K = exp(logKvc(subjind));
                        eta = 0*randn(20000,1);
                        sigmasubj = data.sigma;
                        F = sigmasubj^2/30;
                        H0samps = betarnd(r0*H0,r0*(1-H0),20000,1);
                        Lexp_final = particle_filter_learnH6(K,H0samps,F,x,musgn,partM,1,eta);
                        
                        choice = double(Lexp_final>0);
                        choice = choice(:);
                        
                    case 'delta_rule_normative'
                        H0 = 0.4;
                        alpha = exp(logalphavc(subjind));
                        eta = 0*randn(20000,1);
                        sigmasubj = data.sigma;
                        F = sigmasubj^2/30;
                        
                        Lexp_final = fixedLR_learnH6(alpha,H0,F,x,musgn,1,zeros(length(x),1));
                        
                        choice = double(Lexp_final>0);
                        choice = choice(:);  
                        
                     case 'delta_rule_leakyint'
                        H0 = 0.4;
                        alpha = exp(logalphavc(subjind));
                        eta = 0*randn(20000,1);
                        sigmasubj = data.sigma;
                        F = sigmasubj^2/30;
                        
                        Lexp_final = fixedLR_learnH6_approx(alpha,5,H0,F,x,musgn,1,zeros(length(x),1));
                        
                        choice = double(Lexp_final>0);
                        choice = choice(:);   

%                         mistake = rand(2000,1)<.0024;
%                         choice(mistake) = 1 - choice(mistake);
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


