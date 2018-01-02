function [correctvc,RTvc,subjvc,Hvc,rvc,musgnvc,evc,initialevc,blockindvc,isprevc,sessionindsvc,paramat,alphavc,datacorrectvc,xvc,choicevc,xprevc,rewardvc] = get_modelfree_data(genmod)

% dirnm = pwd;
% data_dir = [dirnm(1:strfind(dirnm,'Dropbox')+19) filesep 'data' filesep];
% 

% Where to find the data
[~,analysis_data_dir, raw_data_dir] = getDataInfo;

if nargin==0
    genmod = 'data';
end

[substrct,sessionNvc,sigmamedvc,Hlowvc,Hhighvc,signaledmedvc,signaledfirstvc] = get_subjstats;
subjindsval = (sessionNvc>1 & signaledfirstvc>0);

Ntot = sum(sessionNvc(subjindsval))*2000;

load([analysis_data_dir '/subj_info.mat'])
load([analysis_data_dir '/task_params.mat'])

k = 1;
correctvc = zeros(Ntot,1);
RTvc = correctvc;
subjvc = correctvc;
Hvc = correctvc;
rvc = Hvc;
evc = Hvc;
musgnvc = Hvc;
alphavc = Hvc;
blockindvc = Hvc;
initialevc = Hvc;
isprevc = Hvc;
datacorrectvc = Hvc;
paramat = zeros(Ntot,8);
sessionindsvc = Hvc;
xvc = Hvc;
choicevc = xvc;
xprevc = xvc;
rewardvc = xvc;

% fit_data_dir = [dirnm(1:strfind(dirnm,'Dropbox')+19) filesep 'fit_data' filesep];


load(fullfile(analysis_data_dir, 'allparam_fixedvarfits_partfilt_sigma_partMfixed20.mat'),'paramstrct_sigma_partMfixed');
load(fullfile(analysis_data_dir, 'allparam_fixedvarfits_bayesultimate_sigmafixed.mat'),'paramstrct_bayesultimate');
load(fullfile(analysis_data_dir, 'allparam_fixedvarfits_regression.mat'),'paramstrct_regress')

for subjind = 1:length(subjids)

    if subjindsval(subjind)
        subjid = subjids{subjind};
        flnm = fullfile(raw_data_dir, ['data_' subjid '.mat']);
        
        if exist(flnm)
            load(flnm)
            
            if session_ind>1
                
                for condind = 1:session_ind-1
                    
                    eval(['data=data' num2str(condind) ';'])
                    
                    midpt = mean(data.muall(:,1));
                    x = data.X(:,1)-midpt;
                    F = data.sigma^2 / 30;
                    
                    LLR = x / F;
                    absLLRprev = LLR;
                    
                    
                    if data.H(1)>1
                        data.H = 1./data.H;
                    end
                    
                    musgn = sign(data.muall(data.muinds,1)-midpt);
                    musgn(data.signaled==0) = 0;
                    musgn = [0;musgn(1:end-1)];
                    
                    
                    musgn2 = sign(data.muall(data.muinds,1)-midpt);
                    absLLRprev = absLLRprev.*musgn2;
                    
                    cpvc = double(data.r==1);
                    
                    if data.Hset(1)>1
                        data.Hset = 1./data.Hset;
                    end
                    
                    
                    switch genmod
                        case 'data'
                            choice = double(data.pred==2);
                            correct = double(data.pred==data.muinds);
                        case 'partfilt'
                            params = paramstrct_sigma_partMfixed(subjind).params;
                            F = data.sigma^2/30;
                            H0samps = betarnd(exp(params(6))*params(7),exp(params(6))*(1-params(7)),20000,1);
                            eta = params(2)*randn(2000,1);
                            [Lexp,Lvar] = particle_filter_learnH6(exp(params(8)),H0samps,F,x,musgn,params(1),1,eta);
                            choice = double(Lexp>0);
                            whoops = logical(rand(2000,1)<params(3));
                            choice(whoops) = 1 - choice(whoops);
                            correctchoice = double(data.muinds==2);
                            correct = double(choice==correctchoice);
                            
                        case 'bayesultimate'
                            r0 = exp(paramstrct_bayesultimate(subjind).params(5));
                            H0 = paramstrct_bayesultimate(subjind).params(6);
                            K = exp(paramstrct_bayesultimate(subjind).params(7));
                            sigma0 = data.sigma;
       
                            xhat = x + paramstrct_bayesultimate(subjind).params(1)*randn(2000,1);
                            
                            l1 = normpdf(xhat,15,sigma0);
                            l2 = normpdf(xhat,-15,sigma0);
                       
                            q1samp = Bayesian_mixed6(K,r0,H0,numel(x),l1,l2,musgn,.0001);
                            choice = double(q1samp>.5);
                            
                            whoopstmp = rand(2000,1)< paramstrct_bayesultimate(subjind).params(2);
                            choice(whoopstmp) = 1 - choice(whoopstmp);
                            correctchoice = double(data.muinds==2);
                            correct = double(choice==correctchoice);
                    end
                    
                    Hnum = size(data.Hset,2);
                    HN = 2000/Hnum;
                    
                    blockinds = repmat([1:HN],Hnum,1);
                    
                    sessioninds = [1:2000]';
                    
                    datamuinds = double(LLR>0)+1;
                    
                    datacorrect = double(data.muinds==datamuinds);
                    
                    kinds = k:k+2000-1;
                    k = k + 2000;
                    correctvc(kinds) = correct;
                    paramat(kinds,:) = repmat(paramstrct_sigma_partMfixed(subjind).params(:)',2000,1);
                    rvc(kinds) = data.r;
                    Hvc(kinds) = data.H(:,1);
                    subjvc(kinds) = subjind;
                    RTvc(kinds) = data.rt - median(data.rt(sessioninds>200 & data.rt<3000));
                    musgnvc(kinds) = musgn;
                    alphavc(kinds) = paramstrct_regress(subjind).params(2);
                    sessionindsvc(kinds) = sessioninds;
                    datacorrectvc(kinds) = datacorrect;
                    
                    blockindvc(kinds) = blockinds;
                    
                    choicevc(kinds) = choice;
                    
                    ksub = 1;
                    cpinds = find(cpvc);
                    absLLRprevrun = zeros(2000,1);
                    
                    cpinds = [cpinds;2000];
                    
                    for cpind = 2:numel(cpinds)-1
                        runinds = cpinds(cpind):cpinds(cpind+1)-1;
                        absLLprevtmp = absLLRprev(cpinds(cpind)-1);
                        absLLRprevrun(runinds) = absLLprevtmp;
                    end
                    
                    rnext = [data.r(2:end);0];
                    
                    xprev = [0;x(1:end-1)];
                    
                    initialevc(kinds) = absLLRprevrun;
                    evc(kinds) = abs(LLR);
                    xvc(kinds) = x;
                    xprevc(kinds) = xprev;
                    isprevc(kinds) = double(rnext==1);
                    
                    if ~isfield(session_params(condind),'reward') || isempty(session_params(condind).reward)
                        reward = 2;
                        penalty = 20;
                    else
                        reward = session_params(condind).reward;
                        penalty = session_params(condind).penalty;
                    end
                    
                    rewardvc(kinds) = correct*reward - (1-correct)*penalty;
                    
                end
                
            end
            
        end
        
    end
    
end