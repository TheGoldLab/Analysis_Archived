clear
clc
close all

% THIS ROUTINE GENERATES DATA AND SIMULATIONS NOT IN THE MS IN CURRENT FORM 
% (AS OF 12/16), BUT WAS COMMUNICATED TO JOSH VIA EMAIL AND IS INCLUDED
% IN CASE A FUTURE VERSION OR RESEARCH WANTS TO USE THIS. SEE
% analyze_everyone_2AFC_regression_crossval.m FOR COMMENTS THAT CLARIFY
% WHAT'S HAPPENING BELOW.

% THE ROUTINE ASKS WHETHER ADAPTABILITY DEPENDS ON SESSION PARAMETERS,
% SPECIFICALLY SIGNALED VS. UNSIGNALED AND MINIMUM GENERATE VARIANCE CONDITION.
% THE ROUTINE CAN ESTIMATE ADAPTABILITY SEPARATELY FOR EACH SUBJECT BY
% CONDITION FOR REAL DATA, FIT SAMPLING MODEL PARAMETERS AND FIT IDEAL
% OBSERVER. CHRIS GLAZE FOUND THAT BOTH FIT PARAMETER SETS COULD REPRODUCE
% DEPENDENCE ON SIGNALED VS. UNSIGNALED (SESSION PARAMETER THAT DOMINATES
% LEARNING THE MOST), BUT THE SAMPLING MODEL FITS COULD DO IT MORE EASILY
% AND ROBUSTLY.


dirnm = pwd;
data_dir = [dirnm(1:strfind(dirnm,'Dropbox')+19) filesep 'data' filesep];
fit_data_dir = [dirnm(1:strfind(dirnm,'Dropbox')+19) filesep 'fit_data' filesep];

load([data_dir 'subj_info.mat'])
load([data_dir 'task_params.mat'])

all_k = 0;

subjN = 56;

genmod = 'data';
genmod = 'partfilt';
% genmod = 'bayesultimate';

whoops = .0025;
mudist = 30;
sigmamn = zeros(subjN,1);

[substrct,sessionNvc,sigmamedvc,Hlowvc,Hhighvc,signaledmedvc,signaledfirstvc] = get_subjstats;
indsval = (sessionNvc>1 & signaledfirstvc>0);

load([fit_data_dir 'allparam_fixedvarfits_regression.mat'],'paramstrct_regress','fixedparamstrct');
load([fit_data_dir 'random_session_partitions.mat']);

if strcmp(genmod,'partfilt')
    load([fit_data_dir 'allparam_fixedvarfits_partfilt_sigma_partMfixed20.mat'],'paramstrct_sigma_partMfixed');
end

if strcmp(genmod,'bayesultimate')
    load([fit_data_dir 'allparam_fixedvarfits_bayesultimate_sigmafixed.mat'],'paramstrct_bayesultimate');
end

mudist = 30;
sigmaminvc = zeros(subjN,1);

paramat = zeros(subjN,4,2);
sigmaminvc = zeros(subjN,1);

for subjind = 1:length(subjids)
    
    subjind
    
    if indsval(subjind)
        subjid = subjids{subjind};
        flnm = [data_dir 'data_' subjid '.mat'];
        
        choicearr = cell(10,1);
        xarr = choicearr;
        musgnarr = choicearr;
        logJarr = choicearr;
        sigmavc = zeros(10,1);
        tpvc = sigmavc;
        
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
                        case 'bayesultimate'
                            r0 = exp(paramstrct_bayesultimate(subjind).params(5));
                            H0 = paramstrct_bayesultimate(subjind).params(6);
                            K = exp(paramstrct_bayesultimate(subjind).params(7));
                            sigma0 = data.sigma;
       
                            xhat = x + paramstrct_bayesultimate(subjind).params(1)*randn(2000,1);
                            
                            l1 = normpdf(xhat,15,sigma0);
                            l2 = normpdf(xhat,-15,sigma0);
                            
                            Hspc = linspace(.01,.99,100);
                            Hspc = Hspc(:);
                            pH0 = betapdf(Hspc,r0*H0,r0*(1-H0));
                            pH0 = pH0 / sum(pH0);
                            
                            q1exp = ...
                                HHMM_mixed_c(l1,l2,Hspc,pH0,K,musgn);
                            
                            choice = double(q1exp>0);
                            choice = choice(:);
                            
                            whoopstmp = rand(2000,1)< paramstrct_bayesultimate(subjind).params(2);
                            choice(whoopstmp) = 1 - choice(whoopstmp);
                            
                        case 'partfilt'
                            r0 = exp(paramstrct_sigma_partMfixed(subjind).params(6));
                            H0 = paramstrct_sigma_partMfixed(subjind).params(7);
                            K = exp(paramstrct_sigma_partMfixed(subjind).params(8));
                            sigma0 = data.sigma;
                            
                            eta = paramstrct_sigma_partMfixed(subjind).params(2)*randn(2000,1);
                            F = sigma0^2/30;
                            H0samps = betarnd(H0*r0,(1-H0)*r0,20000,1);
                            %   [Hsamp,rsamp,q1samp,q2samp] = Bayesian_mixed_learnsigma_mat(K,r0,H0,numel(x),x,musgn,.0001,sigma0,N0,30);
                            Lexp = particle_filter_learnH6(K,H0samps,F,x,musgn,20,1,eta);
                            choice = double(Lexp>0);
                            
                            whoopstmp = rand(2000,1)<paramstrct_sigma_partMfixed(subjind).params(3);
                            choice(whoopstmp) = 1 - choice(whoopstmp);
                            
                    end
                    
                    session_k = session_k + 1;
                    
                    choicearr{condind} = choice;
                    xarr{condind} = x;
                    
                    musgnarr{condind} = musgn;
                    logJarr{condind} = log(data.H(:,1)./(1-data.H(:,1)));
                    sigmavc(condind) = data.sigma;
                    tpvc(condind) = mean(data.signaled);
                    
                end
                
            end
            
        end
        
        choicearr = choicearr(1:session_k);
        musgnarr = musgnarr(1:session_k);
        xarr = xarr(1:session_k);
        logJarr = logJarr(1:session_k);
        sigmavc = sigmavc(1:session_k);
        tpvc = tpvc(1:session_k);
        
    end
    
    if indsval(subjind)
        
        subjind
        
        if ~strcmp(genmod,'data')
            paramstrct_regress(subjind) =  fit_Hregress_someknown_multsession(xarr,choicearr,musgnarr,logJarr,sigmavc,mudist,fixedparamstrct);
        end
        
        sigmaminvc(subjind) = min(sigmavc);
        
        sessionarr{1} = find(tpvc==0);
        sessionarr{2} = find(tpvc>0);
        
        l1arr = cell(session_k,1);
        l2arr = l1arr;
        
        for arrind = 1:session_k
            l1arr{arrind} = normpdf(xarr{arrind},mudist/2,sigmavc(arrind));
            l2arr{arrind} = normpdf(xarr{arrind},-1*mudist/2,sigmavc(arrind));
        end
        
        fixedparamstrct_subj = fixedparamstrct;
        fixedparamstrct_subj.params0 = paramstrct_regress(subjind).params;
        
        
        for blockind = 1:numel(sessionarr)
            
            traininds = sessionarr{blockind};
            
            free0 = [0 0 0 0];
            
            for paramind = 1:4
                free = free0;
                free(paramind) = 1;
                fixedparamstrct_subj.free = logical(free);
                paramstrct = ...
                    fit_Hregress_someknown_multsession(xarr(traininds),choicearr(traininds),musgnarr(traininds),logJarr(traininds),sigmavc,mudist,fixedparamstrct_subj);
                
                paramat(subjind,paramind,blockind) = paramstrct.params(paramstrct.free);
            end
            
        end
        
    end
    
end

H0mat = squeeze(paramat(:,1,:));
H0mat = 1./(1+exp(-H0mat));
alphamat = squeeze(paramat(:,2,:));
psimat = squeeze(paramat(:,3,:));
whoopsmat = squeeze(paramat(:,4,:));

H0mat = H0mat(indsval,:);
alphamat = alphamat(indsval,:);
psimat = psimat(indsval,:);
whoopsmat = whoopsmat(indsval,:);

