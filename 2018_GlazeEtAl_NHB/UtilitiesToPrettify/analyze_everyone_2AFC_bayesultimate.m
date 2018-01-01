clear
clc
close all

dirnm = pwd;
data_dir = [dirnm(1:strfind(dirnm,'Dropbox')+19) filesep 'data' filesep];

[subjstrct,sessionNvc,sigmamedvc,Hlowvc,Hhighvc,signaledmedvc,signaledfirstvc] = get_subjstats;
indsval = sessionNvc>1 & signaledfirstvc>0;

load([data_dir 'subj_info.mat'])
load([data_dir 'task_params.mat'])

all_k = 0;

subjN = 56;

genmod = 'data';

mudist = 30;

fixedparamstrct.varnms = {'psi','whoops','sigma0','gamma','logr0','H0','logK'};
fixedparamstrct.params0 = [1 .0025 10 1 log(20) .4 log(.01)];
fixedparamstrct.ub = [5 .02 20 2 log(20000) .99 log(.1)];
fixedparamstrct.lb = [.0001 .0001 2 0 log(.1) .01 log(.0001)];
fixedparamstrct.free = logical([1 1 0 0 1 1 1]);

Hspc = linspace(.01,.99,100);
Hspc = Hspc(:);

for subjind = 1:length(subjids)
    
    subjind
    subjid = subjids{subjind};
    flnm = [data_dir 'data_' subjid '.mat'];
    
    choicearr = cell(10,1);
    xarr = choicearr;
    musgnarr = choicearr;
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
                
                choice = double(data.pred==2);
                
                session_k = session_k + 1;
                
                choicearr{condind} = choice;
                xarr{condind} = x;
                musgnarr{condind} = musgn;
                sigmavc(condind) = data.sigma;
                
            end
            
        end
        
    end
    
    choicearr = choicearr(1:session_k);
    musgnarr = musgnarr(1:session_k);
    xarr = xarr(1:session_k);
    sigmavc = sigmavc(1:session_k);

    if session_k>=1 & indsval(subjind)
        
        subjind
        paramstrct_bayesultimate(subjind) = fit_bayesultimate_someknown_multsession(xarr,choicearr,musgnarr,sigmavc,mudist,fixedparamstrct,Hspc);

    end
    
end



