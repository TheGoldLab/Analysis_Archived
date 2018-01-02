% get file information
[~, analysis_data_dir, raw_data_dir] = getDataInfo;

[subjstrct,sessionNvc,sigmamedvc,Hlowvc,Hhighvc,signaledmedvc,signaledfirstvc] = get_subjstats;
indsval = sessionNvc>1 & signaledfirstvc>0;

load(fullfile(analysis_data_dir, 'subj_info.mat'));
load(fullfile(analysis_data_dir, 'task_params.mat'));

all_k = 0;

subjN = 56;

genmod = 'data';

mudist = 30;

fixedparamstrct.varnms = {'partM','psi','whoops','sigma0','gamma','logr0','H0','logK'};
fixedparamstrct.params0 = [20 1.5 .0025 10 1 2 .4 -5];
fixedparamstrct.ub = [50 5 .02 20 2 log(2000) .999 log(.1)];
fixedparamstrct.lb = [2 .0001 .0001 2 0 log(.1) .001 log(.0001)];
fixedparamstrct.free = logical(zeros(1,numel(fixedparamstrct.varnms)));

load(fullfile(analysis_data_dir, 'random_session_partitions.mat'), 'sessionindsarr');
load(fullfile(analysis_data_dir, 'allparam_fixedvarfits_partfilt_sigma_partMfixed20.mat'),'paramstrct_sigma_partMfixed');

for subjind = 1:length(subjids)
    
    subjind
    subjid = subjids{subjind};
    flnm = fullfile(raw_data_dir, ['data_' subjid '.mat']);
    
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
        
        fixedparamstrct_subj = fixedparamstrct;
        
        paramvc = [6 8];
        
        for sessionind = 1:2
            sessioninds = sessionindsarr{subjind,sessionind};
            
            for paramind = 1:numel(paramvc)
                fixedparamstrct_subj.params0 = paramstrct_sigma_partMfixed(subjind).params;
                fixedparamstrct_subj.free(paramvc) = 0;
                fixedparamstrct_subj.free(paramvc(paramind)) = 1;
                fixedparamstrct_subj.params0(paramvc(paramind)) = fixedparamstrct.params0(paramvc(paramind));
                paramstrct_sigma_partMfixed_crossval(subjind,paramind,sessionind) = ...
                    fit_partfilt_someknown_multsession(xarr(sessioninds),choicearr(sessioninds),musgnarr(sessioninds),sigmavc(sessioninds),mudist,fixedparamstrct_subj);
                
            end
        end

    end
    
end


