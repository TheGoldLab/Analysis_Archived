% get file information
[~, analysis_data_dir, raw_data_dir] = getDataInfo;

[subjstrct,sessionNvc,sigmamedvc,Hlowvc,Hhighvc,signaledmedvc,signaledfirstvc] = get_subjstats;
indsval = sessionNvc>1 & signaledfirstvc>0;

load(fullfile(analysis_data_dir, 'subj_info.mat'));
load(fullfile(analysis_data_dir, 'task_params.mat'));

all_k = 0;

subjN = 56;

genmod = 'data';
% genmod = 'ideal';

whoops = .0025;
mudist = 30;
sigmamn = zeros(subjN,1);

[substrct,sessionNvc,sigmamedvc,Hlowvc,Hhighvc,signaledmedvc,signaledfirstvc] = get_subjstats;
indsval = (sessionNvc>1 & signaledfirstvc>0);

load(fullfile(analysis_data_dir, 'random_session_partitions.mat'));
load(fullfile(analysis_data_dir, 'allparam_fixedvarfits_regression.mat'), ...
   'paramstrct_regress','fixedparamstrct');

mudist = 30;
sigmamn = zeros(subjN,1);

paramat = zeros(subjN,4,2);


for subjind = 1:length(subjids)
    
    subjind
    subjid = subjids{subjind};
    flnm = fullfile(raw_data_dir, ['data_' subjid '.mat']);
    
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
                
                choice = double(data.pred==2);
                
                
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
        
        sigmamn(subjind) = mean(sigmavc);
        
        sessionarr{1} = sessionindsarr{subjind,1};
        sessionarr{2} = sessionindsarr{subjind,2};
        
        l1arr = cell(session_k,1);
        l2arr = l1arr;
        
        for arrind = 1:session_k
            l1arr{arrind} = normpdf(xarr{arrind},mudist/2,sigmavc(arrind));
            l2arr{arrind} = normpdf(xarr{arrind},-1*mudist/2,sigmavc(arrind));
        end
        
        
        fixedparamstrct.params0 = paramstrct_regress(subjind).params;
        
        
        for blockind = 1:numel(sessionarr)
            
            traininds = sessionarr{blockind};
            
            free0 = [0 0 0 0];
            
            for paramind = 1:4
                free = free0;
                free(paramind) = 1;
                fixedparamstrct.free = logical(free);
                paramstrct = ...
                    fit_Hregress_someknown_multsession(xarr(traininds),choicearr(traininds),musgnarr(traininds),logJarr(traininds),sigmavc,mudist,fixedparamstrct);
                
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
