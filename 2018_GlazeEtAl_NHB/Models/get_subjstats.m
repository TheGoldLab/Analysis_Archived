function [subjstrct,sessionNvc,sigmamedvc,Hlowvc,Hhighvc,signaledmedvc,signaledfirstvc] = get_subjstats

% dirnm = pwd;
% data_dir = [dirnm(1:strfind(dirnm,'Dropbox')+19) filesep 'data' filesep];

% Where to find the data
[~,analysis_data_dir, raw_data_dir] = getDataInfo;

load([analysis_data_dir '/subj_info.mat'])
load([analysis_data_dir '/task_params.mat'])

subjN = length(subjids);

sessionNvc = zeros(subjN,1);
sigmamedvc = sessionNvc;
Hlowvc = sessionNvc;
Hhighvc = sessionNvc;
signaledmedvc = sessionNvc;
signaledfirstvc = sessionNvc;

for subjind = 1:subjN
    
    subjid = subjids{subjind};
    flnm = [raw_data_dir '/data_' subjid '.mat'];
    
    sessionNtmp = 0;
    sigmamedtmp = 0;
    Hlowtmp = zeros(10,1)+nan;
    Hhightmp = zeros(10,1)+nan;
    signaledmedtmp = 0;
    signaledfirst = 0;
    
    if exist(flnm)
        load(flnm)
        
        if session_ind>1
            for condind = 1:session_ind-1
                
                if session_params(condind).signaled<10 & size(session_params(condind).Hset,2)<100
                    
                    eval(['data=data' num2str(condind) ';'])
                    
                    sessionNtmp = sessionNtmp + 1;
                    sigmamedtmp = sigmamedtmp + data.sigma;
                    Hlowtmp(condind) = min(data.Hset(1,:));
                    Hhightmp(condind) = max(data.Hset(1,:));
                    signaledmedtmp = signaledmedtmp + mean(data.signaled);
                    
                    if condind==1
                        signaledfirst = mean(data.signaled);
                    end
                    
                end
                
                
            end
            
            
            sessionNvc(subjind) = sessionNtmp;
            sigmamedvc(subjind) = sigmamedtmp / sessionNtmp;
            Hlowvc(subjind) = min(Hlowtmp(1:condind));
            Hhighvc(subjind) = max(Hhightmp(1:condind));
            signaledmedvc(subjind) = signaledmedtmp / sessionNtmp;
            signaledfirstvc(subjind) = signaledfirst;
            
            subjstrct(subjind).sessionN = sessionNtmp;
            subjstrct(subjind).sigmamed = sigmamedtmp / sessionNtmp;
            subjstrct(subjind).Hlow = min(Hlowtmp(1:condind));
            subjstrct(subjind).Hhigh = max(Hhightmp(1:condind));
            subjstrct(subjind).signaledmed = signaledmedtmp / sessionNtmp;
            subjstrct(subjind).signaledfirst = signaledfirst;
            
        end
    end
    
end
