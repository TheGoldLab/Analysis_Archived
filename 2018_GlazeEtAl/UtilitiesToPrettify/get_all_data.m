function [logJarr,xarr,l1arr,l2arr,musgnarr,sigmavc,subjvc,correctarr] = get_all_data

dirnm = pwd;
data_dir = [dirnm(1:strfind(dirnm,'Dropbox')+19) filesep 'data' filesep];

load([data_dir 'subj_info.mat'])
load([data_dir 'task_params.mat'])

subjN = 56;

logJarr = cell(220,1);
l1arr = logJarr;
l2arr = logJarr;
xarr = logJarr;
musgnarr = logJarr;
sigmavc = zeros(220,1);
subjvc = sigmavc;
correctarr = xarr;

[substrct,sessionNvc,sigmamedvc,Hlowvc,Hhighvc,signaledmedvc,signaledfirstvc] = get_subjstats;
indsval = (sessionNvc>1 & signaledfirstvc>0);

blockinds = sort(repmat([1:10]',200,1));

session_k = 0;

for subjind = 1:length(subjids)
    
    indsval(subjind)
    subjid = subjids{subjind};
    flnm = [data_dir 'data_' subjid '.mat'];
    
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
                
                session_k = session_k + 1;
                
                logJarr{session_k} = grpstats(log(data.H(:,1)./(1-data.H(:,1))),blockinds);
                l1arr{session_k} = normpdf(x,15,data.sigma);
                l2arr{session_k} = normpdf(x,-15,data.sigma);
                musgnarr{session_k} = musgn;
                sigmavc(session_k) = data.sigma;
                subjvc(session_k) = subjind;
                xarr{session_k} = x;
                correctarr{session_k} = double(data.muinds==2);
                
            end
            
        end
        
    end
    
end


logJarr = logJarr(1:session_k);
xarr = xarr(1:session_k);
l1arr = l1arr(1:session_k);
l2arr = l2arr(1:session_k);
musgnarr = musgnarr(1:session_k);
correctarr = correctarr(1:session_k);
sigmavc = sigmavc(1:session_k);


