% get file information
[~, analysis_data_dir, raw_data_dir] = getDataInfo;

load(fullfile(analysis_data_dir, 'subj_info.mat'));
load(fullfile(analysis_data_dir, 'task_params.mat'));


blockinds = sort(repmat([1:10],1,200));
blockinds = blockinds(:);

blockinds_sigma = sort(repmat([1:4],1,500));
blockinds_sigma = blockinds_sigma(:);

all_k = 0;

Hmat = zeros(220,10);
psivc = zeros(220,1);
subjvc = psivc;
sigmamat = zeros(220,numel(unique(blockinds_sigma)));
LLvc = psivc;

Hobjmat = zeros(220,10);
sigmaobjvc = zeros(220,1);
signaledmat = Hobjmat;

sessionindvc = zeros(220,1);

Hnumvc = psivc;
tpvc = psivc;

pcorrectvc = zeros(220,1);

for subjind = 1:length(subjids)
    
    subjind
    subjid = subjids{subjind};
    flnm = fullfile(raw_data_dir, ['data_' subjid '.mat']);
    
    if exist(flnm)
        load(flnm)
        
        for condind = 1:session_ind-1
            
            if session_params(condind).signaled<10 & size(session_params(condind).Hset,2)<100
                
                eval(['data=data' num2str(condind) ';'])
                
                midpt = mean(data.muall(:,1));
                x = data.X(:,1)-midpt;
                F = data.sigma^2 / 15;
                
                LLR = 2 * x / F;
                
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
                
                all_k = all_k + 1
                
                choice = double(data.pred==2);
                
                blockN = size(data.Hset,2);
                
                params = fit_dsprt_blocked(x,choice,musgn,blockinds,blockinds_sigma);
                
                
                psivc(all_k) = params(end);
                Hmat(all_k,:) = params(1:10);
                sigmamat(all_k,:) = params(11:11+numel(unique(blockinds_sigma))-1);
                
                Hobjmat(all_k,:) = grpstats(cpvc,blockinds);
                tpvc(all_k) = session_params(condind).signaled;
                sigmaobjvc(all_k) = data.sigma;
                
                correctchoice = double(data.muinds==2);
                
                pcorrectvc(all_k) = mean(choice==correctchoice);
                
                subjvc(all_k) = subjind;
                
                Hnumvc(all_k) = size(data.Hset,2);
                
            end
            
            
        end
        
    end
    
end
