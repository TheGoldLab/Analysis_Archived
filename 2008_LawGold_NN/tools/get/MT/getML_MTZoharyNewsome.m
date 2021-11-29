function [MTf1, MTs1, MTf2, MTs2, Bef1, Bes1, Bef2, Bes2] = getML_MTZoharyNewsome(Monk, recompute)
%  For each MT cell, get neurometric ROC and threshold for whole session and for the
%  first and second half of the session
savepath = ['/work/mirror_jeff/code/matlab/tmp-mat/getML_MTZoharyNewsome_' Monk '.mat'];

if recompute
%% get ROC

crtf = [0 1];
nvf  = [0 1];
dirf = 0;

a      = getML_txt([Monk, 'TRain_MT.txt']);
fn     = a.data{strcmp(a.name,'dat_fn')};
ses    = a.data{strcmp(a.name,'session')};
usable = a.data{strcmp(a.name,'usable')};
train  = a.data{strcmp(a.name,'train')};
uid    = a.data{strcmp(a.name,'uid')};
d1     = a.data{strcmp(a.name,'ddir')};

% MT roc
rocs   = nans(7,length(fn),3);
rocn   = nans(7,length(fn),3);

% performance
ps     = nans(7,length(fn),3);
pn     = nans(7,length(fn),3);

if strcmp(Monk,'Cy')
    N=200;
else
    N=200;
end

global FIRA
for i=1:length(fn)
    if usable(i)==1 & train(i)==1
        % load data
        openFIRA(fn{i})
        fprintf('%s\n',fn{i})

        %%%     get selection arrays    %%%
        % get roc for each trials, select trials based on
        % correct flag (crtf) and no var flag (nvf)
        Lcrt   =  ismember(getFIRA_ecodesByName('correct', 'id'), crtf);
        trials = ~isnan(FIRA.ecodes.data(:,1)) & Lcrt;

        % task
        Ltsk   =  ismember(getFIRA_ecodesByName('task', 'id'), 6);

        % get selection arrays for dir/coherence
        % remove coh other than the standard cohs (i.e. 0 32 64 128
        % 256 512 999)  (in one pre-training file there's a 9% coh)
        [Lcoh, Ucoh] = selectFIRA_trialsByUniqueID('dot_coh');
        Lcoh = Lcoh(:, ismember(Ucoh, [0 3.2 6.4 12.8 25.6 51.2 99.9]));
        Ucoh = Ucoh(ismember(Ucoh, [0 3.2 6.4 12.8 25.6 51.2 99.9]));

        [Ldir, Udir] = selectFIRA_trialsByUniqueID('dot_dir');
        Idir         = [find(round(Udir)==d1(i)) find(round(Udir)==mod(d1(i)+180,360))];
        Ldir         = Ldir(:,Idir);

        bt           = getFIRA_ecodeTimesByName('dot_on',0);
        et           = getFIRA_ecodeTimesByName('dot_off',0);


        %%%     get rate    %%%
        r           = getFIRA_rate(trials, getFIRA_spikeByID(uid(i)), bt, et);


        %%%     get ROC for each epoch of data      %%%
        c   = getFIRA_ecodesByName('dot_coh');
        crt = getFIRA_ecodesByName('correct');
        coh = [0 3.2 6.4 12.8 25.6 51.2 99.9]';
        %if sum(crt>=0)>=2*N % look at only the first 2*N trials
        if length(crt)>=2*N % look at only the first 2*N trials
            for j = 1:7
                % whole session
                Lc = c==coh(j) & crt>=0 & Ltsk;
                if any(Lc)
                    r1 = r(Lc & Ldir(:, 1));
                    r2 = r(Lc & Ldir(:, 2));

                    rocs(j,i,1) = rocN(r1, r2, 75);
                    rocn(j,i,1) = length(r1) + length(r2);

                    pn(j,i,1)   = sum(crt(Lc)>=0);
                    ps(j,i,1)   = sum(crt(Lc)==1)/pn(j,i,1);
                end


                % first 1/2
%                 ind  = find(crt>=0);
%                 Lind = zeros(size(crt));
%                 Lind(ind(1:N)) = 1;
%                 Lc = c==coh(j) & crt>=0 & Ltsk & Lind;                
                Lc = c==coh(j) & crt>=0 & Ltsk ...
                      & logical([ones(N,1); zeros(length(c)-N,1)]);
                 if any(Lc)
                    r1 = r(Lc & Ldir(:, 1));
                    r2 = r(Lc & Ldir(:, 2));

                    rocs(j,i,2) = rocN(r1, r2, 75);
                    rocn(j,i,2) = length(r1) + length(r2);

                    pn(j,i,2)   = sum(crt(Lc)>=0);
                    ps(j,i,2)   = sum(crt(Lc)==1)/pn(j,i,2);
                 end

                
                % second 1/2
%                 ind  = find(crt>=0);
%                 Lind = zeros(size(crt));
%                 Lind(ind(N+1:2*N)) = 1;
%                 Lc = c==coh(j) & crt>=0 & Ltsk & Lind;                
                Lc = c==coh(j) & crt>=0 & Ltsk ...
                      & logical([zeros(N,1); ones(N,1); zeros(length(c)-2*N,1)]);
                 if any(Lc)
                    r1 = r(Lc & Ldir(:, 1));
                    r2 = r(Lc & Ldir(:, 2));

                    rocs(j,i,3) = rocN(r1, r2, 75);
                    rocn(j,i,3) = length(r1) + length(r2);

                    pn(j,i,3)   = sum(crt(Lc)>=0);
                    ps(j,i,3)   = sum(crt(Lc)==1)/pn(j,i,3);
                end
            end
        end
    end
end


%% compute threshold for each cell
MTf1 = nans(3,length(fn));
MTs1 = nans(3,length(fn));
MTf2 = nans(3,length(fn));
MTs2 = nans(3,length(fn));
Bef1 = nans(3,length(fn));
Bes1 = nans(3,length(fn));
Bef2 = nans(3,length(fn));
Bes2 = nans(3,length(fn));

for i = 1:length(fn)
    fprintf('%s\n',fn{i})
    % MT
    L = ~isnan(rocs(:,i,2));
    if sum(L)>=5
        [MTf1(:,i) MTs1(:,i)] = quick_fit([coh(L) rocs(L,i,2) rocn(L,i,2)]);
    end
        
    L = ~isnan(rocs(:,i,3));
    if sum(L)>=5
        [MTf2(:,i) MTs2(:,i)] = quick_fit([coh(L) rocs(L,i,3) rocn(L,i,3)]);
    end
    
    % Be
    L = ~isnan(ps(:,i,2));
    if sum(L)>=5
        [Bef1(:,i) Bes1(:,i)] = quick_fit([coh(L) ps(L,i,2) pn(L,i,2)]);
    end
        
    L = ~isnan(ps(:,i,3));
    if sum(L)>=5
        [Bef2(:,i) Bes2(:,i)] = quick_fit([coh(L) ps(L,i,3) pn(L,i,3)]);
    end
end


save(savepath, 'MTf1', 'MTs1', 'MTf2', 'MTs2', 'Bef1', 'Bes1', 'Bef2', 'Bes2')

else
    load(savepath);
end