function [TRStat PReStat TRE PReE]= getML_BehResponseStat(monk,recompute)
% get percent of broken fixation and no choice per session
% Output: TRStat  = [pc bf nc bf&nc Time to attain fixation]
%         PReStat = [bf nc bf&nc Time to attain fixation]

savepath = ['/work/mirror_jeff/code/matlab/tmp-mat/getML_BehResponseStat_' monk '.mat'];
if recompute
    % get statistics before training (only for broken fixation and nc)
    a     = getML_txt([monk 'PRe_MT.txt']);
    fname = a.data{find(strcmp(a.name, 'dat_fn'))};
    use   = a.data{find(strcmp(a.name, 'usable'))};
    PReStat = nans(length(fname),4);
    PReE    = nans(length(fname),4);
    
    % get data
    global FIRA
    for i = 1:length(fname)
        if use(i)==1
            openFIRA(fname{i});
            fprintf('%d: %s\n', i, fname{i})

            crt   = getFIRA_ecodesByName('correct');
            [PReStat(i,1) CI] = binofit(nansum(crt==-2),length(crt),1-0.68);
            PReE(i,1)         = CI(2)-PReStat(i,1);
            [PReStat(i,2) CI] = binofit(nansum(crt==-1),length(crt),1-0.68);
            PReE(i,2)         = CI(2)-PReStat(i,2);
            [PReStat(i,3) CI] = binofit(nansum(crt==-1|crt==-2),length(crt),1-0.68);
            PReE(i,3)         = CI(2)-PReStat(i,3);
            PReStat(i,4)      = nanmean(FIRA.ecodes.data(crt==1,getFIRA_ecodeColumnByName('fix_time')))/1000;
            PReE(i,4)         = nanstd(FIRA.ecodes.data(crt==1,getFIRA_ecodeColumnByName('fix_time')))/1000;
            
        end
    end

    
    % get statistics during training
    a     = getML_txt([monk 'TRain_psy.txt']);
    fname = a.data{find(strcmp(a.name, 'dat_fn'))};

    TRStat = nans(length(fname),5);
    TRE    = nans(length(fname),5);
    
    % get data
    global FIRA
    for i = 1:length(fname)
        openFIRA(fname{i});
        fprintf('%d: %s\n', i, fname{i})    
        
        crt   = getFIRA_ecodesByName('correct');
        [TRStat(i,1) CI] = binofit(nansum(crt==1),nansum(crt==0|crt==1),1-0.68);
        TRE(i,1)         = CI(2)-TRStat(i,1);
        [TRStat(i,2) CI] = binofit(nansum(crt==-2),length(crt),1-0.68);
        TRE(i,2)         = CI(2)-TRStat(i,2);
        [TRStat(i,3) CI] = binofit(nansum(crt==-1),length(crt),1-0.68);
        TRE(i,3)         = CI(2)-TRStat(i,3);
        [TRStat(i,4) CI] = binofit(nansum(crt==-1|crt==-2),length(crt),1-0.68);
        TRE(i,4)         = CI(2)-TRStat(i,4);
        TRStat(i,5)      = nanmean(FIRA.ecodes.data(crt==1,getFIRA_ecodeColumnByName('fix_time')))/1000;
        TRE(i,5)         = nanstd(FIRA.ecodes.data(crt==1,getFIRA_ecodeColumnByName('fix_time')))/1000;
            
    end
   
    save(savepath, 'PReStat', 'TRStat', 'PReE', 'TRE')
else
    load(savepath)
    
end



















