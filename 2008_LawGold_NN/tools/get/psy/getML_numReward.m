function [rc rcd rcb rcbd] = getML_numReward(fn, recompute)

if nargin < 1
    return;
elseif nargin < 2
    recompute = 1;
end

savepath = ['/work/mirror_jeff/code/matlab/tmp-mat/getML_numReward' fn(1:2) '.mat'];

if recompute

    a     = getML_txt(fn);
    fname = a.data{strcmp(a.name,'dat_fn')};
    rc     = nans(length(fname),1);
    rcd    = nans(length(fname),1);

    % get fits
    global FIRA
    for i = 1:length(fname)
        fprintf('%d: %s\n', i, fname{i})
        openFIRA(fname{i})
        re   = getFIRA_ecodesByName('num_rewards');
        if length(re)>201
            re(1:end-200) = nan;
        end

        c    = getFIRA_ecodesByName('correct');
        rc(i)  = nanmean(re(c>=0));
        rcd(i) = nanse(re(c>=0));
    end

    % get average pre-training
    a      = getML_txt([fn(1:2) 'PRe_MT.txt']);
    fname  = a.data{strcmp(a.name,'dat_fn')};
    rcb     = nans(length(fname),1);
    rcbd    = nans(length(fname),1);

    % get fits
    global FIRA
    for i = 1:length(fname)
        fprintf('%d: %s\n', i, fname{i})
        openFIRA(fname{i})
        re   = getFIRA_ecodesByName('num_rewards');
        c    = getFIRA_ecodesByName('correct');      
        rcb(i)  = nanmean(re(c>=0));
        rcbd(i) = nanse(re(c>=0));
    end
    
    save(savepath, 'rc', 'rcd', 'rcb', 'rcbd')
else
    load(savepath)
end

