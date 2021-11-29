function [lat, dur, vmax, vavg] = getML_meanSac(fn, recompute)
% get means of 'sac_lat' 'sac_dur' 'sac_vmax' 'sac_vavg' per session

if nargin < 1
    return;
elseif nargin < 2
    recompute = 1;
end


savepath = ['/work/mirror_jeff/code/matlab/tmp-mat/getML_meanSac_' fn(1:end-4) '.mat'];

if recompute
    a     = getML_txt(fn);
    fname = a.data{strcmp(a.name,'dat_fn')};
    
    lat   = nans(length(fname),1);
    dur   = nans(length(fname),1);
    vmax  = nans(length(fname),1);
    vavg = nans(length(fname),1);
    
    
    global FIRA
    for i = 1:length(fname)
        fprintf('%d: %s\n', i, fname{i})
        openFIRA(fname{i})

        Lgd     = getFIRA_ecodesByName('correct')>=0;
        
        lat(i)  = nanmean(getFIRA_ecodesByName('sac_lat', 'value', Lgd));
        dur(i)  = nanmean(getFIRA_ecodesByName('sac_dur', 'value', Lgd));
        vmax(i) = nanmean(getFIRA_ecodesByName('sac_vmax', 'value', Lgd));
        vavg(i) = nanmean(getFIRA_ecodesByName('sac_vavg', 'value', Lgd));
    end

    % save
    save(savepath, 'lat', 'dur', 'vmax', 'vavg')
else   
    % load
    load(savepath)
end
    