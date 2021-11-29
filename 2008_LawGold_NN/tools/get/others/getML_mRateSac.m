function [rM, rSD, rN, rMn] = getML_mRateSac(fname, bins, crtf, nvf, dirf, recompute)
% get mean rate for each good session
% INPUTS:
%   fname - 'txt' filename
%   dtype - 'MT' or 'LIP'
%   bins  - time bins for which spikes rate will be calculated
%   crtf  - [0 1]=both correct and incorrect, 1=correct only, 0=incorrect only
%   nvf   - [0 1]=both nv and not nv, 1=nv only, 0=not nv only
%   dirf  - 0=sorted by dot/trg direction, 1=sorted by choice

if ~isempty(findstr(fname, 'LIP'))
    dtype = 'LIP';

elseif ~isempty(findstr(fname, 'MT'))
        dtype = 'MT';
else
    return;
end


bin_size = bins(:,2)-bins(:,1);
bin_step = bins(1,1)-bins(2,1);
savepath = ['/work/mirror_jeff/code/matlab/tmp-mat/getML_mRateSac_' dtype '_' fname(1:2) '_' num2str(nanmean(bin_size)) '_' num2str(nanmean(bin_step)) '.mat'];


if recompute
    if strcmp(dtype, 'MT') | strcmp(dtype, 'PReMT')
        d1name = 'ddir';
    elseif strcmp(dtype, 'LIP')
        d1name = 'trg_dir';
    end

    % get average response before dots on
    a      = getML_txt(fname);
    fn     = a.data{strcmp(a.name,'dat_fn')};
    usable = a.data{strcmp(a.name,'usable')};
    uid    = a.data{strcmp(a.name,'uid')};
    d1     = a.data{strcmp(a.name,d1name)};
    if strcmp(dtype, 'MT')
        train  = a.data{strcmp(a.name,'train')};
    else
        train = ones(size(d1));
    end
   
    
    rM     = nans(14,size(bins,1),length(fn));
    rMn    = nans(14,size(bins,1),length(fn));
    rSD    = nans(14,size(bins,1),length(fn));
    rN     = nans(14,size(bins,1),length(fn));

    for i = 1:length(fn)
        if usable(i) == 1 & train(i) == 1
            % get mean rate
            [rM(:,:,i), rSD(:,:,i), rN(:,:,i)] = getFIRA_mRateSac(fn{i}, dtype, uid(i), d1(i), bins, crtf, nvf, dirf);
            rMn(:,:,i)                         = getML_normalizedRate(mean(bins,2), rM(:,:,i));
        end
    end

    save(savepath, 'rM', 'rSD', 'rN', 'rMn')
else

    load(savepath)
end