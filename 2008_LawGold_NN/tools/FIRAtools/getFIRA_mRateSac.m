function [rateM, rateSD, rateN] = getFIRA_mRateSac(fn, dtype, uid, d1, bins, crtf, nvf, dirf)
% get mean rate for each bins with respect to saccade onset
% INPUTS:
%   fn   - FIRA data filename
%   dtype- 'MT' or 'LIP'
%   uid  - unit id
%   d1   - dots preferred direction
%   bins - time bins for which spikes rate will be calculated
%   crtf - [0 1]=both correct and incorrect, 1=correct only, 0=incorrect only
%   nvf  - [0 1]=both nv and not nv, 1=nv only, 0=not nv only
%   dirf - 0=sorted by dot/trg direction, 1=sorted by choice

if strcmp(dtype, 'MT') | strcmp(dtype, 'PReMT')
    d1name = 'dot_dir';
elseif strcmp(dtype, 'LIP')
    d1name = 'trg1_dir';
end


%% default
rateM  = nans(14, size(bins,1));      % rate mean for each bin and condition
rateSD = nans(14, size(bins,1));      % SD
rateN  = nans(14, size(bins,1));      % N, number of trials in each bin


%% compute spike rate for every good trial, select for dir/coh/correct later
%
% load data
global FIRA
openFIRA(fn)
fprintf('%s\n',fn)


% get firing rate for each trials, only get rate for trial selected by
% correct flag crtf and no var flag nvf
% *******************************************************
% ******ignore nv for now, there's a bug in spm724j******
% *******************************************************  
Ltr = ismember(getFIRA_ecodesByName('correct', 'id'), crtf);
r   = getFIRA_binRate(Ltr, getFIRA_spikeByID(uid),...
                        getFIRA_ecodeTimesByName('sac_on', min(bins(:,1))),...
                        getFIRA_ecodeTimesByName('all_off', 0),...
                        getFIRA_ecodeTimesByName('sac_on', 0),...
                        bins);


% selection arrays for coherences
[Lcoh, Ucoh] = selectFIRA_trialsByUniqueID('dot_coh');
Lcoh = Lcoh(:, ismember(Ucoh, [0 3.2 6.4 12.8 25.6 51.2 99.9])); % remove strange coh other than the standard cohs
Ucoh = Ucoh(ismember(Ucoh, [0 3.2 6.4 12.8 25.6 51.2 99.9]));    % in one pre-training file there's a 9% coh condition

% selection arrays for direction
[Ldir, Udir] = selectFIRA_trialsByUniqueID(d1name);
Idir         = [find(round(Udir)==d1) find(round(Udir)==mod(d1+180,360))];
if dirf  % if neccessary, convert to by choice
    for j = 1:size(Idir)
        L = ~isnan(FIRA.ecodes.data(:,getFIRA_ecodeColumnByName('correct')));
        Ldir(Lt,Idir(j)) = ~xor(Ldir(Lt,Idir(j)), FIRA.ecodes.data(Lt,getFIRA_ecodeColumnByName('correct')));
    end
end

% selection arrays for task
[Ltk, Utk]  = selectFIRA_trialsByUniqueID('task');
Ltk = ismember(getFIRA_ecodesByName('task','id'), [3,6]);

% get selection array for mean rates             
L = zeros(length(Ldir), 14);
c = [0 3.2 6.4 12.8 25.6 51.2 99.9];
for i = 1:7
    if ismember(c(i),Ucoh)
        L(:,i)   = Ldir(:,Idir(1)) & Lcoh(:,Ucoh==c(i)) & Ltk;
        L(:,i+7) = Ldir(:,Idir(2)) & Lcoh(:,Ucoh==c(i)) & Ltk;
    end 
end

% get stats for each conditions
for j = 1:14
    rateM(j,:)  = nanmean(r(logical(L(:,j)),:),1);
    rateSD(j,:) = nanstd(r(logical(L(:,j)),:),1);
    rateN(j,:)  = sum(~isnan(r(logical(L(:,j)),:)),1);
end

