function [rocs rocn roce] = getML_ROCs(fn, time_bins, crtf, nvf, dirf, marf, recompute)
% get roc for each cell


if nargin < 2
    error('getTR_ROCs needs at least two inputs.')
    return;
end


if ~isempty(findstr(fn, 'LIP'))
    fn_type = 'LIP';

elseif ~isempty(findstr(fn, 'MT'))
    if ~isempty(findstr(fn, 'PRe'))
        fn_type = 'PReMT';
    else
        fn_type = 'MT';
    end    
else
    return;
end


if strcmp(fn_type, 'MT') | strcmp(fn_type, 'PReMT')
    utxtcol = 'ddir';
    firacol = 'dot_dir';
elseif strcmp(fn_type, 'LIP')
    utxtcol = 'trg_dir';
    firacol = 'trg1_dir';
end


if nargin < 3
    time_bins = [0 500];
end

if nargin < 4
    crtf = [0 1];
end

if nargin < 5
    nvf = [0 1];
end

if nargin < 6
    dirf = 0;       % 0 - by direction, 1 - by choice   
end

if nargin < 7
    marf = [];   % marginal flag used during the by choice option, (dirf=1)
                 % if it's = [], both preferred and null responses are used in
                 % by choice, if marf = 1, only preferred response is used and so
                 % on
end


bin_size = nanmean(time_bins(:,2)-time_bins(:,1));
[home_dir_, lab_dir_, current_dir_, tmat_dir_] = dirnames;
savepath = [tmat_dir_ '/getML_ROCs_' fn_type '_' fn(1:2) '_' num2str(bin_size) '.mat'];

if recompute
    global FIRA
    global utxt
    utxt = getML_txt(fn);

    fname  = utxt.data{strcmp(utxt.name,'dat_fn')};
    usable = utxt.data{strcmp(utxt.name,'usable')};
    uid    = utxt.data{strcmp(utxt.name,'uid')};
    d1     = utxt.data{strcmp(utxt.name, utxtcol)};

    rocs   = nans(size(time_bins,1), 7, length(fname));
    roce   = nans(size(time_bins,1), 7, length(fname));
    rocn   = nans(size(time_bins,1), 7, length(fname));
    rocp   = nans(size(time_bins,1), 7, length(fname));

    % compute spike rate for every good trial, select for dir/coh/correct later
    for i=1:length(fname)
        if usable(i)==1
            % load data
            openFIRA(fname{i})
            fprintf('%s\n',fname{i})

            % get roc for each trials, select trials by based on
            % correct flag (crtf) and no var flag (nvf)
            Lcrt =  ismember(getFIRA_ecodesByName('correct', 'id'), crtf);
            trials = find(~isnan(FIRA.ecodes.data(:,1))&Lcrt);


            % get selection arrays for dir/coherence
            % remove strange coh other than the standard cohs (i.e. 0 32 64 128
            % 256 512 999)  (in one pre-training file there's a 9% coh)
            [Lcoh, Ucoh] = selectFIRA_trialsByUniqueID('dot_coh');
            Lcoh = Lcoh(:, ismember(Ucoh, [0 3.2 6.4 12.8 25.6 51.2 99.9]));
            Ucoh = Ucoh(ismember(Ucoh, [0 3.2 6.4 12.8 25.6 51.2 99.9]));


            if dirf % choice by direction or by choice
                [Ldir, Udir] = selectFIRA_trialsByUniqueID('choice');
            else
                [Ldir, Udir] = selectFIRA_trialsByUniqueID(firacol);
            end
            Idir         = [find(round(Udir)==d1(i)) find(round(Udir)==mod(d1(i)+180,360))];
            Ldir         = Ldir(:,Idir);
            

            if ~isempty(marf)   % marginal ROC, ROC for pref responses / null responses only
                if marf==1
                    Ldir   = [Ldir(:,1) & (FIRA.ecodes.data(:,getFIRA_ecodeColumnByName('correct'))==1)...
                        Ldir(:,1) & (FIRA.ecodes.data(:,getFIRA_ecodeColumnByName('correct'))==0)];
                elseif marf==0
                    Ldir   = [Ldir(:,2) & (FIRA.ecodes.data(:,getFIRA_ecodeColumnByName('correct'))==1)...
                        Ldir(:,2) & (FIRA.ecodes.data(:,getFIRA_ecodeColumnByName('correct'))==0)];
                end
            end


            begin_times = getFIRA_ecodeTimesByName('dot_on',0);
            end_times   = getFIRA_ecodeTimesByName('dot_off',0);


            if nargout==4
                [rocs_, rocn_, roce_] = getFIRA_neurometricROCT(trials, Ldir(trials, :), Lcoh(trials, :),...
                    getFIRA_spikeByID(uid(i)), ...
                    begin_times(trials),...
                    end_times(trials),...
                    time_bins, 3);
            else
                [rocs_, rocn_, roce_] = getFIRA_neurometricROCT(trials, Ldir(trials, :), Lcoh(trials, :),...
                    getFIRA_spikeByID(uid(i)), ...
                    begin_times(trials),...
                    end_times(trials),...
                    time_bins, 3);
                rocp_ = nans(size(rocs_));
            end

            %%% added on 1/10/07, try to standardize columns of rocs, some sessions
            %%% didn't use all 7 cohs, but i try to force rocs_ to have 7 columns,
            %%% adding values NaNs for the missing coh
            c = [0 3.2 6.4 12.8 25.6 51.2 99.9];
            for j = 1:7
                if ismember(c(j),Ucoh)
                    rocs(:,j,i) = rocs_(:,Ucoh==c(j));
                    rocn(:,j,i) = rocn_(:,Ucoh==c(j));
                    roce(:,j,i) = roce_(:,Ucoh==c(j));
                    rocp(:,j,i) = rocp_(:,Ucoh==c(j));
                end
            end
        end
    end
    save(savepath, 'rocs', 'rocn', 'roce', 'rocp')

else
    load(savepath)
end
