function [mTR, sdTR, nTR, mPR, sdPR, nPR] = getML_MTRate(monk, recompute)
% DESCRIPTION:
%  Get average rate for MT neurons for each stimulus conditions.
%  There are in total 13 stimulus conditions: 6 for each coh at the pref dir,
%  6 for each coh at the null dir, and 0% coh for both dir.
%
% OUTPUT:
%  Each of the following output are mxn matrix, where m is the number of
%  coh and n is the number of cell. The m rows contains rate at
%  [99.9 51.2 25.6 12.8 6.4 3.2 0 -3.2 -6.4 -12.8 -25.6 -51.2 -99.9]
%  percent coh motion, where -ve values mean motion to the neuron's null direction.
%
%  mTR  - mean for neurons recorded during training
%  sdTR - sd for neurons recorded during training
%  nTR  - number of trials
%  mPR  - mean for neurons recorded before training
%  sdPR - sd for neurons recorded before training
%  nPR  - number of trials


savepath = ['/work/mirror_jeff/code/matlab/tmp-mat/getML_MTRate_' monk '.mat'];

if recompute
    % pre-training
    a      = getML_txt([monk 'PRe_MT.txt']);
    fn     = a.data{strcmp(a.name,'dat_fn')};
    usable = a.data{strcmp(a.name,'usable')};
    uid    = a.data{strcmp(a.name,'uid')};
    d1     = a.data{strcmp(a.name,'ddir')};

    mPR    = nans(13, length(fn));
    sdPR   = nans(13, length(fn));
    nPR    = nans(13, length(fn));

    warning off
    global FIRA
    for i = 1:length(fn)
        if usable(i) == 1;
            openFIRA(fn{i})
            fprintf('%s\n',fn{i})


            % get mean rate from dots_on to dots_off
            % selection arrays for coherences
            [Lcoh, Ucoh] = selectFIRA_trialsByUniqueID('dot_coh');
            Lcoh = Lcoh(:, ismember(Ucoh, [0 3.2 6.4 12.8 25.6 51.2 99.9])); % remove strange coh other than the standard cohs
            Ucoh = Ucoh(ismember(Ucoh, [0 3.2 6.4 12.8 25.6 51.2 99.9]));    % in one pre-training file there's a 9% coh condition

            % selection arrays for direction
            [Ldir, Udir] = selectFIRA_trialsByUniqueID('dot_dir');
            Idir         = [find(round(Udir)==d1(i)) find(round(Udir)==mod(d1(i)+180,360))];

            % selection arrays for task
            [Ltk, Utk]  = selectFIRA_trialsByUniqueID('task');
            Ltk = ismember(getFIRA_ecodesByName('task','id'), [3,6]);

            % selection array for correct/incorrect
            Lcrt =  ismember(getFIRA_ecodesByName('correct', 'id'), [0 1]);

            % get selection array for mean rates
            L = zeros(length(Ldir), 13);
            c = [0 3.2 6.4 12.8 25.6 51.2 99.9];
            for j = 7:-1:2          % pref motion
                if ismember(c(j),Ucoh)
                    L(:,7-j+1) = Ldir(:,Idir(1)) & Lcoh(:,Ucoh==c(j)) & Ltk & Lcrt;
                end
            end

            if ismember(c(1),Ucoh)  % 0 %coh
                L(:,7) = (Ldir(:,Idir(1)) | Ldir(:,Idir(2))) & Lcoh(:,Ucoh==c(1)) & Ltk & Lcrt;
            end

            for j = 2:7             % null motion
                if ismember(c(j),Ucoh)
                    L(:,j+6) = Ldir(:,Idir(2)) & Lcoh(:,Ucoh==c(j)) & Ltk & Lcrt;
                end
            end


            % get rate
            for j = 1:13
                r = getFIRA_rate(L(:,j), getFIRA_spikeByID(uid(i)), ...
                    getFIRA_ecodeTimesByName('dot_on', 0), ...
                    getFIRA_ecodeTimesByName('dot_off', 0));
                mPR(j,i)  = nanmean(r);
                nPR(j,i)  = sum(~isnan(r));
                sdPR(j,i) = nanstd(r);
            end
        end
    end
    warning on

    
    
    % training
    a      = getML_txt([monk 'TRain_MT.txt']);
    fn     = a.data{strcmp(a.name,'dat_fn')};
    usable = a.data{strcmp(a.name,'usable')};
    uid    = a.data{strcmp(a.name,'uid')};
    d1     = a.data{strcmp(a.name,'ddir')};

    mTR    = nans(13, length(fn));
    sdTR   = nans(13, length(fn));
    nTR    = nans(13, length(fn));

    warning off
    global FIRA
    for i = 1:length(fn)
        if usable(i) == 1;
            openFIRA(fn{i})
            fprintf('%s\n',fn{i})


            % get mean rate from dots_on to dots_off
            % selection arrays for coherences
            [Lcoh, Ucoh] = selectFIRA_trialsByUniqueID('dot_coh');
            Lcoh = Lcoh(:, ismember(Ucoh, [0 3.2 6.4 12.8 25.6 51.2 99.9])); % remove strange coh other than the standard cohs
            Ucoh = Ucoh(ismember(Ucoh, [0 3.2 6.4 12.8 25.6 51.2 99.9]));    % in one pre-training file there's a 9% coh condition

            % selection arrays for direction
            [Ldir, Udir] = selectFIRA_trialsByUniqueID('dot_dir');
            Idir         = [find(round(Udir)==d1(i)) find(round(Udir)==mod(d1(i)+180,360))];

            % selection arrays for task
            [Ltk, Utk]  = selectFIRA_trialsByUniqueID('task');
            Ltk = ismember(getFIRA_ecodesByName('task','id'), [3,6]);

            % selection array for correct/incorrect
            Lcrt =  ismember(getFIRA_ecodesByName('correct', 'id'), [0 1]);

            % get selection array for mean rates
            L = zeros(length(Ldir), 13);
            c = [0 3.2 6.4 12.8 25.6 51.2 99.9];
            for j = 7:-1:2          % pref motion
                if ismember(c(j),Ucoh)
                    L(:,7-j+1) = Ldir(:,Idir(1)) & Lcoh(:,Ucoh==c(j)) & Ltk & Lcrt;
                end
            end

            if ismember(c(1),Ucoh)  % 0 %coh
                L(:,7) = (Ldir(:,Idir(1)) | Ldir(:,Idir(2))) & Lcoh(:,Ucoh==c(1)) & Ltk & Lcrt;
            end

            for j = 2:7             % null motion
                if ismember(c(j),Ucoh)
                    L(:,j+6) = Ldir(:,Idir(2)) & Lcoh(:,Ucoh==c(j)) & Ltk & Lcrt;
                end
            end


            % get rate
            for j = 1:13
                r = getFIRA_rate(L(:,j), getFIRA_spikeByID(uid(i)), ...
                    getFIRA_ecodeTimesByName('dot_on', 0), ...
                    getFIRA_ecodeTimesByName('dot_off', 0));
                mTR(j,i)  = nanmean(r);
                nTR(j,i)  = sum(~isnan(r));
                sdTR(j,i) = nanstd(r);
            end
        end
    end
    warning on
    
    
    save(savepath, 'mTR', 'sdTR', 'nTR', 'mPR', 'sdPR', 'nPR')
else

    load(savepath)
end



