function [wavg, wse, wpavg, wpse, wnavg, wnse] = getML_RLFigS6AvgWSortedByTh(Monk, NSEN, NTUNE, ELN, DIRS, SIMNUM, TRIALS, THBINS, TAG, recompute)
% Get average weights for each trial in TRIALS bined by THBINS

%%
[h] = dirnames;
fn  = ['getML_RLFigS6AvgWSortedByTh_' TAG '_' Monk '.mat'];
savepath = [h, filesep, 'code', filesep, 'matlab', filesep, 'current', ...
    filesep, '01MTLIP', filesep, 'reports', filesep, '3-ReinforcementLearning', ...
    filesep, 'Preparation', filesep, 'mat', filesep, fn];


if recompute
    % get weights
    if strcmp(Monk, 'Cy')
        N = 115985;
    else
        N = 77447;
    end

    trialnum  = floor(TRIALS/1000)+1;
  
    
    % get threshold of neurons
    bb   = 100;
    be   = 1000;
    bw   = 100;
    bs   = 25;
    bins = [zeros(size([bb+bw:bs:be]')) [bb+bw:bs:be]'];

    [fitsp, semsp] = getML_neurometricROCT([Monk 'PRe_MT.txt'], bins, 1, [0 1], 0, 0);
    [fits, sems]   = getML_neurometricROCT([Monk 'TRain_MT.txt'], bins, 1, [0 1], 0, 0);
    [pp, np, fp  ] = getML_SigNoise2([Monk 'PRe_MT.txt'], 1, [0 1], 1, 0);
    [p, n, f]      = getML_SigNoise2([Monk 'TRain_MT.txt'], 1, [0 1], 1, 0);

    % combine pre- and during training data
    PREF = [pp p];
    NULL = [np n];
    FANO = [fp f];
    TH   = getML_RLSimulatedThForMT(Monk, 10000, 0);
    if strcmp(Monk, 'Cy')
        Lgd      = ~isnan(PREF(1,:)) & ~isnan(NULL(1,:)) & ~isnan(FANO(1,:));
        Lgd(178) = logical(0);  % this cell is weird, with very small gain to both pref and null
    else
        Lgd      =~isnan(PREF(1,:)) & ~isnan(NULL(1,:)) & ~isnan(FANO(1,:));
    end
    PREF = PREF(:,Lgd);
    NULL = NULL(:,Lgd);
    FANO = FANO(:,Lgd);
    TH   = TH(Lgd);

    % sort data by threshold computed from d-prime
    [TH,xxx] = sort(TH);
    PREF     = PREF(:,xxx);
    NULL     = NULL(:,xxx);
    FANO     = FANO(:,xxx);
    
    
    % get weights for early/mid/late trial, for each simulation
    w   = nans(NSEN, NTUNE, length(SIMNUM), length(trialnum));
    wp  = nans(NSEN, NTUNE, length(SIMNUM), length(trialnum));
    wn  = nans(NSEN, NTUNE, length(SIMNUM), length(trialnum));
    ind = nans(NSEN, NTUNE, length(SIMNUM), length(trialnum));
    for i = 1:length(SIMNUM)
        [hdir, ldir, cdir, tdir] = dirnames;
        if strcmp(TAG,'ONEPOOL')
            fname = ['/getML_RL_' Monk '_' int2str(NSEN) '_' int2str(NTUNE) '_' num2str(ELN) '_v' num2str(SIMNUM(i)) '.mat'];
            load([tdir fname])
            fprintf([fname '\n']);

            for j = 1:length(trialnum)
                w(:,:,i,j)   = W(:,:,trialnum(j));
                ind(:,:,i,j) = repmat(MTindex',1,NTUNE);
            end
        elseif strcmp(TAG,'TWOPOOLS')
            fname = ['/getML_RLTwoLIPPools_' Monk '_' int2str(NSEN) '_' int2str(NTUNE) '_' num2str(ELN) '_v' num2str(SIMNUM(i)) '.mat'];
            load([tdir fname])
            fprintf([fname '\n']);

            
            for j = 1:length(trialnum)
%                 % renormalize weights    
%                 wss  = sqrt(sum(sum(WP(:,:,trialnum(j)).^2)));
%                 WP(:,:,trialnum(j))   = sqrt(0.01)*WP(:,:,trialnum(j))./wss;
%                 wss  = sqrt(sum(sum(WN(:,:,trialnum(j)).^2)));
%                 WN(:,:,trialnum(j))   = sqrt(0.01)*WN(:,:,trialnum(j))./wss;
                wp(:,:,i,j)   = WP(:,:,trialnum(j));
                wn(:,:,i,j)   = WN(:,:,trialnum(j));
                w(:,:,i,j)    = WP(:,:,trialnum(j))-WN(:,:,trialnum(j));
                ind(:,:,i,j)  = repmat(MTindex',1,NTUNE);
            end
        end    
        
    end
    th = TH(ind);


    % compute average weights binned by thresholds
    wavg  = nans(length(THBINS), NTUNE, length(trialnum));
    wse   = nans(length(THBINS), NTUNE, length(trialnum));
    wpavg = nans(length(THBINS), NTUNE, length(trialnum));
    wpse  = nans(length(THBINS), NTUNE, length(trialnum));
    wnavg = nans(length(THBINS), NTUNE, length(trialnum));
    wnse  = nans(length(THBINS), NTUNE, length(trialnum));

    % average weight for each threshold bin
    for i = 1:length(trialnum)
        for j = 1:NTUNE
            if strcmp(TAG, 'TWOPOOLS')
                wp_   = wp(:,j,:,i);
                wp_   = wp_(:);
                wn_   = wn(:,j,:,i);
                wn_   = wn_(:);
                th_  = 100*th(:,j,:,i);
                th_  = th_(:);

                for k = 1:length(THBINS)
                    L            = th_>=THBINS(k,1) & th_<THBINS(k,2);
                    wpavg(k,j,i) = nanmean(wp_(L));
                    wpse(k,j,i)  = nanse(wp_(L));
                    wnavg(k,j,i) = nanmean(wn_(L));
                    wnse(k,j,i)  = nanse(wn_(L));
                end
            end

            w_   = w(:,j,:,i);
            w_   = w_(:);
            th_  = 100*th(:,j,:,i);
            th_  = th_(:);

            for k = 1:length(THBINS)
                L           = th_>=THBINS(k,1) & th_<THBINS(k,2);
                wavg(k,j,i) = nanmean(w_(L));
                wse(k,j,i)  = nanse(w_(L));
            end
        end
    end
    save(savepath, 'wavg', 'wse', 'wpavg', 'wpse', 'wnavg', 'wnse')

else
    load(savepath)

end








