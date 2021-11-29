[wavg, wse] = getML_RLFig5PooledSNR(Monk, NSEN, NTUNE, ELN, DIRS, SIMNUM, recompute)
% get SNR of pooled signals during training averaged across simulations

%%
[h] = dirnames;
fn  = ['getML_RLFig5PooledSNR_' Monk '.mat'];
savepath = [h, filesep, 'code', filesep, 'matlab', filesep, 'current', ...
    filesep, '01MTLIP', filesep, 'reports', filesep, '3-ReinforcementLearning', ...
    filesep, 'Preparation', filesep, 'mat', filesep, fn];




if recompute
    for i = 1:length(SIMNUM)
        % get time constant for lapse
        [hdir, ldir, cdir, tdir] = dirnames;
        fname = ['/getML_RL_' Monk '_' int2str(NSEN) '_' int2str(NTUNE) '_' num2str(ELN) '_v' num2str(SIMNUM(i)) '.mat'];
        load([tdir fname])
        fprintf([fname '\n']);

        N      = max(FIRA(:,1));
        coh    = [coh, FIRA(:,2)];
        dir    = [dir, FIRA(:,3)];
        vt     = [vt, FIRA(:,4)];
        crt    = [crt, FIRA(:,9)];
        x      = [x, [1:N]'];
    end

    trialnum  = floor(TRIALS/1000)+1;
  
    % get threshold of neurons
    [fitsp, semsp] = getML_neurometricROCT([Monk 'PRe_MT.txt'], [], 1, [0 1], 0, 0);
    [fits, sems]   = getML_neurometricROCT([Monk 'TRain_MT.txt'], [], 1, [0 1], 0, 0);
    [pp, np, fp]   = getML_SigNoise2([Monk 'PRe_MT.txt'], 1, [0 1], 1, 0);
    [p, n, f]      = getML_SigNoise2([Monk 'TRain_MT.txt'], 1, [0 1], 1, 0);

    % combine pre- and during training data
    TH   = [fitsp(1,:) fits(1,:)];
    Lgd  = ~isnan(TH);
    TH   = TH(Lgd);
    PREF = [pp p];
    PREF = PREF(:,Lgd);
    NULL = [np n];
    NULL = NULL(:,Lgd);
    FANO = [fp f];
    FANO = FANO(:,Lgd);

    % sort data by threshold coherence
    [TH,xxx] = sort(TH);
    PREF     = PREF(:,xxx);
    NULL     = NULL(:,xxx);
    FANO     = FANO(:,xxx);

    % get weights for early/mid/late trial, for each simulation
    w   = nans(NSEN, NTUNE, length(SIMNUM), length(trialnum));
    ind = nans(NSEN, NTUNE, length(SIMNUM), length(trialnum));
    for i = 1:length(SIMNUM)
        [hdir, ldir, cdir, tdir] = dirnames;
        fname = ['/getML_RL_' Monk '_' int2str(NSEN) '_' int2str(NTUNE) '_' num2str(ELN) '_v' num2str(SIMNUM(i)) '.mat'];
        load([tdir fname])
        fprintf([fname '\n']);

        for j = 1:length(trialnum)
            w(:,:,i,j)   = W(:,:,trialnum(j));
            ind(:,:,i,j) = repmat(MTindex',1,NTUNE);
        end
    end
    th = TH(ind);


    % compute average weights binned by thresholds
    wavg = nans(length(THBINS), NTUNE, length(trialnum));
    wse  = nans(length(THBINS), NTUNE, length(trialnum));

    % average weight for each threshold bin
    for i = 1:length(trialnum)
        for j = 1:NTUNE
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
    save(savepath, 'wavg', 'wse')

else
    load(savepath)

end












