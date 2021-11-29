function [mpool] = getML_RLFigS1WNormalization(NORM, recompute)
% mean of pooled signal with optimal weight as a function of normalization
% factor
[h] = dirnames;
fn  = ['getML_RLFigS1WNormalization' '.mat'];
savepath = [h, filesep, 'code', filesep, 'matlab', filesep, 'current', ...
    filesep, '01MTLIP', filesep, 'reports', filesep, '3-ReinforcementLearning', ...
    filesep, 'Preparation', filesep, 'mat', filesep, fn];



if recompute
    % get optimal weight
    recompute = 0;

    Monk      = 'Cy';
    NSEN      = 200;
    NTUNE     = 36;
    ELN       = [7]*1e-8;
    DIRS      = [-90,90];
    SIMNUM    = [51:60];

    bb     = 2.32;
    be     = 6.53;
    bw     = (6.53-2.32)/50;
    THBINS = [2.^[bb:bw:be-bw]' 2.^[bb+bw:bw:be]'];

    
    % optimal weights
    NORMSCALE = sqrt(0.02);
    [woavg_c, wod, woc_c, wocd, wo, wocorr, ind] = getML_RLFig3OptimalWSortedByTh(Monk, NSEN, NTUNE, ELN, DIRS, SIMNUM, THBINS, NORMSCALE, recompute);
    clear woavg_c wod woc_c wocd wocorr
    wo  = wo(:,:,1)/NORMSCALE;
    ind = ind(:,1,1);




    %% load data
    %
    mpool  = nans(length(NORM),1);

    %%%
    %   Get statistics of MT responses for each neuron
    %   (0% response, pref gain, null gain, and fano factor)
    %%%
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



    th      = TH(ind);
    pref    = PREF(:,ind);
    null    = NULL(:,ind);
    fano    = FANO(:,ind);


    %%%
    %   create 2d arrays, describing pref/null intercepts and slopes, fano, and
    %   direction tunings for each neuron
    %%%
    DTUNE = 360/NTUNE;
    k0    = repmat(nanmean([pref(1,:); null(1,:)],1)', 1, NTUNE);    % pref and null intercepts are highly correlated (r>0.99), so treat them as equal
    kps   = repmat(pref(2,:)', 1, NTUNE);
    kns   = repmat(null(2,:)', 1, NTUNE);
    kf    = repmat(fano', 1, NTUNE);
    kdir  = repmat(linspace(-180,180-DTUNE,NTUNE), NSEN, 1);



    %% get mean responses as a function of for optimal weight, at
    %% different normalization factor
    atune = normpdf(kdir, 90, 40);
    atune = atune/max(max(atune));          % tuning modulation function
    acoh  = (kps-kns).*atune+kns;           % acoh will be a gaussian function, with peak at pref-slope, and baseline equal to null-slope (which can be a negative number)
    a     = k0 + 0.999*acoh;                  % mean firing of MT cell for this trial


    for i = 1:length(NORM)
        mpool(i)  = sqrt(NORM(i))*wo(:)'*a(:);
    end

    save(savepath, 'mpool')

else
    load(savepath)

end


