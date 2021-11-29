function [wo_sorted, wod_sorted, wocorr_sorted, wocorrd_sorted, wo, wocorr, ind] = getML_RLFigS9OptimalWSortedByTh3(Monk, NSEN, NTUNE, ELN, DIRS, SIMNUM, THBINS, NORMSCALE, recompute)
% Get optimal weights, sorted by thresholds
% the optimal weights is computed as the normalized fisher linear discriminant,
% with and without taking into account of the correlations
%
% Correlation independent of sensitivity
%%
[h] = dirnames;
fn  = ['getML_RLFigS9OptimalWSortedByTh3_' Monk '.mat'];
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
    wo     = nans(NSEN, NTUNE, length(SIMNUM));
    wocorr = nans(NSEN, NTUNE, length(SIMNUM));
    ind    = nans(NSEN, NTUNE, length(SIMNUM));
    for nn = 1:length(SIMNUM)
        [hdir, ldir, cdir, tdir] = dirnames;
        fname = ['/getML_RL_' Monk '_' int2str(NSEN) '_' int2str(NTUNE) '_' num2str(ELN) '_v' num2str(SIMNUM(nn)) '.mat'];
        load([tdir fname])
        fprintf([fname '\n']);

        
        %%%
        %   pick randomly 200 neurons from the pool with replacement
        %%%
        th      = TH(MTindex);
        pref    = PREF(:,MTindex);
        null    = NULL(:,MTindex);
        fano    = FANO(:,MTindex);



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



        %%%
        %   create correlation matrix
        %%%
        DTUNE = 360/NTUNE;
        xx    = linspace(-180,180-DTUNE,NTUNE);
        yy    = exp(-(xx+180)/30);
        yy    = [yy(1:NTUNE/2), fliplr(yy(2:NTUNE/2+1))];


        rsen  = 0.5;       % max correlation between neurons with different sensitivity
        rslp  = -rsen/20;  % number of neighbouring neurons that it is correlated to
        R0_   = [];
        th    = th(:);
        for j = 1:NSEN
            %R0_ = [R0_, rsen*exp(-abs(MTindex-MTindex(i))/rtau)];
            %R0_ = [R0_, rsen*exp(-abs([1:NSEN]'-i)/rtau)];
            R0_ = [R0_, rsen*((NSEN-j+1)/NSEN)*ones(size([1:NSEN]'))];
        end
        R0_(R0_<0) = 0;



        R_    = [];
        for j = 1:NTUNE
            R_ = [R_, yy(j)*R0_];
        end
        for j = 1:NSEN      % force diagonal to 1
            R_(j,j) = 1;
        end

        R  = [];
        for j = 1:NTUNE
            R = [R; circshift(R_,[0,(j-1)*NSEN])];
        end
        
        %%%
        %   compute mean and variance for each cell when DIRS(1) and
        %   DIRS(2) are shown
        %%%%
        load([Monk 'Combined.mat'])
        crt  = data(:,3);
        Lgd  = crt>=0 & ~isnan(data(:,6)) & data(:,5)>0;
        coh  = data(Lgd,5)/100;
        COH   = nangeomean(coh);
        clear data coh crt
        
        atune = WrappedGaussian(kdir, 1, 40.^2, DIRS(1));
        atune = atune/max(max(atune));          % tuning modulation function
        acoh  = (kps-kns).*atune+kns;           % acoh will be a gaussian function, with peak at pref-slope, and baseline equal to null-slope (which can be a negative number)
        m1    = k0 + COH*acoh;                      % mean firing of MT cell for this trial
        v1    = kf.*abs(m1);

        atune = WrappedGaussian(kdir, 1, 40.^2, DIRS(2));
        atune = atune/max(max(atune));          % tuning modulation function
        acoh  = (kps-kns).*atune+kns;           % acoh will be a gaussian function, with peak at pref-slope, and baseline equal to null-slope (which can be a negative number)
        m2    = k0 + COH*acoh;                      % mean firing of MT cell for this trial
        v2    = kf.*abs(m2);

        m     = m2-m1;
        s     = sqrt(v1+v2);

        clear acoh atune R_ R0_ k0 kdir kf kns kps kps m1 m2 v1 v2 U
        

    
        % get fisher linear discriminant taking into account of correlations
        R             = s(:)*s(:)'.*R;
        Rinv          = inv(R);
        w_            = Rinv*m(:);
        w_            = NORMSCALE*w_./sqrt(sum(sum(w_.^2)));
        wocorr(:,:,nn) = reshape(w_,NSEN,NTUNE);
        
        
%         % get fisher linear discriminant without taking into account of correlations
%         R             = diag(diag(R));
%         Rinv          = inv(R);
%         w_            = Rinv*m(:);
%         w_            = NORMSCALE*w_./sqrt(sum(sum(w_.^2)));
%         wo(:,:,i)     = reshape(w_,NSEN,NTUNE);
         
        ind(:,:,nn)    = repmat(MTindex',1,NTUNE);
    end
    th = TH(ind);

    
    
    % compute average weights binned by thresholds
    wo_sorted       = nans(length(THBINS), NTUNE);
    wod_sorted      = nans(length(THBINS), NTUNE);
    wocorr_sorted   = nans(length(THBINS), NTUNE);
    wocorrd_sorted  = nans(length(THBINS), NTUNE);

    for i = 1:NTUNE
        wo_   = wo(:,i,:);
        wo_   = wo_(:);
        woc_  = wocorr(:,i,:);
        woc_  = woc_(:);
        th_   = 100*th(:,i,:);
        th_   = th_(:);

        for j = 1:length(THBINS)
            L                   = th_>=THBINS(j,1) & th_<THBINS(j,2);
            wo_sorted(j,i)      = nanmean(wo_(L));
            wod_sorted(j,i)     = nanse(wo_(L));
            wocorr_sorted(j,i)  = nanmean(woc_(L));
            wocorrd_sorted(j,i) = nanse(woc_(L));
        end
    end
    save(savepath, 'wo_sorted', 'wod_sorted', 'wocorr_sorted', 'wocorrd_sorted', 'wo', 'wocorr', 'ind')

else
    load(savepath)

end

