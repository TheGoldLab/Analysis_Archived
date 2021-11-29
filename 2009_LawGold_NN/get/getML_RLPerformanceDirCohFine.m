function [p] = getML_RLPerformanceDirCohFine(Monk, NSEN, NTUNE, DIRS, COHS, w, MTindex, N, TAG, RECOMPUTE)
% Simulate percent correct for a particular weight profile for each
% direction in DIRS and coherence in COHS

[hdir, ldir, cdir, tdir] = dirnames;
savepath = [tdir '/getML_RLPerformanceDirCohFine_' Monk '_'  TAG '.mat'];

if RECOMPUTE

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
    
    
    DTUNE = 360/NTUNE;
    k0    = repmat(nanmean([PREF(1,:); NULL(1,:)],1)', 1, NTUNE);    % pref and null intercepts are highly correlated (r>0.99), so treat them as equal
    kps   = repmat(PREF(2,:)', 1, NTUNE);
    kns   = repmat(NULL(2,:)', 1, NTUNE);
    kf    = repmat(FANO', 1, NTUNE);
    kdir  = repmat(linspace(-180,180-DTUNE,NTUNE), length(TH), 1);
    
    
    % sort data by threshold computed from d-prime
    [TH,xxx] = sort(TH);
    PREF     = PREF(:,xxx);
    NULL     = NULL(:,xxx);
    FANO     = FANO(:,xxx);



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
    for i = 1:NSEN
        %R0_ = [R0_, rsen*exp(-abs(MTindex-MTindex(i))/rtau)];
        %R0_ = [R0_, rsen*exp(-abs([1:NSEN]'-i)/rtau)];
        R0_ = [R0_, rsen+rslp*abs([1:NSEN]'-i)];
    end
    R0_(R0_<0) = 0;
  
    
    
    R_    = [];
    for i = 1:NTUNE
        R_ = [R_, yy(i)*R0_];
    end
    for i = 1:NSEN      % force diagonal to 1
        R_(i,i) = 1;
    end
    
    R  = [];
    for i = 1:NTUNE
        R = [R; circshift(R_,[0,(i-1)*NSEN])];
    end

    % NOTE: the matlab CHOL function decompose R into U*U', where U is a lower triangular matrix.
    %       We thus can use U to "correlate" the normal random deviates
    %       responses of MT neurons (this is analogous to computing the matrix
    %       square root of R and use the square root to correlate the
    %       responses)
    U      = chol(R, 'lower');
    clear R xx yy;



    %%%
    %   Define variables
    %%%
    lambda  = 1;                                             % 'penalty' factor
    p       = nans(size(DIRS,2),length(COHS));
    
    

    %%%
    %   Simulate performance for each dir and coh
    %%%
    for i = 1:size(DIRS,2)
        for j = 1:length(COHS)
            rew = [];
            [i j]
            for nn = 1:N
                
                if ~mod(nn,50)
                    fprintf([int2str(DIRS(1,i)) '_' num2str(COHS(j)) '_' int2str(nn) '\n'])
                end
                
                %%%   [1] plan action:
                %           - pool signals
                %           - form direction decision
                %%%
                % compute a, the mean respose of neuron given trial
                % condition
                dir_ = DIRS(unidrnd(2),i);
                T     = 1;
                atune = WrappedGaussian(kdir, 1, 40.^2, dir_);
                atune = atune/max(max(atune));          % tuning modulation function
                acoh  = (kps-kns).*atune+kns;           % acoh will be a gaussian function, with peak at pref-slope, and baseline equal to null-slope (which can be a negative number)
                a     = T*(k0 + COHS(j)*acoh);          % mean firing of MT cell for this trial

                
                % compute neural responses (normal deviates) for current trial
                r   = normrnd(zeros(size(a)), ones(size(a)));

                % multiple normal deviates with correlation matrix to generate
                % correlated responses, then scale it with the appropriate mean and
                % variance
                r      = reshape(U*r(:),NSEN,NTUNE);
                r      = a+sqrt(kf.*abs(a)).*r;
                r(r<0) = 0;

                % pool signals (add constant pooling noise)
                pool   = sum(sum(w.*r));
                pool   = pool+normrnd(0,sqrt(2*abs(pool))+5);


                % decision
                choice  = pool>0;
                rew(nn) = choice==(sign(dir_-nanmean(DIRS(:,i)))+1)/2;
            end
            p(i,j) = nansum(rew)/N;
        end
    end
       
    
    % save data
    save(savepath, 'p');

else
    load(savepath)
end

