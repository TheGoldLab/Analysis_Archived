function [cp cpd cpraw] = getML_RLFigS3CorrCP(Monk, DIRS, rmax, bdir, bsen, N, recompute)
% Test how correlation structure affects choice probability


[h] = dirnames;
fn  = ['getML_RLFigS3CorrCP_' Monk '_' 'rmax', int2str(length(rmax))  ...
    '_' 'bdir', int2str(length(bdir))  '_' 'bsen', int2str(length(bsen))  '.mat'];
savepath = [h, filesep, 'code', filesep, 'matlab', filesep, 'current', ...
    filesep, '01MTLIP', filesep, 'reports', filesep, '3-ReinforcementLearning', ...
    filesep, 'Preparation', filesep, 'mat', filesep, fn];





if recompute

    %%%
    %   Load task info (the actual coh/dir/viewing-time used during training)
    %%%
    load([Monk 'Combined.mat'])
    crt  = data(:,3);
    Lgd  = crt>=0 & ~isnan(data(:,6));
    crt  = crt(Lgd);
    coh  = data(Lgd,5)/100;
    clear data

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
    TH   = [fitsp(1,:) fits(1,:)];
    PREF = [pp p];
    NULL = [np n];
    FANO = [fp f];

    if strcmp(Monk, 'Cy')
        Lgd      = ~isnan(PREF(1,:)) & ~isnan(NULL(1,:)) & ~isnan(FANO(1,:));
        Lgd(178) = logical(0);  % this cell is weird, with very small gain to both pref and null
    else
        Lgd      =~isnan(PREF(1,:)) & ~isnan(NULL(1,:)) & ~isnan(FANO(1,:));
    end
    TH   = TH(Lgd);
    PREF = PREF(:,Lgd);
    NULL = NULL(:,Lgd);
    FANO = FANO(:,Lgd);

    % sort data by threshold coherence
    [TH,xxx] = sort(TH);
    PREF     = PREF(:,xxx);
    NULL     = NULL(:,xxx);
    FANO     = FANO(:,xxx);


    %%%
    %   create array of neuron base on MTINDEX
    %%%
    NSEN    = 200;
    NTUNE   = 36;
    MTindex = sort(unidrnd(length(TH),1,NSEN));
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




    %% sim
    cp    = nans(length(rmax), length(bdir), length(bsen));
    cpd   = nans(length(rmax), length(bdir), length(bsen));
    cpraw = cell(length(rmax), length(bdir), length(bsen));
    rALL  = nans(NSEN, NTUNE, N);

    for rr1 = 1:length(bsen)
        for rr2 = 1:length(bdir)
            for rr3 = 1:length(rmax)


                %%%
                %   create correlation matrix
                %%%
                DTUNE = 360/NTUNE;
                xx    = linspace(-180,180-DTUNE,NTUNE);
                yy    = exp(-(xx+180)/bdir(rr2));
                yy    = [yy(1:NTUNE/2), fliplr(yy(2:NTUNE/2+1))];


                rsen  = rmax(rr3);       % max correlation between neurons with different sensitivity
                rslp  = -rsen/bsen(rr1);  % number of neighbouring neurons that it is correlated to
                R0_   = [];
                th    = th(:);
                for i = 1:NSEN
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



                if rr1==1   % only compute optimal weight once, because the 'covariance' matrix i used to compute optimal weight ignores correlations
                    % compute optimal weights
                    %%%
                    %   compute mean and variance for each cell when DIRS(1) and
                    %   DIRS(2) are shown
                    %%%%
                    COHM = nangeomean(coh(coh~=0));

                    atune = WrappedGaussian(kdir, 1, 40.^2, DIRS+180);
                    atune = atune/max(max(atune));          % tuning modulation function
                    acoh  = (kps-kns).*atune+kns;           % acoh will be a gaussian function, with peak at pref-slope, and baseline equal to null-slope (which can be a negative number)
                    m1    = k0 + COHM*acoh;                      % mean firing of MT cell for this trial
                    v1    = kf.*abs(m1);

                    atune = WrappedGaussian(kdir, 1, 40.^2, DIRS);
                    atune = atune/max(max(atune));          % tuning modulation function
                    acoh  = (kps-kns).*atune+kns;           % acoh will be a gaussian function, with peak at pref-slope, and baseline equal to null-slope (which can be a negative number)
                    m2    = k0 + COHM*acoh;                      % mean firing of MT cell for this trial
                    v2    = kf.*abs(m2);

                    m     = m2-m1;
                    s     = sqrt(v1+v2);

                    clear m1 m2 v1 v2

                    % get fisher linear discriminant without taking into account of correlations
                    NORMSCALE     = sqrt(0.02);
                    R_            = diag(diag(R));
                    Rinv          = inv(R_);
                    w             = Rinv*m(:);
                    w             = NORMSCALE*w./sqrt(sum(sum(w.^2)));
                    w             = reshape(w,NSEN,NTUNE);
                    clear R_
                end


                % NOTE: the matlab CHOL function decompose R into U*U', where U is a lower triangular matrix.
                %       We thus can use U to "correlate" the normal random deviates
                %       responses of MT neurons (this is analogous to computing the matrix
                %       square root of R and use the square root to correlate the
                %       responses)
                U      = chol(R, 'lower');
                clear R xx yy;




                [rr1 rr2 rr3]
                choice = nans(N,1);
                for k = 1:N
                    if ~mod(k,50)
                        fprintf('%d\n',k)
                    end

                    %%%   [1] plan action:
                    %           - pool signals
                    %           - form direction decision
                    %%%
                    % Stimulus condition
                    % always 0% coh, as a result, dir is meaningless
                    % so average response of the neuron will just be the response at 0%
                    % coh, or K0.
                    a     = k0;       % mean firing of MT cell for this trial


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
                    pool   = pool+normrnd(0,sqrt(2*abs(pool(i)))+5);


                    % decision
                    choice(k) = pool>0;


                    % store response
                    rALL(:,:,k) = r;
                end

                % compute choice probability for neurons tuned to DIR
                cp_ = nans(NSEN,1);
                fprintf('Computing choice probability\n');

                LR = choice==1;
                LL = choice==0;
                for cc = 1:NSEN
                    r       = shiftdim(rALL(cc,linspace(-180,180-DTUNE,NTUNE)==0,:),2);
                    cp_(cc) = rocN(r(LR), r(LL), 100);
                end

                cp(rr3,rr2,rr1)    = nanmean(cp_);
                cpd(rr3,rr2,rr1)   = nanstd(cp_);
                cpraw{rr3,rr2,rr1} = cp_;
            end
        end
    end



    save(savepath, 'cp', 'cpd', 'cpraw')

else
    load(savepath)

end



