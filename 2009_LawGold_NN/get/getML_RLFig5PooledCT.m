function [mpool] = getML_RLFig5PooledCT(Monk, NSEN, NTUNE, ELN, SIMNUM, COHS, TIMES, DIRS_, recompute)
% compute pooled responses as a function of coh and time 
[h] = dirnames;
fn  = ['getML_RLFig5PooledCT_' Monk '.mat'];
savepath = [h, filesep, 'code', filesep, 'matlab', filesep, 'current', ...
    filesep, '01MTLIP', filesep, 'reports', filesep, '3-ReinforcementLearning', ...
    filesep, 'Preparation', filesep, 'mat', filesep, fn];



if recompute
    %% load data
    %
    mpool  = nans(length(COHS), length(TIMES), length(DIRS_), 3, length(SIMNUM));
    sdpool = nans(length(COHS), length(TIMES), length(DIRS_), 3, length(SIMNUM));
    npool  = nans(length(COHS), length(TIMES), length(DIRS_), 3, length(SIMNUM));
    
    for i = 1:length(SIMNUM)
        %% load data
        %
        [hdir, ldir, cdir, tdir] = dirnames;
        fname = ['/getML_RL_' Monk '_' int2str(NSEN) '_' int2str(NTUNE) '_' num2str(ELN) '_v' num2str(SIMNUM(i)) '.mat'];
        load([tdir fname])
        fprintf([fname '\n']);

        TR       = [4000 40000 length(FIRA)];
        trialnum = floor(TR/1000)+1;
  
        
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


        rsen  = 0.4;       % max correlation between neurons with different sensitivity
        rslp  = -rsen/20;  % number of neighbouring neurons that it is correlated to
        R0_   = [];
        th    = th(:);
        for j = 1:NSEN
            %R0_ = [R0_, rsen*exp(-abs(MTindex-MTindex(i))/rtau)];
            %R0_ = [R0_, rsen*exp(-abs([1:NSEN]'-i)/rtau)];
            R0_ = [R0_, rsen+rslp*abs([1:NSEN]'-j)];
        end
        R0_(R0_<0) = 0;



        R_    = [];
        for j = 1:NTUNE
            R_ = [R_, yy(j)*R0_];
        end
        for j = 1:NSEN      % force diagonal to 1
            R_(i,j) = 1;
        end

        R  = [];
        for j = 1:NTUNE
            R = [R; circshift(R_,[0,(j-1)*NSEN])];
        end
         
        % NOTE: the matlab CHOL function decompose R into U*U', where U is a lower triangular matrix.
        %       We thus can use U to "correlate" the normal random deviates
        %       responses of MT neurons (this is analogous to computing the matrix
        %       square root of R and use the square root to correlate the
        %       responses)
        U      = chol(R, 'lower');
        clear R xx yy;



        
        
        %% get mean responses as a function of coh, time and dir
        %
        % compute a, the mean respose of neuron given trial condition
        for j = 1:length(TR)
            for dd = 1:length(DIRS)
                atune = WrappedGaussian(kdir, 1, 40.^2, DIRS_(dd));
                atune = atune/max(max(atune));          % tuning modulation function
                acoh  = (kps-kns).*atune+kns;           % acoh will be a gaussian function, with peak at pref-slope, and baseline equal to null-slope (which can be a negative number)
                for cc = 1:length(COHS)
                    for tt = 1:length(TIMES)
                        [j dd cc tt]
                
                        a                   = TIMES(tt)*(k0 + COHS(cc)*acoh);                  % mean firing of MT cell for this trial
                        w                   = W(:,:,trialnum(j));

                        xx = [];
                        for nn = 1:100
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
                            pool   = pool+normrnd(0,2*sqrt(abs(pool))+5);

                            % decision
                            if pool*sign(cos(pi/180*DIRS_(dd)))>0
                                xx = [xx; pool];
                            end
                        end
                        %mpool(cc,tt,dd,j,i)  = w(:)'*a(:);
                        mpool(cc,tt,dd,j,i)  = nanmean(xx);
                        sdpool(cc,tt,dd,j,i) = nanstd(xx);
                        npool(cc,tt,dd,j,i)  = length(xx);
                        
                        
                    end
                end
            end
        end
    end
    save(savepath, 'mpool', 'sdpool', 'npool')

else
    load(savepath)

end


