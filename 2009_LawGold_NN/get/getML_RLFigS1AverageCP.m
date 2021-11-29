function [cpm cpsd] = getML_RLFigS1AverageCP(Monk, NSEN, NTUNE, N, W, DNOISE, RECOMPUTE)
% Simulate choice probability for random initial weights, at different
% decision noise levels

[h] = dirnames;
fn  = ['/getML_RLFigS1AverageCP_' Monk '.mat'];
savepath = [h, filesep, 'code', filesep, 'matlab', filesep, 'current', ...
    filesep, '01MTLIP', filesep, 'reports', filesep, '3-ReinforcementLearning', ...
    filesep, 'Preparation', filesep, 'mat', filesep, fn];


if RECOMPUTE
    cpm  = nans(1,length(DNOISE));
    cpsd = nans(1,length(DNOISE));


    for nn = 1:length(DNOISE)
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
        %   Define constants
        %%%
        lambda  = 1;                                % 'penalty' factor
        choice  = nans(N,1);
        rALL    = nans(NSEN, NTUNE, N); 
 
%         W       = randn(size(k0));   % random initial condition for weight
%         wss     = sqrt(sum(sum(W.^2)));
%         W       = sqrt(0.02)*W./wss;


        %%%
        %   begin simulation
        %%%
        for i = 1:N
            if ~mod(i,100)
                fprintf([int2str(DNOISE(nn)) int2str(i) '\n'])
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

            % pool signals
            pool(i)   = sum(sum(W.*r));
            pool(i)   = pool(i)+normrnd(0,DNOISE(nn));

            % decision
            choice(i) = pool(i)>0;

            % store response
            rALL(:,:,i) = r;
        end

        nanmean(pool(choice==1))
        nanmean(pool(choice==0))

        % compute choice probability for each neuron
        cp = nans(NSEN,NTUNE);
        fprintf('Computing choice probability\n');

        LR = choice==1;
        LL = choice==0;
        for i = 1:NSEN
            for j = 1:NTUNE
                r       = shiftdim(rALL(i,j,:),2);
                cp(i,j) = rocN(r(LR), r(LL), 100);
            end
        end

        pdir     = -180:10:170;
        cpm(nn)  = nanmean(cp(:,find(pdir==90)));
        cpsd(nn) = nanstd(cp(:,find(pdir==90)));
    end

    
    % save data
    save(savepath, 'cpm', 'cpsd');

else
    load(savepath)
end




