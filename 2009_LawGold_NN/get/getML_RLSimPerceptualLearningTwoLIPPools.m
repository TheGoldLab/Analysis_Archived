function [FIRA, WP, WN, MTindex] = getML_RLSimPerceptualLearningTwoLIPPools(Monk, NSEN, NTUNE, eln, DIRS, MTDIRS, SIMNUM, RECOMPUTE)
% Simulation learning of the motion discrimination task with a
% reinforcement learning rule (A-reward-penalty rule) using response
% statistics of real MT neurons
% 
% % two pools of LIP neurons
% % add hybird pooling noise (poisson+additive)
% % with weight normalization
% % not fixed directions

[hdir, ldir, cdir, tdir] = dirnames;
savepath = [tdir '/getML_RLTwoLIPPools_' Monk '_' int2str(NSEN) '_' int2str(NTUNE) '_' num2str(eln) '_v' num2str(SIMNUM) '.mat'];

if RECOMPUTE

    %%%
    %   Load task info (the actual coh/dir/viewing-time used during training)
    %%%
    load([Monk 'Combined.mat'])
    crt  = data(:,3);
    Lgd  = crt>=0 & ~isnan(data(:,6));
    crt  = crt(Lgd);
    dir  = data(Lgd,4);
    dirsgn = sign(cos(pi/180.*dir));
    coh  = data(Lgd,5)/100;
    vt   = data(Lgd,6)/1000;
    N    = length(vt);
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
    FIRA    = nans(N,9);                        % [trial #, coh, dir, vt, pooled response, reward prob, ...
                                                %   beta, choice, reward, MT responses (optional)]
    poolp   = nans(N,1);
    pooln   = nans(N,1);
    
    w       = randn(size(k0));   % random initial condition for weight
    wss     = sqrt(sum(sum(w.^2)));
    wp      = sqrt(0.01)*w./wss;
    
    w       = randn(size(k0));   % random initial condition for weight
    wss     = sqrt(sum(sum(w.^2)));
    wn      = sqrt(0.01)*w./wss;
    
    WN      = nans(size(w,1), size(w,2), floor(N/1000)+2);   % store weights during training
    lambda  = 1;                                             % 'penalty' factor



    %%%
    %   begin training
    %%%
    WP(:,:,1) = wp;   % store initial weights
    WN(:,:,1) = wn;   % store initial weights

    for i = 1:N
        % store weights every 1000 trials
        if ~mod(i,1000)
            fprintf([int2str(i) '\n'])
            WP(:,:,(i/1000)+1) = wp;
            WN(:,:,(i/1000)+1) = wn;
        end


        %%%   [1] plan action:
        %           - pool signals
        %           - form direction decision
        %%%
        % compute a, the mean respose of neuron given trial condition
        %atune = normpdf(kdir, dir(i), 40);
        atune = WrappedGaussian(kdir, 1, 40.^2, dir(i));
        atune = atune/max(max(atune));          % tuning modulation function
        acoh  = (kps-kns).*atune+kns;           % acoh will be a gaussian function, with peak at pref-slope, and baseline equal to null-slope (which can be a negative number)
        a     = vt(i)*(k0 + coh(i)*acoh);       % mean firing of MT cell for this trial
        
        
        % compute neural responses (normal deviates) for current trial
        r   = normrnd(zeros(size(a)), ones(size(a)));

        % multiple normal deviates with correlation matrix to generate
        % correlated responses, then scale it with the appropriate mean and
        % variance
        r      = reshape(U*r(:),NSEN,NTUNE);
        r      = a+sqrt(kf.*abs(a)).*r;
        r(r<0) = 0;

        % pool signals (add constant pooling noise)
        poolp(i)   = sum(sum(wp.*r));
        poolp(i)   = poolp(i)+normrnd(0,sqrt(2*abs(poolp(i)))+5);
        pooln(i)   = sum(sum(wn.*r));
        pooln(i)   = pooln(i)+normrnd(0,sqrt(2*abs(pooln(i)))+5);
       
        
        % decision
        choice = poolp(i)>pooln(i);
        rew    = choice==(dirsgn(i)+1)/2;

    

        %%%   [2] evaluate action:
        %           - update internal representation of reward probability
        %           - compute prediction error of reward
        %           - change weight according to Ar-p rule
        %%%
        
        %%%             old way, fit logistic to previous 300 trials             %%%
        %         % assume relationship between pooled signal and the probability of
        %         % reward for choosing rightward direction is a logistic, estimate beta
        %         % for the logistic function based on the previous 1000 choices
        %         NBHIST = 300;
        %         if i<NBHIST
        %             beta = 1;
        %         elseif ~mod(i,NBHIST)
        %             fprintf([int2str(i) '\n'])
        %             fits = logist_fit([pool(i-NBHIST+1:i) (sign(dir(i-NBHIST+1:i))+1)/2], 0, 0);
        %             beta = fits(3);
        %         else
        %             beta = FIRA(i-1,7);
        %         end
        %%%
        
        %%%             new way, use sequential update rule                         %%%
        %               laplace approximation (Spiegelhalter and Lauritzen 1990)      %    
        if i == 1
            betap  = 1;  % prior of mean
            dbeta = 10; % prior of variability
        else
            gg    = 1/(1+exp(-betap*poolp(i)));    % probability of reward for choosing right given pool signal
            yy    = (sign(dirsgn(i))+1)/2;          % actual reward for choosing right
            betap = betap+(yy-gg)*(dbeta^-1)*poolp(i);
        end
       
        if i == 1
            betan  = 1;  % prior of mean
            dbeta = 10; % prior of variability
        else
            gg    = 1/(1+exp(-betan*pooln(i)));    % probability of reward for choosing left given pool signal
            yy    = (-sign(dirsgn(i))+1)/2;          % actual reward for choosing left
            betan  = betan+(yy-gg)*(dbeta^-1)*pooln(i);
        end
       
        

        %%%
        %   update weights according to reinforcement learning rule
        %%%
%         p   = 1/(1+exp(-betap*poolp(i)));
%         dw  = eln*rew*(choice-p)*r + lambda*eln*(1-rew)*(1-choice-p)*r;
%         
%         wp  = wp + dw;
%         wss = sqrt(sum(sum(wp.^2)));
%         wp   = sqrt(0.01)*wp./wss;
%         
%         
%         
%         p   = 1/(1+exp(-betan*pooln(i)));
%         dw  = eln*rew*(~choice-p)*r + lambda*eln*(1-rew)*(1-~choice-p)*r;
%         
%         wn  = wn + dw;
%         wss = sqrt(sum(sum(wn.^2)));
%         wn  = sqrt(0.01)*wn./wss;

        if choice
            if rew
                p   = 1/(1+exp(-betap*poolp(i)));
                dw  = eln*rew*(choice-p)*r;

                wp  = wp + dw;
                wss = sqrt(sum(sum(wp.^2)));
                wp   = sqrt(0.02)*wp./wss;
            else
                p   = 1/(1+exp(-betap*poolp(i)));
                dw  = lambda*eln*(1-rew)*(1-choice-p)*r;

                wp  = wp + dw;
                wss = sqrt(sum(sum(wp.^2)));
                wp   = sqrt(0.02)*wp./wss;

            end

        else
            if rew
                p   = 1/(1+exp(-betan*pooln(i)));
                dw  = eln*rew*(~choice-p)*r;

                wn  = wn + dw;
                wss = sqrt(sum(sum(wn.^2)));
                wn  = sqrt(0.02)*wn./wss;
            else
                p   = 1/(1+exp(-betan*pooln(i)));
                dw  = lambda*eln*(1-rew)*(1-~choice-p)*r;

                wn  = wn + dw;
                wss = sqrt(sum(sum(wn.^2)));
                wn  = sqrt(0.02)*wn./wss;

            end

        end

        FIRA(i,:) = [i, coh(i), dir(i), vt(i), poolp(i)-pooln(i), p, betap, choice, rew];
    end

    % store weights for the last trial
    WP(:,:,end) = wp;
    WN(:,:,end) = wn;
    
    
    % save data
    save(savepath, 'FIRA', 'WP', 'WN', 'MTindex', 'NSEN', 'NTUNE', 'eln', 'DIRS', 'MTDIRS');

else
    load(savepath)
end




