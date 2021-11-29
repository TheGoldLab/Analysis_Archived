function [FIRA, W, MTindex] = getML_RLSimPerceptualLearning2(Monk, NSEN, NTUNE, eln, DIRS, MTDIRS, SIMNUM, RECOMPUTE)
% Simulation learning of the motion discrimination task with a
% reinforcement learning rule (A-reward-penalty rule) using response
% statistics of real MT neurons
% 
% get asymtotic thresholds

[hdir, ldir, cdir, tdir] = dirnames;
savepath = [tdir '/getML_RL_' Monk '_' int2str(NSEN) '_' int2str(NTUNE) '_' num2str(eln) '_v' num2str(SIMNUM) '.mat'];

if RECOMPUTE
    KK = 5;
    
    %%%
    %   Load task info (the actual coh/dir/viewing-time used during training)
    %%%
    load([Monk 'Combined.mat'])
    crt  = data(:,3);
    Lgd  = crt>=0 & ~isnan(data(:,6));
    crt  = crt(Lgd);
    dir  = data(Lgd,4);
    dir  = sign(cos(pi/180.*dir));
    dir(dir<0) = DIRS(DIRS<0);
    dir(dir>0) = DIRS(DIRS>0);
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
    if strcmp(Monk, 'Cy')
        Lgd      = ~isnan(PREF(1,:)) & ~isnan(NULL(1,:)) & ~isnan(FANO(1,:));
        Lgd(137) = logical(0);  % this cell is weird, with very small gain to both pref and null
    else
        Lgd  =~isnan(PREF(1,:)) & ~isnan(NULL(1,:)) & ~isnan(FANO(1,:));
    end    
    PREF = PREF(:,Lgd);
    NULL = NULL(:,Lgd);
    FANO = FANO(:,Lgd);

    
    % compute thresholds based on pref and null gain, baseline firing and fano factor 
    % DA  = ((PREF(1,:)+PREF(2,:))-(NULL(1,:)+NULL(2,:)))./sqrt(abs(FANO.*(PREF(1,:)+PREF(2,:))-FANO.*(NULL(1,:)+NULL(2,:))));
    DTH   = 0.9002; % D-prime value at which 81.6% correct is achieved
    rsum  = (PREF(2,:)+NULL(2,:)).^2;   % parameters to compute thresholds
    rdiff = (PREF(2,:)-NULL(2,:)).^2;
    bl    = nanmean([PREF(1,:); NULL(1,:)]);
    f     = FANO;
    TH    = (sqrt(rsum)+sqrt(rsum+8.*bl.*rdiff./f./DTH^2))./(2*rdiff./f./DTH^2);
    clear rsum rdiff bl f DTH
    
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

%     rsen   = 0.3;       % max correlation between neurons with different sensitivity
%     tausen = 0.15;      % time constant for rsen (drop to 0.32*rsen when the threshold is 20%coh larger or smaller)  
%     R0_    = [];
%     th     = th(:);
%     thmin  = min(th);
%     rthmin = rsen*exp(-abs(th-thmin)/tausen);   % correlation of most sensitive to others
%     R0_    = rthmin*rthmin';
     
    rsen  = 0.4;       % max correlation between neurons with different sensitivity
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
    %   Define constants
    %%%
    FIRA    = nans(KK*N,9);          % [trial #, coh, dir, vt, pooled response, reward prob, ...
                                                %   beta, choice, reward, MT responses (optional)]
    pool    = nans(KK*N,1);
    w       = 1e4*eln*(randn(size(k0)));   % random initial condition for weight
    W       = nans(size(w,1), size(w,2), floor(KK*N/1000)+2);   % store weights during training
    lambda  = 1;                                % 'penalty' factor



    %%%
    %   begin training
    %%%
    W(:,:,1) = w;   % store initial weights

    for i = 1:KK*N
        % store weights every 1000 trials
        if ~mod(i,1000)
            fprintf([int2str(i) '\n'])
            W(:,:,(i/1000)+1) = w;
        end


        %%%   [1] plan action:
        %           - pool signals
        %           - form direction decision
        %%%
        % compute a, the mean respose of neuron given trial condition
        if i<=N
            atune = normpdf(kdir, dir(i), 40);
            atune = atune/max(max(atune));          % tuning modulation function
            acoh  = (kps-kns).*atune+kns;           % acoh will be a gaussian function, with peak at pref-slope, and baseline equal to null-slope (which can be a negative number)
            a     = vt(i)*(k0 + coh(i)*acoh);       % mean firing of MT cell for this trial
            dir_  = dir(i);
            coh_  = coh(i);
            vt_   = vt(i);
            
        else
            dir_  = DIRS(unidrnd(length(DIRS)));
            COHS  = [0 3.2 6.4 12.8 25.6 51.2 99.9]/100;
            coh_  = COHS(unidrnd(length(COHS)));
            VTS   = [0:0.1:1];
            vt_   = VTS(unidrnd(length(VTS)));
            
            atune = normpdf(kdir, dir_, 40);
            atune = atune/max(max(atune));          % tuning modulation function
            acoh  = (kps-kns).*atune+kns;           % acoh will be a gaussian function, with peak at pref-slope, and baseline equal to null-slope (which can be a negative number)
            a     = vt_*(k0 + coh_*acoh);       % mean firing of MT cell for this trial
        end    
            
            
            
        
        % compute neural responses (normal deviates) for current trial
        r   = normrnd(zeros(size(a)), ones(size(a)));

        % multiple normal deviates with correlation matrix to generate
        % correlated responses, then scale it with the appropriate mean and
        % variance
        r      = reshape(U*r(:),NSEN,NTUNE);
        r      = a+sqrt(kf.*abs(a)).*r;
        r(r<0) = 0;

        % pool signals (try 0.3 noise, but cp is too high (about 0.55 early
        % in training, so try 1.5 noise)
        pool(i)   = sum(sum(w.*r));
        pool(i)   = pool(i)+normrnd(0,sqrt(0.3*abs(sum(sum(w.*a)))));
     
        
        % decision
        choice = pool(i)>0;
        rew    = choice==(sign(dir_)+1)/2;

    

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
            beta  = 1;  % prior of mean
            dbeta = 10; % prior of variability
        else
            gg    = 1/(1+exp(-beta*pool(i)));    % probability of reward for choosing right given pool signal
            yy    = (sign(dir_)+1)/2;          % actual reward for choosing right
            beta  = beta+(yy-gg)*(dbeta^-1)*pool(i);
        end
       
        
        

        %%%
        %   update weights based on reinforcement learning rule
        %%%
        p  = 1/(1+exp(-beta*pool(i)));
        dw = eln*rew*(choice-p)*r + lambda*eln*(1-rew)*(1-choice-p)*r;
        w  = w + dw;


  
        %  MT responses from cell MTind(1) to MTind(end) for neurons tuned to +ve dir]
        FIRA(i,:) = [i, coh_, dir_, vt_, pool(i), p, beta, choice, rew];
    end

    % store weights for the last trial
    W(:,:,end) = w;
    
    
    % save data
    save(savepath, 'FIRA', 'W', 'MTindex', 'NSEN', 'NTUNE', 'eln', 'DIRS', 'MTDIRS');

else
    load(savepath)
end




