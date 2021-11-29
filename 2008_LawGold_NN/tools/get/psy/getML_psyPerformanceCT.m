function [fits, sems, th] = getML_psyPerformanceCT(fn, recompute)
% use  [fits(:,i),sems(:,i),stats] = quickTs_fit(fit_dat,afun,bfun,[],[],time_bins);
% to get behavioral threshold

if nargin < 1
    return;
elseif nargin < 2
    recompute = 1;
end

[home_dir, lab_dir, current_dir, tmat_dir] = dirnames;
savepath = [tmat_dir '/getML_psyPerformanceCT_' fn(1:2) '.mat'];


if recompute
    a      = getML_txt(fn);
    fname  = a.data{strcmp(a.name,'dat_fn')};
    
    fits   = nans(4,length(fname));
    sems   = nans(4, 2, length(fname));
    th     = nans(length(fname),1);     % threshold

    % get fits
    global FIRA
    for i = length(fname):-1:1
        fprintf('%d: %s\n', i, fname{i})
        openFIRA(fname{i})

        %% get data for fit
        % data(1) = coh  (0 .. 1)
        % data(2) = time (fractional seconds)
        % data(3) = dot dir: left (-1) / right (1)
        % data(4) = pct or correct (1) / error (0)
        % data(5) = (optinal) n
        d = [getFIRA_ecodesByName('dot_coh')./100 ...
            (getFIRA_ecodesByName('dot_off')-getFIRA_ecodesByName('dot_on'))./1000 ...
            sign(cos(pi/180.*getFIRA_ecodesByName('dot_dir'))) ...
            getFIRA_ecodesByName('correct') ...
            getFIRA_ecodesByName('choice')];

        % select trials
        Lgood = d(:,4)>=0 & d(:,2)<=1.5;
        d     = d(Lgood,:);

        % set initial conditions for these cells to improve the fit
        if strcmp(fn(1:2), 'Cy')
            if i == 19
                % this session only have 139 trials, set inital condition
                % to get better fit
                init = [0.4; -0.5; 1.9; 0.2];
            
                 
            elseif i == 29
                % this session has no lapse, but high coh depends on time.
                % However, the resulting fit from QuickTs has about 10%
                % lapse and no time dependence. Try to set initial value to
                % get a better fit to the function.
                init = [0.28; -5; 1.5; 0];
             
            elseif i == 51
                % this session doesn't fit well
                init = [0.418; -1; 1; 0];

            elseif i == 59
                % this session doesn't fit well
                init = [0.1; -5; 1; 0.1];

            elseif i == 62
                % this session doesn't fit well
                init = [0.2; 0; 1; 0];
             
            elseif i == 109
                % this session doesn't fit well
                init = [0.2; -5; 1; 0];
                
                % ignore 99% for this session
                LL      = d(:,1)>0.7;
                d(LL,:) = [];
            
            elseif i == 117
                % this is a good session with many trials, but quickTs
                % doesn't fit it unless I send in the initial condition,
                % which is estimated by fiting quick function to behavior,
                % ignoring time
                init = [0.1; 0; 1; 0];
            
            elseif i == 127
                % this session doesn't fit well
                init = [0.15; -1; 1; 0];

            elseif i == 132
                % this session doesn't have enough trials
                init = [0.15; -1; 1; 0];
                
            elseif i == 139
                init = [0.1; -1; 0.1; 0];
         
                
            elseif i == 137
                % this session doesn't have enough trials
                init = [0.10; -1; 1; 0];
                
            elseif i == 160
                % this session doesn't fit well
                init = [0.15; -1; 1; 0];

            else
                init = [];

            end
            
        elseif strcmp(fn(1:2), 'ZZ')
            if     i == 3
                init = [0.6; 0; 1; 0.2];
            
            elseif     i == 5
                init = [0.6; 0; 10; 0];
            
            elseif     i == 7
                init = [0.6; 0; 10; 0];
                LL   = d(:,2)>1.2;  % remove a few long vt trials, which doesn't have much data
                d(LL,:) = [];
            
            elseif     i == 8
                init = [0.6; 0; 10; 0];
            
            elseif     i == 10
                init = [0.6; 0; 10; 0];
            
            elseif     i == 13
                init = [0.6; -5; 10; 0];
            
            elseif     i == 14
                init = [0.6; 0; 20; 0];
            
            elseif     i == 17
                init = [0.6; 0; 10; 0];
            
            elseif     i == 21
                init = [0.6; -5; 20; 0];
            
            elseif     i == 23
                init = [0.3; -5; 10; 0];
            
            elseif i == 24
                init = [0.8; -5; 5; 0];
                LL   = d(:,2)>1.3;  % remove a few long vt trials, which doesn't have much data
                d(LL,:) = [];
            
            elseif     i == 27
                init = [0.3; -1; 5; 0];
                
            elseif i == 29
                init = [0.2; -1; 10; 0];
            
            elseif i == 30
                init = [0.8; -1; 10; 0];
            
            elseif i == 31
                init = [0.5; -5; 5; 0];
            
            elseif i == 32
                init = [0.8; -3; 10; 0];
                
            elseif i == 34
                init = [0.6; -5; 20; 0];
                
            elseif i==33 | (i>=35 & i<=37) 
                % early in training, zsa zsa's performance is not very
                % good, so forcing the parameter search to start from bad
                % performance (threshold around 80%) helps to improve the
                % fit.
                init = [0.8; -3; 5; 0];
            
            elseif i==39
                init = [0.5; -5; 20; 0];
            
            elseif i==41
                init = [0.6; -5; 10; 0];
            
            elseif i==43
                init = [0.6; -5; 10; 0];
            
            elseif i==48
                init = [0.5; 0; 10; 0];
            
            elseif i==50
                init = [0.5; 0; 10; 0];
            
            elseif i==51
                init = [0.3; -5; 10; 0];
            
                
            elseif i==52
                init = [0.3; 0; 10; 0];
            
            elseif i==58
                init = [0.5; -5; 10; 0];
            
            elseif (i >= 40 & i <= 42) | i == 44 ...
                    | (i >= 53 & i <= 54) | i == 59
                init = [0.8; -5; 5; 0];
            
            elseif i == 55
                % this session is bad, need to set initial condition
                % manually, threshold is still high but at least it makes
                % the se smaller.
                init = [0.5; -5; 10; 0];
                LL   = d(:,2)>1.2;  % remove a few long vt trials, which doesn't have much data
                d(LL,:) = [];
     
            elseif i == 60 
                init = [0.5; -5; 10; 0];
            
            elseif i == 62 
                init = [0.5; -5; 10; 0];
                
            elseif i == 63 
                init = [0.5; -5; 3; 0];
                
            elseif i == 64 | i == 66
                init = [0.2; -5; 3; 0];
                
            elseif i == 71
                init = [0.2; -5; 3; 0];
            
            elseif i == 75
                init = [0.3; -5; 3; 0];
                
            elseif i == 79
                init = [0.4; -5; 5; 0];
            
            elseif i == 82
                init = [0.3; -5; 10; 0];
            
            elseif i == 87
                init = [0.8; 0; 20; 0];
                
            elseif i == 90 
                init = [0.8; -1; 5; 0];
            
            elseif i == 94
                init = [0.4; -5; 2; 0.05];
            
            elseif i == 102
                init = [0.8; -1; 5; 0];
                
            elseif i == 103 | i == 111 | i == 112 ...
                    | i == 115 | i == 121 
                % these sessions fit much better with this initial condition 
                init = [0.4; -5; 5; 0];
           
            elseif i == 113
                init = [0.2; -1; 1; 0];
           
            elseif i == 117
                init = [0.2; -1; 1; 0];
           
            elseif i == 119
                init = [0.2; 0; 1; 0];
           
            elseif i == 120
                init = [0.8; -1; 1; 0];
                
            elseif i == 122
                init = [0.4; 0; 20; 0];
            
            elseif i == 124
                init = [0.4; 0; 20; 0];
                
            elseif i == 125 | i == 126
                % this session only have 187 trials, set inital condition
                % to get better fits, this initial condition works for
                % session 126 too
                init = [0.2; 0; 20; 0];
            
            else
                init = [];
            end
        end
        
        
        % get fits
        citype = {80, 68, 81.61};
        
        %[fits(:,i),sems(:,i)] = ctPsych_fit(@quickTs, d(:,1:3), d(:,4), citype, init, [1,0,1,0]);
        %[a,b] = ctPsych_fit(@quickTs, d(:,[1:2]), d(:,4), [], citype, init, [1,0,1,0]);
        
        [fits(:,i), sems(:,:,i)] = ctPsych_fit(@quickTs, d(:,1:2), d(:,4), [], citype, init, [1,0,1,0]);
        
        
        th(i) = fits(1,i);
        
        fits(:,i)
        sems(:,:,i)
        if 0
            warning off
            % plot and check fit
            bb   = 0;
            be   = 1500;
            bw   = 200;
            bs   = 50;
            bins = [[bb:bs:be-bw]' [bb+bw:bs:be]'];
            mbin = mean(bins,2)./1000;
            coh  = [0 3.2 6.4 12.8 25.6 51.2 99.9]';

            % get performance for each coh and time
            p = nans(length(coh),size(bins,1));
            N = nans(length(coh),size(bins,1));
            for j = 1:size(bins,1)
                Lt = FIRA.ecodes.data(:,getFIRA_ecodeColumnByName('dot_dur'))>=bins(j,1) ...
                    & FIRA.ecodes.data(:,getFIRA_ecodeColumnByName('dot_dur'))<bins(j,2);
                for k = 1:length(coh)
                    Lcoh = FIRA.ecodes.data(:,getFIRA_ecodeColumnByName('dot_coh'))==coh(k);
                    N(k,j) = nansum(FIRA.ecodes.data(Lcoh&Lt,getFIRA_ecodeColumnByName('correct'))>=0);
                    p(k,j) = nansum(FIRA.ecodes.data(Lcoh&Lt,getFIRA_ecodeColumnByName('correct'))==1)/N(k,j);
                end
            end
            
            % prepare plot parameters
            lc  = [];
            for j = 1:7
                lc = [lc; [0.9 0.9 0.9]-(j-1)/7];
            end

            clf
            hold on
            for j = 1:7
                plot(mbin, p(j,:)', 'o', 'MarkerSize', 8, 'MarkerFaceColor', lc(j,:), 'MarkerEdgeColor', 'none')
                ym = 0.5 + (1- 0.5 - fits(4,i)).*(1-exp(-((coh(j)/100)./(fits(1,i).*mbin.^fits(2,i))).^fits(3,i)));
                plot(mbin, ym, 'LineWidth', 2, 'Color', lc(j,:));
            end
            title(sprintf('%.2f, %.2f, %.2f, %.2f', 100*fits(1,i), fits(2,i), fits(3,i), 100*fits(4,i)))
            hold off
            pause
            warning on
        end       
    end

    fits = real(fits);
    sems = real(sems);
    th   = real(th);
    
    clear FIRA
    pack
        
        
    % save
    save(savepath, 'fits', 'sems', 'th')
else
    % load
    load(savepath)
end





