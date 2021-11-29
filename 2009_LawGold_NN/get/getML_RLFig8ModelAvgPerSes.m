function [fits, sems, lp, lpd, lN] = getML_RLFig8ModelAvgPerSes(Monk, NSEN, NTUNE, ELN, DIRS, SIMNUM, recompute)
% Get average data for model simulations and fits for plot

%%
[h] = dirnames;
fn  = ['getML_RLFig8ModelAvgPerSes_' Monk '.mat'];
savepath = [h, filesep, 'code', filesep, 'matlab', filesep, 'current', ...
    filesep, '01MTLIP', filesep, 'reports', filesep, '3-ReinforcementLearning', ...
    filesep, 'Preparation', filesep, 'mat', filesep, fn];


if recompute
    warning off
    % get thresholds
    load([Monk 'Combined.mat'])
    Lgd  = data(:,3)>=0 & ~isnan(data(:,6));
    ses  = data(Lgd,1);
    clear data
    
    % load data
    coh = [];
    dir = [];
    vt  = [];
    crt = [];
    x   = [];
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
        x      = [x, ses];
    end

    

    
    
    fits   = nans(3,max(ses));
    sems   = nans(3,max(ses));

    for k = 1:max(ses)
        % select trials
        LTR = x==k;

        % compute thresholds using time-dependent weibull function (don't
        % fix lapse)
        coh_ = coh(LTR);
        vt_  = vt(LTR);
        dir_ = dir(LTR);
        crt_ = crt(LTR);
        d = [coh_(:) vt_(:) sign(dir_(:)) crt_(:) ones(size(crt_(:)))];
        
        % compute lapse
        Llp = coh>=0.9 & vt>=1 & crt>=0;
        lapse  = sum(sum(crt(LTR&Llp)))/sum(sum(LTR&Llp));
       
        if isnan(lapse)
            lapse = 1;
        elseif lapse<0.8
            lapse = 0.8;
        end

        citype  = [];
        [fits(:,k), sems(:,k)] = ctPsych_fit(@quickTsFixLapse, d(:,1:2), d(:,4), [], citype, [], [1,0,0,0], 1-lapse);
        
    end
  
    
          
    % get lapses
    lp     = nans(max(ses),1);
    lpd    = nans(max(ses),1);
    lN     = nans(max(ses),1);
    for k = 1:max(ses)
        % select trials
        LTR = x==k & coh>0.9 & vt>=1;

        coh_ = coh(LTR);
        vt_  = vt(LTR);
        crt_ = crt(LTR);

        % compute lapse
        lp(k)  = binofit(sum(crt_(:)),sum(LTR(:)));
        lN(k)  = sum(LTR(:));
        lpd(k) = lp(k)*(1-lp(k))/sqrt(lN(k)); 
    end

    
    save(savepath, 'fits', 'sems', 'lp', 'lpd', 'lN')

    warning on
    
else
    load(savepath)

end

