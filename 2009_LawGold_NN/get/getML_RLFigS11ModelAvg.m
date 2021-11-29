function [p, n, fits, sems, mbins, cbinsm] = getML_RLFigS11ModelAvg(Monk, NSEN, NTUNE, ELN, DIRS, SIMNUM, DFLAG, recompute)
% Get average data for model simulations and fits for plot

%%
[h] = dirnames;
fn  = ['getML_RLFigS11ModelAvg_' int2str(DFLAG), '_' 'Monk' '.mat'];
savepath = [h, filesep, 'code', filesep, 'matlab', filesep, 'current', ...
    filesep, '01MTLIP', filesep, 'reports', filesep, '3-ReinforcementLearning', ...
    filesep, 'Preparation', filesep, 'mat', filesep, fn];


if recompute
    warning off
    
    % load data
    coh = [];
    dir = [];
    vt  = [];
    crt = [];
    x   = [];
    for i = 1:length(SIMNUM)
        % get time constant for lapse
        [hdir, ldir, cdir, tdir] = dirnames;
        fname = ['/getML_RLEasyDiff_' Monk '_' int2str(DFLAG) '_' int2str(NSEN) '_' int2str(NTUNE) '_' num2str(ELN) '_v' num2str(SIMNUM(i)) '.mat'];
        load([tdir fname])
        fprintf([fname '\n']);

        N      = max(FIRA(:,1));
        coh    = [coh, FIRA(:,2)];
        dir    = [dir, FIRA(:,3)];
        vt     = [vt, FIRA(:,4)];
        crt    = [crt, FIRA(:,9)];
        x      = [x, [1:N]'];
    end

    

    % get thresholds
    bb     = 0;
    be     = N;
    bw     = 1000;
    bs     = bw/2;
    bins   = [[bb:bs:be-bw]' [bb+bw:bs:be]'];
    
    bb     = 0;
    be     = 1;
    bw     = 0.1;
    bs     = bw;
    cbins  = [[bb:bs:be-bw]' [bb+bw:bs:be]'];
    cbinsm = nanmean(cbins,2);
    
    fits   = nans(3,length(bins));
    sems   = nans(3,length(bins));
    
    p      = nans(length(cbins),length(bins));
    n      = nans(length(cbins),length(bins));
    for k = 1:length(bins)
        % select trials
        LTR = x>bins(k,1) & x<=bins(k,2);

        % compute thresholds using time-dependent weibull function (don't
        % fix lapse)
        coh_ = coh(LTR);
        vt_  = vt(LTR);
        dir_ = dir(LTR);
        crt_ = crt(LTR);
        d = [coh_(:) vt_(:) sign(dir_(:)) crt_(:) ones(size(crt_(:)))];
        
        for ii = 1:length(cbins)
            L = coh_>=cbins(ii,1) & coh_<cbins(ii,2);
            p(ii,k) = sum(crt_(L))/sum(L);  
            n(ii,k) = sum(L);  
        end
         
        Lgd = ~isnan(p(:,k)) & n(:,k)>20;
        [fits(:,k) sems(:,k)] = quick_fit([cbinsm(Lgd) p(Lgd,k) sqrt(n(Lgd,k))]);
        
 
%         citype  = [];
%         [fits(:,k), sems(:,k)] = ctPsych_fit(@quick2, cbinsm(Lgd), p(Lgd,k), [], citype, []);
       
    end
    mbins = nanmean(bins,2);

     
    save(savepath, 'p', 'n', 'fits', 'sems', 'mbins', 'cbinsm')

    warning on
    
else
    load(savepath)

end

