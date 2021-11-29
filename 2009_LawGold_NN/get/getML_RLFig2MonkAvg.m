function [fits, sems, tbins, lp, lpd, lN, lbins] = getML_RLFig2MonkAvg(Monk, recompute)
% Get behavioral thresholds and lapse for running trial bins

%%
home_dir = fullfile(filesep, 'Users', 'jigold', 'Desktop', 'FinalSubmission');
fn  = ['getML_RLFig2MonkAvg_' Monk '.mat'];
savepath = fullfile(home_dir, 'mat', fn);


if recompute
    warning off
      
    %% plot the same figure for monkeys
    COHS = [0 3.2 6.4 12.8 25.6 51.2 99.9]/100;
    load(fullfile(home_dir, 'mat', [Monk 'Combined.mat']));

    crt  = data(:,3);
    Lgd  = crt>=0;
    crt  = crt(Lgd);
    ddir = data(Lgd,4);
    coh  = data(Lgd,5);
    vt   = data(Lgd,6);
    N    = length(crt);

   
    % get running lapse
    bb     = 0;
    be     = N;
    bw     = 250;
    bs     = bw/2;
    lbins  = [[bb:bs:be-bw]' [bb+bw:bs:be]'];
    lp     = nans(length(lbins),1);
    lpd    = nans(length(lbins),1);
    lN     = nans(length(lbins),1);
    for k = 1:length(lbins)
        % select trials
        LTR  = [1:N]'>lbins(k,1) & [1:N]'<=lbins(k,2) & coh/100>0.9 & vt/1000>=1;
        if sum(LTR)>15
            coh_ = coh(LTR);
            vt_  = vt(LTR);
            crt_ = crt(LTR);
            lp(k)  = binofit(sum(crt_(:)),sum(LTR(:)));
            lN(k)  = sum(LTR(:));
            lpd(k) = lp(k)*(1-lp(k))/sqrt(lN(k)); 
        else
            lp(k)    = nan; 
            lpd(k,:) = nan;
            lN(k)    = sum(LTR(:));
        end
    end
    lbins = nanmean(lbins,2);
    
    
    
    % get running thresholds
    bb     = 0;
    be     = N;
    bw     = 1000;
    bs     = bw;
    bins   = [[bb:bs:be-bw]' [bb+bw:bs:be]'];
    fits   = nans(3, length(bins));
    sems   = nans(3, 2, length(bins));

    for i = 1:length(bins)
        % select for trials
        LTR = [1:N]'>bins(i,1) & [1:N]'<=bins(i,2);

        % use time-dependent weibull to get thresholds and lapse

        % prepare data
        % data(1) = coh  (0 .. 1)
        % data(2) = time (fractional seconds)
        % data(3) = dot dir: left (-1) / right (1)
        % data(4) = pct or correct (1) / error (0)
        % data(5) = (optinal) n
        d = [coh(LTR)./100 vt(LTR)./1000 sign(cos(pi/180.*ddir(LTR))) ...
            crt(LTR) ones(size(crt(LTR)))];

        % compute lapse
        Llp = coh>=90 & vt>=1000 & crt>=0;
        lapse  = sum(crt(LTR&Llp))/sum(LTR&Llp);
      
        if isnan(lapse)
            lapse = 1;
        elseif lapse<0.8
            lapse = 0.8;
        end

        citype  = {100, 68, 82};
        [fits(:,i), sems(:,:,i)] = ctPsych_fit(@quickTsFixLapse, d(:,1:2), d(:,4), [], citype, [], [1,0,0,0], 1-lapse);

   end

    fits  = real(fits);
    sems  = real(sems);
    tbins = nanmean(bins,2); 
    
    save(savepath, 'fits', 'sems', 'tbins', 'lp', 'lpd', 'lN', 'lbins')
    
    warning on

else
    load(savepath)

end

