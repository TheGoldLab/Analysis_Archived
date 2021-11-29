function [bt, btd, bl, bld] = getML_RLFig2Tau(Monk, NSEN, NTUNE, ELN, DIRS, SIMNUM, BINS, recompute)
% Fit an exponential function to lapse and thresholds data

%%
[h] = dirnames;
fn  = ['getML_RLFig2Tau_' Monk '.mat'];
savepath = [h, filesep, 'code', filesep, 'matlab', filesep, 'current', ...
    filesep, '01MTLIP', filesep, 'reports', filesep, '3-ReinforcementLearning', ...
    filesep, 'Preparation', filesep, 'mat', filesep, fn];


if recompute
    warning off
    [fits, fitsd] = getML_RLFig2Th(Monk, NSEN, NTUNE, ELN, DIRS, SIMNUM, BINS, 0);
    [lp, lpd, lN] = getML_RLFig2Lapse(Monk, NSEN, NTUNE, ELN, DIRS, SIMNUM, BINS, 0);

    th  = squeeze(fits(1,:,:));
    thd = squeeze(fitsd(1,:,:));

    bl  = nans(3,length(ELN));
    bt  = nans(3,length(ELN));
    bld = nans(3,2,length(ELN));
    btd = nans(3,2,length(ELN));

    % get time constant for lapse
    [hdir, ldir, cdir, tdir] = dirnames;
    fname = ['/getML_RL_' Monk '_' int2str(NSEN) '_' int2str(NTUNE) '_' num2str(ELN(1)) '_v' num2str(SIMNUM(1)) '.mat'];
    load([tdir fname])
    fprintf([fname '\n']);
    N = max(FIRA(:,1));

    for i = 1:length(ELN)
        [ELN(i)]

        % get time constant for lapse
        bb     = 0;
        be     = N;
        bw     = 250;
        bs     = bw;
        bins   = [[bb:bs:be-bw]' [bb+bw:bs:be]'];
        lbins  = nanmean(bins,2);
        Lgd    = lN(:,i)>0;
        lpdd   = lpd;
        lpdd(lpdd==0) = nanmean(lpdd(lpdd~=0));
        [bl(:,i), bld(:,:,i)] = exp_fitWd(lbins(Lgd), 1-lp(Lgd,i), lpdd(Lgd), [], [0 0; 0.5 0.5; 125 N], {100; 68});


        % get time constant for threshold
        bb     = 0;
        be     = N;
        bw     = 1000;
        bs     = bw;
        bins   = [[bb:bs:be-bw]' [bb+bw:bs:be]'];
        tbins  = nanmean(bins,2);
        Lbd    = th(:,i)<0.01;    % remove bad fits (usually only 1-2 sessions early in training)
        if strcmp(Monk, 'Cy')
            [bt(:,i), btd(:,:,i)] = exp_fitWd(tbins(~Lbd), th(~Lbd,i), 0.05*ones(size(thd(~Lbd,i))), ...
                [], [0.08 0.13; 0.5 1; 0 100000], {100; 68});
        else
            [bt(:,i), btd(:,:,i)] = exp_fitWd(tbins(~Lbd), th(~Lbd,i), 0.05*ones(size(thd(~Lbd,i))), ...
                [], [0.12 0.2; 0.7 1; 0 100000], {100; 68});
        end
        
        

        % plot
        if 0
            clf
            subplot(1,2,1)
            plot(lbins,1-lp(:,i), '.k')
            hold on
            yyy = bl(1,i)+(bl(2,i)-bl(1,i)).*exp(-1/bl(3,i).*lbins);
            plot(lbins, yyy, 'r')
            hold off
            ylim([0 0.5])
            xlim([1 N])
            title(num2str(bl(3,i)))


            % plot thresholds
            subplot(1,2,2)
            plot(tbins,th(:,i), '.k')
            hold on
            yyy = bt(1,i)+(bt(2,i)-bt(1,i)).*exp(-1/bt(3,i).*tbins);
            plot(tbins, yyy, 'r')
            hold off
            ylim([0.05 1])
            xlim([1 N])
            set(gca, 'yscale', 'log')
            title(num2str(bt(3,i)))
            pause
        end
    end

    

    

    warning on
    save(savepath, 'bt', 'bl', 'btd', 'bld')

else
    load(savepath)

end

