function [bt, btd, bl, bld] = getML_RLFig2Tau(Monk, NSEN, NTUNE, ELN, DIRS, SIMNUM, BINS, SIMNAME, recompute, POW)
% Fit an exponential function to lapse and thresholds data


%%
if nargin<10
    POW = 1;
end

if ~strcmp(SIMNAME, 'NonLinPool')
    [h] = dirnames;
    fn  = ['getML_RLSuppDocTau_' SIMNAME '_' Monk '.mat'];
    savepath = [h, filesep, 'code', filesep, 'matlab', filesep, 'current', ...
        filesep, '01MTLIP', filesep, 'reports', filesep, '3-ReinforcementLearning', ...
        filesep, 'Preparation', filesep, 'mat', filesep, fn];
else
    [h] = dirnames;
    fn  = ['getML_RLSuppDocTau_' SIMNAME '_' num2str(POW) '_' Monk '.mat'];
    savepath = [h, filesep, 'code', filesep, 'matlab', filesep, 'current', ...
        filesep, '01MTLIP', filesep, 'reports', filesep, '3-ReinforcementLearning', ...
        filesep, 'Preparation', filesep, 'mat', filesep, fn];
end



if recompute
    warning off
    [fits, fitsd] = getML_RLSuppDocTh(Monk, NSEN, NTUNE, ELN, DIRS, SIMNUM, BINS, SIMNAME, 0, POW);
    [lp, lpd, lN] = getML_RLSuppDocLapse(Monk, NSEN, NTUNE, ELN, DIRS, SIMNUM, BINS, SIMNAME, 0, POW);

    th  = squeeze(fits(1,:,:));
    thd = squeeze(fitsd(1,:,:));

    bl  = nans(3,length(ELN));
    bt  = nans(3,length(ELN));
    bld = nans(3,2,length(ELN));
    btd = nans(3,2,length(ELN));

    % get time constant for lapse
    [hdir, ldir, cdir, tdir] = dirnames;
    if ~strcmp(SIMNAME, 'NonLinPool')
        fname = ['/getML_RL' SIMNAME '_' Monk '_' int2str(NSEN) '_' int2str(NTUNE) '_' num2str(ELN(1)) '_v' num2str(SIMNUM(1)) '.mat'];
    else
        fname = ['/getML_RL' SIMNAME '_' Monk '_' num2str(POW) '_' int2str(NSEN) '_' int2str(NTUNE) '_' num2str(ELN(1)) '_v' num2str(SIMNUM(1)) '.mat'];
    end
    load([tdir fname])
    fprintf([fname '\n']);
    N = max(FIRA(:,1));

    for i = 1:length(ELN)
        [ELN(i)]

%         % get time constant for lapse
%         bb     = 0;
%         be     = N;
%         bw     = 250;
%         bs     = bw;
%         bins   = [[bb:bs:be-bw]' [bb+bw:bs:be]'];
%         lbins  = nanmean(bins,2);
%         Lgd    = lN(:,i)>0;
%         lpdd   = lpd;
%         lpdd(lpdd==0) = nanmean(lpdd(lpdd~=0));
%         [bl(:,i), bld(:,:,i)] = exp_fitWd(lbins(Lgd), 1-lp(Lgd,i), lpdd(Lgd), [], [0 0; 0.5 0.5; 125 N], {100; 68});

        
        
        
%         % get time constant for lapse
        bb     = 0;
        be     = N;
        bw     = 250;
        bs     = bw;
        bins   = [[bb:bs:be-bw]' [bb+bw:bs:be]'];
        lbins  = nanmean(bins,2);
%         Lgd    = lN(:,i)>0;
%         lpdd   = lpd;
%         lpdd(lpdd==0) = nanmean(lpdd(lpdd~=0));
%         [bl(:,i), bld(:,:,i)] = exp_fitWd(lbins(Lgd), 1-lp(Lgd,i), lpdd(Lgd), [], [0 0; 0.5 0.5; 125 N], {100; 68});


        
        
        % get time constant for lapse
        [hdir, ldir, cdir, tdir] = dirnames;
        if ~strcmp(SIMNAME, 'NonLinPool')
            fname = ['/getML_RL' SIMNAME '_' Monk '_' int2str(NSEN) '_' int2str(NTUNE) '_' num2str(ELN(i)) '_v' num2str(SIMNUM(1)) '.mat'];
        else
            fname = ['/getML_RL' SIMNAME '_' Monk '_' num2str(POW) '_' int2str(NSEN) '_' int2str(NTUNE) '_' num2str(ELN(i)) '_v' num2str(SIMNUM(1)) '.mat'];
        end
        load([tdir fname])
        fprintf([fname '\n']);

        N      = max(FIRA(:,1));
        coh    = FIRA(:,2);
        dir    = FIRA(:,3);
        vt     = FIRA(:,4);
        crt    = FIRA(:,9);
        L      = coh>0.9 & vt>=0.9;
        x      = [1:N]';
    
        if strcmp(SIMNAME, 'TwoLIPPools') & ELN(i)==3.5e-8  % this simulation doesn't fit well, set initial condition
            [bl(:,i), bld(:,:,i)] = exp_fitd_bino(x(L), 1-crt(L), [0; 0.5; 1000], [0 0; 0.5 0.5; 0 N], {100; 68});
        elseif strcmp(SIMNAME, 'NonLinPool') & POW==1.41 & ELN(i)==0.6e-8  % this simulation doesn't fit well, set initial condition
            [bl(:,i), bld(:,:,i)] = exp_fitd_bino(x(L), 1-crt(L), [0; 0.5; 2000], [0 0; 0.5 0.5; 0 N], {100; 68});
        elseif strcmp(SIMNAME, 'NonLinPool') & POW==1 & ELN(i)==5e-8  % this simulation doesn't fit well, set initial condition
            [bl(:,i), bld(:,:,i)] = exp_fitd_bino(x(L), 1-crt(L), [0; 0.5; 3000], [0 0; 0.5 0.5; 0 N], {100; 68});
        else
            [bl(:,i), bld(:,:,i)] = exp_fitd_bino(x(L), 1-crt(L), [], [0 0; 0.5 0.5; 0 N], {100; 68});
        end    


        % get time constant for threshold
        bb     = 0;
        be     = N;
        bw     = 1000;
        bs     = bw;
        bins   = [[bb:bs:be-bw]' [bb+bw:bs:be]'];
        tbins  = nanmean(bins,2);
        Lbd    = th(:,i)<0.01;    % remove bad fits (usually only 1-2 sessions early in training)

        if strcmp(SIMNAME,'TwoLIPPools')
            bcon = [0.08 0.13; 0.5 1; 0 100000];
        elseif strcmp(SIMNAME, 'AddNoise')
            bcon = [0.08 0.13; 0.6 1; 0 100000];
        elseif strcmp(SIMNAME, 'MulNoise')
            bcon = [0.03 0.13; 0.4 1; 0 100000];
        elseif strcmp(SIMNAME, 'NoNorm') 
            bcon = [0.08 0.15; 0.5 1; 0 100000];
        elseif strcmp(SIMNAME, 'SubNorm') 
            bcon = [0.08 0.13; 0.5 0.8; 0 100000];
        elseif strcmp(SIMNAME, 'Oja') 
            bcon = [0.08 0.11; 0.4 1; 0 100000];
        elseif strcmp(SIMNAME, 'NonLinPool') & POW == 1.19
            bcon = [0.11 0.11; 0.45 0.45; 0 100000];
        elseif strcmp(SIMNAME, 'NonLinPool') & POW == 1.41
            bcon = [0.11 0.11; 0.4 1; 0 N];
        elseif strcmp(SIMNAME, 'NonLinPool') & POW == 2
            bcon = [0.10 0.15;  0.7 0.7; 0 N];
        elseif strcmp(SIMNAME, 'FitSig') 
            bcon = [0.08 0.13; 0.5 1; 0 100000];
        elseif strncmp(SIMNAME, 'GenRL0', 6)
            bcon = [0.08 0.13; 0.4 0.6; 0 N];
        elseif strncmp(SIMNAME, 'GenRL1', 6)
            bcon = [0.15 0.15; 0.6 1; 0 100000];
        else
            bcon = [0.08 0.13; 0.4 1; 0 100000];
        end    
        
        
        if strcmp(SIMNAME, 'NonLinPool') & POW==1.19 & ELN(i)==0.55e-8  % this simulation doesn't fit well, set initial condition
            [bt(:,i), btd(:,:,i)] = exp_fitWd(tbins(~Lbd), th(~Lbd,i), 0.05*ones(size(thd(~Lbd,i))), ...
                                                    [], [0.11 0.11; 0.45 1; 0 100000], {100; 68});
        elseif strcmp(SIMNAME, 'NonLinPool') & POW==1.19 & ELN(i)==0.5e-8  % this simulation doesn't fit well, set initial condition
            [bt(:,i), btd(:,:,i)] = exp_fitWd(tbins(~Lbd), th(~Lbd,i), 0.05*ones(size(thd(~Lbd,i))), ...
                                                    [], [0.11 0.11; 0.7 0.7; 0 100000], {100; 68});
        elseif strcmp(SIMNAME, 'NonLinPool') & POW==1.41 & ELN(i)==0.2e-8  % this simulation doesn't fit well, set initial condition
            [bt(:,i), btd(:,:,i)] = exp_fitWd(tbins(~Lbd), th(~Lbd,i), 0.05*ones(size(thd(~Lbd,i))), ...
                                                    [], [0.11 0.11; 0.8 0.8; 0 100000], {100; 68});
        elseif strcmp(SIMNAME, 'NonLinPool') & POW==2 & ELN(i)==2e-8  % this simulation doesn't fit well, set initial condition
            [bt(:,i), btd(:,:,i)] = exp_fitWd(tbins(~Lbd), th(~Lbd,i), 0.05*ones(size(thd(~Lbd,i))), ...
                                                    [], [0.03 0.05; 0.7 0.7; 0 100000], {100; 68});
        elseif strcmp(SIMNAME, 'GenRL11') & ELN(i)==9e-8  % this simulation doesn't fit well, set initial condition
            [bt(:,i), btd(:,:,i)] = exp_fitWd(tbins(~Lbd), th(~Lbd,i), 0.05*ones(size(thd(~Lbd,i))), ...
                                                    [], [0.11 0.11; 0.5 1; 0 100000], {100; 68});
        else
            [bt(:,i), btd(:,:,i)] = exp_fitWd(tbins(~Lbd), th(~Lbd,i), 0.05*ones(size(thd(~Lbd,i))), ...
                                                    [], bcon, {100; 68});
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

