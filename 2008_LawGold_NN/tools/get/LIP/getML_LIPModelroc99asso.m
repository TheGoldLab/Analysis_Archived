function [a, ae, rs, rsd, rn] = getML_LIPModelroc99asso(Monk, recompute)
% fit LIP roc areas at 99.9% coh to a piecewise-linear model
% (early training sessions with only 99% coh)

savepath = ['/work/mirror_jeff/code/matlab/tmp-mat/getML_LIPModelroc99asso_' Monk '.mat'];

if recompute
    a      = getML_txt([Monk 'TRain_asso.txt']);
    fn     = a.data{strcmp(a.name,'neuro_fn')};
    uid    = a.data{strcmp(a.name,'uid')};
    d1     = a.data{strcmp(a.name,'trg1_dir')};

    
    bb     = -200;
    be     = 1500;
    bs     = 50;
    bw     = 100;
    bins   = [[bb:bs:be-bw]' [bb+bw:bs:be]'];
    mbins  = mean(bins,2)/1000;

    
    a   = nans(3,length(fn));
    ae  = nans(3,2,length(fn));
    rs  = nans(length(mbins), length(fn));
    rsd = nans(length(mbins), length(fn));
    rn  = nans(length(mbins), length(fn));

    global FIRA
    for i = 1:length(fn)
        if ~isempty(fn{i})
            %%%     get ROC areas   %%%
            % define parameters
            crtf = [1];
            nvf  = [0 1];
            dirf = 0;

            % use dot dir to sort pref/null LIP responses because early
            % in training, target directions of pref and null targets are
            % not 180 degrees away, but dot dir is.
            utxtcol = 'ddir';
            firacol = 'dot_dir';

      
            % get rocs
            openFIRA(fn{i})
            fprintf('%.d: %s\n',i, fn{i})

            % get roc for each trials, select trials by based on
            % correct flag (crtf) and no var flag (nvf)
            Lcrt   = ismember(getFIRA_ecodesByName('correct', 'id'), crtf);
            trials = find(~isnan(FIRA.ecodes.data(:,1))&Lcrt);

            % get selection arrays for 99% coh
            Lcoh   = getFIRA_ecodesByName('dot_coh')==99.9;
            [Ldir, Udir] = selectFIRA_trialsByUniqueID(firacol);
            Idir         = [find(round(Udir)==d1(i)) find(round(Udir)==mod(d1(i)+180,360))];
            Ldir         = Ldir(:,Idir);

            begin_times = getFIRA_ecodeTimesByName('dot_on',0);
            end_times   = getFIRA_ecodeTimesByName('dot_off',0);

            [rs(:,i), rn(:,i), rsd(:,i)] = getFIRA_neurometricROCT(trials, Ldir(trials, :), Lcoh(trials, :),...
                getFIRA_spikeByID(uid(i)), ...
                begin_times(trials),...
                end_times(trials),...
                bins, 2);


            % get fit parameters
            if strcmp(Monk, 'ZZ')
                if i==1         L = mbins>=0   & mbins<=0.7; bl = 0.5;
                elseif i==2     L = mbins>=0   & mbins<=0.5; bl = 0.65;
                elseif i==3     L = mbins>=0   & mbins<=1; bl = 0.6;
                elseif i==18    L = mbins>=0.5 & mbins<=1; bl = 0.5;
                elseif i==20    L = mbins>=0   & mbins<=1; bl = 0.65;
                elseif i==24    L = mbins>=0   & mbins<=1; bl = 0.6;
                elseif i==25    L = mbins>=0   & mbins<=0.5; bl = 0.55;
                else
                    L = mbins>=0 & mbins<=1; bl=0.5;
                end
            else
                L = mbins>=0 & mbins<=1; bl=0.5;
            end
            m   = rs(:,i);                  % select between normalized and unnormalized spikes
            sd  = rsd(:,i)./(rn(:,i).^0.5); % standard error

            bcon = [min(mbins(L)) 0 0; max(mbins(L)) 6 1];

            mbP   = mbins(L);
            mbP   = mbP(:);
            c     = repmat(99.9,1,sum(L));
            c     = c(:);
            r     = m;
            rd    = sd;
            r     = r(L);
            rd    = rd(L);
            rd(rd==0) = nanmean(rd(rd~=0));

            if ~all(isnan(r))
                % fit data to piecewise-linear model to get plateau
                [a(:,i) xx gof ym] = sramp_fitW2(mbP, r, rd, bl, [], bcon);
                a(3,i) = ym(end);
                
                % estimate se
                % Check for input etype, which indicates that we will use
                %    Monte Carlo resampling (parametric bootstrap).
                %   See Wichmann & Hill (2001)
                %       The psychometric function: II. Bootstrap-
                %           based confidence intervals and sampling
                % citype{1} is number of simulated data sets to use
                mcn = 100;
                CI  = 68;
                mcfits = zeros(size(a,1),mcn);
                for ii = 1:mcn
                    if ~mod(ii,10)
                        disp(sprintf('Bootstrap CIs, set %d', ii))
                    end
                    % compute fit on simulated data set
                    ysim = normrnd(ym,rd);
                    [mcfits(:,ii) xx gof mcm] = sramp_fitW2(mbP, ysim, rd, bl, [], bcon);
                    mcfits(3,ii) = mcm(end);
                end
                
                ae(:,:,i) = [myprctile(mcfits',50-CI/2)' myprctile(mcfits',50+CI/2)'];
                a(:,i)
                ae(:,:,i)
                % plot and check fit
                if 0
                    clf
                    hold on
                    plot(mbins, m, 'o', 'MarkerSize', 6, ...
                        'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none')
                    plot(mbP, ym, '-', 'LineWidth', 3, 'Color', 'k')
                    xlim([0 1.5]);
                    ylim([0.4 1]);
                    title(sprintf('%.2f (%.2f), %.2f (%.2f)', a(2,i), ae(2,i), a(3,i), ae(3,i)))
                    pause
                end
            end

        end
    end


    save(savepath, 'a', 'ae', 'rs', 'rsd', 'rn')
else

    load(savepath)
end

