function getBIAS_fitsForIlenna

% get list of monkeys
[monks,monkn,~] = getBIAS_monks;

pmfdat = cell(monkn, 5);
fcs    = getBIAS_fcs;
funs   = {@ddExp2_L; @ddExp3z_L; @ddExp41B_LF; @ddExp5RL_LF; @ddExp3fz_LF; @ddExp4fz_LF; @ddExp6RL_LF};
nfuns  = length(funs);
acsz   = 100;
    
    for mm = 1:monkn

        dat          = FS_getDotsTrainingData(monks{mm});
        Lgood        = dat(:,2) <= 2 & dat(:,3)>=0 & isfinite(fcs{mm});% & dat(:,6)<0.8;
        L51          = dat(:,5) == 0.5120;
        L99          = dat(:,5) >  0.9;
        Lcor         = dat(:,3) == 1;
        sessions     = unique(dat(:,1));
        num_sessions = length(sessions);
        chc          = dat(:,[9 9]);
        chc(~Lcor,1) = 0;
        chc( Lcor,2) = 0;

        pmfdat{mm,1} = nans(num_sessions, 3);           % corrcoeff, CIs
        pmfdat{mm,2} = nans(num_sessions, 7, 2, nfuns); % fits/sems; last arg is thresh
        pmfdat{mm,3} = nans(num_sessions, 3, nfuns);    % bs mean/sem, bs/dv
        pmfdat{mm,4} = nans(size(dat,1),  3, nfuns);    % preds, resids, biases
        pmfdat{mm,5} = nans(num_sessions, 3, 2, nfuns); % mean/median/sem; resid/chc/dir; DC/AC

        for ss = 1:num_sessions

            Lses  = Lgood & dat(:,1) == sessions(ss);
            disp([ss sessions(ss) sum(Lses)])

            % check 99 & 51 coh
            Ltm = Lses & dat(:,6) > nanmedian(dat(Lses,6));
            if sum(Ltm&L99&Lcor)/sum(Ltm&L99) < sum(Ltm&L51&Lcor)/sum(Ltm&L51)
                Lses = Lses & ~L99;
            end
    
            if sum(Lses) > 100

                % Get lapse rate
                lapse  = getSLOPE_lapse(dat(Lses, [3 5 6]));
                ltofit = min(0.2, max(lapse, 0.001));
                
                % fit to simple model
                [f1,s1,t1,p1,r1] = ctPsych_fit(funs{1}, dat(Lses, [5 6 8 9]), ...
                    dat(Lses, 3), [], [], [], [], ltofit);
                % ctPsych_plot(fun1, f1, dat(Lses, [5 6 8 9]), dat(Lses, 3))
                
                pmfdat{mm,2}(ss,1:7,1,1) = [f1' nans(1,3) lapse ctPsych_thresh(funs{1}, f1)];
                pmfdat{mm,2}(ss,1:2,2,1) = s1;
                pmfdat{mm,4}(Lses,1:2,1) = [p1, r1(:,4)];

                pmfdat{mm,5}(ss,:,1,1) = [nanmean(r1(:,4)) nanmean(r1(:,4).^2) nanse(r1(:,4))];
                ac = xcorr(r1(:,4), acsz, 'coeff');
                pmfdat{mm,5}(ss,:,2,1) = [ac(acsz+2) abs(mean(ac(acsz+(2:51)))) ...
                        abs(mean(ac(acsz+(2:101))))];

                % correlation coef between wk-filtered choices and resids
                [R, P, RLO, RUP]    = corrcoef(fcs{mm}(Lses), r1(:,4));
                pmfdat{mm,1}(ss, :) = [R(2,1) RLO(2,1) RUP(2,1)];
                
                % fit to funs, using A and lapse from ddExp3
                for ff = 2:nfuns
                    if ff == 2
                        va = [];
                    elseif ff==3 || ff==4 || ff==7
                        va = chc(Lses,:);
                    elseif ff==5 || ff==6
                        va = fcs{mm}(Lses);
                    end
                    [fits,sems,ts,ps,rs,d]  = ctPsych_fit(funs{ff}, dat(Lses,[5 6 8 9]), ...
                        dat(Lses, 3), [], [], f1(1), [], ltofit, va);
                    bs                      = feval(funs{ff}, fits, d, nan, va);
                    pmfdat{mm,2}(ss,:,1,ff) = [fits' nans(1,7-length(fits)-1) ctPsych_thresh(funs{ff}, fits)];
                    pmfdat{mm,2}(ss,:,2,ff) = [sems' nans(1,7-length(fits))];
                    nma                     = nanmean(abs(bs(:,1)));
                    pmfdat{mm,3}(ss,:,ff)   = [nma, nma./nanmean(abs(bs(:,2))), nanse(abs(bs(:,1)))];
                    pmfdat{mm,4}(Lses,:,ff) = [ps rs(:,4) bs(:,1)];
                    
                    pmfdat{mm,5}(ss,:,1,ff) = [nanmean(rs(:,4)) nanmean(rs(:,4).^2) nanse(rs(:,4))];
                    ac = xcorr(rs(:,4), acsz, 'coeff');
                    pmfdat{mm,5}(ss,:,2,ff) = [ac(acsz+2) abs(mean(ac(acsz+(2:51)))) ...
                        abs(mean(ac(acsz+(2:101))))];
                end                    
            end
        end
    end

    % save it to disk
    FS_saveProjectFile('2008_Bias', 'figBIAS_pmfCompare', pmfdat);
end
