function [bt, btd, bl, bld] = getML_RLFigS3Tau(Monk, NSEN, NTUNE, ELN, DIRS, SIMNUM, rmax, bdir, bsen, recompute)
% Fit an exponential function to lapse and thresholds data for simulations
% with different correlation structures

%%
[h] = dirnames;
fn  = ['getML_RLFigS3Tau_' Monk '_' num2str(length(rmax))  '_' num2str(length(bdir))  '_' num2str(length(bsen)) '.mat'];
savepath = [h, filesep, 'code', filesep, 'matlab', filesep, 'current', ...
    filesep, '01MTLIP', filesep, 'reports', filesep, '3-ReinforcementLearning', ...
    filesep, 'Preparation', filesep, 'mat', filesep, fn];


if recompute
    warning off
    
    bl  = nans(3,length(rmax),length(bdir),length(bsen));
    bt  = nans(3,length(rmax),length(bdir),length(bsen));
    bld = nans(3,2,length(rmax),length(bdir),length(bsen));
    btd = nans(3,2,length(rmax),length(bdir),length(bsen));
    
    
    for rr1 = 1:length(rmax)
        for rr2 = 1:length(bdir)
            for rr3 = 1:length(bsen)

                coh = [];
                dir = [];
                vt  = [];
                crt = [];
                x   = [];
                for i = 1:length(SIMNUM)
                    % get time constant for lapse
                    [hdir, ldir, cdir, tdir] = dirnames;
                    fname = ['/getML_RLCorr_' Monk '_' int2str(NSEN) '_' int2str(NTUNE) '_' num2str(ELN) '_' num2str(rmax(rr1)) '_' num2str(bdir(rr2)) '_' num2str(bsen(rr3)) '_v' num2str(SIMNUM) '.mat'];
                    load([tdir fname])
                    fprintf([fname '\n']);

                    N      = max(FIRA(:,1));
                    coh    = [coh, FIRA(:,2)];
                    dir    = [dir, FIRA(:,3)];
                    vt     = [vt, FIRA(:,4)];
                    crt    = [crt, FIRA(:,end)];
                    x      = [x, [1:N]'];
                end



 

               % get lapse time constants
                N      = max(FIRA(:,1));
                coh    = FIRA(:,2);
                dir    = FIRA(:,3);
                vt     = FIRA(:,4);
                crt    = FIRA(:,9);
                L      = coh>0.9 & vt>=1;
                x      = [1:N]';
                [bl(:,rr1,rr2,rr3), bld(:,:,rr1,rr2,rr3)] = exp_fitd_bino(x(L), 1-crt(L), [], [0 0; 0.5 0.5; 0 N], {100, 68});
            
                
                
                
                % get threshold time constants
                bb     = 0;
                be     = N;
                bw     = 1000;
                bs     = bw;
                bins   = [[bb:bs:be-bw]' [bb+bw:bs:be]'];
                tmb    = nanmean(bins,2);
                th     = nans(length(bins),1);
                thd    = nans(length(bins),1);
                tbins  = nanmean(bins,2);            
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
                    
                    % compute lapse
                    Llp = coh_>=0.9 & vt_>=1 & crt_>=0;
                    lapse  = sum(crt_(Llp))/sum(Llp);
                    if isnan(lapse)
                        lapse = 1;
                    elseif lapse<0.8
                        lapse = 0.8;
                    end
                    
                    
                    citype      = [];
                    [fits sems] = ctPsych_fit(@quickTsFixLapse, d(:,1:3), d(:,4), [], citype, [], [1,0,0,0], 1-lapse);
                    th(k)  = real(fits(1));
                    thd(k) = real(sems(1));
                end

                Lbd    = th<0.01 | isnan(th) | thd==0;    % remove bad fits (usually only 1-2 sessions early in training)
                [bt(:,rr1,rr2,rr3), btd(:,:,rr1,rr2,rr3)] = exp_fitWd(tbins(~Lbd), th(~Lbd), thd(~Lbd), ...
                                                                    [], [0.03 0.4; 0.4 1; bw N], {1000, 68});
                
                
            end
        end
    end

        
    warning on
    save(savepath, 'bt', 'bl', 'btd', 'bld')

else
    load(savepath)

end

