function [bl, bld, bt, btd] = getML_RLFig2ModelTau(Monk, NSEN, NTUNE, ELN, DIRS, SIMNUM, recompute)
% Get average data for model simulations and fits for plot

%%
[h] = dirnames;
fn  = ['getML_RLFig2ModelTau_' Monk '.mat'];
savepath = [h, filesep, 'code', filesep, 'matlab', filesep, 'current', ...
    filesep, '01MTLIP', filesep, 'reports', filesep, '3-ReinforcementLearning', ...
    filesep, 'Preparation', filesep, 'mat', filesep, fn];


if recompute
    warning off

    bl  = nans(3,length(ELN));
    bt  = nans(3,length(ELN));
    bld = nans(3,2,length(ELN));
    btd = nans(3,2,length(ELN));


    for i = 1:length(ELN)
    
            % load data
            coh = [];
            dir = [];
            vt  = [];
            crt = [];
            x   = [];
            for j = 1:length(SIMNUM)
                % get time constant for lapse
                [hdir, ldir, cdir, tdir] = dirnames;
                fname = ['/getML_RL_' Monk '_' int2str(NSEN) '_' int2str(NTUNE) '_' num2str(ELN(i)) '_v' num2str(SIMNUM(j)) '.mat'];
                load([tdir fname])
                fprintf([fname '\n']);

                N      = max(FIRA(:,1));
                coh    = [coh, FIRA(:,2)];
                dir    = [dir, FIRA(:,3)];
                vt     = [vt, FIRA(:,4)];
                crt    = [crt, FIRA(:,9)];
                x      = [x, [1:N]'];
            end



%             % get thresholds
%             bb     = 0;
%             be     = N;
%             bw     = 1000;
%             bs     = bw;
%             bins   = [[bb:bs:be-bw]' [bb+bw:bs:be]'];
%             fits   = nans(3,length(bins));
%             sems   = nans(3,length(bins));
% 
%             for k = 1:length(bins)
%                 % select trials
%                 LTR = x>bins(k,1) & x<=bins(k,2);
% 
%                 % compute thresholds using time-dependent weibull function (don't
%                 % fix lapse)
%                 coh_ = coh(LTR);
%                 vt_  = vt(LTR);
%                 dir_ = dir(LTR);
%                 crt_ = crt(LTR);
%                 d = [coh_(:) vt_(:) sign(dir_(:)) crt_(:) ones(size(crt_(:)))];
% 
%                 % compute lapse
%                 Llp = coh>=0.9 & vt>=1 & crt>=0;
%                 lapse  = sum(sum(crt(LTR&Llp)))/sum(sum(LTR&Llp));
% 
%                 if isnan(lapse)
%                     lapse = 1;
%                 elseif lapse<0.8
%                     lapse = 0.8;
%                 end
% 
%                 citype  = [];
%                 [fits(:,k), sems(:,k)] = ctPsych_fit(@quickTsFixLapse, d(:,1:2), d(:,4), [], citype, [], [1,0,0,0], 1-lapse);
% 
%             end
%             mbins = nanmean(bins,2);

            % get time constants for thresholds
            
            
            

            % get lapses
            bb     = 0;
            be     = N;
            bw     = 250;
            bs     = bw;
            lbins  = [[bb:bs:be-bw]' [bb+bw:bs:be]'];
            lp     = nans(length(lbins),1);
            lpd    = nans(length(lbins),1);
            lN     = nans(length(lbins),1);
            for k = 1:length(lbins)
                % select trials
                LTR = x>lbins(k,1) & x<=lbins(k,2) & coh>0.9 & vt>=1;

                coh_ = coh(LTR);
                vt_  = vt(LTR);
                crt_ = crt(LTR);

                % compute lapse
                lp(k)  = binofit(sum(crt_(:)),sum(LTR(:)));
                lN(k)  = sum(LTR(:));
                lpd(k) = lp(k)*(1-lp(k))/sqrt(lN(k));
            end
            lbins = nanmean(lbins,2);

            % get time constants for lapse
            Lgd = ~isnan(lp);
            [bl, bld] = exp_fitWd(lbins(Lgd), 1-lp(Lgd), lpd(Lgd), [], [0 0; 0.5 0.5; 10 100000], {100,68});

 
    end

    save(savepath, 'bl', 'bld', 'bt', 'btd')

    warning on

else
    load(savepath)

end

