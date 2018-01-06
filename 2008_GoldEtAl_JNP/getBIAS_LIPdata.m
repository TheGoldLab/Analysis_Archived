function [wrd_,pd_,exdat_,wk_]   = getBIAS_LIPdata(sid, slopes, exd)
% function [wrd_,pd_,exdat_,wk_] = getBIAS_LIPdata(sid, slopes, exd)

%global FIRA

% make output arrays
wrd_ = nans(1, 2, 7);
pd_  = nans(1, 3, 7, 6); % b/sem/n 7 selectors; trg/pre/peri/post/sac; 3 Rs

% get tin
tin = getBIAS_tinFromFIRA(sid);

% get wk value
wk = getBIAS_wkFromFIRA(tin);

times = getFIRA_ecodesByName({'trg_on', 'dot_on', 'fp_off', 'sac_lat', 'dot_off'});
cohs  = getFIRA_ecodesByName('dot_coh')./100;
cors  = getFIRA_ecodesByName('correct');
stin  = sign(cos(tin*pi/180));
chcs  = stin.*sign(cos(getFIRA_ecodesByName('choice')*pi/180));
dirs  = stin.*sign(cos(getFIRA_ecodesByName('dot_dir')*pi/180));
Lgood = getFIRA_ecodesByName('task')==6 & cors>=0;
Lcor  = cohs==0 | cors==1;
Ltin  = cohs==0 | chcs==1;

% fit to PMF, with wk
ctdat = [cohs (times(:,5)-times(:,2))./1000 dirs chcs cors];
fits  = ctPsych_fit(@ddExp5fz_F, ctdat(Lgood, 1:4), ctdat(Lgood, 5), ...
    [], [], [], [], wk(Lgood));
wk_   = wk.*fits(4);

% get data
rdat                 = nans(size(Lgood,1), 6);
rdat(:,4)            = slopes;
[r,zs,rdat(Lgood,1)] = getFIRA_rateByBin(Lgood, sid, times(Lgood,1)-300, times(Lgood,1));
[r,zs,rdat(Lgood,2)] = getFIRA_rateByBin(Lgood, sid, times(Lgood,1), times(Lgood,1)+300);
[r,zs,rdat(Lgood,3)] = getFIRA_rateByBin(Lgood, sid, times(Lgood,2)-300, times(Lgood,2));
[r,zs,rdat(Lgood,5)] = getFIRA_rateByBin(Lgood, sid, times(Lgood,3)-300, times(Lgood,3));
[r,zs,rdat(Lgood,6)] = getFIRA_rateByBin(Lgood, sid, times(Lgood,3)+times(Lgood,4)-200, ...
    times(Lgood,3)+times(Lgood,4));
Lr = isfinite(rdat);

% 7 data sets: all, tin, tout, tin/cor, tin/err, tout/cor, tout/err
Lds = [Lgood Lgood&Ltin Lgood&~Ltin Lgood&Lcor&Ltin ...
    Lgood&~Lcor&Ltin Lgood&Lcor&~Ltin Lgood&~Lcor&~Ltin];
for ff = 1:size(Lds,2) % by selector

    % save mean,sem abs wk
    wrd_(1,:,ff) = [nanmean(abs(wk_(Lds(:,ff)))) nanmedian(abs(wk_(Lds(:,ff))))];
    
    for pp = 1:size(rdat,2) % by epoch
        Ly = Lds(:,ff) & Lr(:,pp);
        if sum(Ly) > 10
            if pp<4
                [R,P] = corr(rdat(Ly,pp), wk(Ly), 'type', 'Spearman');
            else
                [R,P] = partialcorr([rdat(Ly,pp) wk(Ly)], cohs(Ly), 'type', 'Spearman');
                R = R(2,1);
                P = P(2,1);
            end
            V     = (1+R.^2/2)*(1-R.^2).^2./(sum(Ly)-3);
            pd_(1,:,ff,pp) = [R, P, V];
        end
    end
end

% save example data
if nargin > 2 && exd
    
    % get regression inputs
    regs = {};
    for pp = 1:size(rdat,2)
        Ly = Lds(:,2) & Lr(:,pp);
        regs = cat(2,regs,[wk(Ly) rdat(Ly,pp)]);
    end

    % get PSTH, correct trials only
    ucohs = nonanunique(cohs);
    bsz = [800 1500 1000];
    tms = [1 -400 1 500; 2 -400 5 0; 3 -600 3 400];
    bs  = cell(length(bsz),1);
    es  = cell(length(bsz),1);
    for bb = 1:length(bsz)
        bs{bb} = 0:10:bsz(bb);
        bs{bb} = [bs{bb}' bs{bb}'+100];
        es{bb} = nans(size(bs{bb},1), length(ucohs), 2);
    end

    % correct only, tin/tout, all cohs
    Lps   = [Lgood&Lcor&Ltin Lgood&Lcor&~Ltin];
    for pp = 1:2
        for cc = 1:length(ucohs)
            Lall = Lps(:,pp) & cohs==ucohs(cc);
            if any(Lall)
                for bb = 1:length(bsz)
                    [rb,zbs,rbs] = getFIRA_rateByBin(Lall, sid, ...
                        times(Lall,tms(bb,1))+tms(bb,2), ...
                        times(Lall,tms(bb,3))+tms(bb,4), ...
                        times(Lall,tms(bb,1))+tms(bb,2), bs{bb});
                    es{bb}(:,cc,pp) = nanmean(rbs);
                end
            end
        end
    end
    exdat_ = {regs bs es};
else
    exdat_ = [];
end

%                 rmsz = 5;
%                 for gg=1:5
%                     subplot(5,1,gg);cla reset;hold on;
%                     plot([0 sum(Lgood)], [0 0], 'k:');
%                     Lg = Lgood & isfinite(zdat(:,gg));
%                     plot(zdat(Lgood,gg), 'k-');
%                     plot(nanrunmean(zdat(Lgood,gg),rmsz), 'k-', 'LineWidth', 2);
%                     plot(wk(Lgood).*10, 'r-');
%                     plot(nanrunmean(wk(Lgood),rmsz)*10, 'r-', 'LineWidth', 2);
%                     if sum(Lg) > 1
%                         [r,p] = corrcoef(zdat(Lg,gg), wk(Lg));
%                         title(sprintf('%.3f, p=%.3f, r=%.3f, a=%.3f', ...
%                             r(2,1), p(2,1), wrdat(ss,1), nanmean(abs(wk(Lgood)))))
%                     end
%                     ylim([-3 3]);
%                 end
%                 f=input('next')
