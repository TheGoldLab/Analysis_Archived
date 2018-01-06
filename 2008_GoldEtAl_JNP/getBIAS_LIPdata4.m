function [bs_,rocs_,ds_]   = getBIAS_LIPdata4(sid, slopes)
% function [bs_,rocs_,ds_] = getBIAS_LIPdata4(sid, slopes)

global FIRA

% get tin, wk
tin   = getBIAS_tinFromFIRA(sid);
wk    = getBIAS_wkFromFIRA(tin);

times = getFIRA_ecodesByName({'trg_on', 'dot_on', 'dot_off', 'fp_off', 'sac_lat'});
cohs  = getFIRA_ecodesByName('dot_coh');
cors  = getFIRA_ecodesByName('correct');
stin  = sign(cos(tin*pi/180));
chcs  = stin.*sign(cos(getFIRA_ecodesByName('choice')*pi/180));
dirs  = stin.*sign(cos(getFIRA_ecodesByName('dot_dir')*pi/180));
Lgood = getFIRA_ecodesByName('task')==6 & cors>=0;
Lcor  = cohs==0 | cors==1;
Ltin  = chcs==1;

% fit to PMF, with wk
ctdat = [cohs./100 (times(:,3)-times(:,2))./1000 dirs chcs cors];
fits  = ctPsych_fit(@ddExp5fz_F, ctdat(Lgood, 1:4), ctdat(Lgood, 5), ...
    [], [], [], [], wk(Lgood));
wk    = wk.*fits(4) + fits(5);

% save mean abs wk
bs_ = [nanmean(abs(wk(Lgood))) nanse(abs(wk(Lgood))) nanmean(wk(Lgood))];

% get data
rs = nans(size(Lgood,1), 4);
[r,zs,rs(Lgood,1)] = getFIRA_rateByBin(Lgood, sid, times(Lgood,1)-300, times(Lgood,1));
[r,zs,rs(Lgood,2)] = getFIRA_rateByBin(Lgood, sid, times(Lgood,1),     times(Lgood,2));
[r,zs,rs(Lgood,3)] = getFIRA_rateByBin(Lgood, sid, times(Lgood,2),     times(Lgood,3));
[r,zs,rs(Lgood,4)] = getFIRA_rateByBin(Lgood, sid, times(Lgood,4)-300, times(Lgood,4));
[r,zs,rs(Lgood,5)] = getFIRA_rateByBin(Lgood, sid, times(Lgood,4)+times(Lgood,5)-200, ...
    times(Lgood,4)+times(Lgood,5));

% pres -- roc tin vs tout, sep by bias
rocs_ = nans(1,8,size(rs,2));
for rr = 1:size(rs,2)
    rocs_(1,:,rr) = [ ...
        rocN(rs(Lgood&Lcor&Ltin,rr),      rs(Lgood&Lcor&~Ltin,rr)), ...
        rocN(rs(Lgood&Lcor&Ltin&wk>0,rr), rs(Lgood&Lcor&~Ltin&wk<0,rr)), ...
        rocN(rs(Lgood&Lcor&Ltin&wk<0,rr), rs(Lgood&Lcor&~Ltin&wk>0,rr)), ...
        rocN(rs(Lgood&Ltin,rr),           rs(Lgood&~Ltin,rr)), ...
        rocN(rs(Lgood&Ltin&wk>0,rr),      rs(Lgood&~Ltin&wk<0,rr)), ...
        rocN(rs(Lgood&Ltin&wk<0,rr),      rs(Lgood&~Ltin&wk>0,rr)), ...
        rocN(rs(Lgood&Ltin&wk>0,rr),      rs(Lgood&Ltin&wk<0,rr)), ...
        rocN(rs(Lgood&~Ltin&wk>0,rr),     rs(Lgood&~Ltin&wk<0,rr))];
end

% fit the dots stuff
bins = [(0:50:1400)', (100:50:1500)'];
nb   = size(bins,1);
ds   = nans(size(Lgood,1), nb);
[r,zs,ds(Lgood,:)] = getFIRA_rateByBin(Lgood, sid, times(Lgood,2), times(Lgood,3), [], bins);
ucoh = nonanunique(cohs);
nc   = length(ucoh);
Lwk  = [wk<0 wk>0];
drocs = nans(nb,nc,2);
for cc = 1:nc
    Lcoh = Lgood & Lcor & cohs == ucoh(cc);
    if sum(Lcoh) > 5
        for bb = 1:nb
            Lb =  Lcoh & isfinite(ds(:,bb));
            for ww = 1:2
                Lw = Lwk(:,ww) & Lb;
                if sum(Lw&Ltin) > 4 && sum(Lw&~Ltin) > 4
                    drocs(bb,cc,ww) = rocN(ds(Lw&Ltin,bb), ds(Lw&~Ltin,bb));
                end
            end
        end
    end
end

diffs = reshape(drocs(3:17,:,2)-drocs(3:17,:,1),[],1);
[H,P] = ttest(diffs);
ds_   = [nanmean(diffs) nanmedian(diffs) nanse(diffs) P signtest(diffs) sum(isfinite(diffs))];

return

co = {'y' 'c' 'm' 'r' 'g' 'b' 'k'};
for ww = 1:2
    subplot(3,1,ww); cla reset; hold on;
    plot([0 1000], [0.5 0.5], 'k:');
    for cc = 1:nc
        plot(mean(bins,2), drocs(:,cc,ww), co{cc}, 'LineStyle', '-');
    end
    axis([0 1000 0.2 0.9]);
end
subplot(3,1,3); cla reset; hold on;
plot([0 1], [0 1], 'k:');
plot([0 1], [0.5 0.5], 'k:');
plot([0.5 0.5], [0 1], 'k:');
for cc = 1:nc
    plot(drocs(:,cc,1), drocs(:,cc,2), [co{cc} '.']);
end
axis([0.2 1 0.2 1]);

r = input('next')
