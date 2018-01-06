function [wrd_, pd_, exdat_] = getBIAS_MTdata2(sid, exd)
% 8 data sets of pd_ are: all, pref, null, pref/cor, pref/err, null/cor, null/err, 0% coh

global FIRA

% make output arrays
wrd_ = nans(1, 2, 8);
pd_  = nans(1, 3, 8, 3);

% dumb fix
dd = getFIRA_ecc('dot_dur');
co = getFIRA_ecc('correct');
FIRA.ecodes.data(~isfinite(FIRA.ecodes.data(:,dd)), co) = -2;

%size(find(FIRA.ecodes.data(:,co)==1&isnan(FIRA.ecodes.data(:,dd))))

% get pref, 
pref = getBIAS_tinFromFIRA(sid);
wk   = getBIAS_wkFromFIRA(pref);

% selectors
times = getFIRA_ecodesByName({'dot_on', 'dot_off', 'fp_off'});
cohs  = getFIRA_ecodesByName('dot_coh');
cors  = getFIRA_ecodesByName('correct');
spref = sign(cos(pref*pi/180));
chcs  = spref.*sign(cos(getFIRA_ecodesByName('choice')*pi/180));
dirs  = spref.*sign(cos(getFIRA_ecodesByName('dot_dir')*pi/180));
Lgood = isfinite(wk) & getFIRA_ecodesByName('task')==6 & cors>=0;
Lcor  = cohs==0 | cors==1;
Lpref = cohs==0 | chcs==1;
wk    = wk.*chcs;

% fit to PMF, with wk
ctdat = [cohs./100 (times(:,2)-times(:,1))./1000 dirs chcs cors];
fits  = ctPsych_fit(@ddExp5fz_F, ctdat(Lgood, 1:4), ctdat(Lgood, 5), ...
    [], [], [], [], wk(Lgood));
wkf   = wk.*fits(4);

% get data
rdat = nans(size(Lgood,1),3);
[r,zs,rdat(Lgood,1)] = getFIRA_rateByBin(Lgood, sid, times(Lgood,1)-300, times(Lgood,1));
[r,zs,rdat(Lgood,2)] = getFIRA_rateByBin(Lgood, sid, times(Lgood,1)+75,  times(Lgood,1)+1000);
[r,zs,rdat(Lgood,3)] = getFIRA_rateByBin(Lgood, sid, times(Lgood,3)-300, times(Lgood,3));
Lr                   = isfinite(rdat);

% 8 data sets: all, pref, null, pref/cor, pref/err, null/cor, null/err, 0% coh
Lds = [Lgood Lgood&Lpref Lgood&~Lpref Lgood&Lpref&Lcor Lgood&Lpref&~Lcor ...
    Lgood&~Lpref&Lcor Lgood&~Lpref&~Lcor Lgood&cohs==0];
for ff = 1:size(Lds,2)

    % save mean,sem abs wk
    wrd_(1,:,ff) = [nanmean(abs(wkf(Lds(:,ff)))) nanmedian(abs(wkf(Lds(:,ff))))];

    % compute correlations
    for pp = 1:size(rdat,2) % pre/peri/post
        Ly = Lds(:,ff) & Lr(:,pp);
        if sum(Ly) > 10
            if ff<8
                [R,P] = partialcorr([rdat(Ly,pp) wk(Ly)], cohs(Ly), 'type', 'Spearman');
                R = R(2,1);
                P = P(2,1);
            else
                [R,P] = corr(rdat(Ly,pp), wk(Ly), 'type', 'Spearman');
            end
            V     = (1+R.^2/2)*(1-R.^2).^2./(sum(Ly)-3);
            pd_(1,:,ff,pp) = [R, P, V];
        end
    end
end

% save example data
if nargin > 1 && exd
    % get regression inputs
    exdat_ = {};
    for pp = 1:size(rdat,2)
        Ly = Lds(:,1) & Lr(:,pp);
        exdat_ = cat(2,exdat_,[wk(Ly) rdat(Ly,pp)]);
    end
    % get PSTH, correct trials only
    b1s = 0:10:1500;
    b1s = [b1s' b1s'+100];
    b2s = -600:10:100;
    b2s = [b2s' b2s'+100];

    % correct only, pref/null, all cohs
    ucohs = nonanunique(cohs);
    e1s  = nans(size(b1s,1), length(ucohs), 2);
    e2s  = nans(size(b2s,1), length(ucohs), 2);
    Lps  = [Lgood&Lcor&Lpref Lgood&Lcor&~Lpref];
    for pp = 1:2
        for cc = 1:length(ucohs)
            Lall = Lps(:,pp) & cohs==ucohs(cc);
            if any(Lall)
                [r1,z1s,r1s] = getFIRA_rateByBin(Lall, sid, times(Lall,1)-400, ...
                    times(Lall,2), times(Lall,1)-400, b1s);
                e1s(:,cc,pp) = nanmean(r1s)';
                [r2,z2s,r2s] = getFIRA_rateByBin(Lall, sid, times(Lall,3)-600, ...
                    times(Lall,3)+200, times(Lall,3), b2s);
                e2s(:,cc,pp) = nanmean(r2s)';
            end
        end
    end
    exdat_ = cat(2, exdat_, {b1s e1s b2s e2s});
else
    exdat_ = {};
end