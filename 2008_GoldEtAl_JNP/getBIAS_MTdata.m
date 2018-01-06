function [wrd_, pd_, exdat_, lag_, r4_] = getBIAS_MTdata(sid, randomize, exd, lag, r4)
% 8 data sets of pd_ are: all, pref, null, pref/cor, pref/err, null/cor, null/err, 0% coh

global FIRA

% make output arrays
wrd_ = nans(1, 4, 2);
pd_  = nans(1, 3, 8, 3);

% check for randomize
if nargin > 2 && randomize
    Frnd             = randperm(size(FIRA.ecodes.data, 1));
    FIRA.ecodes.data = FIRA.ecodes.data(Frnd, :);
    FIRA.spikes.data = FIRA.spikes.data(Frnd, :);
    if nargin > 4 && ~isempty(r4)
        [a,b] = ismember(Frnd, r4(:,1));
        r4    = r4(b(b>0),2);
    end
else
    randomize = false;
end

% dumb fix
dd = getFIRA_ecc('dot_dur');
co = getFIRA_ecc('correct');
FIRA.ecodes.data(~isfinite(FIRA.ecodes.data(:,dd)), co) = -2;

%size(find(FIRA.ecodes.data(:,co)==1&isnan(FIRA.ecodes.data(:,dd))))

% get pref, 
[pref,wrd_(1,4)] = getBIAS_tinFromFIRA(sid);

% get wk value
if randomize && nargin > 4 && ~isempty(lag) && ~isempty(r4)
    [wk,wrd_(1,1:3,1)]          = getBIAS_wkFromFIRA(pref, lag, r4);
else
    [wk,wrd_(1,1:3,1),lag_,r4_] = getBIAS_wkFromFIRA(pref);
end
wk = wk./nanmean(abs(wk));

% selectors
times = getFIRA_ecodesByName({'dot_on', 'dot_off', 'fp_off'});
cohs  = getFIRA_ecodesByName('dot_coh');
cors  = getFIRA_ecodesByName('correct');
spref = sign(cos(pref*pi/180));
chcs  = spref.*sign(cos(getFIRA_ecodesByName('choice')*pi/180));
%dirs  = spref.*sign(cos(getFIRA_ecodesByName('dot_dir')*pi/180));
Lgood = getFIRA_ecodesByName('task')==6 & cors>=0;
Lcor  = cohs==0 | cors==1;
Lpref = cohs==0 | chcs==1;

% fit to PMF, with wk
%ctdat = [cohs./100 (times(:,2)-times(:,1))./1000 dirs chcs cors];
%fits  = ctPsych_fit(@ddExp5fz_F, ctdat(Lgood, 1:4), ctdat(Lgood, 5), ...
%    [], [], [], [], wk(Lgood));
%wk    = wk.*fits(4);

% save mean,sem abs wk
wrd_(1,1:2,2) = [nanmean(abs(wk(Lgood))) nanse(abs(wk(Lgood)))];

% get data
[r,zs,pre_rs] = getFIRA_rateByBin(Lgood, sid, times(Lgood,1)-300, times(Lgood,1));
[r,zs,per_rs] = getFIRA_rateByBin(Lgood, sid, times(Lgood,1)+75,  times(Lgood,1)+1000);
[r,zs,pos_rs] = getFIRA_rateByBin(Lgood, sid, times(Lgood,3)-300, times(Lgood,3));
zdat          = nans(size(Lgood,1), 3);
zdat(Lgood,:) = [pre_rs per_rs-pre_rs pos_rs-pre_rs];
Lz            = isfinite(zdat);

% z-score all pre
zdat(Lz(:,1),1) = zscore(zdat(Lz(:,1),1));

% z-score peri,post by dir/coh
Ltc   = [Lgood&Lpref&Lcor Lgood&~Lpref&Lcor Lgood&Lpref&~Lcor Lgood&~Lpref&~Lcor];
ucohs = nonanunique(cohs);
for tt = 1:size(Ltc,2)
    for cc = 1:length(ucohs)
        for pp = 2:3
            Lcz = Lz(:,pp) & Ltc(:,tt) & cohs==ucohs(cc);
            if any(Lcz)
                zdat(Lcz,pp) = zscore(zdat(Lcz,pp));
            end
        end
    end
end

% 8 data sets: all, pref, null, pref/cor, pref/err, null/cor, null/err, 0% coh
Lds = [Lgood Lgood&Lpref Lgood&~Lpref Lgood&Lpref&Lcor Lgood&Lpref&~Lcor ...
    Lgood&~Lpref&Lcor Lgood&~Lpref&~Lcor Lgood&cohs==0];
if randomize
    num = 1;
else
    num = size(Lds,2);
end
for ff = 1:num
    for pp = 1:size(zdat,2) % pre/peri/post
        Ly = Lds(:,ff) & Lz(:,pp);
        if sum(Ly) > 10
            [b,sem] = lscov([ones(sum(Ly),1) wk(Ly)], zdat(Ly,pp));
            pd_(1,:,ff,pp) = [b(2) sem(2) sum(Ly)];
%             [R,P,RLO,RUP]      = corrcoef(wk(Ly), zdat(Ly,pp));
%             pd_(1,:,ff,pp,1) = [R(2,1) RLO(2,1) RUP(2,1) P(2,1)];
        end
    end
end

% save example data
if nargin > 2 && exd
    % get regression inputs
    exdat_ = {};
    for pp = 1:size(zdat,2)
        Ly = Lds(:,1) & Lz(:,pp);
        exdat_ = cat(2,exdat_,[wk(Ly) zdat(Ly,pp)]);
    end
    % get PSTH, correct trials only
    b1s = 0:10:1500;
    b1s = [b1s' b1s'+100];
    b2s = -600:10:100;
    b2s = [b2s' b2s'+100];

    % correct only, pref/null, all cohs
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