function [fits, fsem, fstat, N] = getML_ROCT(uid, d1name, d1, bins, lambda, dirf)
% [fits, fsem, fstat, N] = getML_ROCT(uid, d1name, d1, bins)
%  get ROC for spikes in time bins BINS=[t0 t1; t1 t2; ... ; tn-1 tn]
%   INPUTS:
%       uid    - unit id in nex (eg 3001, 9001, 9002)
%       d1name - name for direction 1 (preferred direction) in
%                     FIRA.ecodes.name. (e.g. 'trg1_dir', 'dot_dir')
%       d1     - preferred direction in degree
%       bins   - time bins for which ROC will be computed
%                BINS=[t0 t1; t1 t2; ... ; tn-1 tn]
%
%   OUTPUTS:
%       fits - 3xn fits parameters for each time bins
%       fsem - 3xn sem of the fitted parameters
%       fstat- 3xn LLR and p-value of the fit
%       N    - 1xn number of trials in each time bin
%
% assumed that data is loaded in global variable FIRA

% created by jcl on 9/17/05

% check inputs
if nargin<3
    error('getML_ROCT needs at least three inputs')
elseif nargin<4
    bins = [0 1500];
elseif nargin<5
    lambda = [];
    dirf = 0;
elseif nargin<6
    dirf = 0;
end


global FIRA

% default outs
fits  = nans(3,size(bins,1));
fsem  = nans(3,size(bins,1));
fstat = nans(3,size(bins,1));
N     = nans(1,size(bins,1));

% get selection arrays for plotFIRA_neurometricROC
trials = find(~isnan(FIRA.ecodes.data(:,1)));
[Ltsk, Utsk] = selectFIRA_trialsByUniqueID('task', trials); % select only task 3 (fixation) and task 6 (dots)

if any(ismember(Utsk,3) | ismember(Utsk,6))
    if ismember(Utsk,6) & ismember(Utsk,3)
        Ltsk = Ltsk(:,Utsk==3) | Ltsk(:,Utsk==6);
    elseif any(ismember(Utsk,3))
        Ltsk = Ltsk(:,Utsk==3);
    else
        Ltsk = Ltsk(:,Utsk==6);
    end
else
    fprintf('FIRA contains neither task 3 (fixation) nor task 6 (dots)\n')
    return;
end

[Lcrt, Ucrt] = selectFIRA_trialsByUniqueID('correct', trials); % select for correct and incorrect trials
L1 = Lcrt(:,Ucrt==1);
if ~isempty(Lcrt(:,Ucrt==0))
    L1 = Lcrt(:,Ucrt==0)|L1;
end

[Lcoh, Ucoh] = selectFIRA_trialsByUniqueID('dot_coh', trials);
[Ldir, Udir] = selectFIRA_trialsByUniqueID(d1name, trials);
Ldir = [(Ldir(:,round(Udir)==d1) & L1 & Ltsk)...
         (Ldir(:, round(Udir)==mod(d1+180,360)) & L1 & Ltsk)];
Lt = L1 & Ltsk;

if dirf  % transfor dot/trg dir to choice dir
    for j = 1:2
        L = ~isnan(FIRA.ecodes.data(Lt,getFIRA_ecodeColumnByName('correct')));
        Ldir(Lt,j) = ~xor(Ldir(Lt,j), FIRA.ecodes.data(Lt,getFIRA_ecodeColumnByName('correct')));
    end
end

durs = getFIRA_ecodeTimesByName('dot_off', 0, trials)...
        - getFIRA_ecodeTimesByName('dot_on', 0, trials);
bt    = getFIRA_ecodeTimesByName('dot_on', 0, trials);

for i = 1:size(bins,1)
    % remove short trials
    Ldt = durs > max(0, bins(i,2));
    if any(Ldt)
        trials_ = trials(Ldt);
        Lcoh_   = Lcoh(Ldt, :);
        Ldir_   = Ldir(Ldt, :);
        bt_     = bt(Ldt);        
        [fits_, stats_] = plotFIRA_neurometricROC(trials_, Ldir_, Lcoh_, Ucoh,...
                                                getFIRA_spikeByID(uid),...
                                                bt_+bins(i,1), bt(Ldt)+bins(i,2), lambda);
    fits(:,i)  = fits_(:,1);
    fsem(:,i)  = fits_(:,2);
    fstat(:,i) = stats_';
    N(:,i)     = sum(Ldt);
    end
end