function [fits, fsem, fstat] = getML_ROC(uid, d1name, d1, lambda, dirf)
% [fits, fsem, fstat, N] = getML_ROCT(uid, d1name, d1, bins)
%  get ROC for spikes between dot_on dot_off
%   INPUTS:
%       uid    - unit id in nex (eg 3001, 9001, 9002)
%       d1name - name for direction 1 (preferred direction) in
%                     FIRA.ecodes.name. (e.g. 'trg1_dir', 'dot_dir')
%       d1     - preferred direction in degree
%
%   OUTPUTS:
%       fits - 3 fits parameters
%       fsem - 3 sem of the fitted parameters
%       fstat- 3 LLR and p-value of the fit
%
% assumed that data is loaded in global variable FIRA

% created by jcl on 9/17/05

% check inputs
if nargin<3
    error('getML_ROC needs at least three inputs')
elseif nargin<4
    lambda = [];
    dirf = 0;
elseif nargin<5
    dirf = 0;
end


global FIRA

% default outs
fits  = nans(3);
fsem  = nans(3);
fstat = nans(3);

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

bt    = getFIRA_ecodeTimesByName('dot_on', 0, trials);
et    = getFIRA_ecodeTimesByName('dot_off', 0, trials);
[fits_, stats_] = plotFIRA_neurometricROC(trials, Ldir, Lcoh, Ucoh,...
    getFIRA_spikeByID(uid), bt, et, lambda);


fits  = fits_(:,1);
fsem  = fits_(:,2);
fstat = stats_';
