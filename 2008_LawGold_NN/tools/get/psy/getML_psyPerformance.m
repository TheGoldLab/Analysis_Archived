function [fits, sems] = getML_psyPerformance(fn, recompute)

if nargin < 1
    return;
elseif nargin < 2
    recompute = 1;
end

if recompute
    a     = getML_txt(fn);
    fname = a.data{strcmp(a.name,'dat_fn')};
    fits = nans(3, length(fname));
    sems = nans(3, length(fname));

    % get threshold
    global FIRA
    for i = 1:length(fname)
        fprintf('%d: %s\n', i, fname{i})
        openFIRA(fname{i})

        [Lcoh, Ucoh] = selectFIRA_trialsByUniqueID('dot_coh');
        [Lcrt, Ucrt] = selectFIRA_trialsByUniqueID('correct');
        Lgd          = Lcrt(:,ismember(Ucrt,0)) | Lcrt(:,ismember(Ucrt,1));     % remove no choice, broken fixation
        [fits(:,i), sems(:,i)]  = quick_fit(quick_formatData2(Ucoh,Lcoh(Lgd, :),Lcrt(Lgd,ismember(Ucrt,1))));
    end

    % save
    s = which('getML_psyPerformance.m');
    save([s(1:end-2) fn(1:end-4) '.mat'], 'fits', 'sems')
else   
    % load
    s = which('getML_psyPerformance.m');
    load([s(1:end-2) fn(1:end-4) '.mat'])
end
    