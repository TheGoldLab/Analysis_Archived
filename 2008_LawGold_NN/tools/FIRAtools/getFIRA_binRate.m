function rate = getFIRA_binRate(trials, spikei, rate_begin, rate_end, rate_wrt, bins)
% Compute spike rate for each trial TRIALS for each bin in BINS.
%  function rate = getFIRA_rateBined(trials, spikei, rate_begin, rate_end, rate_wrt, bins)%
% 
% INPUTS:
%   trials      : selection array for trials
%   spikei      : spike id
%   rate_begin  : rate begin time
%   rate_end    : rate end time
%   rate_wrt    : rate with respect to time
%   bins        : time bins for which spike rate is computed for each trials.
%                 BINS = [t0 t1; t1 t2;... tn-1 tn], bin edges can be overlapped.
% OUTPUTS:
%   rate        : mxn matrix, where m=total number of trials, including
%                 those that are not asserted (0) in selection array TRIALS,
%                 and n=total number of bins. Each entry contains the
%                 instantaneous spike rate in the m-th trial and the n-th
%                 bin.
%

% modified from plotFIRA_raterPSTH.m (written by jig) by ctl on 8/10/05

global FIRA


% easy check for now
if nargin < 6 | isempty(FIRA) | isempty(FIRA.spikes.data) | ...
        isempty(trials) | isempty(spikei) | spikei < 1 | ...
        spikei > size(FIRA.spikes.data,2)
    return
end


num_bins = size(bins,1);


% default out
rate = nans(length(trials), num_bins);


%%%
% DATA LOOP
%%%
%
% loop through the trials, get spike rate for each bins
% get the trial indices
trs = find(trials);
for i = 1:length(trs)
    % get spikes from this trial
    t  = trs(i);
    sp = FIRA.spikes.data{t, spikei};

    if ~isempty(sp)
        % loop through bins
        for j = 1:num_bins
            if rate_end(t) > rate_begin(t) ...
                    & bins(j,1)>=rate_begin(t)-rate_wrt(t) ...
                    & bins(j,2)<=rate_end(t)-rate_wrt(t)

                rate(t,j) = sum(sp>=(bins(j,1)+rate_wrt(t)) & sp<=(bins(j,2)+rate_wrt(t)))...
                             / ((bins(j,2)-bins(j,1))/1000);
            end
        end
    else
        % check if it's a good trial, it is either a correct or an
        % incorrect trial, give it a zero
        if ismember(FIRA.ecodes.data(t,getFIRA_ecodeColumnByName('correct')), [0 1])
            % loop through bins
            for j = 1:num_bins
                if rate_end(t) > rate_begin(t) ...
                        & bins(j,1)>=rate_begin(t)-rate_wrt(t) ...
                        & bins(j,2)<=rate_end(t)-rate_wrt(t)

                    rate(t,j) = 0;
                end
            end
        end
    end
end


