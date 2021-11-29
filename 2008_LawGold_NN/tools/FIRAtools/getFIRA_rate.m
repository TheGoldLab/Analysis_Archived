function rate = getFIRA_rate(trials, spikei, rate_begin, rate_end)
% Compute spike rate for each trial TRIALS between rate begin and rate end
%  function rate = getFIRA_rate(trials, spikei, rate_begin, rate_end)%
% 
% INPUTS:
%   trials      : selection array for trials
%   spikei      : spike id
%   rate_begin  : rate begin time
%   rate_end    : rate end time
% OUTPUTS:
%   rate        : mx1 matrix, where m=total number of trials, including
%                 those that are not asserted (0) in selection array
%                 TRIALS.
%                 Each entry contains the
%                 instantaneous spike rate in the m-th trial.
%

% modified from plotFIRA_raterPSTH.m (written by jig) by ctl on 8/10/05

global FIRA
warning off
% easy check for now
if nargin < 4 | isempty(FIRA) | isempty(FIRA.spikes.data) | ...
        isempty(trials) | isempty(spikei) | spikei < 1 | ...
        spikei > size(FIRA.spikes.data,2)
    return
end


% default out
rate = nans(length(trials),1);


%%%
% DATA LOOP
%%%
% loop through trials and get instantaneous spike rate
trs = find(trials);
for i = 1:length(trs)
    t  = trs(i);
    sp = FIRA.spikes.data{t, spikei};

    if ~isempty(sp)
        rate(t) = sum(sp>=(rate_begin(t)) & sp<=(rate_end(t)))...
            / ((rate_end(t)-rate_begin(t))/1000);
    else
        rate(t) = 0;
    end
end
warning on