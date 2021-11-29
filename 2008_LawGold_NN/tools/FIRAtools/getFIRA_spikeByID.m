function index = getFIRA_spikeByID(id)
% function ind = getFIRA_spikeByID(id)
%
% Returns the index of the spike id in FIRA.spikes.id
% Usage:
%      index = getFIRA_spikeByID(id);

% created 8/10/05 by ctl

global FIRA

% check
if isempty(FIRA.spikes.id) | nargin<1 | length(id)>1
    return
end

% default out
index = [];

% get spikes ind
index = find(FIRA.spikes.id==id);