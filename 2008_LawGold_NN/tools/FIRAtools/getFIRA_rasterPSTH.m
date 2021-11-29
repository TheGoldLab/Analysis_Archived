function [ras, psth_x, psth_y, psth_ts] = getFIRA_rasterPSTH(trials, Lsb, spikei, ...
    raster_begin, raster_end, raster_wrt, bin_size, sortras)
% function [raster, xs, ys, yd] = getFIRA_rasterPSTH(fname, trials, Lsb, spikei, ...
%    raster_begin, raster_end, bin_size)
%
% Returns:
%   ras:     mx1 cell array, raster, one per Lsb column
%   psth_y:  nxm matrix. The Nth row represents spike rate in the nth time bin.
%            The mth column represents spike rate for trials selected by the mth Lsb column
%   psth_ts: nxm martix. Number of trials in each conditions. Organized
%            similarly to psth_y.
%
%

% modified from plotFIRA_raterPSTH.m (written by jig) by ctl on 8/10/05

global FIRA

% easy check for now
if nargin < 7 | isempty(FIRA) | isempty(FIRA.spikes.data) | ...
        isempty(trials) | isempty(Lsb) | isempty(spikei) | ...
        spikei < 1 | spikei > size(FIRA.spikes.data,2)
    return
end

if isempty(sortras)
    sortras = 0;    % default don't sort
end



% length of ras_hs & psth_hs are either 0 or same
% as the number of columns of Lsb
num_axes = size(Lsb, 2);
num_bins = ceil(max(raster_end-raster_begin)/bin_size);

% default out
ras = cell(num_axes,1);
psth_x  = nans(num_bins,num_axes);
psth_y  = nans(num_bins,num_axes);
psth_ts = nans(num_bins,num_axes);

%% check times...
FIX_TIME = 0;
% raster_begin determines start time of the plotted raster
if isempty(raster_begin)
    raster_begin = -FIX_TIME*ones(size(trials));
else
    raster_begin(isnan(raster_begin)) = -FIX_TIME;
end

% raster_end determines end time of the plotted raster
if isempty(raster_end)
    raster_end = FIX_TIME*ones(size(trials));
else
    raster_end(isnan(raster_end)) = FIX_TIME;
end

% fix wrt to raster_begin for now
raster_wrt(isnan(raster_wrt)) = 0;
wrt                               = raster_wrt;

%%%
% DATA LOOP
%%%
%
% loop through the axes (angles)
for i = 1:num_axes

    % get the trial indices
    trs    = find(Lsb(:,i));
    num_tr = length(trs);
    raster = [];
        
    % sort, if neccesary
    if sortras % sort by raster time
        [I,Y] = sort(raster_end(trs) - raster_begin(trs));
        trs = trs(Y);
    end

    % (hopefully) faster way of keeping track of the number of
    %   trials per bin, for the psth
    mint = floor(min(raster_begin(trs)-wrt(trs)));
    maxt = ceil(max(raster_end(trs)-wrt(trs)));
    lent = maxt - mint;
    
    if  lent < 10000 && ...
            all(~isnan(raster_begin(trs))) && all(~isnan(raster_end(trs))) && ...
            all(~isnan(wrt(trs)))
        psth_times = zeros(1, lent+1);
        psth_bts   = floor(raster_begin(trs) - wrt(trs)) - mint + 1;
        psth_ets   = ceil(raster_end(trs) - wrt(trs)) - mint + 1;
    else
        lent       = [];
        psth_times = [];
    end
    
    % loop through the trials
    for ti = 1:num_tr

        t = trs(ti);    
        
        % get spikes from this trial
        sp = FIRA.spikes.data{trials(t), spikei};

        if ~isempty(sp)

            % raster
            msp = sp(sp>=raster_begin(t)&sp<=raster_end(t));
            if ~isempty(msp)
                raster = [raster; msp-wrt(t) ti*ones(size(msp))];
            end

       
            % psth
            if ~isempty(lent)
                
                % this is the way we like to compute it -- increment counts
                % for the given time interval
                
                psth_times(psth_bts(ti):psth_ets(ti)) = ...
                    psth_times(psth_bts(ti):psth_ets(ti)) + 1;
                
            elseif ~isinf(raster_begin(t)) && ~isinf(raster_end(t))
                
                % this is a much slower way -- keep track of all the times
                psth_times = [psth_times; [raster_begin(t):raster_end(t)]'-wrt(t)];

            else
                
                % this way sucks -- we don't really even know when the
                % trial begins & ends
                psth_times = [psth_times; [min(sp):max(sp)]'-wrt(t)];
                
            end
        end
    end

    if ~isempty(raster)
        % plot psth
        if ~isempty(lent)
            xs = [mint:bin_size:maxt]';
            ys = histc(raster(:,1), xs)*1000; % ctl & jig removed bin_size scale term, 6/27/05
            ys = ys(:);
            yd  = zeros(bin_size, length(xs));
            len = min(length(psth_times), length(xs)*bin_size);
            yd(1:len) = psth_times(1:len);
            yd = sum(yd,1)';
        else
            xs = [min(psth_times):bin_size:max(psth_times)]';
            ys = histc(raster(:,1), xs)*1000/bin_size;
            ys = ys(:);
            yd = histc(psth_times, xs);
        end
        Lp      = yd>0;
        ys(Lp)  = ys(Lp)./yd(Lp);
        ys(~Lp) = 0;

    end

    % assign output
    ras{i,1}     = raster;
    psth_x(1:length(xs),i)  = xs;
    psth_y(1:length(ys),i)  = ys;
    psth_ts(1:length(yd),i) = yd;

end