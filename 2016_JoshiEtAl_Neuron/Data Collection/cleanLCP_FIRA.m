 function cleanLCP_FIRA(clean_name, pulseOxFlag)
% function cleanLCP_FIRA(clean_name, pulseOxFlag)
%
% Cleans up and saves siteData (in file "clean_name"):
%  siteData{1}: trialsxcols matrix, cols are:
%   1 ... fix start time wrt fixation on
%   2 ... fix end time wrt fix start time (fix duration)
%   3 ... reported correct
%   4 ... beep on time (when appropriate), wrt to fix start time
%   5 ... trial begin time, wrt to fix start time
%   6 ... trial end time, wrt to fix start time
%   7 ... trial wrt time (cpu clock)
%   8 ... LFP index corresponding to fix start time (coded above)
%   9 ... ELESTM on time (when appropriate), wrt fix start time
%
%  siteData{2}: Analog:
%   dim1: trial
%   dim2: sample
%   dim3: 1 = x, 2 = y, 3 = z-pupil, 4 = corrected z-pupil, 5 = pupil slope
%   [remember first sample is considered time=0 (i.e., wrt fix start time)]
%     = eyedat(Lgood,:,:);
%
%  siteData{3}: spikes, re-coded wrt fix start time
%
%  siteData{4}: LFP
%
%  siteData{5}: pupil events
%   1. trial number
%   2. start time of event (wrt fix start time)
%   3. end time of event (wrt fix start time)
%   4. magnitude at start of event (raw z-score)
%   5. magnitude at end of event (raw z-score)
%   6. magnitude at start of event (corrected z-score)
%   7. magnitude at end of event (corrected z-score)
%   8. time of subsequent max slope
%   9. magnitude of subsequent max slope (corrected z/sample)
%
%  siteData{6}: microsaccades
%   1. trial number
%   2. start time of event (wrt fix start time)
%   3. duration of event (wrt fix start time)
%   4. maximum velocity (deg/ms)
%   5. magnitude of microsaccade event (deg)
%   6. onset time wrt phase of associated pupil event (fraction)
%   7. magnitude of associated pupil event
%
%  siteData{7}: spike/LFP channels
%   1. spike channels
%   2. LFP names

% acts on existing FIRA
global FIRA

% make filter for running linear slope
FILTER_HW = 75;
filter_c  = (3 / (2*FILTER_HW^3+3*FILTER_HW^2+FILTER_HW))*(FILTER_HW:-1:-FILTER_HW);

% create cell array
siteData  = cell(1,7); % ecodes, eye, spikes, LFP, pupil events, microsaccades

% Analog channel indices
ai_p = find(strncmp('pup', FIRA.analog.name, 3), 1);
ai_h = find(strncmp('hor', FIRA.analog.name, 3), 1);
ai_v = find(strncmp('ver', FIRA.analog.name, 3), 1);
if isempty(ai_p) || isempty(ai_h) || isempty(ai_v)
    disp('BAD ANALOG CHANNELS')
end
ai_l = find(strncmp('LFP', FIRA.analog.name, 3) | ...
    strncmp('unk', FIRA.analog.name, 3));

% selection stuff
fp_times   = getFIRA_ec([], {'fp_on', 'fp_off'});
num_trials = size(FIRA.ecodes.data, 1);
correct    = getFIRA_ec([], 'correct');
beep_on    = getFIRA_ec([], 'beep_on');
stim_on    = getFIRA_ec([], 'stim_time');
disp(sprintf('num stim = %d', sum(isfinite(stim_on))))
maxx       = 0;

% array of start/end times of good fixation
fix_times  = nans(num_trials, 2);

% for eye data:1-Horiz, 2-Vert, 3-Z-pupil,4-Corrected z-pupil, 5-pupil slope
eyedat     = nans(num_trials, 5500, 5); %

% special case of getting pulse-ox data
if nargin < 2
    pulseOxFlag = false;
    poxdat = nans(num_trials, 5500);
end

% loop through the trials to find good fixation epochs
for tt = 1:num_trials
    
    % find beginning, end times
    fp_start = round(fp_times(tt,1) + 500 - FIRA.analog.data(tt,ai_p).start_time);
    fp_end   = min(FIRA.analog.data(tt,ai_p).length-20, ...
        round(fp_times(tt,2) - FIRA.analog.data(tt,ai_p).start_time));
    if fp_end > 10000
        fp_start = fp_end - 5000;
    end
    
    % if there is something there
    if fp_end > fp_start + 100 && fp_start > 0
        
        % get all saccades between fix on and fix off
        sacs = findSaccadesA( ...
            FIRA.analog.data(tt,ai_h).values(fp_start:fp_end), ...
            FIRA.analog.data(tt,ai_v).values(fp_start:fp_end), ...
            FIRA.analog.store_rate(ai_p), 7, false);
        if isempty(sacs)
            sacs = cat(1, zeros(1,2), ...
                [fp_end-fp_start 0]);
        else
            sacs = cat(1, zeros(1,2), ...
                sacs(isfinite(sacs(:,1)),[1 2]), ...
                [fp_end-fp_start 0]);
        end
        
        %% SPECIAL CASE, ADDED BY JIG 1/2/15
        % on ustim trials, keep first saccade after fixation
        if isfinite(stim_on(tt)) && size(sacs,1)>2 && ...
                any(sacs(2:end-1,1)>stim_on(tt))
            sacs(find(sacs(2:end-1)>stim_on(tt),1),:) = [];
        end
        nsacs = size(sacs,1);
        
        % find max time between saccades while maintaining fixation
        gi = zeros(1,3);
        for ss = 1:nsacs-1
            
            % begin/end indices of this interval
            i1 = fp_start + round(sacs(ss,1) + sacs(ss,2)) + 1;
            i2 = fp_start + round(sacs(ss+1,1)) - 1;
            
            % get median values of x,y positions during this interval
            medx = median(FIRA.analog.data(tt,ai_h).values(i1:i2));
            medy = median(FIRA.analog.data(tt,ai_v).values(i1:i2));
            
            % make sure x/y medians are reasonable
            if abs(medx) < 5 && abs(medy) < 5
                
                % subtract out median
                FIRA.analog.data(tt,ai_h).values(i1:i2) = ...
                    FIRA.analog.data(tt,ai_h).values(i1:i2) - medx;
                FIRA.analog.data(tt,ai_v).values(i1:i2) = ...
                    FIRA.analog.data(tt,ai_v).values(i1:i2) - medy;
                
                % find lag to fixation
                offset = find(abs(FIRA.analog.data(tt,ai_h).values(i1:i2))<2.5 & ...
                    abs(FIRA.analog.data(tt,ai_v).values(i1:i2))<2.5, 1);
                i1 = i1 - 1 + offset;
                if i1 + 200 < i2 % at least 200 ms
                    
                    % check good fixation
                    Lfix = abs(FIRA.analog.data(tt,ai_h).values(i1:i2))<2.5 & ...
                        abs(FIRA.analog.data(tt,ai_v).values(i1:i2))<2.5 & ...
                        FIRA.analog.data(tt,ai_p).values(i1:i2)>0;
                    nfix = sum(Lfix);
                    
                    % got enough good samples & this is longest
                    % interval
                    if nfix > 0.8*length(Lfix) && gi(1) < nfix
                        if any(~Lfix)
                            % get rid o' crap
                            FnoFix = find(~Lfix);
                            for aa = [ai_p ai_h ai_v]
                                FIRA.analog.data(tt,aa).values(i1+FnoFix-1) = nan;
                            end
                        end
                        gi = [nfix i1 i2];
                    end
                end
            end
        end
        
        % found something, save associated data
        if gi(1) > 200
            if gi(3)-gi(2)+1 > 5500
                gi(3) = gi(2) + 5499;
            end
            len = gi(3)-gi(2)+1;
            if len > maxx
                maxx = len;
            end
            
            % save start,end times
            fix_times(tt,:) = FIRA.analog.data(tt,ai_p).start_time + gi([2 3]);
            % eyedat(tt,1:len,1) = FIRA.analog.data(tt,ai_h).values(gi(2):gi(3));
            % eyedat(tt,1:len,2) = FIRA.analog.data(tt,ai_v).values(gi(2):gi(3));
            % pd = FIRA.analog.data(tt,ai_p).values(gi(2):gi(3))';            
            
            % possible linear interpolation of nans in x/y/pd
            inds = [ai_h ai_v ai_p];
            for xx = 1:length(inds)
                vals = FIRA.analog.data(tt,inds(xx)).values(gi(2):gi(3))';
                i1   = find(isfinite(vals),1);
                i2   = find(isfinite(vals),1,'last');
                
                % fill in nans by linear interpolation
                if ~isempty(i1) && ~isempty(i2) && any(~isfinite(vals(i1:i2)))
                    valsn       = vals(i1:i2);
                    Lnan        = ~isfinite(valsn);
                    valsn(Lnan) = linterp(find(~Lnan), valsn(~Lnan), find(Lnan));
                    vals(i1:i2) = valsn;
                end
                eyedat(tt,1:len,xx) = vals;%nanrunmean(pd,50);
            end
            
            % re-code spikes wrt fixation start
            if isfield(FIRA, 'spikes')
                for ss = 1:size(FIRA.spikes.data,2)
                    FIRA.spikes.data{tt,ss} = FIRA.spikes.data{tt,ss} - fix_times(tt,1);
                end
            end
            
            % re-code LFP start time OR ADD PULSE-OX
            if pulseOxFlag
                poxdat(tt,1:len) = ...
                    FIRA.analog.data(tt,ai_l).values(gi(2):gi(3));
            else
                for ll = 1:length(ai_l)
                    FIRA.analog.data(tt,ai_l(ll)).start_time = gi(2);
                end
            end
        end
        
        %% Maybe plot stuff
        if false
            % check out trial data
            tax = FIRA.analog.data(tt,ai_p).start_time + ...
                (1:FIRA.analog.data(tt,ai_p).length)';
            
            % eye position
            subplot(3,1,1); cla reset; hold on;
            plot(tax, FIRA.analog.data(tt,ai_h).values, 'k-');
            plot(tax, FIRA.analog.data(tt,ai_v).values, 'k-');
            plot(tax([1 end]), [0 0], 'k:');
            plot(fp_times([tt tt], 1), [-5 5], 'g-');
            plot(fp_times([tt tt], 2), [-5 5], 'r-');
            for ss = 1:size(sacs,1)
                plot(FIRA.analog.data(tt,ai_p).start_time+fp_start+sacs(ss,[1 1]), [-5 5], 'b--');
            end
            if isfinite(fix_times(tt,1))
                plot(tax(gi(2):gi(3)), FIRA.analog.data(tt,ai_h).values(gi(2):gi(3)), 'c-');
                plot(tax(gi(2):gi(3)), FIRA.analog.data(tt,ai_v).values(gi(2):gi(3)), 'c-');
            end
            title(sprintf('Correct=%d, Beep=%d', ...
                correct(tt), beep_on(tt)));
            axis([-1000 5000 -20 20])
            
            % pupil diameter
            subplot(3,1,2); cla reset; hold on;
            plot(tax, FIRA.analog.data(tt,2).values, 'k-');
            if isfinite(fix_times(tt,1))
                plot(tax(gi(2):gi(3)), FIRA.analog.data(tt,ai_p).values(gi(2):gi(3)), 'c-');
            end
            axis([-1000 5000 10 30])
            
            % spikes
            subplot(3,1,3); cla reset; hold on;
            for ss = 1:length(FIRA.spikes.channel)
                plot(repmat(FIRA.spikes.data{tt,ss}',2,1), ...
                    repmat([ss-1 ss-0.5]', 1, size(FIRA.spikes.data{tt,ss},1)), 'k-');
            end
            axis([-1000 5000 -0.5 ss]);
            pause(0.1); %
            r = input('next')
        end
    end
end

% for ii = 1:size(eyedat,1)
%     plot(eyedat(ii,:,3));
%     title(sprintf('%d', ii))
%     r = input('next')
%     axis([0 5500 0 25]);
% end

% zscore pupil data ... first simply
eyedat = eyedat(:,1:maxx,:);
pdat   = reshape(eyedat(:,:,3),[],1);
pdat(isfinite(pdat)) = zscore(pdat(isfinite(pdat)));
eyedat(:,:,3) = reshape(pdat,size(eyedat,1),size(eyedat,2));

% now wrt trends across trials -- subtract mean trend (all zscored)
Lpg      = isfinite(eyedat(:,:,3));
Lgood    = any(Lpg,2);
LfixOnly = Lgood & ~isfinite(beep_on) & ~isfinite(stim_on);
pd       = eyedat(:,:,3) - ...
    repmat(nanmean(eyedat(LfixOnly,:,3)), [num_trials,1]);
% Standardize variance??
% pd     = pd ./ repmat(nanstd(eyedat(LnoBeep,:,3)), [num_trials,1]);
pd       = nanrunmean(pd',75)';
pd(~Lpg) = nan;
pd(:, sum(isfinite(eyedat(:,:,3))) < 15) = nan; % need at least 15 trials
eyedat(:,:,4) = pd;

% now loop through and compute slope
for tt = find(Lgood)'
    
    % Compute on non-nan (already interpolated and smoothed) RAW pupil
    % data
    i1  = find(isfinite(eyedat(tt,:,4)),1);
    i2  = find(isfinite(eyedat(tt,:,4)),1,'last');
    pds = filter(filter_c, 1, eyedat(tt,i1:i2,4));
    eyedat(tt,(i1+FILTER_HW):(i2-FILTER_HW-1),5) = ...
        pds(length(filter_c)+1:end);
    %eyedate(tt,i1:i2,5)=movingslope(eyedat(tt,i1:i2,4),150);
end

% ecodes:
%   1 ... fix start time wrt fixation on
%   2 ... fix end time wrt fix start time (fix duration)
%   3 ... reported correcdt
%   4 ... beep on time, wrt to fix start time
%   5 ... trial begin time, wrt to fix start time
%   6 ... trial end time, wrt to fix start time
%   7 ... trial wrt time (cpu clock)
%   8 ... LFP index corresponding to fix start time (coded above)
%   9 ... ELESTM on time (when appropriate), wrt fix start time
siteData{1} = [ ...
    fix_times(Lgood,1), ...
    diff(fix_times(Lgood,:),[],2), ...
    correct(Lgood), ...
    beep_on(Lgood) - fix_times(Lgood,1), ...
    getFIRA_ec(Lgood, 'trial_begin') - fix_times(Lgood,1), ...
    getFIRA_ec(Lgood, 'trial_end') - fix_times(Lgood,1), ...
    getFIRA_ec(Lgood, 'trial_wrt'), ...
    [FIRA.analog.data(Lgood,1).start_time]', ...
    stim_on(Lgood) - fix_times(Lgood,1)];

% analog: trial/sample/1:x, 2:y, 3:pupil ... remember first sample
%    is considered time=0 (i.e., wrt fix start time)
siteData{2} = eyedat(Lgood,:,:);

% spikes, re-coded wrt fix start time
if isfield(FIRA, 'spikes') && ~isempty(FIRA.spikes.data)
    siteData{3} = FIRA.spikes.data(Lgood, :);
end

% LFP OR PULSE-OX
if pulseOxFlag && ~isempty(ai_l)
    siteData{4} = poxdat(:,1:maxx);
elseif ~isempty(ai_l)
    for ll = 1:length(ai_l)
        siteData{4} = cat(2, siteData{4}, ...
            {FIRA.analog.data(Lgood,1).values}');
    end
end

% Pupil events
% pdat_ is mxn matrix of PEAK data
%   rows are events, columns are:
%   1. trial number
%   2. start time of event (wrt fix start time)
%   3. end time of event (wrt fix start time)
%   4. magnitude at start of event (raw z-score)
%   5. magnitude at end of event (raw z-score)
%   6. magnitude at start of event (corrected z-score)
%   7. magnitude at end of event (corrected z-score)
%   8. time of subsequent max slope
%   9. magnitude of subsequent max slope (corrected z/sample)
siteData{5} = getLCP_pupilEvents(eyedat(Lgood,:,:));

% Microsaccades
%   1. trial number
%   2. start time of event (wrt fix start time)
%   3. duration of event (wrt fix start time)
%   4. maximum velocity (deg/ms)
%   5. magnitude of microsaccade event (deg)
%   6. onset time wrt phase of associated pupil event (fraction)
%   7. magnitude of associated pupil event
if ~isempty(siteData{5})
    mdat = [];
    for tt = 1:size(siteData{2},1)
        
        % find microsaccades
        Ltt  = isfinite(siteData{2}(tt,:,1));
        sacs = findMicroSaccades(siteData{2}(tt,Ltt,1)', ...
            siteData{2}(tt,Ltt,2)', find(Ltt,1)-1, false);
        
        if ~isempty(sacs)
            
            % recalibrate start time to t=0
            sacs(:,1) = sacs(:,1)+find(Ltt,1)-1;
            
            % get pupil events from the same trial
            Fp = find(ismember(siteData{5}(:,1), tt));
            
            for ss = 1:size(sacs,1)
                pv = nan;
                pm = nan;
                if ~isempty(Fp)
                    Lpi = sacs(ss,1)>=siteData{5}(Fp,2) & ...
                        sacs(ss,1)<=siteData{5}(Fp,3);
                    if any(Lpi)
                        Fpi = find(Lpi,1);
                        pv = (sacs(ss,1)-siteData{5}(Fp(Fpi),2))./ ...
                            (siteData{5}(Fp(Fpi),3)-siteData{5}(Fp(Fpi),2));
                        pm = siteData{5}(Fp(Fpi),7)-siteData{5}(Fp(Fpi),6);
                    end
                end
                mdat = cat(1, mdat, [tt sacs(ss,:) pv pm]);
            end
        end
    end
    siteData{6} = mdat;
end

%     % add residuals from regressions
%     %   8. Residual magnitude
%     %   9. Residual max slope
%     pdat(:,8:9) = nan;
%     Lp = pdat(:,7)>0;
%     A  = [pdat(:,4) ones(size(pdat,1),1)];
%     Bs = [diff(pdat(:,4:5),[],2) pdat(:,7)];
%     fdat = nans(2,2);
%
%     % save residuals
%     for bb = 1:2
%         Xp = lscov(A(Lp,:), Bs(Lp,bb));
%         pdat(Lp,7+bb) = Bs(Lp,bb) - A(Lp,:)*Xp;
%         fdat(bb,1) = Xp(1);
%
%         Xn = lscov(A(~Lp,:), Bs(~Lp,bb));
%         pdat(~Lp,7+bb) = Bs(~Lp,bb) - A(~Lp,:)*Xn;
%         fdat(bb,2) = Xn(1);
%     end
% siteData{6} = fdat;
 
% save spike/LFP indices
siteData{7} = {FIRA.spikes.id FIRA.analog.name(ai_l)};

% save it
if nargin >= 1 && ~isempty(clean_name)
    save(clean_name, 'siteData');
end

%
%      Ld = siteData{5}(:,7)>0;
%     for tt = 1:size(siteData{1},1)
%         cla reset; hold on
%         plot(siteData{2}(tt,:,3))
%         plot(siteData{2}(tt,:,4),'r')
%         Fed = find( Ld&siteData{5}(:,1)==tt);
%         if any(Fed)
%             plot(repmat(siteData{5}(Fed,2),1,2)', repmat([-1 1]',1,length(Fed)), 'g-')
%         end
%         Fec = find(~Ld&siteData{5}(:,1)==tt);
%         if any(Fec)
%             plot(repmat(siteData{5}(Fec,2),1,2)', repmat([-1 1]',1,length(Fec)), 'k-')
%         end
%         r = input('rr')
%     end

% subplot(2,1,1); cla reset; hold on; set(gca,'FontSize',14);
% plot([-6 6],[0 0],'k:')
% plot(pdat(Lp,4), Bs(Lp,1), 'r.');
% plot(pdat(~Lp,4), Bs(~Lp,1), 'b.');
% ylabel('Pupil change (z-score)')
% xlabel('Initial pupil size (z-score)')
%
% subplot(2,1,2); cla reset; hold on; set(gca,'FontSize',14);
% plot([-6 6],[0 0],'k:')
% plot(pdat(Lp,4), pdat(Lp,8), 'r.');
% plot(pdat(~Lp,4), pdat(~Lp,8), 'b.');
% ylabel('Residual pupil change (z-score)')
% xlabel('Initial pupil size (z-score)')
%
%


