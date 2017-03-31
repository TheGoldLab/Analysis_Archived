% script for collecting pupil/pulse-ox data from Yin
% pulse-ox should be analog channel 13 in nex file

base_dir =  fullfile(dirnames('local'), 'Data', 'Projects', '2013_LCPupil', ...
    'Data', 'Recording', 'Cheetah', 'PCC');

file_list = { ...
    'CHfix20150226pccs', 'unknown2'; ...
    'CHfix20150302pccs', 'unknown2'; ...
    'CHfix20150304pccm', 'unknown2'; ...
    'CHfix20150305pccs', 'unknown2'; ...
    };

global FIRA

%% load mat file
for ff = 1:size(file_list,1)

    % make FIRA, without spikes
    bNex(fullfile(base_dir, 'nex', [file_list{ff,1} '-vf.nex']), ...
        'spmRKLCPupil', [], [], 'all', 0, 1, [], []); % [4,14:16]
    
    % keep only these as analog channels
    Lgood = ismember(FIRA.analog.name, ...
        cat(2, file_list{ff,2}, {'pupil_diameter', 'horiz_eye', 'vert_eye'}));
    FIRA.analog.name         = FIRA.analog.name(Lgood);
    FIRA.analog.acquire_rate = FIRA.analog.acquire_rate(Lgood);
    FIRA.analog.store_rate   = FIRA.analog.store_rate(Lgood);
    FIRA.analog.data         = FIRA.analog.data(:,Lgood);
    
    % save it!
    saveFIRA(fullfile(base_dir, 'pulseOx', [file_list{ff,1} '_OUT']));

    % Make the clean file -- remember the pulseOx is in the LFP channel
    %   (siteData{4})
    cleanLCP_FIRA(fullfile(base_dir, 'pulseOx', [file_list{ff,1} '_clean']), true);    
end

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

%% compare distributions of pupil/pulse-ox events
pdat  = [];
hdat  = [];
for ff = 1:size(file_list,1)
    load(fullfile(base_dir, 'pulseOx', [file_list{ff,1} '_clean']))
    disp(ff)
    Fgood = find(isnan(siteData{1}(:,4)) & isnan(siteData{1}(:,9)));
    
    for tt = 1:length(Fgood)        
        tr   = Fgood(tt);
        
        % distrubtion of full pupil event durations (pos+neg)
        Le = siteData{5}(:,1)==tr;
        if sum(Le) > 1
            evts = siteData{5}(Le,[2 3 9]); % start/end/mag
            while size(evts,1)>1
                if evts(1,2) == evts(2,1) && sign(evts(1,3)) == -sign(evts(2,3))
                    pdat = cat(1, pdat, evts(2,2)-evts(1,1));
                    evts(1:2,:) = [];
                else
                    evts(1,:) = [];
                end
            end
        end
        
        % pulse periods        
        Lg = siteData{4}(tr,:)~=0;
        if sum(Lg) > 200
            vals = sgolayfilt(siteData{4}(tr,Lg), 1, 151);
            [peaks,locs] = findpeaks(vals, 'MinPeakDistance', 300, 'MinPeakProminence', 2);
            if length(locs)>1
                hdat = cat(1, hdat, diff(locs)');
            end
        end
    end
end

wid        = 8.5; % total width
ht         = 3;
cols       = {1,1,2,2};
[axs,fig_] = getPLOT_axes(0, wid, ht, cols, 2, 2, [], 'Joshi et al', true);
set(axs,'Units','normalized');

%% Example PD, HR
ex1 = 24;
load(fullfile(base_dir, 'pulseOx', [file_list{1,1} '_clean']))
Fgood = find(isnan(siteData{1}(:,4)) & isnan(siteData{1}(:,9)));
x  = siteData{2}(Fgood(ex1),:,4);
y  = siteData{4}(Fgood(ex1),:);
Lg = isfinite(x) & isfinite(y) & y~=0;
axes(axs(1)); cla reset; hold on;
plot(x(Lg)-nanmean(x));
ylabel('z-score PD')
axes(axs(2)); cla reset; hold on;
plot(y(Lg)-nanmean(y(Lg)));
ylabel('Pulse-ox')
xlabel('Time (ms)')

%% period histograms

bins = 0:50:2000;
Np = hist(pdat(pdat<2000),bins);
Nh = hist(hdat(hdat<2000),bins);
axes(axs(3)); cla reset; hold on;
H = bar(bins,Np);
set(H,'FaceColor','k')
axis([0 2000 0 450])
ylabel('Count')
xlabel('Duration (ms)');

axes(axs(5)); cla reset; hold on;
H = bar(bins,Nh);
set(H,'FaceColor','k')
axis([0 2000 0 5000])

disp(sprintf('95 pct CI hr period = %d ? %d', ...
    prctile(hdat(hdat<900),2.5), prctile(hdat(hdat<900),97.5)))

%% find phase of pulse-ox peaks wrt pupil events
pbins = 0:0.1:1;
npbins = length(pbins)-1;
pdat   = zeros(npbins, 1);
for ff = 1:size(file_list,1)
    load(fullfile(base_dir, 'pulseOx', [file_list{ff,1} '_clean']))
    disp(ff)
    Fgood  = find(isnan(siteData{1}(:,4)) & isnan(siteData{1}(:,9)));
    Fevts  = find(ismember(siteData{5}(:,1), Fgood));
    for ee = 1:length(Fevts)
        tr = siteData{5}(Fevts(ee),1);
        Lg   = siteData{4}(tr,:)~=0;
        if sum(Lg) > 200
            vals = sgolayfilt(siteData{4}(tr,Lg), 1, 151);
            [peaks,locs] = findpeaks(vals, 'MinPeakDistance', 300, 'MinPeakProminence', 2);
            if ~isempty(peaks)
                estart = round(siteData{5}(Fevts(ee),2));
                eend   = round(siteData{5}(Fevts(ee),3));
                efrac  = (locs-estart)./(eend-estart);
                Le     = efrac>=0 & efrac<=1;
                if any(Le)
                    for pp = find(Le)
                        ind = efrac(pp)>=pbins(1:end-1)&efrac(pp)<=pbins(2:end);
                        pdat(ind) = pdat(ind) + 1;
                    end
                end
            end
        end
    end
end
axes(axs(4)); cla reset; hold on;
plot(pbins(1:end-1), pdat./sum(pdat));
axis([0 1 0 1])
xlabel('Phase (fraction)')
ylabel('Fraction of events')

MAXLAG = 500;
ccdat   = [];
cdat    = [];
for ff = 1:size(file_list,1)
    load(fullfile(base_dir, 'pulseOx', [file_list{ff,1} '_clean']))
    
    Fgood  = find(isnan(siteData{1}(:,4)) & isnan(siteData{1}(:,9)));
    ntr    = length(Fgood);
    tccdat = nans(ntr, 2*MAXLAG+1);
    tcdat  = nans(ntr, 1);
    for tt = 1:ntr
        x  = siteData{2}(Fgood(tt),:,4);
        y  = siteData{4}(Fgood(tt),:);
        Lg = isfinite(x) & isfinite(y) & y~=0;
        if sum(Lg) > 1024
            tccdat(tt,:) = xcorr(x(Lg)-nanmean(x), y(Lg)-nanmean(y(Lg)), MAXLAG, 'coeff');
            tcdat(tt) = corr((x(Lg)-nanmean(x))', (y(Lg)-nanmean(y(Lg)))', 'type', 'Spearman');
        end
    end
    ccdat = cat(1, ccdat, tccdat);
    cdat  = cat(1, cdat, tcdat);
end

cdat = cdat(isfinite(cdat));
disp(sprintf('median=%.2f [%.2f %.2f], p=%.2f', median(cdat), ...
    prctile(cdat,25), prctile(cdat,75), signrank(cdat)))

clf; hold on;
tax = -500:500;
plot(tax, nanmean(cdat))
for tt = 1:size(cdat,2)
    if signrank(cdat(:,tt))<0.05
        plot(tax(tt), nanmean(cdat(:,tt)), 'r.')
    end
end

    r = input('next')
end
 

MAXLAG = 500;
%CSIZE   = 4096;
for ff = 1:size(file_list,1)
    load(fullfile(base_dir, 'pulseOx', [file_list{ff,1} '_clean']))
    
    Fgood = find(isnan(siteData{1}(:,4)) & isnan(siteData{1}(:,9)));
    ntr   = length(Fgood);
    cdat  = nans(ntr, 2*MAXLAG+1);
    %cdat  = nans(ntr, CSIZE/2+1);
    for tt = 1:ntr
        x  = siteData{2}(Fgood(tt),:,4);
        y  = siteData{4}(Fgood(tt),:);
        Lg = isfinite(x) & isfinite(y) & y~=0;
        if sum(Lg) > 1024
            %[cdat(tt,:),F] = mscohere(x(Lg)-nanmean(x),y(Lg)-nanmean(y(Lg)),...
            %    hanning(1024),[],CSIZE,1000);
            cdat(tt,:) = xcorr(x(Lg)-nanmean(x), y(Lg)-nanmean(y(Lg)), MAXLAG, 'coeff');
            subplot(3,1,1); cla reset; hold on;
            plot(x(Lg)-nanmean(x));
            subplot(3,1,2); cla reset; hold on;
            vals = y(Lg)-nanmean(y(Lg));
            valf = sgolayfilt(vals, 1, 151);
            [peaks,locs] = findpeaks(valf, 'MinPeakDistance', 300, 'MinPeakProminence', 2);
            plot([0 length(vals)], [0 0], 'k:');
            plot(vals);
            plot(valf,'r');
            for pp = 1:length(locs)
                plot(locs(pp), peaks(pp), 'go');
            end
            subplot(3,1,3); cla reset; hold on;
            plot(-500:500, cdat(tt,:))
            r = input('next')
        end
    end
    plot(F, nanmean(cdat));
    %plot(-500:500, nanmean(cdat))
    r = input('next')
end
 
          
          