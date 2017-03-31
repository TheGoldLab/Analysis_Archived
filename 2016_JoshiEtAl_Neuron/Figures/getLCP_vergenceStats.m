function getLCP_vergenceStats
% function getLCP_vergenceStats
%
% Check if eyes are moving slowly (or at all) during fixation
%
% Uses siteData:
%  siteData{1}: trialsxcols matrix, cols are:
%   1 ... fix start time wrt fixation on
%   2 ... fix end time wrt fix start time (fix duration)
%   3 ... reported correcdt
%   4 ... beep on time, wrt to fix start time
%   5 ... trial begin time, wrt to fix start time
%   6 ... trial end time, wrt to fix start time
%   7 ... trial wrt time (cpu clock)
%   8 ... LFP index corresponding to fix start time (coded above)
%   9 ... ELESTM on time (when appropriate), wrt fix start time
%
%  siteData{2}: Analog:
%   dim1: trial
%   dim2: sample
%   dim3: 1=x, 2=y, 3=z-pupil, 4=corrected z-pupil, 5=slope
%   [remember first sample is considered time=0 (i.e., wrt fix start time)]
%     = eyedat(Lgood,:,:);
%
%  siteData{3}: spikes, re-coded wrt fix start time
%       first is multi-unit
%       rest are single-unit
%
%  siteData{4}: LFP, re-coded wrt fix start time
%
% SiteData{5}: pupil events... rows are events, columns are:
%   1. trial number
%   2. start time of event (wrt fix start time)
%   3. end time of event (wrt fix start time)
%   4. magnitude at start of event (raw z-score)
%   5. magnitude at end of event (raw z-score)
%   6. magnitude at start of event (corrected z-score)
%   7. magnitude at end of event (corrected z-score)
%   8. time of subsequent max slope
%   9. magnitude of subsequent max slope (corrected z/sample)

%% NON-USTIM
monkeys = {...
    'Oz'      {'LC' 'IC' 'SC'}; ...
    'Cicero'  {'LC' 'IC' 'SC'}; ...
    'Sprout'  {'ACC' 'PCC'}; ...
    'Atticus' {'ACC'}; ...
    'Cheetah' {'PCC'}; ...
    };
nm   =  size(monkeys,1);
vdat = cell(nm,1);
pdat = [];
for mm = 1:nm
    ns = size(monkeys{mm,2},2);
    mdat = [];
    for ss = 1:ns
        [base_dir, fnames] = getLCP_cleanDataDir(monkeys{mm,1}, monkeys{mm,2}{ss});
        for ff = 1:length(fnames)
            disp(sprintf('%s: %s, file %d/%d', ...
                monkeys{mm,1}, monkeys{mm,2}{ss}, ff, length(fnames)))
            load(fullfile(base_dir, fnames{ff}));
            Lnb  = ~isfinite(siteData{1}(:,4)) & ~isfinite(siteData{1}(:,9));
            Fnb  = find(Lnb);
            fdat = nans(length(Fnb),1);
            
            % collect difference in start/end points of fixation per trial
            for tt = 1:length(Fnb)
                ei = find(isfinite(siteData{2}(Fnb(tt),:,1)),1,'last');
                fdat(tt) = nanmean(siteData{2}(Fnb(tt),ei-50:ei,1)) - ...
                    nanmean(siteData{2}(Fnb(tt),1:50,4));
            end
            pdat = cat(1, pdat, signrank(fdat));
            mdat = cat(1, mdat, fdat);
        end
    end
    vdat{mm} = mdat;
end

for vv = 1:size(vdat,1)
    Lg = isfinite(vdat{vv});
    disp([median(vdat{vv}(Lg)) signrank(vdat{vv}(Lg))])
end

%% USTIM
sites = { ...
    'LC', {'Cicero', {'short'}; 'Oz', {'short' 'long'}}, [2 1 6]; ...
    'IC', {'Cicero', {'short'}; 'Oz', {'short' 'long'}}, [2 1 6]; ...
    'SC', {'Cicero', {'short'}; 'Oz', {'short'}},        [2 1 4]};
num_sites = size(sites,1);
vdat = cell(nm,2);
pdat = [];

for ss = 1:num_sites
    for mm = 1:size(sites{ss,2},1)
        for tt = 1:size(sites{ss,2}{mm,2},2)
                
            % get directory
            cdir = fullfile(dirnames('local'), 'Data', 'Projects', ...
                '2013_LCPupil', 'Data', 'uStim', ...
                sites{ss,2}{mm,1}, sites{ss,1}, sites{ss,2}{mm,2}{tt}, 'clean');
            files = dir(fullfile(cdir, '*.mat'));
            for ff = 1:size(files,1)
                    
                % load file
                disp(sprintf('site = %s, monkey = %s, type = %s, %d/%d files', ...
                    sites{ss,1}, sites{ss,2}{mm,1}, ...
                    sites{ss,2}{mm,2}{tt}, ff, size(files,1)))
                load(fullfile(cdir, files(ff).name));
                
                % all uStim trials with at least 1 sec baseline
                Fstim = find(siteData{1}(:,9)>=1000);
                
                % collect difference in start/end points of fixation per trial
                for ii = 1:length(Fstim)
                    ei = find(isfinite(siteData{2}(Fstim(ii),:,1)),1,'last');
                    fdat(ii) = nanmean(siteData{2}(Fstim(ii),ei-50:ei,1)) - ...
                        nanmean(siteData{2}(Fstim(ii),1:50,4));
                end
                pdat = cat(1, pdat, signrank(fdat));
                mdat = cat(1, mdat, fdat);
            end
        end
        vdat{ss,mm} = mdat;
    end
end

for ss = 1:size(vdat,1)
    for mm = 1:size(vdat,2)
        Lg = isfinite(vdat{ss,mm});
    disp([median(vdat{ss,mm}(Lg)) signrank(vdat{ss,mm}(Lg))])
    end    
end

