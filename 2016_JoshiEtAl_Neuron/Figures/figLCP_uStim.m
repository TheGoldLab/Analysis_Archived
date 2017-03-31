function figLCP_uStim(num, collectData)
% function figLCP_uStim(num, collectData)
%
%  siteData{1}: trialsxcols matrix, cols are:
%   1 ... fix start time wrt fixation on
%   2 ... fix end time wrt fix start time (fix duration)
%   3 ... reported correcdt
%   4 ... beep on time, wrt to fix start time
%   5 ... trial begin time, wrt to fix start time
%   6 ... trial end time, wrt to fix start time
%   7 ... trial wrt time (cpu clock)
%   8 ... LFP index corresponding to fix start time (coded above)
%
%  siteData{2}: Analog:
%   dim1: trial
%   dim2: sample
%   dim3: 1 = x, 2 = y, 3 = pupil, 4 = corrected z-pupil, 5 = pupil slope
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

%% parse args
if nargin < 1 || isempty(num)
    num = 7;
end

if nargin < 2 || isempty(collectData)
    collectData = false;
end

%% Set up figure
wid  = 11.6; % total width
ht   = 3;
cols = {3,3};
axs  = getPLOT_axes(num, wid, ht, cols, 1.5, .5, [], 'Joshi et al', true);
set(axs,'Units','normalized');

%% Stuff we need
% site, {monkey, type}, example_data (monkey type file)
sites = { ...
    'LC', {'Cicero', {'short'}; 'Oz', {'short' 'long'}}, [2 1 6]; ...
    'IC', {'Cicero', {'short'}; 'Oz', {'short' 'long'}}, [2 1 6]; ...
    'SC', {'Cicero', {'short'}; 'Oz', {'short'}},        [2 1 4]};
num_sites = size(sites,1);

% time axis for uStim-aligned PD
tmin   = -1000;
tmax   = 1000;
pax    = (tmin:tmax)';
npax   = length(pax);

%% Possibly get data
if collectData
    
    udat = cell(num_sites, 2); % pupil/sites
    
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
                    
                    %                     tax = 1:size(siteData{2},2);
                    %                     subplot(2,1,1); cla reset; hold on;
                    %                     subplot(2,1,2); cla reset; hold on;
                    %                     for rr = 1:length(Fstim)
                    %                         ntax = tax-siteData{1}(Fstim(rr),9);
                    %                         for ii = 1:2
                    %                             if ii == 1
                    %                                 yy = siteData{2}(Fstim(rr),:,4);
                    %                             else
                    %                                 yy = sqrt( ...
                    %                                     siteData{2}(Fstim(rr),:,1).^2 + ...
                    %                                     siteData{2}(Fstim(rr),:,2).^2);
                    %                             end
                    %                             subplot(2,1,ii);
                    %                             plot(ntax, yy-nanmean(yy(ntax>-100&ntax<0)), 'k-');
                    %                         end
                    %                     end
                    %                     for ii = 1:2
                    %                         subplot(2,1,ii);
                    %                         plot([0 0], [-4 4], 'r:');
                    %                         plot([1000 1000], [-4 4], 'b:');
                    %                         axis([-1500 1500 -4 4]);
                    %                     end
                    %                     r = input('next');
                    %                 end
                    
                    % fudge factor to account for the fact that the ustim
                    % ecode was dropped at the END of ustim
                    if strcmp(sites{ss,2}{mm,2}{tt}, 'long')
                        uoff = 400;
                    else
                        if strcmp(sites{ss,2}{mm,1}, 'Cicero') || ...
                                ~strcmp(files(ff).name, 'Oz_061013_d6_435_FixnStim.mat')
                            uoff = 100;
                        else
                            uoff = 50;
                        end
                    end
                    
                    if ~isempty(Fstim)
                        ntr   = length(Fstim);
                        lpd   = size(siteData{2},2);
                        pdat  = nans(ntr, npax, 2);
                        pmdat = nans(ntr, 2);
                        
                        for rr = 1:ntr
                            
                            % PD wrt ustim
                            stm   = siteData{1}(Fstim(rr),9)-uoff;
                            lastp = find(isfinite(siteData{2}(Fstim(rr),:,5)),1,'last');
                            bpax  = pax+round(stm);
                            Lp    = bpax>=0&bpax<lpd;
                            if any(Lp) && (lastp >= stm+700)
                                
                                pdat(rr,Lp,1) = siteData{2}(Fstim(rr),bpax(Lp)+1,4);
                                pdat(rr,Lp,2) = siteData{2}(Fstim(rr),bpax(Lp)+1,5);
                                
                                % PD stats:
                                %   1. baseline (previous sec)
                                %   2. response (subsequent sec)
                                pmdat(rr,1) = nanmean(pdat(rr,pax<0,2));
                                pmdat(rr,2) = nanmean(pdat(rr,pax>0,2));
                            end
                        end
                        
                        % save PD data
                        udat{ss,1} = cat(1, udat{ss,1}, pdat);
                        udat{ss,2} = cat(1, udat{ss,2}, ...
                            cat(2, repmat([mm tt ff], ntr, 1), pmdat));
                    end
                end
            end
        end
    end
    
    % save data to file
    FS_saveProjectFile('2013_LCPupil', 'uStim', udat);
else
    % load data from file
    udat = FS_loadProjectFile('2013_LCPupil', 'uStim');
end

Lmp = pax>0;
pv  = 0.05;
for ss = 1:num_sites
    
    %% Example site, mean±sem PD response
    axes(axs(ss)); cla reset; hold on;
    [co, sy] = getLCP_colorsAndSymbol(sites{ss,2}(sites{ss,3}(1)), ...
        sites{ss,1});
    plot(pax([1 end]), [0 0], 'k:');
    plot([0 0], [-4 4], 'k:');
    Lmtf = udat{ss,2}(:,1)==sites{ss,3}(1) & ... % monkey
        udat{ss,2}(:,2)==sites{ss,3}(2) & ...    % type
        udat{ss,2}(:,3)==sites{ss,3}(3);         % file
    nms = nanmean(udat{ss,1}(Lmtf,:,1))';
    ses = nanse(udat{ss,1}(Lmtf,:,1))';
    Lg  = isfinite(nms) & isfinite(ses);% & sum(isfinite(udat{ss,1}))'>5;
    h   = patch([pax(Lg); flipud(pax(Lg))], [nms(Lg)+ses(Lg); ...
        flipud(nms(Lg)-ses(Lg))], co{2});
    set(h,'LineStyle','none');
    plot(pax(Lg), nms(Lg), '-', 'Color', co{1});
    axis([tmin tmax -0.7 3.5]);
    title(sites{ss,1})
    if ss == 1
        ylabel('Pupil diameter (z-score)')
    end
    
    %% Scatters of peak amplitude vs time, per site
    axes(axs(ss+num_sites)); cla reset; hold on;
    plot([-5 5], [0 0],  'k:');
    plot([0 0],  [-5 5], 'k:');
    plot([-5 5], [-5 5], 'k:');
    Lg = isfinite(udat{ss,2}(:,4)) & isfinite(udat{ss,2}(:,5));
    sdat = [];
    for mm = 1:size(sites{ss,2},1) % each monkey
        
        % get color, symbol
        [co, sy] = getLCP_colorsAndSymbol(sites{ss,2}(mm), sites{ss,1});
        
        Lm = udat{ss,2}(:,1)==mm;
        
        % show fraction of evoked saccades
        disp(sprintf('%s %s: %d/%d (%.2f pct) evoked', ...
            sites{ss,1}, sites{ss,2}{mm,1}, ...
            sum(Lm&~Lg), sum(Lm), sum(Lm&~Lg)./sum(Lm)));
        
        for tt = 1:size(sites{ss,2}{mm,2},2) % each type
            Lmt = Lg & Lm & udat{ss,2}(:,2)==tt;
            for ff = 1:max(udat{ss,2}(Lmt,3))
                Lf = Lmt & udat{ss,2}(:,3)==ff;
               
                if sum(Lf)>4
                    
                    % get PD at mean peak, plot index vs mean value
                    meanp = nanmean(udat{ss,1}(Lf,:,2));
                    mvs   = meanp(Lmp);
                    maxi  = find(mvs==max(mvs),1);
                    pis   = (-10:10)+maxi-tmin+1;
                    pis   = pis(pis>=1&pis<=npax);
                    fdat  = nans(sum(Lf), 1);
                    Ff    = find(Lf);
                    for rr = 1:sum(Lf)
                        fdat(rr) = nanmean(udat{ss,1}(Ff(rr),pis,2)-...
                            udat{ss,2}(Ff(rr),4));
                    end
                    mn = nanmean(fdat);
                    %disp([ss mm tt ff maxi])
                    se = nanse(fdat);
                    plot([maxi maxi], [mn-se mn+se], 'k-');
                    h=plot(maxi,mn,sy,'Color',co{1});
                    %                     if tt==2
                    %                         set(h, 'MarkerSize', 12);
                    %                     end
                    
                    if signrank(fdat)<pv
                        sdat = cat(1, sdat, [maxi mn 1]);
                        set(h, 'MarkerFaceColor', co{1});
                    else
                        sdat = cat(1, sdat, [maxi mn 0]);
                    end
                end
            end
        end
        axis([0 tmax 0 .02])
    end
    disp(sprintf('%d sessions, %d sig, %d-%d range', ...
        size(sdat,1), sum(sdat(:,3)), min(sdat(sdat(:,3)==1,1)), ...
        max(sdat(sdat(:,3)==1,1))))
    if ss == 1
        ylabel('Change in pupil diameter (z-score/ms)')
    elseif ss == 2
        xlabel('Time of peak change (ms)')
    end
end

setPLOT_panelLabel(axs(1), 1);
setPLOT_panelLabel(axs(4), 2);


